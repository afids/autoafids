# === DATASET-LEVEL REGISTRATION QC SUMMARY ===
"""
Aggregates per-subject regqc CSVs into a single interactive HTML report.

Columns expected in each CSV:
    AFID, dx (mm), dy (mm), dz (mm), ED (mm)

Output: one self-contained HTML file with:
  - Dataset stats tiles
  - Subject summary table with hyperlinks to individual reports
  - Heatmap: subjects x AFIDs coloured by ED
  - Boxplot: ED distribution per AFID across subjects
  - Bar chart: mean ED per subject (sorted, worst first)
  - 3D scatter: all subjects' AFIDs in MNI space vs ground truth
"""

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Template

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_bids_entities(filename: str) -> dict:
    """Extract BIDS key-value entities from a filename."""
    entities = {}
    for match in re.finditer(r"([a-zA-Z]+)-([a-zA-Z0-9]+)", filename):
        entities[match.group(1)] = match.group(2)
    return entities


def _subject_label(entities: dict) -> str:
    """Build a compact human-readable subject label from BIDS entities."""
    parts = []
    for key in ("sub", "ses", "acq", "run"):
        if key in entities:
            parts.append(f"{key}-{entities[key]}")
    return "_".join(parts) if parts else "unknown"


def load_all_csvs(csv_paths: list, html_paths: list) -> pd.DataFrame:
    """
    Read every per-subject CSV and concatenate into a long-form DataFrame.

    Returns a DataFrame with columns:
        subject_label, AFID, dx (mm), dy (mm), dz (mm), ED (mm),
        html_relpath  (relative path to per-subject HTML)
    """
    records = []
    # Build a mapping from subject_label → html_path
    html_map = {}
    for hp in html_paths:
        ents = _parse_bids_entities(Path(hp).name)
        label = _subject_label(ents)
        html_map[label] = hp

    for cp in csv_paths:
        ents = _parse_bids_entities(Path(cp).name)
        label = _subject_label(ents)
        try:
            df = pd.read_csv(cp)
        except Exception:
            continue
        df["subject_label"] = label
        df["html_path"] = html_map.get(label, "")
        records.append(df)

    if not records:
        raise ValueError("No valid regqc CSVs found.")
    return pd.concat(records, ignore_index=True)


def load_fcsv(path: str) -> tuple:
    """Load a Slicer FCSV file (skip 2 header lines).

    Returns
    -------
    coords : np.ndarray of shape (N, 3)
    labels : list of AFID labels (ints or strings)
    """
    df = pd.read_csv(path, skiprows=2)
    coords = df[["x", "y", "z"]].to_numpy(dtype=float)
    labels = df["label"].tolist()
    return coords, labels


def load_all_fcsvs(
    fcsv_paths: list,
) -> list:
    """Load all per-subject FCSVs.

    Returns list of dicts: {subject_label, coords, labels}
    """
    result = []
    for fp in fcsv_paths:
        ents = _parse_bids_entities(Path(fp).name)
        label = _subject_label(ents)
        try:
            coords, afid_labels = load_fcsv(fp)
        except Exception:
            continue
        result.append(
            {
                "subject_label": label,
                "coords": coords,
                "labels": afid_labels,
            }
        )
    return result


def compute_error_decomposition(
    subject_fcsvs: list,
    gt_coords: np.ndarray,
    gt_labels: list,
) -> pd.DataFrame:
    """Decompose per-AFID error into systematic bias and random scatter.

    For each AFID, across all subjects:
      - Signed error vector = subject_coord - gt_coord
      - Bias = mean signed error vector (systematic component)
      - Scatter = std of signed error vectors (random component)
      - Bias magnitude |bias|, scatter magnitude (RMS of stds)
      - Bias-to-total ratio = |bias| / (|bias| + scatter_rms)

    Parameters
    ----------
    subject_fcsvs : list of dicts with {subject_label, coords (N,3), labels}
    gt_coords : (N, 3) ground-truth MNI coordinates
    gt_labels : list of AFID labels

    Returns
    -------
    DataFrame with one row per AFID and columns:
        label, bias_x, bias_y, bias_z, bias_mag,
        scatter_x, scatter_y, scatter_z, scatter_rms,
        total_error, bias_ratio, n_subjects,
        gt_x, gt_y, gt_z
    """
    n_afids = len(gt_labels)
    # Collect signed error vectors: shape (n_subjects, n_afids, 3)
    all_errors = []
    for sub in subject_fcsvs:
        if sub["coords"].shape[0] != n_afids:
            continue  # skip if AFID count mismatch
        signed_err = sub["coords"] - gt_coords  # (N, 3)
        all_errors.append(signed_err)

    if not all_errors:
        return pd.DataFrame()

    errors = np.stack(all_errors, axis=0)  # (S, N, 3)

    # Per-AFID statistics
    bias = errors.mean(axis=0)  # (N, 3) — systematic
    scatter = (
        errors.std(axis=0, ddof=1) if errors.shape[0] > 1 else np.zeros_like(bias)
    )  # (N, 3) — random

    bias_mag = np.linalg.norm(bias, axis=1)  # (N,)
    scatter_rms = np.sqrt((scatter**2).mean(axis=1))  # (N,)
    total_error = np.linalg.norm(errors, axis=2).mean(axis=0)  # (N,)
    bias_ratio = bias_mag / (bias_mag + scatter_rms + 1e-9)

    return pd.DataFrame(
        {
            "label": gt_labels,
            "bias_x": bias[:, 0],
            "bias_y": bias[:, 1],
            "bias_z": bias[:, 2],
            "bias_mag": bias_mag,
            "scatter_x": scatter[:, 0],
            "scatter_y": scatter[:, 1],
            "scatter_z": scatter[:, 2],
            "scatter_rms": scatter_rms,
            "total_error": total_error,
            "bias_ratio": bias_ratio,
            "n_subjects": errors.shape[0],
            "gt_x": gt_coords[:, 0],
            "gt_y": gt_coords[:, 1],
            "gt_z": gt_coords[:, 2],
        }
    )


# ---------------------------------------------------------------------------
# Plotly chart builders
# ---------------------------------------------------------------------------


def make_heatmap(long_df: pd.DataFrame) -> str:
    """
    Subjects (rows) x AFIDs (cols), colour = ED (mm).
    Rows sorted by mean ED descending (worst on top).
    """
    pivot = long_df.pivot_table(index="subject_label", columns="AFID", values="ED (mm)")
    row_order = pivot.mean(axis=1).sort_values(ascending=False).index
    pivot = pivot.loc[row_order]

    fig = go.Figure(
        data=go.Heatmap(
            z=pivot.values,
            x=[f"AFID {c}" for c in pivot.columns],
            y=list(pivot.index),
            colorscale=[[0.0, "rgb(255,255,255)"], [1.0, "rgb(220,0,0)"]],
            colorbar=dict(title="ED (mm)"),
            hovertemplate=("%{y}<br>AFID %{x}<br>ED = %{z:.2f} mm<extra></extra>"),
        )
    )
    fig.update_layout(
        title="Registration Error (ED mm) — Subjects x AFIDs",
        height=max(300, 30 * len(pivot)),
        margin=dict(t=50, b=40, l=160, r=20),
        template="plotly_white",
        xaxis=dict(tickangle=-45),
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=True)


def make_afid_boxplot(long_df: pd.DataFrame) -> str:
    """Box-plot of ED distribution across subjects, one box per AFID."""
    afid_ids = sorted(long_df["AFID"].unique())
    n_subjects = long_df["subject_label"].nunique()
    # Show all points for small datasets, outliers-only for large
    show_points = "all" if n_subjects <= 30 else "outliers"
    fig = go.Figure()
    for afid in afid_ids:
        vals = long_df.loc[long_df["AFID"] == afid, "ED (mm)"]
        fig.add_trace(
            go.Box(
                y=vals,
                name=f"AFID {afid}",
                boxpoints=show_points,
                jitter=0.4,
                pointpos=0,
                marker=dict(size=4, opacity=0.5),
                line=dict(width=1.5),
            )
        )
    fig.update_layout(
        title="ED Distribution per AFID (across subjects)",
        yaxis_title="ED (mm)",
        height=500,
        margin=dict(t=50, b=80, l=60, r=20),
        template="plotly_white",
        showlegend=False,
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


def make_subject_bar(subject_stats: pd.DataFrame) -> str:
    """
    Horizontal bar chart of mean ED per subject, sorted worst → best.
    Bars coloured by value (green → red).
    """
    df = subject_stats.sort_values("mean_ED", ascending=True)
    colorscale_vals = df["mean_ED"].values
    norm = (colorscale_vals - colorscale_vals.min()) / (np.ptp(colorscale_vals) + 1e-9)
    colours = [f"rgb({int(255*v)},{int(255*(1-v))},80)" for v in norm]
    fig = go.Figure(
        go.Bar(
            x=df["mean_ED"],
            y=df["subject_label"],
            orientation="h",
            marker_color=colours,
            hovertemplate="%{y}<br>Mean ED = %{x:.2f} mm<extra></extra>",
        )
    )
    fig.update_layout(
        title="Mean ED per Subject (sorted best → worst)",
        xaxis_title="Mean ED (mm)",
        height=max(300, 28 * len(df)),
        margin=dict(t=50, b=40, l=160, r=20),
        template="plotly_white",
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


def make_3d_landmark_scatter(
    subject_fcsvs: list,
    gt_coords: np.ndarray = None,
    gt_labels: list | None = None,
) -> str:
    """Interactive 3D scatter of all subjects' AFIDs in MNI space.

    Parameters
    ----------
    subject_fcsvs : list of dicts with keys subject_label, coords, labels
    gt_coords : optional (N, 3) array of ground-truth MNI coordinates
    gt_labels : optional list of AFID labels for ground truth
    """
    fig = go.Figure()

    # Ground truth (template) as large green diamonds
    if gt_coords is not None:
        hover = [f"AFID {i}" for i in (gt_labels or range(1, len(gt_coords) + 1))]
        fig.add_trace(
            go.Scatter3d(
                x=gt_coords[:, 0],
                y=gt_coords[:, 1],
                z=gt_coords[:, 2],
                mode="markers",
                name="MNI Ground Truth",
                marker=dict(
                    size=6,
                    color="green",
                    opacity=0.9,
                    symbol="diamond",
                ),
                text=hover,
                hovertemplate="%{text}<br>x=%{x:.1f} y=%{y:.1f} z=%{z:.1f}"
                "<extra>MNI Ground Truth</extra>",
            )
        )

    # Per-subject data — merge into a single trace for scalability
    # (avoids creating N traces which kills performance at N>30)
    all_x, all_y, all_z, all_hover, all_colour = [], [], [], [], []
    colours = [
        "#e6194b",
        "#3cb44b",
        "#4363d8",
        "#f58231",
        "#911eb4",
        "#42d4f4",
        "#f032e6",
        "#bfef45",
        "#fabed4",
        "#469990",
        "#dcbeff",
        "#9A6324",
        "#fffac8",
        "#800000",
        "#aaffc3",
        "#808000",
        "#ffd8b1",
        "#000075",
        "#a9a9a9",
        "#000000",
    ]
    for i, sub in enumerate(subject_fcsvs):
        colour = colours[i % len(colours)]
        n = sub["coords"].shape[0]
        all_x.extend(sub["coords"][:, 0].tolist())
        all_y.extend(sub["coords"][:, 1].tolist())
        all_z.extend(sub["coords"][:, 2].tolist())
        all_hover.extend(f"{sub['subject_label']} — AFID {i}" for i in sub["labels"])
        all_colour.extend([colour] * n)

    fig.add_trace(
        go.Scatter3d(
            x=all_x,
            y=all_y,
            z=all_z,
            mode="markers",
            name=f"Subjects (n={len(subject_fcsvs)})",
            marker=dict(
                size=3,
                color=all_colour,
                opacity=0.5,
            ),
            text=all_hover,
            hovertemplate="%{text}<br>x=%{x:.1f} y=%{y:.1f} z=%{z:.1f}"
            "<extra></extra>",
        )
    )

    fig.update_layout(
        scene=dict(
            xaxis_title="X (mm)",
            yaxis_title="Y (mm)",
            zaxis_title="Z (mm)",
            bgcolor="rgb(250,250,250)",
        ),
        height=700,
        width=900,
        margin=dict(t=50, b=30, l=20, r=20),
        template="plotly_white",
        title="AFID Landmarks in MNI Space (all subjects)",
        legend=dict(
            x=0.82,
            y=0.95,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="#ccc",
            borderwidth=1,
        ),
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


def make_bias_vector_3d(
    decomp_df: pd.DataFrame,
) -> str:
    """3D cone plot showing bias direction and magnitude at each AFID.

    Cones are anchored at the ground-truth AFID location, pointing in the
    direction of the mean signed error (bias).  Colour encodes the
    bias-to-total ratio (red = mostly systematic, blue = mostly random).
    """
    df = decomp_df.copy()
    # Scale factor for cone size visibility
    max_bias = df["bias_mag"].max() + 1e-9

    fig = go.Figure()

    # Ground-truth markers
    fig.add_trace(
        go.Scatter3d(
            x=df["gt_x"],
            y=df["gt_y"],
            z=df["gt_z"],
            mode="markers",
            name="AFID Location (MNI)",
            marker=dict(size=4, color="#888", opacity=0.6, symbol="diamond"),
            text=[f"AFID {i}" for i in df["label"]],
            hovertemplate="%{text}<extra>Ground Truth</extra>",
        )
    )

    # Bias cones — Plotly Cone trace
    fig.add_trace(
        go.Cone(
            x=df["gt_x"],
            y=df["gt_y"],
            z=df["gt_z"],
            u=df["bias_x"],
            v=df["bias_y"],
            w=df["bias_z"],
            sizemode="absolute",
            sizeref=max_bias * 0.8,
            colorscale=[
                [0.0, "rgb(49,130,189)"],  # blue = mostly random
                [0.5, "rgb(255,255,100)"],  # yellow
                [1.0, "rgb(215,48,39)"],  # red = mostly systematic
            ],
            cmin=0,
            cmax=1,
            colorbar=dict(
                title="Bias Ratio",
                tickvals=[0, 0.5, 1],
                ticktext=["Random", "Mixed", "Systematic"],
                len=0.6,
                y=0.5,
            ),
            name="Bias Vector",
            text=[
                f"AFID {row.label}<br>"
                f"Bias: {row.bias_mag:.2f} mm<br>"
                f"Scatter: {row.scatter_rms:.2f} mm<br>"
                f"Ratio: {row.bias_ratio:.0%} systematic"
                for _, row in df.iterrows()
            ],
            hoverinfo="text",
        )
    )

    # Endpoint markers (bias tip) for error lines
    tip_x = df["gt_x"] + df["bias_x"]
    tip_y = df["gt_y"] + df["bias_y"]
    tip_z = df["gt_z"] + df["bias_z"]
    fig.add_trace(
        go.Scatter3d(
            x=tip_x,
            y=tip_y,
            z=tip_z,
            mode="markers",
            name="Bias Endpoint",
            marker=dict(size=2, color="red", opacity=0.5),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    fig.update_layout(
        scene=dict(
            xaxis_title="X (mm)",
            yaxis_title="Y (mm)",
            zaxis_title="Z (mm)",
            bgcolor="rgb(250,250,250)",
            aspectmode="data",
        ),
        height=700,
        width=950,
        margin=dict(t=50, b=30, l=20, r=20),
        template="plotly_white",
        title="Systematic Bias Vectors at Each AFID",
        legend=dict(x=0.01, y=0.99),
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


def make_decomposition_bar(decomp_df: pd.DataFrame) -> str:
    """Grouped bar chart: bias magnitude vs scatter (RMS) per AFID."""
    df = decomp_df.sort_values("label")
    labels = [f"AFID {i}" for i in df["label"]]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=labels,
            y=df["bias_mag"],
            name="Systematic Bias |μ|",
            marker_color="rgb(215,48,39)",
            hovertemplate="%{x}<br>Bias = %{y:.2f} mm<extra></extra>",
        )
    )
    fig.add_trace(
        go.Bar(
            x=labels,
            y=df["scatter_rms"],
            name="Random Scatter o",
            marker_color="rgb(49,130,189)",
            hovertemplate="%{x}<br>Scatter = %{y:.2f} mm<extra></extra>",
        )
    )
    fig.update_layout(
        barmode="group",
        title="Error Decomposition per AFID: Bias (red) vs Scatter (blue)",
        yaxis_title="mm",
        height=450,
        margin=dict(t=50, b=80, l=60, r=20),
        template="plotly_white",
        legend=dict(x=0.65, y=0.95),
        xaxis=dict(tickangle=-45),
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

SUMMARY_TEMPLATE = Template(
    """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <title>AutoAFIDs — Dataset Registration QC</title>

  <style>
    body {
      font-family: 'Segoe UI', sans-serif;
      background: #f4f6f8;
      margin: 0; padding: 0;
    }
    header {
      background: linear-gradient(
            135deg,
            #1a1a2e 0%,
            #16213e 60%,
            #0f3460 100%
          );
      color: white;
      padding: 36px 48px;
      box-shadow: 0 4px 16px rgba(0,0,0,0.2);
    }
    header h1 { margin: 0;
                font-size: 28px;
                font-weight: 700;
                letter-spacing: 0.5px;
              }
    header p  { margin: 6px 0 0;
                font-size: 14px;
                color: rgba(255,255,255,0.7);
              }
    .container { padding: 36px 48px;
                 max-width: 1500px;
                 margin: auto;
                }
    .card {
      background: white;
      border-radius: 14px;
      box-shadow: 0 4px 16px rgba(0,0,0,0.05);
      padding: 36px;
      margin-bottom: 36px;
    }
    .section-title {
      font-size: 20px; font-weight: 600;
      margin-bottom: 18px;
      border-bottom: 2px solid #eaeaea;
      padding-bottom: 8px;
    }
    .tile-row { display: flex; flex-wrap: wrap; gap: 14px; }
    .tile {
      background: #f0f2f5; border-radius: 10px;
      padding: 16px 20px; flex: 1 1 140px; text-align: center;
      transition: background 0.2s;
    }
    .tile:hover { background: #e2e6ea; }
    .tile-label { font-size: 13px; font-weight: 500; color: #666; }
    .tile-value {
                  font-size: 22px;
                  font-weight: 700;
                  color: #111;
                  margin-top: 4px;
                }
    table {
      width: 100%; border-collapse: collapse; font-size: 14px;
    }
    thead th {
      background: #f0f2f5; padding: 10px 14px;
      text-align: left; font-weight: 600; cursor: pointer;
      white-space: nowrap;
    }
    thead th:hover { background: #e0e4ea; }
    tbody tr:nth-child(even) { background: #fafafa; }
    tbody tr:hover { background: #eef2ff; }
    tbody td { padding: 9px 14px; }
    .badge {
      display: inline-block; padding: 2px 10px;
      border-radius: 12px; font-size: 12px; font-weight: 600;
    }
    .badge-ok   { background: #d4edda; color: #155724; }
    .badge-warn { background: #fff3cd; color: #856404; }
    .badge-bad  { background: #f8d7da; color: #721c24; }
    a.report-link {
      color: #3d5af1; text-decoration: none; font-weight: 500;
    }
    a.report-link:hover { text-decoration: underline; }
    .plot-wrap { overflow-x: auto; }
    @media (max-width: 900px) {
      .container { padding: 18px; }
      .card { padding: 18px; }
    }
  </style>
</head>
<body>
<header>
  <h1>AutoAFIDs — Dataset Registration QC</h1>
  <p>
    {{ n_subjects }} subjects &nbsp;|&nbsp;
    {{ n_afids }} AFIDs per subject
  </p>
</header>
<div class="container">

  <!-- ── Dataset Stats ── -->
  <div class="card">
    <div class="section-title">Dataset Summary</div>
    <div class="tile-row">
      <div class="tile">
        <div class="tile-label">Subjects</div>
        <div class="tile-value">{{ n_subjects }}</div>
      </div>
      <div class="tile">
        <div class="tile-label">Dataset Mean ED</div>
        <div class="tile-value">{{ "%.2f"|format(dataset_mean_ed) }} mm</div>
      </div>
      <div class="tile">
        <div class="tile-label">Dataset Median ED</div>
        <div class="tile-value">{{ "%.2f"|format(dataset_median_ed) }} mm</div>
      </div>
      <div class="tile">
        <div class="tile-label">Dataset Std Dev</div>
        <div class="tile-value">{{ "%.2f"|format(dataset_std_ed) }} mm</div>
      </div>
      <div class="tile">
        <div class="tile-label">Best Subject (mean ED)</div>
        <div class="tile-value">{{ "%.2f"|format(best_ed) }} mm</div>
      </div>
      <div class="tile">
        <div class="tile-label">Worst Subject (mean ED)</div>
        <div class="tile-value">{{ "%.2f"|format(worst_ed) }} mm</div>
      </div>
    </div>
  </div>

  <!-- ── Subject Table ── -->
  <div class="card">
    <div class="section-title">Per-Subject Summary</div>
    <table id="subjectTable">
      <thead>
        <tr>
          <th onclick="sortTable(0)">Subject ▾</th>
          <th onclick="sortTable(1)">Mean ED (mm) ▾</th>
          <th onclick="sortTable(2)">Median ED (mm) ▾</th>
          <th onclick="sortTable(3)">Max ED (mm) ▾</th>
          <th onclick="sortTable(4)">Worst AFID ▾</th>
          <th>Report</th>
        </tr>
      </thead>
      <tbody>
        {% for row in subject_rows %}
        <tr>
          <td><code>{{ row.subject_label }}</code></td>
          <td>
            <span class="badge {{ row.badge }}">
                            {{ "%.2f"|format(row.mean_ED) }}
            </span>
          </td>
          <td>{{ "%.2f"|format(row.median_ED) }}</td>
          <td>{{ "%.2f"|format(row.max_ED) }}</td>
          <td>AFID {{ row.worst_afid }}</td>
          <td>
            {% if row.html_relpath %}
            <a class="report-link"
              href="{{ row.html_relpath }}" target="_blank">
              Open ↗
            </a>
            {% else %}
            —
            {% endif %}
          </td>
        </tr>
        {% endfor %}
      </tbody>
    </table>
  </div>

  <!-- ── Heatmap ── -->
  <div class="card">
    <div class="section-title">
    Registration Error Heatmap (Subjects x AFIDs)
    </div>
    <div class="plot-wrap">{{ heatmap_html | safe }}</div>
  </div>

  <!-- ── AFID Boxplots ── -->
  <div class="card">
    <div class="section-title">
      ED Distribution per AFID (across subjects)
    </div>
    <div class="plot-wrap">{{ boxplot_html | safe }}</div>
  </div>

  <!-- ── Subject Bar Chart ── -->
  <div class="card">
    <div class="section-title">Mean ED per Subject</div>
    <div class="plot-wrap">{{ bar_html | safe }}</div>
  </div>

  <!-- ── 3D Landmark Scatter ── -->
  {% if scatter3d_html %}
  <div class="card">
    <div class="section-title">AFID Landmarks in MNI Space</div>
    <p style="color:#666;font-size:13px;margin-top:-10px;margin-bottom:14px;">
      Each subject's warped AFIDs are plotted alongside the MNI ground-truth
      template (green diamonds). Hover to identify individual AFIDs.
    </p>
    <div class="plot-wrap">{{ scatter3d_html | safe }}</div>
  </div>
  {% endif %}

  <!-- ── Error Decomposition ── -->
  {% if decomp_html %}
  <div class="card">
    <div class="section-title">
    Error Decomposition: Systematic Bias vs Random Scatter
    </div>
    <p style="color:#666;font-size:13px;margin-top:-10px;margin-bottom:14px;">
      Registration error at each AFID is decomposed into a <strong>systematic
      bias</strong> (mean signed error vector, consistent across subjects —
      shown in <span style="color:#d73027;font-weight:600">red</span>) and
      <strong>random scatter</strong> (standard deviation, varies between
      subjects — shown in
      <span style="color:#3182bd;font-weight:600">
                            blue
      </span>).
      Systematic bias is potentially correctable; random scatter represents
      the precision floor of the registration.
    </p>

    <!-- Decomposition summary tiles -->
    <div class="tile-row" style="margin-bottom:24px;">
      <div class="tile">
        <div class="tile-label">Overall Bias (mean)</div>
        <div class="tile-value">{{ "%.2f"|format(decomp_bias_mean) }} mm</div>
      </div>
      <div class="tile">
        <div class="tile-label">Overall Scatter (mean)</div>
        <div class="tile-value">
          {{ "%.2f"|format(decomp_scatter_mean) }} mm
        </div>
      </div>
      <div class="tile">
        <div class="tile-label">Systematic Fraction</div>
        <div class="tile-value">{{ "%.0f"|format(decomp_bias_pct) }}%</div>
      </div>
      <div class="tile">
        <div class="tile-label">Most Biased AFID</div>
        <div class="tile-value">
        AFID {{ decomp_worst_bias_afid }}
        ({{ "%.2f"|format(decomp_worst_bias_val) }} mm)
        </div>
      </div>
      <div class="tile">
        <div class="tile-label">Most Variable AFID</div>
        <div class="tile-value">
        AFID {{ decomp_worst_scatter_afid }}
        ({{ "%.2f"|format(decomp_worst_scatter_val) }} mm)
        </div>
      </div>
    </div>

    <!-- Decomposition bar chart -->
    <div class="plot-wrap">{{ decomp_bar_html | safe }}</div>
  </div>

  <!-- Bias vectors in 3D -->
  <div class="card">
    <div class="section-title">Systematic Bias Vectors in MNI Space</div>
    <p style="color:#666;font-size:13px;margin-top:-10px;margin-bottom:14px;">
      Cones show the direction and magnitude of systematic bias at each AFID.
      Colour encodes the bias-to-total ratio:
      <span style="color:#d73027;font-weight:600">red</span> = mostly
      systematic (correctable),
      <span style="color:#3182bd;font-weight:600">blue</span> = mostly
      random (precision limit). Hover for details.
    </p>
    <div class="plot-wrap">{{ decomp_3d_html | safe }}</div>
  </div>
  {% endif %}

</div>

<script>
/* Simple client-side table sort */
let sortDir = {};
function sortTable(col) {
  const table = document.getElementById("subjectTable");
  const tbody = table.querySelector("tbody");
  const rows  = Array.from(tbody.querySelectorAll("tr"));
  const asc   = (sortDir[col] = !sortDir[col]);
  rows.sort((a, b) => {
    const av = a.cells[col].innerText.trim();
    const bv = b.cells[col].innerText.trim();
    const an = parseFloat(av), bn = parseFloat(bv);
    if (!isNaN(an) && !isNaN(bn)) return asc ? an - bn : bn - an;
    return asc ? av.localeCompare(bv) : bv.localeCompare(av);
  });
  rows.forEach(r => tbody.appendChild(r));
}
</script>
</body>
</html>
"""
)


def _ed_badge(mean_ed: float) -> str:
    if mean_ed < 1.5:
        return "badge-ok"
    elif mean_ed < 3.0:
        return "badge-warn"
    return "badge-bad"


def render_summary_html(
    output_path: str,
    long_df: pd.DataFrame,
    output_dir: str,
    subject_fcsvs: list | None = None,
    gt_coords: np.ndarray = None,
    gt_labels: list | None = None,
) -> str:
    """Build and write the dataset summary HTML file."""

    # Per-subject stats
    grp = long_df.groupby("subject_label")
    stats = (
        grp["ED (mm)"]
        .agg(mean_ED="mean", median_ED="median", max_ED="max")
        .reset_index()
    )
    worst_afid = long_df.loc[long_df.groupby("subject_label")["ED (mm)"].idxmax()][
        ["subject_label", "AFID"]
    ].rename(columns={"AFID": "worst_afid"})
    html_paths_ser = long_df.groupby("subject_label")["html_path"].first().reset_index()
    subject_stats = (
        stats.merge(worst_afid, on="subject_label")
        .merge(html_paths_ser, on="subject_label")
        .sort_values("mean_ED", ascending=False)
    )

    # Make relative paths from the output HTML's directory
    out_dir = Path(output_path).parent

    def _relpath(p):
        if not p:
            return ""
        try:
            return os.path.relpath(p, out_dir)
        except ValueError:
            return p

    subject_stats["html_relpath"] = subject_stats["html_path"].apply(_relpath)
    subject_stats["badge"] = subject_stats["mean_ED"].apply(_ed_badge)

    # Dataset-level stats (all subjects x all AFIDs)
    all_ed = long_df["ED (mm)"]
    dataset_mean = float(all_ed.mean())
    dataset_median = float(all_ed.median())
    dataset_std = float(all_ed.std())
    best_ed = float(subject_stats["mean_ED"].min())
    worst_ed = float(subject_stats["mean_ED"].max())
    n_subjects = int(subject_stats["subject_label"].nunique())
    n_afids = int(long_df["AFID"].nunique())

    # Build charts
    heatmap_html = make_heatmap(long_df)
    boxplot_html = make_afid_boxplot(long_df)
    bar_html = make_subject_bar(subject_stats)

    # 3D landmark scatter (optional — only if fcsv data provided)
    scatter3d_html = ""
    if subject_fcsvs:
        scatter3d_html = make_3d_landmark_scatter(subject_fcsvs, gt_coords, gt_labels)

    # Error decomposition (requires ≥1 subject fcsv + gt)
    decomp_html = False
    decomp_bar_html = ""
    decomp_3d_html = ""
    decomp_bias_mean = 0.0
    decomp_scatter_mean = 0.0
    decomp_bias_pct = 0.0
    decomp_worst_bias_afid = "-"
    decomp_worst_bias_val = 0.0
    decomp_worst_scatter_afid = "-"
    decomp_worst_scatter_val = 0.0
    if subject_fcsvs and gt_coords is not None:
        decomp_df = compute_error_decomposition(subject_fcsvs, gt_coords, gt_labels)
        if not decomp_df.empty:
            decomp_html = True
            decomp_bar_html = make_decomposition_bar(decomp_df)
            decomp_3d_html = make_bias_vector_3d(decomp_df)
            decomp_bias_mean = float(decomp_df["bias_mag"].mean())
            decomp_scatter_mean = float(decomp_df["scatter_rms"].mean())
            decomp_bias_pct = (
                decomp_bias_mean / (decomp_bias_mean + decomp_scatter_mean + 1e-9) * 100
            )
            worst_bias_idx = decomp_df["bias_mag"].idxmax()
            decomp_worst_bias_afid = decomp_df.loc[worst_bias_idx, "label"]
            decomp_worst_bias_val = decomp_df.loc[worst_bias_idx, "bias_mag"]
            worst_scatter_idx = decomp_df["scatter_rms"].idxmax()
            decomp_worst_scatter_afid = decomp_df.loc[worst_scatter_idx, "label"]
            decomp_worst_scatter_val = decomp_df.loc[worst_scatter_idx, "scatter_rms"]

    # Render template
    html = SUMMARY_TEMPLATE.render(
        n_subjects=n_subjects,
        n_afids=n_afids,
        dataset_mean_ed=dataset_mean,
        dataset_median_ed=dataset_median,
        dataset_std_ed=dataset_std,
        best_ed=best_ed,
        worst_ed=worst_ed,
        subject_rows=subject_stats.to_dict("records"),
        heatmap_html=heatmap_html,
        boxplot_html=boxplot_html,
        bar_html=bar_html,
        scatter3d_html=scatter3d_html,
        decomp_html=decomp_html,
        decomp_bar_html=decomp_bar_html,
        decomp_3d_html=decomp_3d_html,
        decomp_bias_mean=decomp_bias_mean,
        decomp_scatter_mean=decomp_scatter_mean,
        decomp_bias_pct=decomp_bias_pct,
        decomp_worst_bias_afid=decomp_worst_bias_afid,
        decomp_worst_bias_val=decomp_worst_bias_val,
        decomp_worst_scatter_afid=decomp_worst_scatter_afid,
        decomp_worst_scatter_val=decomp_worst_scatter_val,
    )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html, encoding="utf-8")
    return output_path


# ---------------------------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------------------------


def main(csv_paths, html_paths, output_html, fcsv_paths=None, gt_fcsv_path=None):
    long_df = load_all_csvs(csv_paths, html_paths)
    out_dir = str(Path(output_html).parent)

    # Load FCSV landmark data if available
    subject_fcsvs = []
    gt_coords, gt_labels = None, None
    if fcsv_paths:
        subject_fcsvs = load_all_fcsvs(fcsv_paths)
    if gt_fcsv_path:
        try:
            gt_coords, gt_labels = load_fcsv(gt_fcsv_path)
        except Exception:
            pass

    render_summary_html(
        output_html,
        long_df,
        out_dir,
        subject_fcsvs=subject_fcsvs,
        gt_coords=gt_coords,
        gt_labels=gt_labels,
    )
    print(f"[regqc_summary] Written → {output_html}")


# Snakemake script entry point
fcsv_list = list(snakemake.input.get("fcsvs", []))
gt_path = snakemake.params.get("gt_fcsv", None)
main(
    csv_paths=list(snakemake.input["csvs"]),
    html_paths=list(snakemake.input["htmls"]),
    output_html=snakemake.output["summary_html"],
    fcsv_paths=fcsv_list,
    gt_fcsv_path=gt_path,
)
