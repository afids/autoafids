# === REGISTRATION QC SCRIPT WITH MRI VIEWER ===
import base64
import csv
from io import BytesIO
from pathlib import Path

import ants
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import scipy.io
import SimpleITK as SimpleITK
from jinja2 import Template

FCSV_TEMPLATE = (
    Path(__file__).parent
    / ".."
    / ".."
    / "resources"
    / "tpl-MNI152NLin2009cAsym_res-01_T1w.fcsv"
)


def afids_to_fcsv(
        transformed_coords, output_path, template_path=FCSV_TEMPLATE
    ):
    """
    Writes transformed AFID coordinates
    to a new FCSV file using a template FCSV.

    Parameters
    ----------
    transformed_coords : np.ndarray
        Transformed coordinates of shape (N, 3)
    template_path : str or Path
        Path to the template FCSV file (Slicer-compatible)
    output_path : str or Path
        Output path for the new FCSV file
    """
    with open(template_path, encoding="utf-8", newline="") as template_file:
        # Read the first three header lines
        header = [template_file.readline() for _ in range(3)]
        reader = csv.DictReader(
            template_file,
            fieldnames=[
                "id",
                "x",
                "y",
                "z",
                "ow",
                "ox",
                "oy",
                "oz",
                "vis",
                "sel",
                "lock",
                "label",
                "desc",
                "associatedNodeID"
                ]
            )
        rows = list(reader)

    # Update coordinates
    for i, row in enumerate(rows):
        row['x'], row['y'], row['z'] = map(str, transformed_coords[i])

    # Write updated FCSV
    with open(output_path, "w", newline="") as out_file:
        out_file.writelines(header)
        writer = csv.DictWriter(out_file, fieldnames=reader.fieldnames)
        writer.writerows(rows)

# --- LOAD FCSV ---
def load_fcsv(path):
    """
    Loads an FCSV file and returns both coordinates and full DataFrame.

    Returns:
    - coords: np.ndarray of shape (N, 3)
    - labels: list of coordinate labels
    """
    df = pd.read_csv(path, skiprows=2)
    coords = df[['x', 'y', 'z']].to_numpy()
    labels = df['label'].tolist()
    return coords, labels

# --- APPLY 4X4 ---
def apply_affine_transform(mat_path, coords):
    """
    Applies a transformation matrix to 3D coordinates.

    Parameters:
    - mat_path: str, path to the .mat file with the transformation matrix
    - coords: np.ndarray of shape (N, 3), input coordinates

    Returns:
    - transformed_coords: np.ndarray of shape (N, 3)
    """
    mat_contents = scipy.io.loadmat(mat_path)
    transformation_matrix = mat_contents['tmat']
    homogeneous_coords = np.hstack((coords, np.ones((coords.shape[0], 1))))
    transformed_homogeneous = (
        np.linalg.inv(transformation_matrix) @ homogeneous_coords.T
    ).T
    return transformed_homogeneous[:, :3]

# --- APPLY WARP (HARMONIZED) ---
def apply_warp_deformation(transform_path, coords, flip_ras_lps=True, leaddbs=True):
    """
    Applies a non-linear warp deformation to coordinates
    using a SimpleITK displacement field.

    Parameters:
    - transform_path: str, path to the transformation file (.nii or .h5)
    - coords: np.ndarray of shape (N, 3), coordinates to transform
    - flip_ras_lps: bool, whether to flip x/y axes for RAS ↔ LPS conversion

    Returns:
    - transformed_coords: np.ndarray of shape (N, 3), transformed coordinates
    """
    if flip_ras_lps:
        coords = coords * np.array([-1, -1, 1])

    if leaddbs:

        transform_image = SimpleITK.ReadImage(transform_path)
        transform_image = SimpleITK.Cast(
            transform_image, SimpleITK.sitkVectorFloat64
        )
        transform = SimpleITK.Transform(transform_image)

        transformed_coords = np.array([
            transform.TransformPoint(point.tolist()) for point in coords
        ])
    
    else:
        d = pd.DataFrame(data=coords, columns=['x','y','z'])
        transformed_coords_dataframe = ants.apply_transforms_to_points( 3, d, transform_path)
        transformed_coords = transformed_coords_dataframe.to_numpy()

    if flip_ras_lps:
        transformed_coords = transformed_coords * np.array([-1, -1, 1])

    return transformed_coords

# --- COMPUTE ERROR ---
def compute_error_components(gt, pred):
    dx = np.abs(pred[:, 0] - gt[:, 0])
    dy = np.abs(pred[:, 1] - gt[:, 1])
    dz = np.abs(pred[:, 2] - gt[:, 2])
    ed = np.linalg.norm(pred - gt, axis=1)
    return dx, dy, dz, ed


# --- HEATMAP PLOT ---
def make_toggleable_heatmap(error_components, afid_ids):
    components = ["x", "y", "z", "ED"]
    data_matrix = np.stack(error_components, axis=1)
    y_labels = [f"AFID {i}" for i in afid_ids]
    fig = go.Figure(data=go.Heatmap(
        z=data_matrix, x=components, y=y_labels,
        colorscale=[[0.0, 'rgb(255,255,255)'], [1.0, 'rgb(220,0,0)']],
        colorbar=dict(title="Error (mm)"),
        hovertemplate="AFID %{y}<br>%{x} Error = %{z:.2f} mm<extra></extra>"
    ))
    fig.update_layout(
        height=600,
        width=500,
        margin=dict(t=0, b=0, l=0, r=0),
        template='plotly_white'
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs='cdn')


# --- 3D SCATTER PLOT ---
def make_3d_plot(gt, pred, afid_ids):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=gt[:, 0],
            y=gt[:, 1],
            z=gt[:, 2],
            mode='markers',
            name='Stereotactic Space',
            marker=dict(
                size=6, color='green', opacity=0.7, symbol='diamond'
            ),
            text=afid_ids
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=pred[:, 0],
            y=pred[:, 1],
            z=pred[:, 2],
            mode='markers',
            name='Registered Subject',
            marker=dict(
                size=6, color='red', opacity=0.7, symbol='cross'
            ),
            text=afid_ids
        )
    )
    for i in range(gt.shape[0]):
        fig.add_trace(
            go.Scatter3d(
                x=[gt[i, 0], pred[i, 0]],
                y=[gt[i, 1], pred[i, 1]],
                z=[gt[i, 2], pred[i, 2]],
                mode='lines',
                line=dict(color='orange', width=7),
                showlegend=False
            )
        )
    fig.update_layout(
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            bgcolor='rgb(250,250,250)'
        ),
        height=600,
        width=800,
        margin=dict(t=0, b=0, l=0, r=0),
        template='plotly_white',
        legend=dict(x=0.7, y=0.9)
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs='cdn')

# --- Coordinate Utilities ---
def world_to_voxel_coords(img, world_coords):
    """Convert world (RAS+) coordinates to
    voxel coordinates using inverse affine."""
    affine_inv = np.linalg.inv(img.affine)
    return np.round(affine_inv.dot(np.append(world_coords, 1))).astype(int)

def get_afid_world_and_voxel_coords(img, fcsv_path, label_num):
    """Get both world and voxel coordinates for a given AFID label."""
    df = pd.read_csv(fcsv_path, skiprows=2)
    row = df[df['label'] == label_num]
    if row.empty:
        return None, None
    world_coords = row[['x', 'y', 'z']].to_numpy()[0]
    voxel_coords = world_to_voxel_coords(img, world_coords)
    return world_coords, voxel_coords

# --- Visualization Utilities ---
def render_mri_planes(data, coords, zoom=50, show_crosshairs=True):
    """
    Render sagittal, coronal, and axial slices
    centered at given voxel coordinates.

    Parameters:
    - data: 3D NumPy array of the MRI volume
    - coords: tuple of 3 voxel coordinates (x, y, z)
    - zoom: int, number of voxels to include around center
    - show_crosshairs: bool, whether to draw red crosshairs

    Returns:
    - base64-encoded PNG image with the 3 slices side-by-side
    """
    x, y, z = coords
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    slice_data_list = [data[x, :, :], data[:, y, :], data[:, :, z]]
    cross_coords_list = [(y, z), (x, z), (x, y)]
    titles = ["Sagittal", "Coronal", "Axial"]

    for i, ax in enumerate(axes):
        slice_2d = slice_data_list[i].T  # Transpose for display
        x_cross, y_cross = cross_coords_list[i]

        ax.imshow(slice_2d, cmap="gray", origin="lower")

        if show_crosshairs:
            ax.axvline(x_cross, color="r", lw=1)
            ax.axhline(y_cross, color="r", lw=1)

        ax.set_xlim(x_cross - zoom, x_cross + zoom)
        ax.set_ylim(y_cross - zoom, y_cross + zoom)
        ax.set_title(titles[i])
        ax.axis("off")

    buf = BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()

# --- HTML Utilities ---
def build_afid_viewer_block(label, error_val, subject_image, ref_image):
    return f"""
    <div class="viewer-card" style="text-align: center">

        <!-- Header -->
        <div class="viewer-header"
            style="
                font-weight: bold;
                font-size: 18px;
                margin-bottom: 1px;">
            AFID {label}: {error_val:.2f} mm
        </div>

        <!-- Image container -->
        <div class="fade-image-wrapper"
            style="
                position: relative;
                width: 100%;
                height: 400px;
                border-radius: 10px;
                overflow: hidden;">
            <img src="data:image/png;base64,{subject_image}"
                style="
                    position: absolute;
                    top: 0;
                    left: 0;
                    width: 100%;
                    height: 100%;
                    object-fit: contain;
                    z-index: 1;" />
            <img id="fade_{label}"
                src="data:image/png;base64,{ref_image}"
                style="
                    position: absolute;
                    top: 0;
                    left: 0;
                    width: 100%;
                    height: 100%;
                    object-fit: contain;
                    z-index: 2;
                    opacity: 0;
                    ransition: opacity 1s ease-in-out;" />
        </div>

        <!-- Labels below -->
        <div style="margin-top: 0px; font-size: 14px;">
            <span id="label_subject_{label}"
            style="
                padding: 3px 10px;
                background-color: #ddd;
                border-radius: 6px;
                font-weight: bold;">Subject Space</span>
            &nbsp;←→&nbsp;
            <span id="label_ref_{label}"
            style="
                padding: 3px 10px;
                background-color: #f0f0f0;
                border-radius: 6px;
                color: #999;">Stereotactic Space</span>
        </div>

        <!-- Scoped JS -->
        <script>
        (function() {{
            const fadeImg = document.getElementById("fade_{label}");
            const labelSub = document.getElementById("label_subject_{label}");
            const labelRef = document.getElementById("label_ref_{label}");
            let fadeIn = false;

            setInterval(() => {{
                fadeIn = !fadeIn;
                fadeImg.style.opacity = fadeIn ? 1 : 0;

                labelRef.style.color = fadeIn ? "#000" : "#999";
                labelRef.style.backgroundColor = fadeIn ? "#ddd" : "#f0f0f0";
                labelRef.style.fontWeight = fadeIn ? "bold" : "normal";

                labelSub.style.color = fadeIn ? "#999" : "#000";
                labelSub.style.backgroundColor = fadeIn ? "#f0f0f0" : "#ddd";
                labelSub.style.fontWeight = fadeIn ? "normal" : "bold";
            }}, 2500);
        }})();
        </script>
    </div>
    """

# --- Main Viewer Builder ---
def generate_qualitative_mri_html(
        subject_nii,
        ref_nii,
        ref_fcsv,
        ed_errors,
        labels=list(range(1, 33)),
        zoom=40
    ):
    """Generates HTML for visual comparison of AFID localization
    on subject vs template MRI."""
    subject_img = nib.as_closest_canonical(nib.load(subject_nii))
    ref_img = nib.as_closest_canonical(nib.load(ref_nii))
    subject_data = subject_img.get_fdata()
    ref_data = ref_img.get_fdata()

    viewer_blocks = []
    for i, label in enumerate(labels):
        world_coords, template_voxel = get_afid_world_and_voxel_coords(
            ref_img, ref_fcsv, label
        )
        if world_coords is None:
            continue
        subject_voxel = world_to_voxel_coords(subject_img, world_coords)
        subject_image = render_mri_planes(
            subject_data, subject_voxel[:3], zoom=zoom
        )
        ref_image = render_mri_planes(ref_data, template_voxel[:3], zoom=zoom)
        viewer_blocks.append(
            build_afid_viewer_block(
                label, ed_errors[i], subject_image, ref_image
            )
        )

    return f"""
    <div class='card'>
        <div class='section-title'>Registration Error Qualitative Review</div>
        <div
            style="
                display: flex;
                align-items: center;
                justify-content:
                space-between;">
            <button onclick="scrollAfid(-1)"
                style="font-size: 20px;">&#8592;</button>
            <div style="flex: 1; overflow: hidden; position: relative;">
                <div class='afid-strip'
                    id="afid-strip"
                    onscroll="updateAfidCounter()">
                    {''.join(viewer_blocks)}
                </div>
                <div id="afid-counter"
                    style="position: absolute;
                    bottom: 10px; right: 16px;
                    background: rgba(255,255,255,0.8);
                    padding: 4px 10px;
                    border-radius: 8px;
                    font-weight: 600;">
                    1 / {len(viewer_blocks)}
                </div>
            </div>
            <button onclick="scrollAfid(1)"
                style="font-size: 20px;">&#8594;</button>
        </div>
    </div>
    """

# --- RENDER HTML DASHBOARD ---
def render_dashboard_html(
        output_path, heatmap_html, scatter_html, ed_errors, qualitative_html
    ):
    stats = pd.Series(ed_errors).describe(percentiles=[.25, .5, .75]).to_dict()
    std = np.std(ed_errors)
    template = Template("""
    <html>
    <head>
        <title>AFID Inspector</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            #zoom-wrapper {
                transform: scale(0.8);
                transform-origin: top left;
                width: 125%;
            }
            body {
                font-family: 'Segoe UI', sans-serif;
                background: #f4f6f8;
                margin: 0;
                padding: 0;
            }
            header {
                background: white;
                padding: 30px;
                text-align: center;
                box-shadow: 0 4px 10px rgba(0,0,0,0.05);
            }
            header h1 {
                margin: 0;
                font-size: 32px;
                font-weight: 700;
            }
            .container {
                padding: 40px;
                max-width: 1400px;
                margin: auto;
            }
            .card {
                background: white;
                border-radius: 12px;
                box-shadow: 0 4px 16px rgba(0,0,0,0.05);
                padding: 40px;
                margin-bottom: 40px;
            }
            .section-title {
                font-size: 22px;
                font-weight: 600;
                margin-bottom: 20px;
                border-bottom: 2px solid #eaeaea;
                padding-bottom: 8px;
            }
            .tile-row {
                display: flex;
                flex-wrap: wrap;
                gap: 16px;
            }
            .tile {
                background-color: #f0f2f5;
                border-radius: 10px;
                padding: 10px;
                flex: 1 1 160px;
                text-align: center;
            }
            .tile-label {
                font-size: 14px;
                font-weight: 500;
                color: #666;
            }
            .tile-value {
                font-size: 19px;
                font-weight: 700;
                color: #111;
            }
            .flex-row {
                display: flex;
                flex-direction: row;
                gap: 24px;
                flex-wrap: wrap;
            }
            .plot-box-left {
                flex: 0.4;
                min-width: 200px;
            }
            .plot-box-right {
                flex: 0.6;
                min-width: 200px;
            }
            .tile:hover {
                background-color: #e0e0e0;
                cursor: pointer;
                transition: background-color 0.3s ease;
            }
            .afid-strip {
                display: flex;
                overflow-x: auto;
                gap: 16px;
                padding: 10px 0;
                scroll-snap-type: x mandatory;
                scroll-behavior: smooth
            }

            .afid-strip::-webkit-scrollbar {
                display: none;
            }

            .viewer-card {
                flex: 0 0 auto;
                width: 100%;
                height: 480px;
                background-color: white;
                border-radius: 16px;
                scroll-snap-align: start;
                text-align: center;
                position: center;
                padding: 16px;
            }

            .viewer-header {
                font-weight: 600;
                font-size: 20px;
                background: #E5E5E5;
                padding: 3px 0;
                border-bottom: 1px solid #e0e0e0;
            }
            @media (max-width: 1000px) {
                .flex-row {
                    flex-direction: column;
                }
                .plot-box-left,
                .plot-box-right {
                    flex: 1 1 100%;
                    max-height: 420px;
                }
            }
        </style>
    </head>
    <body>
        <div id="zoom-wrapper">
            <header>
                <h1>AFID Registration QC</h1>
            </header>
            <div class="container">
                <!-- Summary Card -->
                <div class="card">
                    <div class="section-title">Registration Error Summary</div>
                    <div class="tile-row">
                    {% for label, value in stats.items() if label != 'Mean' %}
                        <div class="tile">
                            <div class="tile-label">{{ label }}</div>
                            <div class="tile-value">
                                {{ '%.2f' % value }} mm
                            </div>
                        </div>
                        {% endfor %}
                    </div>
                    <div class="tile-row" style="margin-top: 20px;">
                        <div class="tile">
                            <div class="tile-label">Mean</div>
                            <div class="tile-value">
                                {{ '%.2f' % stats['Mean'] }} mm
                            </div>
                        </div>
                        <div class="tile">
                            <div class="tile-label">Std Dev</div>
                            <div class="tile-value">{{ '%.2f' % std }} mm</div>
                        </div>
                    </div>
                </div>

                <!-- MRI Viewer Card -->
                {{ qualitative|safe }}

                <!-- Heatmap + 3D Plot Card -->
                <div class="card">
                    <div class="section-title">
                        Registration Error Visualization
                    </div>
                    <div class="flex-row">
                        <div class="plot-box-left">
                            {{ heatmap|safe }}
                        </div>
                        <div class="plot-box-right">
                            {{ scatter|safe }}
                        </div>
                    </div>
                </div>
            </div>
        <script>
        function scrollAfid(direction) {
            const container = document.getElementById('afid-strip');
            const card = container.querySelector('.viewer-card');
            if (!card) return;
            const cardWidth = card.offsetWidth + 16;
            container.scrollBy(
                        { left: direction * cardWidth, behavior: 'smooth' });
            setTimeout(updateAfidCounter, 300);
        }

        function updateAfidCounter() {
            const container = document.getElementById('afid-strip');
            const cards = container.querySelectorAll('.viewer-card');
            const scrollLeft = container.scrollLeft;
            const cardWidth = cards[0]?.offsetWidth + 16;
            const activeIndex = Math.round(scrollLeft / cardWidth);
            const counter = document.getElementById('afid-counter');
            if (counter) counter.textContent =
                        `${activeIndex + 1} / ${cards.length}`;
        }

        document.addEventListener('keydown', function(event) {
            if (event.key === "ArrowRight") scrollAfid(1);
            else if (event.key === "ArrowLeft") scrollAfid(-1);
        });
        </script>
        </div>
    </body>
    </html>
    """)
    html = template.render(
        stats={
            "Min": stats["min"],
            "25%": stats["25%"],
            "Median": stats["50%"],
            "75%": stats["75%"],
            "Max": stats["max"],
            "Mean": stats["mean"]
        },
        std=std,
        heatmap=heatmap_html,
        scatter=scatter_html,
        qualitative=qualitative_html
    )
    Path(output_path).write_text(html)
    return output_path


# --- MAIN WRAPPER ---
def generate_afid_qc_dashboard(
        template_name,
        gt_fcsv_path,
        pred_fcsv_path,
        output_html_path,
        output_fcsv_path,
        output_csv_path,subject_nii,
        template_dir,
        matfile_path,
        warpfile_path
    ):
    native_coords, afid_ids = load_fcsv(gt_fcsv_path)
    if matfile_path:
        gt_coords_a= apply_affine_transform(matfile_path,native_coords)
        gt_coords= apply_warp_deformation(warpfile_path,gt_coords_a)
    else:
        gt_coords= apply_warp_deformation(warpfile_path,native_coords, leaddbs=False)
    pred_coords, _ = load_fcsv(pred_fcsv_path)
    afids_to_fcsv(gt_coords, output_fcsv_path)
    dx, dy, dz, ed = compute_error_components(gt_coords, pred_coords)

    # --- EXPORT TO CSV ---
    error_df = pd.DataFrame({
        "AFID": afid_ids,
        "dx (mm)": dx,
        "dy (mm)": dy,
        "dz (mm)": dz,
        "ED (mm)": ed
    })

    if template_name in ['MNI152NLin2009bAsym','MNI152NLin2009bSym']:
        template_nii = str(Path(template_dir) / f'tpl-{template_name}' / f'tpl-{template_name}_res-1_T1w.nii.gz')
    elif template_name in ['MNI305', 'MNIColin27']:
        template_nii = str(Path(template_dir) / f'tpl-{template_name}' / f'tpl-{template_name}_T1w.nii.gz')
    else:
        template_nii = str(Path(template_dir) / f'tpl-{template_name}' / f'tpl-{template_name}_res-01_T1w.nii.gz')
    error_df.to_csv(output_csv_path, index=False)
    heatmap_html = make_toggleable_heatmap([dx, dy, dz, ed], afid_ids)
    scatter_html = make_3d_plot(gt_coords, pred_coords, afid_ids)
    qualitative_html = generate_qualitative_mri_html(
        subject_nii, template_nii, pred_fcsv_path, ed,zoom=50)
    return render_dashboard_html(
        output_html_path, heatmap_html, scatter_html, ed, qualitative_html)


if __name__ == "__main__":
    generate_afid_qc_dashboard(
        template_name=snakemake.params['template'],
        gt_fcsv_path=snakemake.input["afidfcsv"],
        pred_fcsv_path=snakemake.input["refcoord"],
        output_html_path=snakemake.output["html"],
        output_fcsv_path=snakemake.output["fcsv"],
        output_csv_path=snakemake.output["csv"],
        subject_nii=snakemake.input["im"],
        template_dir=snakemake.input["refim_dir"],
        matfile_path=snakemake.input["optional_matrix"],
        warpfile_path=snakemake.input["warp"]
    )
