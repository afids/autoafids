from pathlib import Path
import pandas as pd
import numpy as np
from jinja2 import Template
import plotly.graph_objects as go
import plotly.io as pio


def load_fcsv(path):
    df = pd.read_csv(path, skiprows=2)
    coords = df[['x', 'y', 'z']].to_numpy()
    labels = df['label'].tolist()
    return coords, labels


def compute_error_components(gt, pred):
    dx = np.abs(pred[:, 0] - gt[:, 0])
    dy = np.abs(pred[:, 1] - gt[:, 1])
    dz = np.abs(pred[:, 2] - gt[:, 2])
    ed = np.linalg.norm(pred - gt, axis=1)
    return dx, dy, dz, ed


def make_toggleable_heatmap(error_components, afid_ids):
    components = ["x", "y", "z", "ED"]
    custom_colorscale = [
        [0.0, 'rgb(255,255,255)'],
        [1.0, 'rgb(220,0,0)']
    ]
    traces = []
    for i, (comp, data) in enumerate(zip(components, error_components)):
        trace = go.Heatmap(
            z=np.reshape(data, (-1, 1)),
            x=[comp],
            y=[f"AFID {i}" for i in afid_ids],
            colorscale=custom_colorscale,
            colorbar=dict(title="Error (mm)"),
            hovertemplate=f"%{{y}}<br>{comp} Error = %{{z:.2f}} mm<extra></extra>",
            visible=(i == 0)
        )
        traces.append(trace)

    buttons = [
        dict(label=comp,
             method="update",
             args=[{"visible": [i == j for j in range(4)]},
                   {"title": f"AFID Registration Error – {comp}"}])
        for i, comp in enumerate(components)
    ]

    fig = go.Figure(data=traces)
    fig.update_layout(
        title="AFID Registration Error – x",
        updatemenus=[dict(
            buttons=buttons,
            direction="down",
            showactive=True,
            x=0.5,
            xanchor="center",
            y=1.2,
            yanchor="top"
        )],
        height=600,
        margin=dict(t=80, b=40)
    )

    return pio.to_html(fig, full_html=False, include_plotlyjs='cdn')


def make_3d_plot(gt, pred, afid_ids):
    fig = go.Figure()

    fig.add_trace(go.Scatter3d(
        x=gt[:, 0], y=gt[:, 1], z=gt[:, 2],
        mode='markers',
        name='Ground Truth',
        marker=dict(size=6, color='green', opacity=0.7),
        hovertemplate='GT AFID %{text}<br>(%{x:.2f}, %{y:.2f}, %{z:.2f})<extra></extra>',
        text=afid_ids
    ))

    fig.add_trace(go.Scatter3d(
        x=pred[:, 0], y=pred[:, 1], z=pred[:, 2],
        mode='markers',
        name='Predicted',
        marker=dict(size=6, color='red', opacity=0.7, symbol='diamond'),
        hovertemplate='Pred AFID %{text}<br>(%{x:.2f}, %{y:.2f}, %{z:.2f})<extra></extra>',
        text=afid_ids
    ))

    fig.update_layout(
        title="AFID Locations (GT vs Predicted)",
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            bgcolor='rgb(250,250,250)'
        ),
        height=700,
        margin=dict(t=50),
        template='plotly_white',
        legend=dict(x=0.8, y=0.9)
    )

    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


def make_brain_placeholder_slider():
    placeholder_img_url = "https://upload.wikimedia.org/wikipedia/commons/8/88/MNI152_T1_1mm_brain.png"
    z_positions = list(range(40, 90, 10))
    frames = [
        go.Frame(
            data=[go.Image(source=placeholder_img_url)],
            name=f"Z={z}"
        ) for z in z_positions
    ]
    fig = go.Figure(
        data=[go.Image(source=placeholder_img_url)],
        layout=go.Layout(
            title="Brain Slice Viewer (Placeholder)",
            updatemenus=[{
                "type": "buttons",
                "buttons": [{
                    "label": "Play",
                    "method": "animate",
                    "args": [None, {"frame": {"duration": 800, "redraw": True}, "fromcurrent": True}]
                }]
            }],
            sliders=[{
                "steps": [{
                    "args": [[f"Z={z}"], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                    "label": f"Z={z}",
                    "method": "animate"
                } for z in z_positions],
                "x": 0.1,
                "len": 0.8
            }]
        ),
        frames=frames
    )
    fig.update_layout(height=400, width=400)
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)


def render_dashboard_html(output_path, heatmap_html, scatter_html, brain_slider_html, ed_errors):
    stats = pd.Series(ed_errors).describe(percentiles=[.25, .5, .75]).to_dict()
    std = np.std(ed_errors)

    template = Template("""
    <html>
    <head>
        <title>AFID Registration QC Dashboard</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {
                font-family: 'Segoe UI', sans-serif;
                background-color: #f4f6f8;
                margin: 0;
                padding: 0;
                color: #2e2e2e;
            }
            header {
                background: white;
                padding: 40px 20px 20px;
                box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
                text-align: center;
            }
            header h1 {
                margin: 0;
                font-size: 32px;
                font-weight: 600;
                letter-spacing: -0.5px;
            }
            .container {
                padding: 40px;
                max-width: 1400px;
                margin: 0 auto;
            }
            .card {
                background: white;
                border-radius: 14px;
                box-shadow: 0 4px 16px rgba(0, 0, 0, 0.06);
                padding: 24px;
                margin-bottom: 40px;
                overflow: visible;
            }
            .section-title {
                font-size: 20px;
                font-weight: 600;
                margin-bottom: 16px;
                color: #1e1e1e;
                border-bottom: 2px solid #eaecef;
                padding-bottom: 8px;
            }
            .flex-row {
                display: flex;
                flex-direction: row;
                gap: 32px;
            }
            .tile-row {
                display: flex;
                flex-wrap: wrap;
                gap: 16px;
                margin-top: 16px;
            }
            .tile {
                background-color: #f0f2f5;
                border-radius: 10px;
                padding: 12px;
                text-align: center;
                flex: 1 1 120px;
                min-width: 120px;
            }
            .tile-label {
                font-weight: 500;
                font-size: 14px;
                color: #555;
            }
            .tile-value {
                font-size: 20px;
                font-weight: 600;
                color: #111;
                margin-top: 4px;
            }
            .viewer-container {
                width: 400px;
                flex-shrink: 0;
            }
            .plot-box {
                padding-top: 10px;
            }
            .responsive-stack {
                display: flex;
                flex-direction: column;
                gap: 20px;
            }
            @media (min-width: 1000px) {
                .responsive-stack {
                    flex-direction: row;
                }
            }
        </style>
    </head>
    <body>
        <header>
            <h1>AFID Registration QC Dashboard</h1>
        </header>
        <div class="container">

            <div class="card flex-row">
                <div style="flex: 1;">
                    <div class="section-title">Error Summary (ED)</div>
                    <div class="tile-row">
                        {% for label, value in stats.items() %}
                        <div class="tile">
                            <div class="tile-label">{{ label }}</div>
                            <div class="tile-value">{{ '%.2f' % value }} mm</div>
                        </div>
                        {% endfor %}
                        <div class="tile">
                            <div class="tile-label">Std Dev</div>
                            <div class="tile-value">{{ '%.2f' % std }} mm</div>
                        </div>
                    </div>
                </div>
                <div class="viewer-container plot-box">
                    {{ brain_slider|safe }}
                </div>
            </div>

            <div class="card">
                <div class="section-title">Error Heatmap and 3D Visualization</div>
                <div class="responsive-stack">
                    <div style="flex: 1;" class="plot-box">
                        {{ heatmap|safe }}
                    </div>
                    <div style="flex: 1.4;" class="plot-box">
                        {{ scatter|safe }}
                    </div>
                </div>
            </div>
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
            "Mean": stats["mean"],
        },
        std=std,
        heatmap=heatmap_html,
        scatter=scatter_html,
        brain_slider=brain_slider_html
    )

    Path(output_path).write_text(html)
    return output_path


def generate_afid_qc_dashboard(gt_fcsv_path, pred_fcsv_path, output_html_path):
    gt_coords, afid_ids = load_fcsv(gt_fcsv_path)
    pred_coords, _ = load_fcsv(pred_fcsv_path)

    dx, dy, dz, ed = compute_error_components(gt_coords, pred_coords)
    heatmap_html = make_toggleable_heatmap([dx, dy, dz, ed], afid_ids)
    scatter_html = make_3d_plot(gt_coords, pred_coords, afid_ids)
    brain_slider_html = make_brain_placeholder_slider()

    return render_dashboard_html(output_html_path, heatmap_html, scatter_html, brain_slider_html, ed)

if __name__ == "__main__":
    generate_afid_qc_dashboard(
        gt_fcsv_path = snakemake.input["afidfcsv"],
        pred_fcsv_path = snakemake.params["refcoord"],
        output_html_path=snakemake.output["html"]
    )
