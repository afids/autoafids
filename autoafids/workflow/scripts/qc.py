import base64
import os
from io import BytesIO

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

# Dictionary for AFID labels
afids_labels = {
    1: "AC",
    2: "PC",
    3: "ICS",
    4: "PMJ",
    5: "SIPF",
    6: "RSLMS",
    7: "LSLMS",
    8: "RILMS",
    9: "LILMS",
    10: "CUL",
    11: "IMS",
    12: "RMB",
    13: "LMB",
    14: "PG",
    15: "RLVAC",
    16: "LLVAC",
    17: "RLVPC",
    18: "LLVPC",
    19: "GENU",
    20: "SPLE",
    21: "RALTH",
    22: "LALTH",
    23: "RSAMTH",
    24: "LSAMTH",
    25: "RIAMTH",
    26: "LIAMTH",
    27: "RIGO",
    28: "LIGO",
    29: "RVOH",
    30: "LVOH",
    31: "ROSF",
    32: "LOSF",
}


def generate_html_with_keypress(
    subject_images, reference_images, output_html="mri_viewer.html"
):
    """Generate an interactive HTML viewer with sticky instructions."""
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>MRI Viewer</title>
        <style>
            body { display: flex;
                   font-family: Arial,
                   sans-serif;
                   margin: 0;
                   padding: 0; }
            .instructions {
                width: 20%; padding: 20px; background-color: #f4f4f4;
                border-right: 2px solid #ddd;
                position: sticky; top: 0; height: 100vh; overflow-y: auto;
                box-sizing: border-box;
            }
            .instructions h2 { margin-top: 0; }
            .viewer { flex: 1; padding: 20px; text-align: center; }
            .slider { width: 80%; margin: 20px auto; }
            .image { display:
                     block;
                     margin: 0 auto;
                     max-width: 80%;
                     transition: opacity 0.2s ease-in-out; }
        </style>
    </head>
    <body>
        <div class="instructions">
            <h2>Tips & Tricks</h2>
            <ul>
                <li>Use the <strong>slider</strong> to navigate through MRI.
                Keyboard can also be used after clicking on slider.</li>
                <li>Check the <strong>Red</strong> crosshair which represents
                the placement of a given landmark.</li>
                <li>Press <strong>R</strong> to toggle between the subject
                scan and the protocol-defined placement (if one exists).</li>
                <li>Use the <strong>TAB</strong> key to navigate to
                next slider & MRI.</li>
            </ul>
        </div>
        <div class="viewer">
    """

    for label, images in subject_images.items():
        num_slices = len(images)
        has_reference = reference_images is not None and label in reference_images

        html_content += f"""
        <div class="container">
            <h2>Landmark: {afids_labels[label]}</h2>
            <img id="{label}_image" class="image"
            src="data:image/png;base64,{images[0]}"
            alt="MRI Slice">
            <input type="range" min="0"
            max="{num_slices - 1}" value="0"
            class="slider" id="{label}_slider">
        </div>
        """

    html_content += """
        </div>
        <script>
            const subjects = {};
    """

    for label, images in subject_images.items():
        has_reference = reference_images is not None and label in reference_images
        html_content += f"""
            subjects["{label}"] = {{
                targetImages: {images},
                referenceImages: {
                    reference_images[label] if has_reference else "null"
                    },
                showingReference: false
            }};
        """

    html_content += """
            document.addEventListener('DOMContentLoaded', () => {
                for (const [label, data] of Object.entries(subjects)) {
                    const slider = document.getElementById(`${label}_slider`);
                    const image = document.getElementById(`${label}_image`);

                    slider.addEventListener('input', () => updateImage(label));
                    document.addEventListener('keydown', (event) => {
                        if (event.key.toLowerCase() === 'r'
                            && data.referenceImages) {
                            data.showingReference = !data.showingReference;
                            updateImage(label);
                        }
                    });

                    function updateImage(label) {
                        const sliceIndex = document.getElementById(
                        `${label}_slider`
                        ).value;
                        const imageArray = subjects[label].showingReference
                            ? subjects[label].referenceImages
                            : subjects[label].targetImages;
                        document.getElementById(`${label}_image`).src =
                        'data:image/png;base64,' + imageArray[sliceIndex];
                    }
                }
            });
        </script>
    </body>
    </html>
    """

    with open(output_html, "w") as f:
        f.write(html_content)


def save_single_slice_in_memory(data, x, y, z, offset, zoom_radius, show_crosshairs):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (slice_data, coord, title) in enumerate(
        zip(
            [data[x + offset, :, :], data[:, y + offset, :], data[:, :, z + offset]],
            [(y, z), (x, z), (x, y)],
            ["Sagittal", "Coronal", "Axial"],
        )
    ):
        axes[i].imshow(slice_data.T, origin="lower", cmap="gray")
        if offset == 0 and show_crosshairs:
            axes[i].axhline(y=coord[1], color="r", lw=1)
            axes[i].axvline(x=coord[0], color="r", lw=1)
        axes[i].set_xlim(coord[0] - zoom_radius, coord[0] + zoom_radius)
        axes[i].set_ylim(coord[1] - zoom_radius, coord[1] + zoom_radius)
        axes[i].set_title(title)
        axes[i].axis("off")

    buffer = BytesIO()
    plt.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)

    return base64.b64encode(buffer.read()).decode()


def save_mri_slices_as_images(data, x, y, z, jitter, zoom_radius, show_crosshairs=True):
    return Parallel(n_jobs=-1)(
        delayed(save_single_slice_in_memory)(
            data, x, y, z, offset, zoom_radius, show_crosshairs
        )
        for offset in range(-jitter, jitter + 1)
    )


def extract_coordinates_from_fcsv(file_path, label_description):
    df = pd.read_csv(file_path, comment="#", header=None)
    row = df[df[11] == label_description]
    return tuple(row.iloc[0, 1:4]) if not row.empty else None


def generate_interactive_mri_html(
    nii_path,
    fcsv_path,
    labels,
    ref_nii_path=None,
    ref_fcsv_path=None,
    jitter=2,
    zoom_radius=20,
    out_file_prefix="mri_viewer",
):
    """
    Generates an interactive HTML viewer for
    MRI slices based on fiducial coordinates for a single subject.

    Parameters:
    - nii_path: str
        Full path to the subject's NIfTI (.nii.gz) file.
    - fcsv_path: str
        Full path to the subject's FCSV (.fcsv) file.
    - labels: list of int
        Landmark indices to extract for the subject.
    - ref_nii_path: str or None
        Full path to reference NIfTI file (optional).
    - ref_fcsv_path: str or None
        Full path to reference FCSV file (optional).
    - jitter: int, optional (default=2)
        The number of pixels to expand around
        the coordinate in each slice for visualization.
    - zoom_radius: int, optional (default=20)
        The radius (in pixels) around the coordinate to extract for display.
    - out_file_prefix: str, optional (default="mri_viewer")
        The prefix for the output HTML file name.

    Notes:
    - Extracts specified fiducial coordinates from the subject's .fcsv file.
    - Maps the coordinates from world space to voxel
        space using the NIfTI affine transformation.
    - MRI slices centered around these coordinates are saved as images.
    - If a reference MRI and coordinate file are provided,
        the same process is applied.
    - Generates an interactive HTML viewer allowing keypress navigation.

    Returns:
    - None (outputs an HTML file with the
        interactive MRI viewer for the subject).
    """
    subject_key = os.path.basename(nii_path).replace(".nii.gz", "")
    target_img = nib.load(nii_path, mmap=True)
    affine_inv = np.linalg.inv(target_img.affine)
    target_data = target_img.get_fdata(dtype=np.float32)

    target_images = {}  # Dictionary to store images per landmark

    for label in labels:
        target_coord = extract_coordinates_from_fcsv(fcsv_path, label)
        if not target_coord:
            print(
                f"Coordinates for label '{label}'"
                + f"not found in subject '{subject_key}'."
            )
            continue

        target_voxel = np.round(affine_inv.dot((*target_coord, 1))).astype(int)
        target_images[label] = save_mri_slices_as_images(
            target_data, *target_voxel[:3], jitter, zoom_radius
        )

    reference_images = None
    if ref_nii_path and ref_fcsv_path:
        ref_img = nib.load(ref_nii_path, mmap=True)
        ref_data = ref_img.get_fdata(dtype=np.float32)
        reference_images = {}

        for label in labels:
            ref_coord = extract_coordinates_from_fcsv(ref_fcsv_path, label)
            if ref_coord:
                ref_voxel = np.round(
                    np.linalg.inv(ref_img.affine).dot((*ref_coord, 1))
                ).astype(int)
                reference_images[label] = save_mri_slices_as_images(
                    ref_data, *ref_voxel[:3], jitter, zoom_radius
                )

    out_file = f"{out_file_prefix}"

    # Pass the correct image dictionaries
    generate_html_with_keypress(target_images, reference_images, out_file)


if __name__ == "__main__":
    generate_interactive_mri_html(
        nii_path=snakemake.input["im"],
        fcsv_path=snakemake.input["afidfcsv"],
        labels=[
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
        ],
        ref_nii_path=snakemake.params["refim"],
        ref_fcsv_path=snakemake.params["refcoord"],
        out_file_prefix=snakemake.output["html"],
    )
