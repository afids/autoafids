#!/usr/bin/env python3
import os
import shutil
import tempfile
import subprocess
import sys
from argparse import ArgumentParser
from pathlib import Path

IMAGE_MODALITY = {
    "T1w": {
        "suffix": "T1w",
        "use_derivatives": False
    },
    "T2w": {
        "suffix": "T2w",
        "use_derivatives": False
    },
    "FLAIR": {
        "suffix": "flair",
        "use_derivatives": False
    },
    "CT": {
        "suffix": "ct",
        "use_derivatives": False
    }
}

def check_conda_installation():
    try:
        conda_version = subprocess.run(
            ["conda", "--version"], capture_output=True, text=True, check=True
        )
        print(f"Conda is installed (version: {conda_version.stdout.strip()})")
        return True
    except FileNotFoundError:
        error_message = "Conda is not installed on your system or not found in PATH"
        print(error_message)
        return False
    
def gen_parser():
    parser = ArgumentParser(description="Run autoafids for a single subject.")

    # command line arguments for autoafids-quick
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="File path to your input NIfTI image. Must have .nii.gz extension.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to your desired output folder."
    )
    parser.add_argument(
        "-s", "--subject", required=True, help="SUBJECT ID (e.g., 001)"
    )
    parser.add_argument(
        "-m",
        "--modality",
        required=True,
        choices=list(
            IMAGE_MODALITY.keys()
        ),
        help="Image modality - chose between: " + ", ".join(IMAGE_MODALITY.keys())
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Execute a dry run without actually running the full pipeline."
    )
    
    return parser


def main():

    script_path = Path(__file__).resolve().parent / "run.py"

    if not check_conda_installation():
        print("Please activate conda and install snakebids to continue using autoafids quick.")
        sys.exit(1)

    if "SNAKEMAKE_PROFILE" in os.environ:
        del os.environ["SNAKEMAKE_PROFILE"]

    args = gen_parser().parse_args()

    with tempfile.TemporaryDirectory(prefix=None) as temp_dir:

        # create temporary input directory within temp_dir
        temp_input_dir = Path(temp_dir) / "input"
        temp_input_dir.mkdir(parents=True, exist_ok=True)

        temp_output_dir = Path(temp_dir) / "output"

        # FIX 1: Correct BIDS structure -> sub-<subject>/anat
        subject_folder = temp_input_dir / f"sub-{args.subject}" / "anat"
        subject_folder.mkdir(parents=True, exist_ok=True)

        # FIX 2: Fixed nested double quote syntax in f-string
        input_filename = (
            f"sub-{args.subject}_{IMAGE_MODALITY[args.modality]['suffix']}.nii.gz"
        )

        # add file name to the created subject folder
        temp_input_file = subject_folder / input_filename

        # copy the input file
        shutil.copy(args.input, temp_input_file)

        # run autoafids
        command = [
            str(script_path),
            str(temp_input_dir),
            str(temp_output_dir),
            "participant",
            "-c",
            "all",
            "--force-output",
            "--nolock",
            "--modality",
            args.modality,
            "--detect-with-nnlm"
        ]

        if args.dry_run:
            command.append("-n")

        if IMAGE_MODALITY[args.modality]["use_derivatives"]:
            # need to have desc file in bids dir to use --derivatives
            desc_file = temp_input_dir / "dataset_description.json"
            desc_file.write_text(
                '{"Name": "Generated Derivatives", '
                '"BIDSVersion": "1.0.2", '
                '"GeneratedBy": [{"Name": "autoafids-quick"}]}'
            )
            command.append("--derivatives")
            command.append(str(temp_input_dir))

        try:
            subprocess.run(command, check=True)
            print("autoafids ran successfully.")

            # copy results from temp output to final output
            temp_subject_dir = temp_output_dir / f"sub-{args.subject}"
            final_output_dir = Path(args.output)

            if temp_subject_dir.exists():
                # create the final output structure
                final_subject_dir = final_output_dir / f"sub-{args.subject}"

                # copy the subject directory
                if final_subject_dir.exists():
                    shutil.rmtree(final_subject_dir)
                final_output_dir.mkdir(parents=True, exist_ok=True)
                shutil.copytree(temp_subject_dir, final_subject_dir)
                print(f"Results copied to {final_subject_dir}")
            else:
                print(
                    f"Warning: Expected output directory {temp_subject_dir} not found."
                ) 

        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")

if __name__ == "__main__":
    main()