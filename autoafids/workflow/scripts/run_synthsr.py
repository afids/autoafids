import logging
import os
import subprocess


# run SynthSR
def run_program(
        input_img,
        output_img,
        modality,
        log_file
    ):

    cmd = ["python", "-m", "SynthSR.scripts.predict_command_line", input_img, output_img, "--cpu"]

    if modality == "ct":
        cmd.append["--ct"]

    with open(log_file, "a") as log:
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=log,
                stderr=log
            )
            logging.info(
                f"SynthSR processing completed for {input_img}."
                f"Output saved to {output_img}"
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running SynthSR: {e}")
            raise


log_file=snakemake.log[0]
download_dir=snakemake.params["download_dir"]

run_program(
    input_img=snakemake.input["im"],
    output_img=snakemake.output["SynthSR"],
    modality=snakemake.params["modality"],
    log_file=log_file
)



