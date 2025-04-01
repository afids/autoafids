import logging
import os
import subprocess


# clone SynthSR repo
def clone_repo(download_dir, log_file):
    repo = "https://github.com/mackenziesnyder/SynthSR"
    branch = "mackenzie/scaled-down-repo"

    repo_dir = os.path.join(download_dir, "SynthSR")
    if os.path.exists(repo_dir):
        logging.info(
            f"SynthSR repository already exists at {repo_dir}. Skipping clone."
        )
        return repo_dir

    with open(log_file, "a") as log:
        try:
            subprocess.run(
                ["git", "clone", "--branch", branch, repo, repo_dir],
                check=True,
                stdout=log,
                stderr=log
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Error cloning SynthSR repo: {e}")
            raise

    return repo_dir

# run SynthSR
def run_program(
        repo_dir,
        input_img,
        output_img,
        modality,
        log_file
    ):
    predict_script = os.path.join(
        repo_dir,
        "scripts",
        "predict_command_line.py"
    )
    cmd = [predict_script, input_img, output_img, "--cpu"]

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
repo_dir = clone_repo(download_dir, log_file)
# venv_dir = make_virtual_env(repo_dir, log_file)
# venv_python = install_dependencies(venv_dir, log_file)

run_program(
    repo_dir=repo_dir,
    input_img=snakemake.input["im"],
    output_img=snakemake.output["SynthSR"],
    modality=snakemake.params["modality"],
    log_file=log_file
)



