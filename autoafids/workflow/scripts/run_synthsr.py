import subprocess
import os
import logging

# clone SynthSR repo
def clone_repo(download_dir, log_file):
    repo = "https://github.com/mackenziesnyder/SynthSR"
    branch = "mackenzie/scaled-down-repo"

    repo_dir = os.path.join(download_dir, "SynthSR")
    if os.path.exists(repo_dir):
        logging.info(f"SynthSR repository already exists at {repo_dir}. Skipping clone.")
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

# make virtual environment with python3.8
def make_virtual_env(repo_dir, log_file):
    venv_dir = os.path.join(repo_dir, "venv")
    with open(log_file, "a") as log:
        try:
            subprocess.run(
                ["python3", "-m", "venv", venv_dir], 
                check=True,
                stdout=log,
                stderr=log
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Error making virtual environment for SynthSR repo: {e}")
            raise
    return venv_dir

# install dependencies in the virtual env for SynthSR
def install_dependencies(venv_dir, log_file):
    venv_python = os.path.join(venv_dir, "bin", "python")
    with open(log_file, "a") as log:
        try:
            subprocess.run(
                [
                    venv_python, "-m", "pip", "install",
                    "tensorflow==2.12.0",
                    "keras==2.12.0",
                    "protobuf==3.20.3",
                    "numpy==1.23.5",
                    "nibabel==5.0.1",
                    "matplotlib==3.6.2"
                ], 
                check=True,
                stdout=log,
                stderr=log
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Error installing dependencies for SynthSR repo: {e}")
            raise
    return venv_python

# run SynthSR
def run_program(venv_python, repo_dir, input_img, output_img, modality, log_file):
    predict_script = os.path.join(repo_dir, "scripts", "predict_command_line.py")
    cmd = [venv_python, predict_script, input_img, output_img, "--cpu"]
    
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
            logging.info(f"SynthSR processing completed for {input_img}. Output saved to {output_img}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running SynthSR: {e}")
            raise
            
# main function 
def main():
    log_file=snakemake.log[0]
    download_dir=snakemake.params["download_dir"]
    repo_dir = clone_repo(download_dir, log_file)
    venv_dir = make_virtual_env(repo_dir, log_file)
    venv_python = install_dependencies(venv_dir, log_file)
    
    run_program(
        venv_python=venv_python,
        repo_dir=repo_dir,
        input_img=snakemake.input["im"],
        output_img=snakemake.output["SynthSR"],
        modality=snakemake.params["modality"],
        log_file=log_file
    )

if __name__ == "__main__":
    main()


