import subprocess
import os
import tempfile

# clone SynthSR repo
def clone_repo():
    repo = "https://github.com/mackenziesnyder/SynthSR"
    branch = "main"

    temp_dir = tempfile.mkdtemp()
    repo_dir = os.path.join(temp_dir, "SynthSR")

    subprocess.run(["git", "clone", "--branch", branch, repo, repo_dir], check=True)
    return repo_dir

# make virtual environment with python3.8
def make_virtual_env(repo_dir):
    venv_dir = os.path.join(repo_dir, "venv")
    subprocess.run(["python3.8", "-m", "venv", venv_dir], check=True)
    return venv_dir

# install dependencies in the virtual env for SynthSR
def install_dependencies(venv_dir):
    venv_python = os.path.join(venv_dir, "bin", "python")
    subprocess.run([
        venv_python, "-m", "pip", "install",
        "tensorflow==2.12.0",
        "keras==2.12.0",
        "protobuf==3.20.3",
        "numpy==1.23.5",
        "nibabel==5.0.1",
        "matplotlib==3.6.2"
    ], check=True)
    return venv_python

# run SynthSR
def run_program(venv_python, repo_dir, input_img, output_img):
    predict_script = os.path.join(repo_dir, "scripts", "predict_command_line.py")
    subprocess.run([venv_python, predict_script, input_img, output_img, "--cpu"], check=True)

# main function 
def main():
    repo_dir = clone_repo()
    venv_dir = make_virtual_env(repo_dir)
    venv_python = install_dependencies(venv_dir)
    
    run_program(
        venv_python=venv_python,
        repo_dir=repo_dir,
        input_img=snakemake.input["im"],
        output_img=snakemake.output["SynthSR"]
    )

if __name__ == "__main__":
    main()


