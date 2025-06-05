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

log_file=snakemake.log[0]
download_dir=snakemake.params["download_dir"]
repo_dir = clone_repo(download_dir, log_file)