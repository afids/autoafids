# Running Autoafids with Docker

Brief walkthrough on how to prepare the `autoafids` package in a Docker container.

## Local Build Instructions
To build a docker image of the project locally, do the following:
1. Download the git repo:

```
git clone https://github.com/afids/autoafids.git
```
2. Ensure you have installed [Docker](https://docs.docker.com/engine/install/) and that it's up and running.

3. Build the docker image using the following command; it should be executed within the repository folder (i.e., autoafids/). 

```
docker build -t autoafids:<VERSION> -f ./docker/Dockerfile .
```

4. Once the image is built, you can then run the pipeline as follows:

```
docker run -it --rm \                     
-v local/path/to/dataset:/bids_dataset \ 
-v local/path/to/dataset/derivatives:/output \  
autoafids:<VERSION> \
/bids_dataset /output participant --cores all
```

The first time building this image may take some time because poetry takes a while to install all the required python packages. 