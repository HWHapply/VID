# Use the official Anaconda3 image with Ubuntu 22.04
FROM rapidsai/ci-conda:cuda12.2.2-ubuntu22.04-py3.11-amd64

# Set environment variable to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install git
RUN apt-get update && \
    apt-get install -y git

# Update the package list and install required system dependencies for R
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean

# RUN apt-get update && apt-get install -y \
#     libhdf5-dev \
#     libxml2-dev \
#     libssl-dev \
#     libcurl4-openssl-dev \
#     libglpk-dev \
#     libgmp-dev \
#     libigraph-dev \
#     libzmq3-dev \
#     libboost-all-dev \
#     libfftw3-dev \
#     libpng-dev \
#     libx11-dev \
#     libcairo2-dev \
#     libgit2-dev \
#     git \
#     && rm -rf /var/lib/apt/lists/*

# Set up Conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels bioconda

# Clone the GitHub repository into /opt/VID
RUN git clone https://github.com/HWHapply/VID.git /opt/VID

# Set the working directory to the cloned repository
WORKDIR /opt/VID

# Copy the entrypoint.sh script to /opt/VID
COPY entrypoint.sh /opt/VID/entrypoint.sh

# Create the Conda environment from the env.yml file
RUN conda env create -f ./env/vid_env_Linux.yml -n vid_env 

# Install SeuratDisk
RUN git clone https://github.com/mojaveazure/seurat-disk /opt/seurat-disk
RUN conda run -n vid_env Rscript -e "install.packages('/opt/seurat-disk', repos = NULL, type = 'source')"

# Ensure all files under /opt/VID are executable
RUN chmod -R +x /opt/VID

# Set the default command to run the tool with Conda environment
ENTRYPOINT ["conda", "run", "-n", "vid_env", "bash", "/opt/VID/entrypoint.sh"]

