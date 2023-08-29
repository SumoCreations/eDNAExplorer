# Use an official Ubuntu as a parent image
FROM ubuntu:latest

# Run package updates and install packages
RUN apt-get clean && \
  apt-get update --fix-missing && \
  apt-get upgrade -y && \
  apt-get install -y wget bzip2

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
  bash ~/miniconda.sh -b -p /opt/conda && \
  rm ~/miniconda.sh

# Add conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Copy env.yaml and install.R to the image
COPY env.yml /tmp/env.yml
COPY install.R /tmp/install.R

# Create a Conda environment from env.yml
RUN conda env create -f /tmp/env.yml

# Activate the Conda environment and run install.R
RUN echo "source activate $(head -1 /tmp/env.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/env.yml | cut -d' ' -f2)/bin:$PATH
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate GBIF_env && Rscript /tmp/install.R"