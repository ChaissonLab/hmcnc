FROM continuumio/miniconda3:latest

# Install basic build tools
RUN apt-get update && apt-get install -y build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set up conda channels and create environment in one step
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n hmcnc_env bedtools samtools snakemake boost r-base gxx_linux-64 tabix htslib -y && \
    echo "conda activate hmcnc_env" >> ~/.bashrc

# Set the working directory
WORKDIR /app

# Copy the source code and build files
COPY src/ /app/src/
COPY include/ /app/include/
COPY HMM/ /app/HMM/

# Set the shell to use bash and activate conda environment
SHELL ["/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/hmcnc_env/bin:$PATH


RUN cd /app/src && conda run -n hmcnc_env make all
RUN cd /app/src && ./hmcnc --help

# Default command: interactive bash with environment activated
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "hmcnc_env", "/bin/bash"]