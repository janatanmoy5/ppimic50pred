# Base image with conda
FROM continuumio/miniconda3:latest

# Create conda environment with Python 3.10
RUN conda create -y -n ic50env python=3.10

# Ensure all following commands run inside the conda environment
SHELL ["conda", "run", "-n", "ic50env", "/bin/bash", "-c"]

# Install RDKit from conda-forge
RUN conda install -y -c conda-forge rdkit

# Copy requirements and install pip dependencies
COPY requirements.txt /app/requirements.txt
WORKDIR /app
RUN pip install --no-cache-dir -r requirements.txt

# Copy all project files
COPY . /app

# Expose Streamlit port
EXPOSE 8501

# Start Streamlit inside the conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "ic50env", "streamlit", "run", "run.py", "--server.port=8501", "--server.address=0.0.0.0"]
