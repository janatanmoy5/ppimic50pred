# Use Miniconda base image
FROM continuumio/miniconda3:latest

# Create environment
RUN conda create -y -n ic50env python=3.10
SHELL ["conda", "run", "-n", "ic50env", "/bin/bash", "-c"]

# Install RDKit from conda-forge
RUN conda install -y -c conda-forge rdkit

# Install pip dependencies
COPY requirements.txt /app/requirements.txt
WORKDIR /app
RUN pip install --no-cache-dir -r requirements.txt

# Copy app files
COPY . /app

# Expose Streamlit port
EXPOSE 8501

# Streamlit entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "ic50env", "streamlit", "run", "run.py", "--server.port=8501", "--server.address=0.0.0.0"]
