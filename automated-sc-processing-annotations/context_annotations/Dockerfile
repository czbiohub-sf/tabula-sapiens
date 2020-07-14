FROM continuumio/miniconda3

# Update and install should be on the same line:
# If you add a library, you don't want an old cache of the package list.

# apt-get for python-igraph taken from:
# https://github.com/igraph/python-igraph/blob/master/docker/jupyter/Dockerfile

RUN apt-get -y update && \
    apt-get -y install \
      build-essential \
      libcairo-dev \
      libffi-dev \
      libxml2-dev \
      pkg-config \
      python-dev \
      python-pip \
      zlib1g-dev

COPY requirements-freeze.txt .
RUN pip install  -r ./requirements-freeze.txt


# In development, you may want to pin a single dependency in requirements.txt,
# without throwing away the entire cache layer from requirements-freeze.txt.
# (But once it works, you should check in an updated freeze!)

COPY requirements.txt .
RUN pip install  -r ./requirements.txt

RUN pip freeze  # Copy this to requirements-freeze.txt if there are changes.

COPY . .
CMD [ "python", "onclass_annotations.py", \
      "--input", "/input.h5ad", \
      "--output", "/output.h5ad", \
      "--data_file","/data_file.h5ad", \
      "--cell_type_network_file", "/cell_type_network_file.ontology", \
      "--use_pretrain_emb","/pretrain/tp2emb_500", \
      "--name_mapping_file", "/cl.obo", \
      "--use_pretrain_data", "/pretrain/BilinearNN_50019", \
      "--use_pretrain_data_expression", "/pretrain/BilinearNN_500"]
