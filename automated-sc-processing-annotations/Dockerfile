FROM continuumio/miniconda3

RUN apt-get -y update && \
    apt-get -y install \
	build-essential

COPY environment.yml .
RUN conda env create -f environment.yml
