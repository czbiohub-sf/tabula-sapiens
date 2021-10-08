FROM continuumio/miniconda3

ENTRYPOINT ["/bin/bash", "-c"]

RUN mkdir /usr/src/app
WORKDIR /usr/src/app
COPY environment.yml /usr/src/app
COPY . /usr/src/app

RUN conda env create -f environment.yml

CMD ["source activate backend && flask run"]

