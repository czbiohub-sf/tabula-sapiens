FROM continuumio/miniconda3

ENTRYPOINT ["/bin/bash", "-c"]

RUN apt-get -y update && \
    apt-get -y install \
	build-essential

COPY environment.yml .
RUN conda env create -f environment.yml

RUN [ "/bin/bash", "-c", "source activate tabula_sapiens && python -m ipykernel install --user --name=tabula_sapiens --display-name=tabula_sapiens"]
CMD [ "source activate tabula_sapiens && \
       jupyter lab --port=8889 --no-browser --ip=0.0.0.0 --allow-root"]
