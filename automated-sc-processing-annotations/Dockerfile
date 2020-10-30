FROM continuumio/miniconda3

ENTRYPOINT ["/bin/bash", "-c"]
RUN apt-get -y update && \
    apt-get -y upgrade && \
    apt-get -y install \
	build-essential \
        libcurl4-gnutls-dev \
 	libxml2-dev \
 	libssl-dev \
 	libcairo-dev \
 	libffi-dev \
 	pkg-config \
 	python-dev \
 	python-pip \
 	zlib1g-dev \
	apt-utils \
	unzip \
	nodejs \
	npm

RUN apt-get -y install software-properties-common && \
	apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian buster-cran35/' && \
	apt-get -y update && \
	apt-get install -y r-base

 	#r-base
#r-base-core \
#r-recommended \ 
#r-base-html \


#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && apt update && apt install r-base

#RUN sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/' && \
#	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
#	apt update && \
#	apt install r-base

#RUN pip install jupyterlab
COPY environment.yml .
#COPY tmp.yml environment.yml
RUN conda env create -f environment.yml
RUN R -e "install.packages('devtools')" &&\
    R -e "library(devtools);  devtools::install_github('pcahan1/singlecellnet')"

#SHELL ["conda", "run", "-n", "tabula_sapiens", "/bin/bash", "-c"]
#RUN R -e "install.packages('devtools'); library(devtools); devtools::install_github('pcahan1/singlecellnet')"

#RUN /bin/bash -c "source activate tabula_sapiens &&  conda install r-base r-essentials r-devtools r-dplyr rpy2"

#run /bin/bash -c "source activate tabula_sapiens &&  R -e 'version'"
#run /bin/bash -c "export tar=/bin/tar && source activate tabula_sapiens &&  R -e 'library(devtools); devtools::install_github(\"pcahan1/singlecellnet\")'"
#RUN R -e "library('devtools'); devtools::install_github('pcahan1/singleCellNet')"
#RUN R -e "library('devtools'); devtools::install_github('pcahan1/singleCellNet')"
#RUN R -e "install.packages('devtools');library('devtools'); devtools::install_github('pcahan1/singleCellNet')"
#RUN R -e "install.packages('devtools'); devtools::install_github('pcahan1/singleCellNet')"

RUN [ "/bin/bash", "-c", "source activate tabula_sapiens && python -m ipykernel install --user --name=tabula_sapiens --display-name=tabula_sapiens && jupyter nbextension enable --py widgetsnbextension && jupyter labextension install @jupyter-widgets/jupyterlab-manager"]
#RUN [ "/bin/bash", "-c", "source activate tabula_sapiens && python -m ipykernel install --user --name=tabula_sapiens --display-name=tabula_sapiens"]
CMD [ "cd /home/ && source activate tabula_sapiens && \
       jupyter lab --port=8888 --no-browser --ip=0.0.0.0 --allow-root"]
