FROM jupyter/base-notebook

USER root
RUN apt-get update && \
	apt-get install -y build-essential git && \
	rm -rf /var/lib/apt/lists/*
USER jovyan

RUN pip install --quiet --upgrade \
	"tensorflow>=1.15" \
	git+https://github.com/wangshenguiuc/OnClass@master \
	"git+https://github.com/yoseflab/scvi-tools@master#egg=scvi-tools[tutorials]"

RUN pip install --quiet --upgrade --no-deps scanorama

WORKDIR /home/jovyan/work
