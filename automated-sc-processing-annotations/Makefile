VERBOSE=1
BACKEND_DOCKER_IMAGE = czbiohub/sc-rna-seq-processing
FRONTEND_DOCKER_IMAGE = frontend
ANNOTATION_DOCKER_IMAGE = sapiens-annotation-notebook
NOTEBOOK_PORT = 8888

.PHONY: notebooks

check-env:
ifndef NAME
	$(error You must supply a name for the docker container by setting NAME)
endif

build-annotation-notebook:
	cd notebooks && docker build --tag $(ANNOTATION_DOCKER_IMAGE) .

build-backend:
	cd backend && docker build --tag $(BACKEND_DOCKER_IMAGE) .

build-frontend:
	cd frontend && docker build --tag $(FRONTEND_DOCKER_IMAGE) .

build: build-backend build-frontend

get-data: backend/data

backend/data:
	cd backend && bash download-data.sh

run-annotation-notebook: check-env stop-annotation-notebook build-annotation-notebook get-data
	docker run -t \
		   -u root \
		   --name $(NAME) \
		   --rm -p $(NOTEBOOK_PORT):8888 \
		   -e NB_UID=$$(id -u $$USER) \
		   -e NB_GID=$$(id -g $$USER) \
		   -e JUPYTER_ENABLE_LAB=yes \
		   -e NB_PORT=$(NOTEBOOK_PORT) \
		   -v "$$PWD"/notebooks:/home/jovyan/work:rw \
		   -v "$$PWD/backend/data:/home/jovyan/work/data":rw \
		   $(ANNOTATION_DOCKER_IMAGE) | sed 's/:8888/:$(NOTEBOOK_PORT)/g'

stop-annotation-notebook: check-env
	- docker stop $(NAME)

up: down
	docker-compose up --build -d && \
	echo "Go to: http://localhost:3000"

down:
	docker-compose down
