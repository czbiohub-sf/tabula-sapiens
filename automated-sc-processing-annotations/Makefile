VERBOSE=1
BACKEND_DOCKER_IMAGE = czbiohub/sc-rna-seq-processing
FRONTEND_DOCKER_IMAGE = frontend
ANNOTATION_DOCKER_IMAGE = sapiens-annotation-notebook
NOTEBOOK_PORT = 8888

.PHONY: notebooks

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

run-annotation-notebook: stop-annotation-notebook build-annotation-notebook get-data
	docker run -t \
		   --name $(ANNOTATION_DOCKER_IMAGE) \
		   --rm -p $(NOTEBOOK_PORT):8888 \
		   -e JUPYTER_ENABLE_LAB=yes \
		   -e NB_PORT=$(NOTEBOOK_PORT) \
		   -v "$$PWD"/notebooks:/home/jovyan/work \
		   -v "$$PWD/backend/data:/home/jovyan/work/data" \
		   $(ANNOTATION_DOCKER_IMAGE) | sed 's/:8888/:$(NOTEBOOK_PORT)/g'

stop-annotation-notebook:
	- docker stop $(ANNOTATION_DOCKER_IMAGE)

up: down
	docker-compose up --build -d && \
	echo "Go to: http://localhost:3000"

down:
	docker-compose down
