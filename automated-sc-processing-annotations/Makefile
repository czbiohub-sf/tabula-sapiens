VERBOSE=1
BACKEND_DOCKER_IMAGE = czbiohub/sc-rna-seq-processing
FRONTEND_DOCKER_IMAGE = frontend
ANNOTATION_DOCKER_IMAGE = sapiens-annotation-notebook

.PHONY: notebooks

build-annotation-notebook:
	cd notebooks && docker build --tag $(ANNOTATION_DOCKER_IMAGE) .

build-backend:
	cd backend && docker build --tag $(BACKEND_DOCKER_IMAGE) .

build-frontend:
	cd frontend && docker build --tag $(FRONTEND_DOCKER_IMAGE) .

build: build-backend build-frontend

get-data:
	bash cd backend/download-data.sh

run-annotation-notebook: stop-annotation-notebook build-annotation-notebook
	docker run -it \
		   --name $(ANNOTATION_DOCKER_IMAGE) \
		   --rm -p 8888:8888 \
		   -e JUPYTER_ENABLE_LAB=yes \
		   -v "$$PWD"/notebooks:/home/jovyan/work \
		   -v "$$PWD/backend/data:/home/jovyan/work/data" \
		   $(ANNOTATION_DOCKER_IMAGE)

stop-annotation-notebook:
	- docker stop $(ANNOTATION_DOCKER_IMAGE)

up: down
	docker-compose up --build -d && \
	echo "Go to: http://localhost:3000"

down:
	docker-compose down
