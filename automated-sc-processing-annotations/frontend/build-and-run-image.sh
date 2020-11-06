#!/usr/bin/env bash
set -o errexit
set -o pipefail

NAME=frontend

# to create a new container
docker build --tag $NAME .

# push to dockerhub
echo "if ready to push run 'docker push $NAME' "
# docker push $NAME

#run the container
docker run -p 3000:3000 -v $(pwd):/home/ $NAME 
