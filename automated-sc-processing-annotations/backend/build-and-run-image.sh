#!/usr/bin/env bash
set -o errexit
set -o pipefail

NAME=czbiohub/sc-rna-seq-processing

# to create a new container
docker build --tag $NAME .

# push to dockerhub
echo "if ready to push run 'docker push $NAME' "
# docker push $NAME

#run the container
docker run -p 8888:8888 -v $(pwd):/home/ $NAME 
