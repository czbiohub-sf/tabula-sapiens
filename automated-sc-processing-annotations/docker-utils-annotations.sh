#!/usr/bin/env bash
set -o errexit
set -o pipefail

NAME=czbiohub/sc-rna-seq-annotating:0.0.1

# to create a new container
docker build ./context_annotations/ --tag $NAME

# push to dockerhub
echo "if ready to push run 'docker push $NAME' "
# docker push $NAME

# run inside the container
echo "previous file will be re-written"
touch $PWD/data/output/output_processed_annotated.h5ad;
docker run \
  --mount type=bind,source=$PWD/data/output/output_processed.h5ad,target=/input.h5ad \
  --mount type=bind,source=$PWD/data/output/output_processed_annotated.h5ad,target=/output.h5ad \
  --mount type=bind,source=$PWD/data/OnClass_data/data_used_for_training/tabula-muris-senis-facs_cell_ontology.h5ad,target=/data_file.h5ad \
  --mount type=bind,source=$PWD/data/OnClass_data/OnClass_data_others/cell_ontology/cl.ontology,target=/cell_type_network_file.ontology \
  --mount type=bind,source=$PWD/data/OnClass_data/cell_ontology/cl.obo,target=/cl.obo \
  --mount type=bind,source=$PWD/data/OnClass_data/pretrain/,target=/pretrain \
  $NAME

  # python3 onclass_annotations.py --input "$PWD/data/output/output_processed.h5ad" \
  # --output "$PWD/data/output/output_processed_annotated.h5ad" \
  # --data_file "$PWD/data/OnClass_data/data_used_for_training/tabula-muris-senis-facs_cell_ontology.h5ad" \
  # --cell_type_network_file "$PWD/data/OnClass_data/OnClass_data_others/cell_ontology/cl.ontology" \
  # --name_mapping_file "$PWD/data/OnClass_data/cell_ontology/cl.obo" \
  # --use_pretrain_emb "$PWD/data/OnClass_data/pretrain/tp2emb_500" \
  # --use_pretrain_data "$PWD/data/OnClass_data/pretrain/BilinearNN_50019" \
  # --use_pretrain_data_expression "$PWD/data/OnClass_data/pretrain/BilinearNN_500"


# to make executable: chmod a+x docker-utils-annotations.sh
