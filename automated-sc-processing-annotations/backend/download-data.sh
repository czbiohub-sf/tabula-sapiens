#!/usr/bin/env bash
set -o errexit
set -o pipefail

OTHERS='data/OnClass_data/OnClass_data_others'
[ -d $OTHERS ] || mkdir -p $OTHERS
TRAINING='data/OnClass_data/data_used_for_training'
[ -d $TRAINING ] || mkdir -p $TRAINING

# Train dataset
[ -e data/adata_small_test.h5ad ] \
  || wget https://files.figshare.com/21084318/adata_small_test.h5ad \
      -O ./data/adata_small_test.h5ad


# OnClass pretrained model
[ -e data/OnClass_data/pretrain/tp2emb_500X.npy ] \
  || wget https://files.figshare.com/21009027/pretrain.tar.gz -O - \
  | tar -xz -C ./data/OnClass_data/
# Cell type embeddings, marker genes and 26-datasets
[ -e $OTHERS/cell_ontology/cl.ontology.gz ] \
  || wget https://files.figshare.com/20157284/OnClass_data_others.tar.gz -O - \
  | tar -xz -C ./data/OnClass_data/
# FACS cells used in the study
[ -e $TRAINING/tabula-muris-senis-facs_cell_ontology.h5ad ] \
  || wget https://files.figshare.com/20157215/tabula-muris-senis-facs_cell_ontology.h5ad \
      -O $TRAINING/tabula-muris-senis-facs_cell_ontology.h5ad
