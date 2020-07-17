# Analysis Notebooks
These notebooks should be run successively to reproduce the results used in the paper

## Generating unannotated dataset

## Label Propagation
Each organ expert group had used cellxgene to annotate celltypes from a cellxgene object. 
We read in the csv output by cellxgene (this required a little bit of manual engineering because the column names and numbers are not the same) 
We convert all the non-standard terms to cell ontology standard cell types, except for terms that we deem necessary to be added to the ontology. 
There are unlabelled cells after manual annotation (Smartseq2 cells, or cells with missing expression that cannot be annotated confidently) and all those cells are labelled using the scanvi algorithm.

## Split Compartment 
This notebook uses compartment markers to predict whether each cell is Immune, Epithelial, Endothelial or Stromal. It then combines the compartment and the celltype predictions and generates the final combined object used for most of the downstream analysis. It also generates a seperate object for each compartment. The latent space is computed on all of the data, treating each technology x donor as a batch. Only the 2D UMAP is recomputed for the per-compartment objects. 

## Generate Organ Objects
Recomputing UMAP and generating objects for each organ. 

## Hierarchical DE
Running differential gene expression on the ontology tree 
