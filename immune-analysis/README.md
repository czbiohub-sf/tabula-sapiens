# TabulaSapiens

this repository contains the code used to analyze the immune compartment of Tabula Sapiens. 

I used snakemake to preprocess the files and integrate datasets from smart-seq2 and 10X genomics. 

The resulting summary files are analyzed in the jupyter notebooks. I mostly make use of scanpy and scirpy, two valuable tools for analyzing single-cell RNA sequencing data.

Note that in most notebooks I specify relative paths that would need to be edited by the end user to point to the relevant data. 

explanation of directories:

preprocessing/ contains the snakemake workflow I used to integrate the data

metadata/ contains files used in analyses such as dissociation genes

data/ contains intermediate and final tables used for analysis

noteboks/ contains the ipython notebooks I used analyze and generate figures from the data
