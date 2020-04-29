import h5py
import os
import numpy as np
from anndata import read_h5ad

# Documentation on the structure of the h5 files
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info

rawdata_path = '/data/yosef2/users/chenling/TabulaSapiens/AnnotationsRound1/data/'
tenx = read_h5ad(rawdata_path + 'tabula-sapiens-10X-pilot-filtered.h5ad')
barcodes = tenx.obs.index
barcodes = [x.split('-')[0] for x in barcodes]
unique_barcodes, count = np.unique(barcodes, return_counts=True)
duplicated_barcode = unique_barcodes[count>1]

data_dir = '/data/scratch/users/chenling/TabulaSapiens/10x'
datasets = os.listdir(data_dir)
output = open(data_dir + '/molecule_info.csv', 'w')

for x in datasets:
    filename = '/'.join([data_dir, x, 'molecule_info.h5/molecule_info.h5'])
    f = h5py.File(filename, 'r')
    # for key in f.keys():
    #     print(f[key])

    assert np.unique(f['barcode_info']['pass_filter'][:, 1]) == 0
    assert np.unique(f['barcode_info']['pass_filter'][:, 2]) == 0

    molecules_cb = f['barcode_idx']
    molecules_gene = f['feature_idx']
    genenames = f['features']['name']
    molecules_umi = f['umi']
    all_barcodes = f['barcode_info']['pass_filter'][:, 0]
    barcode_seq = f['barcodes']

    assert len(molecules_umi) == len(molecules_cb)
    assert len(molecules_gene) == len(molecules_cb)

    # for i in range(len(molecules_cb)):
    for i in range(len(molecules_cb)):
        if molecules_cb[i] in all_barcodes:
            a = barcode_seq[molecules_cb[i]]
            if a.astype(str) in duplicated_barcode:
                b = genenames[molecules_gene[i]]
                c = molecules_umi[i]
                d = x
                output.write(','.join([a.astype(str), b.astype(str), c.astype(str), d])+'\n')


output.close()
