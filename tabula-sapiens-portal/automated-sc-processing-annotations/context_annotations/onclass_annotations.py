# Import OnClass and other libs
from OnClass.utils import *
from OnClass.OnClassModel import OnClassModel
from OnClass.other_datasets_utils import my_assemble, data_names_all, load_names

import subprocess
import scanpy as sc
import argparse


def main(
    input_file_path,
    output_file_path,
    data_file,# = '../data/OnClass_data/data_used_for_training/tabula-muris-senis-facs_cell_ontology.h5ad',
    cell_type_network_file,#='../data/OnClass_data/OnClass_data_others/cell_ontology/cl.ontology',
    use_pretrain_emb,#='../data/OnClass_data/pretrain/tp2emb_500',
    name_mapping_file,#='../data/OnClass_data/cell_ontology/cl.obo',
    use_pretrain_data,#='../data/OnClass_data/pretrain/BilinearNN_50019',
    use_pretrain_data_expression#='../data/OnClass_data/pretrain/BilinearNN_500'
):
    print('Embed the cell ontology')

    onclassmodel = OnClassModel()
    tp2emb, tp2i, i2tp = onclassmodel.EmbedCellTypes(dim=500,
        cell_type_network_file=cell_type_network_file,#'../data/OnClass_data/OnClass_data_others/cell_ontology/cl.ontology',
        use_pretrain=use_pretrain_emb)#'../data/OnClass_data/pretrain/tp2emb_500')

    print('Here, we used the pretrain cell type embedding file tp2emb_500')

    data_file = data_file#'../data/OnClass_data/data_used_for_training/tabula-muris-senis-facs_cell_ontology.h5ad' #same as the input
    train_X, train_genes, train_Y = read_data(feature_file=data_file, tp2i = tp2i, AnnData_label='cell_ontology_class_reannotated')

    print('Load the new dataset.')#Scanorama is used autoamtically to correct batch effcts.:
    test_data_file = input_file_path #'../data/output/output_processed.h5ad'
    test_X, test_genes, test_Y, test_AnnData = read_data(feature_file=test_data_file, tp2i = tp2i, return_AnnData=True, AnnData_label='cell_ontology_id')

    print('Predict the labels of cells')# in droplet cells. Scanorama is used autoamtically to correct batch effcts.:
    # currently only running on tensor flow 1; working on the fix for tensor flow 2
    # currently only running with correct_batch=False
    onclassmodel.train(train_X, train_Y, tp2emb, train_genes, nhidden=[500], log_transform = True, use_pretrain = use_pretrain_data, pretrain_expression=use_pretrain_data_expression)


    test_label = onclassmodel.predict(test_X, test_genes,log_transform=True,correct_batch=False)

    # Add the new annotations to the working file
    x = write_anndata_data(test_label, test_AnnData, i2tp, name_mapping_file=name_mapping_file)#'../data/OnClass_data/cell_ontology/cl.obo')#output_file is optional

    #print (x.obs['OnClass_annotation_ontology_ID'])
    #print (x.obs['OnClass_annotation_ontology_name'])

    x.write(output_file_path)#('../data/output/output_processed_annotated.h5ad')
    # Launch cellxgene with experimental_annotations


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Annotates data.
        ''')
    parser.add_argument(
        '--input', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--output', required=True,
        help='directory where arrow files should be written')
    parser.add_argument(
        '--data_file', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--cell_type_network_file', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--use_pretrain_emb', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--name_mapping_file', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--use_pretrain_data', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--use_pretrain_data_expression', required=True,
        help='directory containing h5ad files to read')


    args = parser.parse_args()



    main(input_file_path = args.input,
        output_file_path = args.output,
        data_file = args.data_file,
        cell_type_network_file = args.cell_type_network_file,
        use_pretrain_emb = args.use_pretrain_emb,
        name_mapping_file = args.name_mapping_file,
        use_pretrain_data = args.use_pretrain_data,
        use_pretrain_data_expression = args.use_pretrain_data_expression)
