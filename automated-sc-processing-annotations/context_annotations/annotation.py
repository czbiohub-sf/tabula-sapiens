# Import general packages
import os
import subprocess
import scanpy as sc
import argparse
import numpy as np

# Import OnClass related libs
from OnClass.utils import *
from OnClass.OnClassModel import OnClassModel
from OnClass.other_datasets_utils import my_assemble, data_names_all, load_names

# Import svm related libs
from sklearn.svm import LinearSVC

# Import scVI related libs
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from scvi.dataset import AnnDatasetFromAnnData
from scvi.models import VAE, SCANVI
import torch

def balance_n_labelled(l, labelled, nlabels=30):
    balanced_labelled = []
    for i in np.unique(l[labelled]):
        idx = np.where(l == i)[0]
        if len(idx) >= nlabels:
            subset = np.random.choice(idx, nlabels, replace=False)
        else:
            subset = idx
        balanced_labelled.append(subset)
    balanced_labelled = np.concatenate(balanced_labelled)
    return balanced_labelled

# Import singlecellnet related libs
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
singleCellNet = importr('singleCellNet')
dplyr = importr('dplyr')


def main(
    input_file_path,
    output_file_path,
    data_file,# = '../data/OnClass_data/data_used_for_training/tabula-muris-senis-facs_cell_ontology.h5ad',
    cell_type_network_file,#='../data/OnClass_data/OnClass_data_others/cell_ontology/cl.ontology',
    use_pretrain_emb,#='../data/OnClass_data/pretrain/tp2emb_500',
    name_mapping_file,#='../data/OnClass_data/cell_ontology/cl.obo',
    use_pretrain_data,#='../data/OnClass_data/pretrain/BilinearNN_50019',
    use_pretrain_data_expression,#='../data/OnClass_data/pretrain/BilinearNN_500'
    annotation_method = 'onclass',
    AnnData_label = 'cell_ontology_class_reannotated',
    AnnData_batch='batch',
    use_cuda = False,
    scvi_model=None,
    scanvi_model=None,
):
    if annotation_method == 'onclass':
        print('Embed the cell ontology')

        onclassmodel = OnClassModel()
        tp2emb, tp2i, i2tp = onclassmodel.EmbedCellTypes(dim=500,
            cell_type_network_file=cell_type_network_file,#'../data/OnClass_data/OnClass_data_others/cell_ontology/cl.ontology',
            use_pretrain=use_pretrain_emb)#'../data/OnClass_data/pretrain/tp2emb_500')

        print('Here, we used the pretrain cell type embedding file tp2emb_500')

        data_file = data_file#'../data/OnClass_data/data_used_for_training/tabula-muris-senis-facs_cell_ontology.h5ad' #same as the input
        train_X, train_genes, train_Y = read_data(feature_file=data_file, tp2i = tp2i, AnnData_label=AnnData_label)

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
    elif annotation_method=='scanvi':
        print('data preprocessing')
        train_data = sc.read_h5ad(data_file)
        test_data = sc.read_h5ad(input_file_path)
        data = train_data.concatenate(test_data)

        batch = list(np.unique(data.obs[AnnData_batch]))
        batch_id = [batch.index(x) for x in data.obs[AnnData_batch]]
        data.obs['batch'] = batch_id

        scvi_data = AnnDatasetFromAnnData(data, batch_label = 'batch')
        labels = data.obs[AnnData_label]
        scvi_data.cell_types, scvi_data.labels = np.unique(labels, return_inverse=True)
        scvi_data.labels = scvi_data.labels.reshape(len(scvi_data.labels), 1)
        scvi_data.n_labels = len(scvi_data.cell_types)

        vae = VAE(scvi_data.nb_genes, n_batch=scvi_data.n_batches,
                  n_layers=3, n_latent=50, dispersion='gene-batch',
                  reconstruction_loss='zinb')
        trainer = UnsupervisedTrainer(
            vae,
            scvi_data,
            train_size=0.99,
            use_cuda=use_cuda,
            frequency=5,
            data_loader_kwargs={"pin_memory": False},
            n_epochs_kl_warmup=10
        )
        if os.path.isfile(scvi_model):
            print('loading pre-trained scVI model')
            trainer.model.load_state_dict(torch.load(scvi_model))
            trainer.model.eval()
        else:
            print('train scVI')
            trainer.train(n_epochs=100)
            torch.save(trainer.model.state_dict(), scvi_model)

        if 'unassigned' in scvi_data.cell_types:
            unlabelled_idx = list(scvi_data.cell_types).index('unassigned')
            labelled = np.where(scvi_data.labels.ravel() != unlabelled_idx)[0]
        else:
            labelled = np.arange(len(scvi_data.labels))

            # balance number of labelled cells from each cell type
        labelled = balance_n_labelled(scvi_data.labels.ravel(), labelled, nlabels=50)
        labelled = np.random.choice(labelled, len(labelled), replace=False)

        unlabelled = [x for x in np.arange(len(scvi_data.labels.ravel())) if x not in labelled]
        unlabelled = np.random.choice(unlabelled, len(unlabelled), replace=False)

        scanvi = SCANVI(scvi_data.nb_genes, scvi_data.n_batches, scvi_data.n_labels, n_layers=3, n_latent=50,
                        symmetric_kl=True, dispersion='gene-batch')


        trainer_scanvi = SemiSupervisedTrainer(scanvi, scvi_data,
                                               n_epochs_classifier=100,
                                               lr_classification=5 * 1e-3, seed=1,
                                               n_epochs_kl_warmup=1)

        trainer_scanvi.model.load_state_dict(torch.load(scvi_model), strict=False)
        trainer_scanvi.model.eval()
        trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=labelled)
        trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=unlabelled)

        if os.path.isfile(scanvi_model):
            trainer_scanvi.model.load_state_dict(torch.load(scanvi_model))
            trainer_scanvi.model.eval()
        else:
            trainer_scanvi.train(n_epochs=15)
            torch.save(trainer_scanvi.model.state_dict(), scanvi_model)

        full = trainer_scanvi.create_posterior(trainer_scanvi.model, data, indices=np.arange(len(data)))
        _, pred = full.sequential().compute_predictions()
        data.obs['scanvi_pred'] = pred
        data = data[train_data.obs.index]
        data.write(output_file_path + "scanvi.h5ad")

    elif annotation_method=='svm':
        train_data = sc.read_h5ad(data_file)
        test_data = sc.read_h5ad(input_file_path)
        data = train_data.concatenate(test_data)
        # AnnData_batch must be a column that exist in both the data file and test data file
        batch = list(np.unique(data.obs[AnnData_batch]))
        batch_id = [batch.index(x) for x in data.obs[AnnData_batch]]
        data.obs['batch'] = batch_id

        data = np.log1p(data.X)
        # feature selection
        sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='batch')
        var_select = data.var.highly_variable_nbatches > 1
        var_genes = var_select.index[var_select]
        data = data[:, var_genes]

        train_idx = train_data.obs.index
        test_idx = test_data.obs.index
        train_X = data[train_idx].X
        test_X = data[test_idx].X
        train_Y = data[train_idx].obs[AnnData_label]

        # train model
        Classifier = LinearSVC()
        Classifier.fit(train_X, train_Y)
        pred = Classifier.predict(test_X)
        train_data.obs['SVM_pred'] = pred
        train_data.write(output_file_path + "SVM.h5ad")

    elif annotation_method=='singlecellnet':
        print('running singlecellnet')

        train_data = sc.read_h5ad(data_file)
        test_data = sc.read_h5ad(input_file_path)
        data = train_data.concatenate(test_data)
        # AnnData_batch must be a column that exist in both the data file and test data file
        batch = list(np.unique(data.obs[AnnData_batch]))
        batch_id = [batch.index(x) for x in data.obs[AnnData_batch]]
        data.obs['batch'] = batch_id

        data = np.log1p(data.X)
        # # feature selection
        # sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='batch')
        # var_select = data.var.highly_variable_nbatches > 1
        # var_genes = var_select.index[var_select]
        # data = data[:, var_genes]

        train_idx = train_data.obs.index
        test_idx = test_data.obs.index
        train_X = data[train_idx].X
        train_Y = data[train_idx].obs

        R_train_X = ro.r.matrix(train_X)
        R_train_Y = com.convert_to_r_dataframe(train_Y)

        cgenes2 = singleCellNet.findClassyGenes(R_train_X, R_train_Y, AnnData_label)
        # ro.globalenv['cgenesA'] = cgenes2.rx2('cgenes')
        cgenes = cgenes2.rx2('cgenes')

        train_X = data[train_idx,cgenes].X
        R_train_X = ro.r.matrix(train_X)

        xpairs = singleCellNet.ptGetTop(R_train_X, cgenes2.rx2('grps'), ncores=1)
        pdTrain = singleCellNet.query_transform(R_train_X, xpairs)
        rf = singleCellNet.sc_makeClassifier(pdTrain.rx(xpairs,True), genes=xpairs, groups= cgenes2.rx2('grps'))

        test_X = data[test_idx,cgenes].X
        R_test_X = ro.r.matrix(test_X)
        DataTest = singleCellNet.query_transform(R_test_X, xpairs)
        classRes = singleCellNet.rf_classPredict(rf, DataTest)
        pred = com.load_data(classRes)

        train_data.obsm['singlecellnet'] = pred
        train_data.write(output_file_path + "singlecellnet.h5ad")

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
