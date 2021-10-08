import os
import scanpy as sc
import numpy as np
import torch

from OnClass.utils import *
from OnClass.OnClassModel import OnClassModel
from OnClass.other_datasets_utils import my_assemble, data_names_all, load_names
from sklearn.svm import LinearSVC
from rpy2.robjects.packages import importr


import rpy2.robjects as ro

from rpy2.robjects import pandas2ri

pandas2ri.activate()

dplyr = importr("dplyr")
singleCellNet = importr("singleCellNet")


def setup_dataset(
    input_anndata_path,
    ref_anndata_path,
    tissue=None,
    use_10X_only=False,
    batch_correction_conditions=None,
):
    # TODO check that the input vars are the same
    train_data = sc.read_h5ad(ref_anndata_path)
    test_data = sc.read_h5ad(input_anndata_path)
    
    if tissue is not None:
        train_data = train_data[train_data.obs["tissue"] == tissue]
    if use_10X_only:
        train_data = train_data[train_data.obs["method"] == "10X"]
    if batch_correction_conditions is None:
        full = train_data.concatenate(test_data)
        full.obs["_batch_indices"] = np.zeros(full.n_obs).astype(int)
    else:
        condition = batch_correction_conditions[0]
        train_data.obs["_batch"] = train_data.obs[condition]

        if len(batch_correction_conditions) >= 2:
            for condition in batch_correction_conditions[1:]:
                x1 = train_data.obs["_batch"].to_numpy().astype(str)
                underscore = ["_"] * train_data.n_obs
                x1 = np.char.add(x1, underscore)
                x2 = train_data.obs[condition].to_numpy().astype(str)
                train_data.obs["_batch"] = np.char.add(x1, x2)
        batch_key = "_batch"

        all_batches = np.unique(train_data.obs["_batch"].astype(str))
        batch_indices = np.zeros(train_data.n_obs).astype(int)
        for i, batch in enumerate(all_batches):
            idx = np.where(train_data.obs[batch_key] == batch)
            batch_indices[idx] = i

        train_data.obs["_batch_indices"] = batch_indices
        test_data.obs["_batch_indices"] = np.zeros(test_data.n_obs).astype(int) + len(
            all_batches
        )
        full = train_data.concatenate(test_data)

    full.obs["_dataset"] = ["reference"] * train_data.n_obs + ["user"] * test_data.n_obs
    return full


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


def scanvi_annotation(
    full_dataset,
    batch_key,
    output_folder,
    ts_label_key="manual_cell_ontology_class",
    scvi_model=None,
    scanvi_model=None,
    n_scvi_epochs=400,
    n_scanvi_epochs=15,
    use_gpu=True,
):

    print("data preprocessing")

    train_data = full_dataset[full_dataset.obs["_dataset"] == "reference"]
    test_data = full_dataset[full_dataset.obs["_dataset"] == "user"]

    scvi_data = AnnDatasetFromAnnData(full_dataset, batch_label=batch_key)

    labels = full_dataset.obs[ts_label_key]
    mapping = labels.astype("category").cat.categories
    codes = labels.astype("category").cat.codes
    scvi_data.cell_types = labels.astype("category").unique()
    scvi_data.labels = np.unique(codes, return_inverse=True)[1]
    scvi_data.labels = scvi_data.labels.reshape(len(scvi_data.labels), 1)
    scvi_data.n_labels = len(scvi_data.cell_types)

    vae = VAE(
        scvi_data.nb_genes,
        n_batch=scvi_data.n_batches,
        n_layers=3,
        n_latent=50,
        dispersion="gene-batch",
        reconstruction_loss="zinb",
    )
    trainer = UnsupervisedTrainer(
        vae,
        scvi_data,
        train_size=0.99,
        use_cuda=use_gpu,
        frequency=5,
        data_loader_kwargs={"pin_memory": False},
        n_epochs_kl_warmup=10,
        trainer_kwargs=dict(data_loader_kwargs=dict(pin_memory=False)),
    )

    if scvi_model is None:
        print("Training scVI")
        trainer.train(n_epochs=n_scvi_epochs)
        scvi_model = os.path.join(output_folder, "trained_scvi_model.pt")
        torch.save(trainer.model.state_dict(), scvi_model)
    elif os.path.isfile(scvi_model):
        print("Loading pre-trained scVI model")
        trainer.model.load_state_dict(torch.load(scvi_model))
        trainer.model.eval()
    else:
        raise ValueError("invalid scvi_model arg")

    if "unassigned" in scvi_data.cell_types:
        unlabelled_idx = list(scvi_data.cell_types).index("unassigned")
        labelled = np.where(scvi_data.labels.ravel() != unlabelled_idx)[0]
    else:
        labelled = np.arange(len(scvi_data.labels))

    # balance number of labelled cells from each cell type
    labelled = balance_n_labelled(scvi_data.labels.ravel(), labelled, nlabels=50)
    labelled = np.random.choice(labelled, len(labelled), replace=False)

    unlabelled = [
        x for x in np.arange(len(scvi_data.labels.ravel())) if x not in labelled
    ]
    unlabelled = np.random.choice(unlabelled, len(unlabelled), replace=False)

    scanvi = SCANVI(
        scvi_data.nb_genes,
        scvi_data.n_batches,
        scvi_data.n_labels,
        n_layers=3,
        n_latent=50,
        dispersion="gene-batch",
    )  # symmetric_kl=True
    trainer_scanvi = SemiSupervisedTrainer(
        scanvi,
        scvi_data,
        n_epochs_classifier=100,
        lr_classification=5 * 1e-3,
        seed=1,
        n_epochs_kl_warmup=1,
        trainer_kwargs=dict(data_loader_kwargs=dict(pin_memory=False)),
    )

    trainer_scanvi.model.load_state_dict(torch.load(scvi_model), strict=False)
    trainer_scanvi.model.eval()
    trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=labelled)
    trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(indices=unlabelled)

    if scanvi_model is None:
        print("training scanvi")
        trainer_scanvi.train(n_epochs=n_scanvi_epochs)
        scanvi_model = os.path.join(output_folder, "trained_scanvi_model.pt")
        torch.save(trainer_scanvi.model.state_dict(), scanvi_model)
    elif os.path.isfile(scanvi_model):
        trainer_scanvi.model.load_state_dict(torch.load(scanvi_model))
    else:
        labelled = np.arange(len(scvi_data.labels))

    trainer_scanvi.model.eval()
    user_data = AnnDatasetFromAnnData(test_data, batch_label="batch")

    # should be scvi
    full = trainer_scanvi.create_posterior(
        trainer_scanvi.model, user_data, indices=np.arange(len(user_data))
    )
    _, pred = full.sequential().compute_predictions()
    import pdb

    pdb.set_trace()
    test_data.obs["scanvi_pred"] = [mapping[i] for i in pred]

    scanvi_annotation_path = os.path.join(output_folder, "scanvi_annotations.h5ad")
    test_data.write(scanvi_annotation_path)


def subsample_dataset(train_data, labels_key, n_samples):
    sample_idx = []
    labels, counts = np.unique(train_data.obs[labels_key], return_counts=True)
    for i, label in enumerate(labels):
        label_locs = np.where(train_data.obs[labels_key] == label)[0]
        if counts[i] < n_samples:
            sample_idx.append(label_locs)
        else:
            label_subset = np.random.choice(label_locs, n_samples, replace=False)
            sample_idx.append(label_subset)
    sample_idx = np.concatenate(sample_idx)
    return sample_idx


def svm_annotation(
    full_dataset, batch_key, output_folder, ts_label_key="manual_cell_ontology_class"
):

    train_data = full_dataset[full_dataset.obs["_dataset"] == "reference"]
    train_idx = subsample_dataset(train_data, ts_label_key, 100)
    train_data = train_data[train_idx].copy()
    test_data = full_dataset[full_dataset.obs["_dataset"] == "user"]

    sc.pp.log1p(train_data)
    sc.pp.log1p(test_data)

    train_X = train_data.X
    test_X = test_data.X
    train_Y = train_data.obs[ts_label_key]

    # train model
    Classifier = LinearSVC()
    Classifier.fit(train_X, train_Y)
    pred = Classifier.predict(test_X)

    test_data.obs["SVM_pred"] = pred

    svm_output_path = os.path.join(output_folder, "SVM_annotations.h5ad")
    test_data.write(svm_output_path)


def onclass_annotation(
    input_anndata,
    output_folder,
    ref_anndata_path,
    ref_label_key="manual_cell_ontology_class",
):

    cell_type_network_file = (
        "data/OnClass_data/OnClass_data_others/cell_ontology/cl.ontology"
    )
    use_pretrain_emb = "data/OnClass_data/pretrain/tp2emb_500"
    name_mapping_file = "data/OnClass_data/OnClass_data_others/cell_ontology/cl.obo"
    name_mapping_file = "data/chenling_cl.obo"
    use_pretrain_data = "data/OnClass_data/pretrain/BilinearNN_50019"
    use_pretrain_data_expression = "data/OnClass_data/pretrain/BilinearNN_500"

    print("Embed the cell ontology")

    onclassmodel = OnClassModel()
    tp2emb, tp2i, i2tp = onclassmodel.EmbedCellTypes(
        dim=500,
        cell_type_network_file=cell_type_network_file,
        use_pretrain=use_pretrain_emb,
    )

    print("Here, we used the pretrain cell type embedding file tp2emb_500")

    train_X, train_genes, train_Y = read_data(
        feature_file=ref_anndata_path, tp2i=tp2i, AnnData_label=ref_label_key
    )

    print("Loading the new dataset.")
    test_X, test_genes, test_AnnData = read_data(
        feature_file=input_anndata, tp2i=tp2i, return_AnnData=True
    )

    print(
        "Predicting the labels of cells"
    )  # in droplet cells. Scanorama is used autoamtically to correct batch effcts.:
    # currently only running on tensor flow 1; working on the fix for tensor flow 2
    # currently only running with correct_batch=False
    onclassmodel.train(
        train_X,
        train_Y,
        tp2emb,
        train_genes,
        nhidden=[500],
        log_transform=True,
        use_pretrain=use_pretrain_data,
        pretrain_expression=use_pretrain_data_expression,
    )

    test_label = onclassmodel.predict(
        test_X, test_genes, log_transform=True, correct_batch=False
    )

    # Add the new annotations to the working file
    x = write_anndata_data(
        test_label, test_AnnData, i2tp, name_mapping_file=name_mapping_file
    )  #'../data/OnClass_data/cell_ontology/cl.obo')#output_file is optional

    onclass_output_fn = os.path.join(output_folder, "onclass_annotations.h5ad")
    x.write(onclass_output_fn)


def singlecellnet_annotation(
    full_dataset,
    annotation_key,
    batch_key=None,
):
    print("running singlecellnet")
    sc.pp.log1p(full_dataset)

    train_data = full_dataset[full_dataset.obs["_dataset"] == "reference"]
    test_data = full_dataset[full_dataset.obs["_dataset"] == "user"]
    train_X = pd.DataFrame(data=train_data.X.T.todense())
    nr, nc = train_X.shape

    # fails here
    R_train_X = ro.conversion.py2rpy(train_X)

    train_Y = pd.DataFrame(data=train_data.obs[annotation_key])
    R_train_Y = ro.conversion.py2rpy(train_Y)
    cgenes2 = singleCellNet.findClassyGenes(R_train_X, R_train_Y, annotation_key)

    cgenes = cgenes2.rx2("cgenes")
    cgenes = cgenes.astype(int)
    ctrainX = train_X.iloc[cgenes]
    crtrainX = ro.conversion.py2rpy(ctrainX)

    xpairs = singleCellNet.ptGetTop(
        crtrainX, cgenes2.rx2("grps"), cgenes2.rx2("cgenes_list")
    )

    pdTrain = singleCellNet.query_transform(crtrainX, xpairs)

    tmp_pdTrain = pd.DataFrame(data=pdTrain, index=xpairs)
    rpdTrain = ro.conversion.py2rpy(tmp_pdTrain)
    print("making classifier")
    rf = singleCellNet.sc_makeClassifier(
        rpdTrain, genes=xpairs, groups=cgenes2.rx2("grps")
    )

    test_X = data[test_idx, cgenes].X
    R_test_X = ro.r.matrix(test_X)
    DataTest = singleCellNet.query_transform(R_test_X, xpairs)
    classRes = singleCellNet.rf_classPredict(rf, DataTest)
    pred = com.load_data(classRes)

    test_data.obsm["singlecellnet"] = pred
    test_data.write(output_file_path + "singlecellnet.h5ad")