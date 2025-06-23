import os
from typing import List

import anndata
import numpy as np
import OnClass as oc
import scanorama
import scanpy as sc
from OnClass.OnClassModel import OnClassModel
from OnClass.utils import read_data, write_anndata_data
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

import scvi

def annotate_data(adata, methods, save_path, ref_path='ref.h5ad', query_path='query.h5ad'):
    if "bbknn" in methods:
        run_bbknn(adata, batch_key="_batch")
        run_knn_on_bbknn(adata, labels_key="_labels", result_key="knn_on_bbknn_pred")
        save_results(adata, save_path, 
                     obs_keys=["knn_on_bbknn_pred"], 
                     obsm_keys=['bbknn_umap'])

    if "scvi" in methods:
        training_mode = adata.uns["_training_mode"]
        scvi_obsm_latent_key = 'X_scvi_' + training_mode
        run_scvi(
            adata,
            max_epochs=None,
            n_latent=50,
            dropout_rate=0.1,
            dispersion = "gene-batch",
            obsm_latent_key = scvi_obsm_latent_key,
            pretrained_scvi_path="lung_scvi_reference",
        )
        knn_pred_key = "knn_on_scvi_{}_pred".format(training_mode)
        run_knn_on_scvi(
            adata, obsm_key=scvi_obsm_latent_key,
            result_key= knn_pred_key
        )
        save_results(
            adata,
            save_path,
            obs_keys=[knn_pred_key],
            obsm_keys=[scvi_obsm_latent_key, scvi_obsm_latent_key+'_umap'])
    
    if "scanvi" in methods:
        training_mode = adata.uns['_training_mode']
        obsm_latent_key = 'X_scanvi_{}'.format(training_mode)
        predictions_key = 'scanvi_{}_pred'.format(training_mode)
        run_scanvi(adata, 
                   max_epochs=None,
                   n_latent=100,
                   dropout_rate=0.1,
                   obsm_latent_key=obsm_latent_key,
                   obs_pred_key=predictions_key,
                   pretrained_scanvi_path='lung_scanvi_ts')
        
        save_results(
            adata,
            save_path,
            obs_keys=[predictions_key],
            obsm_keys=[obsm_latent_key],
        )
    if "svm" in methods:
        run_svm_on_hvg(adata)
        save_results(adata, 
             save_path, 
             obs_keys=['svm_pred'])
    
    if 'rf' in methods:
        run_rf_on_hvg(adata)
        save_results(adata, 
             save_path, 
             obs_keys=['rf_pred'])
    
    if 'onclass' in methods:
        # TODO: change the ref adata and query adata paths
        run_onclass(adata=adata, 
            max_iter=2,
            ref_adata_path=ref_path,
            query_adata_path=query_path, 
            cl_obo_file='cl.obo', 
            cl_ontology_file='cl.ontology')
        save_results(adata, 
             save_path, 
             obs_keys=['onclass_pred'])

    if 'scanorama' in methods:
        run_scanorama(adata)
        run_knn_on_scanorama(adata)
        save_results(adata, 
             save_path, 
             obs_keys=['knn_on_scanorama_pred'],
             obsm_keys=['scanorama_umap'])

def save_results(
    adata, results_adata_path, obs_keys=[], obsm_keys=[], compression="gzip"
):
    """
    If results_adata_path exists, will load and save results into it
    Else, will save adata to results_adata_path

    Parameters
    ----------
    adata
          adata with results in it
    results_adata_path
          path to save results. If it already exists, will load and save data to it
    obs_keys
          obs keys to save
    obsm_keys
          obsm keys to save
    compression
          If enabled, will save with compression. Smaller file sizes, but longer save times
    """
    if os.path.exists(results_adata_path):
        results = anndata.read(results_adata_path)
        for key in obs_keys:
            results.obs[key] = adata.obs[key]
        for key in obsm_keys:
            results.obsm[key] = adata.obsm[key]
        results.write(results_adata_path, compression)
    else:
        adata.write(results_adata_path, compression)


def check_genes_is_subset(ref, query):
    ref_genes = ref.var_names.to_numpy().astype("str")
    query_genes = query.var_names.to_numpy().astype("str")

    if len(set(query_genes)) != len(query_genes):
        print("Warning: Your genes are not unique.")

    if set(ref_genes).issubset(set(query_genes)):
        print("All ref genes are in query dataset. Can use pretrained models.")
        is_subset = True
    else:
        print("Not all reference genes are in query dataset. Retraining models.")
        is_subset = False

    return is_subset


def process_query(
    ref_path,
    query_path,
    ref_labels_key,
    ref_batch_key,
    query_labels_key=None,
    query_batch_key=None,
    unknown_celltype_label="unknown",
    join="inner",
    train_from_scratch=True,
    inplace=True,
    training_mode="online",  # online or offline
):
    ref_adata = anndata.read(ref_path)
    query_adata = anndata.read(query_path)
    query_adata.obs[query_labels_key] = "unknown"

    if inplace is not True:
        query_adata = query_adata.copy()

    is_subset = check_genes_is_subset(ref_adata, query_adata)

    if training_mode == "online" and is_subset is False:
        print("Is not subset, training offline.")
        training_mode = "offline"

    # TODO: add this to the datasets themselves
    ref_adata.obs["_labels"] = ref_adata.obs[ref_labels_key]
    ref_adata.obs["_batch"] = ref_adata.obs[ref_batch_key]
    ref_adata.obs["_dataset"] = "ref"

    ref_adata.obs["_ref_subsample"] = False
    ref_subsample_idx = subsample_dataset(ref_adata, ref_labels_key)
    ref_adata.obs["_ref_subsample"][ref_subsample_idx] = True

    query_adata.obs["_batch"] = (
        query_adata.obs[query_batch_key].astype("str") + "_query"
    )
    query_adata.obs["_dataset"] = "query"
    query_adata.obs["_ref_subsample"] = False

    # TODO: delete
    query_adata.obs["_labels"] = "unknown"
    query_labels = query_adata.obs[query_labels_key]
    known_cell_idx = np.where(query_labels != unknown_celltype_label)[0]
    if len(known_cell_idx) != 0:
        query_adata.obs["_labels"][known_cell_idx] = query_labels[known_cell_idx]

    if training_mode == "online":
        query_adata = query_adata[:, ref_adata.var_names].copy()
        adata = anndata.concat((ref_adata, query_adata))
    elif training_mode == "offline":
        adata = anndata.concat((ref_adata, query_adata), join=join)

    adata.uns["_training_mode"] = training_mode

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10, zero_center=False)

    sc.pp.highly_variable_genes(
        adata, n_top_genes=4000, subset=True, layer="counts", flavor="seurat_v3"
    )
    sc.tl.pca(adata, svd_solver="arpack")
    scvi.data.setup_anndata(
        adata, batch_key="_batch", labels_key="_labels", layer="counts"
    )
    return adata


def run_bbknn(adata, batch_key="_batch"):
    print("Running bbknn")
    sc.external.pp.bbknn(
        adata,
        batch_key=batch_key,
        approx=True,
        metric="angular",
        n_pcs=20,
        trim=None,
        n_trees=10,
        use_faiss=True,
        set_op_mix_ratio=1.0,
        local_connectivity=1,
    )
    sc.tl.umap(adata)
    adata.obsm['bbknn_umap'] = adata.obsm['X_umap']
    del adata.obsm['X_umap']
    
    return adata


def run_knn_on_bbknn(adata, labels_key="_labels", result_key="knn_on_bbknn_pred"):
    print("Classifying with knn on bbknn distances")
    distances = adata.obsp["distances"]

    ref_idx = adata.obs["_dataset"] == "ref"
    query_idx = adata.obs["_dataset"] == "query"

    ref_dist_idx = np.where(ref_idx == True)[0]
    query_dist_idx = np.where(query_idx == True)[0]

    train_Y = adata[ref_idx].obs[labels_key].to_numpy()
    train_distances = distances[ref_dist_idx, :][:, ref_dist_idx]

    knn = KNeighborsClassifier(n_neighbors=2, metric="precomputed")
    knn.fit(train_distances, y=train_Y)

    test_distances = distances[query_dist_idx, :][:, ref_dist_idx]
    knn_pred = knn.predict(test_distances)

    # save_results
    adata.obs[result_key] = "na"
    adata.obs[result_key][query_idx] = knn_pred
    print('Saved knn on bbknn results to adata.obs["{}"]'.format(result_key))


def subsample_dataset(
    train_data, labels_key, n_samples_per_label=100, n_total_samples=None
):
    # TODO right now all tissues will subsample to 100.
    # should go through each dataset and set manually
    sample_idx = []
    labels, counts = np.unique(train_data.obs[labels_key], return_counts=True)

    if n_total_samples is not None:
        if n_total_samples > train_data.n_obs:
            print(
                "n_total_samples exceeds size of dataset. "
                + "Setting to input adata size of {} cells".format(train_data.n_obs)
            )
            return train_data.obs_names

        n_samples_per_label = int(n_total_samples / len(labels))
    print("Sampling {} per label".format(n_samples_per_label))

    for i, label in enumerate(labels):
        label_locs = np.where(train_data.obs[labels_key] == label)[0]
        if counts[i] < n_samples_per_label:
            sample_idx.append(label_locs)
        else:
            label_subset = np.random.choice(
                label_locs, n_samples_per_label, replace=False
            )
            sample_idx.append(label_subset)
    sample_idx = np.concatenate(sample_idx)

    if n_total_samples is not None:
        if len(sample_idx) < n_total_samples:
            all_idx = np.arange(train_data.n_obs)
            remaining_idx = list(set(all_idx) - set(sample_idx))
            remaining_samples = n_total_samples - len(sample_idx)
            label_subset = np.random.choice(
                remaining_idx, remaining_samples, replace=False
            )
            sample_idx = np.concatenate((sample_idx, label_subset))
    return train_data.obs_names[sample_idx]


def run_rf_on_hvg(adata, save_key="rf_pred"):
    train_idx = adata.obs["_ref_subsample"]
    test_idx = adata.obs["_dataset"] == "query"

    train_X = adata[train_idx].layers["counts"]
    train_Y = adata[train_idx].obs["_labels"].to_numpy()
    test_X = adata[test_idx].layers["counts"]

    print("Training random forest classifier with {} cells".format(len(train_Y)))
    rf = RandomForestClassifier()
    rf.fit(train_X, train_Y)
    rf_pred = rf.predict(test_X)

    adata.obs[save_key] = "na"
    adata.obs[save_key][test_idx] = rf_pred


def run_onclass(
    adata,
    ref_adata_path,
    query_adata_path,
    cl_obo_file,
    cl_ontology_file,
    cell_ontology_obs_key="cell_ontology_id",
    save_key="onclass_pred",
    n_hidden=500,
    max_iter=20,
    save_model="onclass_model",
    shard_size=10000,
):
    oc = OnClassModel()
    tp2emb, tp2i, i2tp = oc.EmbedCellTypes(
        dim=500, cell_type_network_file=cl_ontology_file, use_pretrain=None
    )
    train_X, train_genes, train_Y = read_data(
        feature_file=ref_adata_path, tp2i=tp2i, AnnData_label=cell_ontology_obs_key
    )
    oc.train(
        train_X,
        train_Y,
        tp2emb,
        train_genes,
        nhidden=[n_hidden],
        max_iter=max_iter,
        use_pretrain=None,
        save_model=save_model,
    )
    test_X, test_genes, test_adata = read_data(
        feature_file=query_adata_path, tp2i=tp2i, return_AnnData=True
    )

    test_adata.obs["onclass_pred"] = "na"
    print(test_adata.n_obs)
    if test_adata.n_obs > shard_size:
        for i in range(0, test_adata.n_obs, shard_size):
            tmp_Anndata = test_adata[i : i + shard_size]
            tmp_X = test_X[i : i + shard_size]
            test_label = oc.predict(tmp_X, test_genes, log_transform=True)
            onclass_pred = write_anndata_data(
                test_label, tmp_Anndata, i2tp, name_mapping_file=cl_obo_file
            )
            onclass_pred = onclass_pred.obs["OnClass_annotation_ontology_name"]
            test_adata.obs[save_key][i : i + shard_size] = onclass_pred
    else:
        test_label = oc.predict(test_X, test_genes, log_transform=True)
        onclass_pred = write_anndata_data(
            test_label, test_adata, i2tp, name_mapping_file=cl_obo_file
        )
        onclass_pred = onclass_pred.obs["OnClass_annotation_ontology_name"]
        test_adata.obs[save_key] = onclass_pred

    adata.obs["onclass_pred"] = "na"
    adata.obs["onclass_pred"][test_adata.obs_names] = test_adata.obs["onclass_pred"]
    return adata


def run_svm_on_hvg(adata):
    train_idx = adata.obs["_ref_subsample"]
    test_idx = adata.obs["_dataset"] == "query"

    train_X = adata[train_idx].layers["counts"]
    train_Y = adata[train_idx].obs["_labels"].to_numpy()
    test_X = adata[test_idx].layers["counts"]

    clf = svm.LinearSVC(max_iter=1000)
    clf.fit(train_X, train_Y)
    svm_pred = clf.predict(test_X)

    # save_results
    adata.obs["svm_pred"] = "na"
    adata.obs["svm_pred"][test_idx] = svm_pred


def run_scvi(
    adata,
    n_latent=50,
    n_layers=3,
    dropout_rate=0.1,
    dispersion="gene",
    max_epochs=None,
    batch_size=1024,
    pretrained_scvi_path=None,
    var_subset_type="inner_join",
    obsm_latent_key="X_scvi",
    save_folder=None,
    overwrite=True,
    save_anndata=False,
):
    training_mode = adata.uns["_training_mode"]
    if training_mode == "online" and pretrained_scvi_path is None:
        raise ValueError("online training but no pretrained_scvi_path passed in.")

    if training_mode == "offline":
        model = scvi.model.SCVI(
            adata,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            dispersion=dispersion,
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
        )
        print("Training scvi offline.")

    elif training_mode == "online":
        if max_epochs is None:
            n_cells = adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 200), 200])

        query = adata[adata.obs["_dataset"] == "query"].copy()
        model = scvi.model.SCVI.load_query_data(query, pretrained_scvi_path)
        print("Training scvi online.")

    model.train(max_epochs=max_epochs, train_size=1.0, batch_size=batch_size)
    adata.obsm[obsm_latent_key] = model.get_latent_representation(adata)
    
    sc.pp.neighbors(adata, use_rep=obsm_latent_key)
    sc.tl.umap(adata)
    adata.obsm[obsm_latent_key+'_umap'] = adata.obsm['X_umap']
    del adata.obsm['X_umap']
    
    if save_folder is not None:
        model.save(save_folder, overwrite=overwrite, save_anndata=save_anndata)


def run_knn_on_scvi(adata, obsm_key="X_scvi", result_key="knn_on_scvi_pred"):
    if obsm_key not in adata.obsm.keys():
        raise ValueError("Please train scVI first or pass in a valid obsm_key.")

    print(
        "Training knn on scvi latent space. "
        + 'Using latent space in adata.obsm["{}"]'.format(obsm_key)
    )

    ref_idx = adata.obs["_dataset"] == "ref"
    query_idx = adata.obs["_dataset"] == "query"

    train_X = adata[ref_idx].obsm[obsm_key]
    train_Y = adata[ref_idx].obs["_labels"].to_numpy()

    test_X = adata[query_idx].obsm[obsm_key]

    knn = KNeighborsClassifier(n_neighbors=15, weights="uniform")
    knn.fit(train_X, train_Y)
    knn_pred = knn.predict(test_X)

    # save_results
    adata.obs[result_key] = "na"
    adata.obs[result_key][query_idx] = knn_pred


def run_knn_on_scanorama(
    adata,
    obsm_key="X_scanorama",
    result_key="knn_on_scanorama_pred",
):
    print('Running knn on scanorama')
    if obsm_key not in adata.obsm.keys():
        print("Please run scanorama first or pass in a valid obsm_key.")
    
    ref_idx = adata.obs["_dataset"] == "ref"
    query_idx = adata.obs["_dataset"] == "query"

    train_X = adata[ref_idx].obsm[obsm_key]
    train_Y = adata[ref_idx].obs["_labels"].to_numpy()

    test_X = adata[query_idx].obsm[obsm_key]
    
    knn = KNeighborsClassifier(n_neighbors=15, weights="uniform")
    knn.fit(train_X, train_Y)
    knn_pred = knn.predict(test_X)

    # save_results
    adata.obs[result_key] = "na"
    adata.obs[result_key][query_idx] = knn_pred


def run_scanorama(adata, batch_key="_batch"):
    # TODO add check if in colab and n_genes > 120000
    # throw warning
    adatas = [adata[adata.obs[batch_key] == i] for i in np.unique(adata.obs[batch_key])]
    scanorama.integrate_scanpy(adatas, dimred=50)
    tmp_adata = anndata.concat(adatas)
    adata.obsm["X_scanorama"] = tmp_adata[adata.obs_names].obsm["X_scanorama"]
    
    print('Computing umap on scanorama')
    sc.pp.neighbors(adata, use_rep='X_scanorama')
    sc.tl.umap(adata)
    adata.obsm['scanorama_umap'] = adata.obsm['X_umap']
    del adata.obsm['X_umap']
    
def run_scanvi(
    adata,
    n_layers=3,
    dropout_rate=0.2,
    n_classifier_layers=1,
    classifier_dropout=0.4,
    max_epochs=None,
    n_latent=100,
    batch_size=1024,
    n_epochs_kl_warmup=20,
    n_samples_per_label=100,
    obsm_latent_key="X_scanvi",
    obs_pred_key="scanvi_pred",
    pretrained_scanvi_path=None,
    save_folder=None,
    save_anndata=False,
    overwrite=True,
):
    training_mode = adata.uns["_training_mode"]
    if training_mode == "online" and pretrained_scanvi_path is None:
        raise ValueError("online training but no pretrained_scvi_path passed in.")

    if training_mode == "offline":
        model_kwargs = dict(
            use_layer_norm="both",
            use_batch_norm="none",
            classifier_parameters={
                "n_layers": n_classifier_layers,
                "dropout_rate": classifier_dropout,
            },
        )
        model = scvi.model.SCANVI(
            adata,
            unlabeled_category="unknown",
            n_layers=n_layers,
            encode_covariates=True,
            dropout_rate=dropout_rate,
            n_latent=n_latent,
            **model_kwargs,
        )

    elif training_mode == "online":
        if max_epochs is None:
            n_cells = adata.n_obs
            max_epochs = np.min([round((20000 / n_cells) * 200), 200])

        query = adata[adata.obs["_dataset"] == "query"].copy()
        model = scvi.model.SCANVI.load_query_data(
            query, pretrained_scanvi_path, freeze_classifier=True
        )

    plan_kwargs = dict(n_epochs_kl_warmup=n_epochs_kl_warmup)
    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=1.0,
        n_samples_per_label=n_samples_per_label,
        plan_kwargs=plan_kwargs,
    )

    adata.obsm[obsm_latent_key] = model.get_latent_representation(adata)
    adata.obs[obs_pred_key] = model.predict(adata)

    if save_folder is not None:
        model.save(save_folder, overwrite=overwrite, save_anndata=save_anndata)
