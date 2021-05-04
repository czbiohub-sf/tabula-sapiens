import os
import obonet

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanorama
import scanpy as sc
import scvi
import seaborn as sns

from OnClass.OnClassModel import OnClassModel

from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.neighbors import KNeighborsClassifier
import matplotlib.backends.backend_pdf


def make_cell_ontology_id(adata, labels_key, celltype_dict, ontology_key=None):
    """
    Convert celltype names to ontology id.

    Parameters
    ----------
    adata
        AnnData object
    labels_key
        Key in adata.obs to convert to ontology id
    celltype_dict
        Dictionary mapping celltype to ontology id
    ontology_key
        Key in adata.obs to save ontology ids to.
        Default will be <labels_key>_cell_ontology_id
    """
    if ontology_key is None:
        ontology_key = labels_key + "_cell_ontology_id"
    ontology_id = []

    for label in adata.obs[labels_key]:
        if label != "unknown":
            if label not in celltype_dict:
                print("Following label not in celltype_dict ", label)
            ontology_id.append(celltype_dict[label])
        else:
            ontology_id.append("unknown")

    adata.obs[ontology_key] = ontology_id
    return ontology_key


def make_celltype_to_cell_ontology_id_dict(obo_file):
    """
    Make celltype to ontology id dict and vice versa.

    Parameters
    ----------
    obo_file
        obofile to read

    Returns
    -------
    name2id
        dictionary of celltype names to ontology id
    id2name
        dictionary of ontology id to celltype names
    """
    with open(obo_file, "r") as f:
        co = obonet.read_obo(f)
        id2name = {id_: data.get("name") for id_, data in co.nodes(data=True)}
        id2name = {k: v.lower() for k, v in id2name.items() if v is not None}
        name2id = {v: k for k, v in id2name.items()}

    return name2id, id2name


def prediction_eval(
    pred,
    labels,
    name,
    x_label="",
    y_label="",
    res_dir="./",
):
    """
    Generate confusion matrix
    """
    x = np.concatenate([labels, pred])
    types, temp = np.unique(x, return_inverse=True)
    prop = np.asarray([np.mean(np.asarray(labels) == i) for i in types])
    prop = pd.DataFrame([types, prop], index=["types", "prop"], columns=types).T
    mtx = confusion_matrix(labels, pred, normalize="true")
    df = pd.DataFrame(mtx, columns=types, index=types)
    df = df.loc[np.unique(labels), np.unique(pred)]
    df = df.rename_axis(
        x_label, axis="columns"
    )  # TODO: double check the axes are correct
    df = df.rename_axis(y_label)
    df.to_csv(res_dir + "/%s_prediction_accuracy.csv" % name)
    plt.figure(figsize=(15, 12))
    sns.heatmap(df, linewidths=0.005, cmap="OrRd")
    plt.tight_layout()
    output_pdf_fn = os.path.join(res_dir, "confusion_matrices.pdf")
    pdf = matplotlib.backends.backend_pdf.PdfPages(output_pdf_fn)
    for fig in range(1, plt.gcf().number + 1):
        pdf.savefig(fig)
    pdf.close()


def make_agreement_plots(adata, methods, save_folder):
    # TODO should this be pulling from resultsadata?

    # clear all existing figures first
    # or else this will interfere with the pdf saving capabilities
    fig_nums = plt.get_fignums()
    for num in fig_nums:
        plt.close(num)

    for method in methods:
        print("Making confusion matrix for {}".format(method))
        x_label = method
        y_label = "consensus_prediction"
        prediction_eval(
            adata.obs[x_label],
            adata.obs[y_label],
            name=method,
            x_label=x_label,
            y_label=y_label,
            res_dir=save_folder,
        )
    plt.close()


def get_pretrained_model_genes(scvi_model_path):
    """
    Get the genes used to train a saved scVI model

    Parameters
    ----------
    scvi_model_path
        Path to saved scvi model

    Returns
    -------
    var_names
        Names of genes used to train the saved scvi model
    """
    varnames_path = os.path.join(scvi_model_path, "var_names.csv")
    var_names = np.genfromtxt(varnames_path, delimiter=",", dtype=str)
    return var_names


def try_method(log_message):
    """
    Decorator which will except an Exception if it failed.
    """

    def try_except(func):
        def wrapper(*args, **kwargs):
            try:
                print("{}.".format(log_message))
                func(*args, **kwargs)
            except Exception as e:
                print("{} failed. Skipping.".format(log_message))
                print(e)

        return wrapper

    return try_except


def subsample_dataset(
    adata,
    labels_key,
    n_samples_per_label=100,
    ignore_label=None,
):
    """
    Subsamples dataset per label to n_samples_per_label.

    If a label has fewer than n_samples_per_label examples, then will use
    all the examples. For labels in ignore_label, they won't be included
    in the resulting subsampled dataset.

    Parameters
    ----------
    adata
        AnnData object
    labels_key
        Key in adata.obs for label information
    n_samples_per_label
        Maximum number of samples to use per label
    ignore_label
        List of labels to ignore (not subsample).

    Returns
    -------
    Returns list of obs_names corresponding to subsampled dataset

    """
    sample_idx = []
    labels, counts = np.unique(adata.obs[labels_key], return_counts=True)

    print("Sampling {} per label".format(n_samples_per_label))

    for label in ignore_label:
        if label in labels:
            idx = np.where(labels == label)
            labels = np.delete(labels, idx)
            counts = np.delete(counts, idx)

    for i, label in enumerate(labels):
        label_locs = np.where(adata.obs[labels_key] == label)[0]
        if counts[i] < n_samples_per_label:
            sample_idx.append(label_locs)
        else:
            label_subset = np.random.choice(
                label_locs, n_samples_per_label, replace=False
            )
            sample_idx.append(label_subset)
    sample_idx = np.concatenate(sample_idx)
    return adata.obs_names[sample_idx]


def check_genes_is_subset(ref_genes, query_genes):
    """
    Check whether query_genes is a subset of ref_genes.

    Parameters
    ----------
    ref_genes
        List of reference genes
    query_genes
        List of query genes

    Returns
    -------
    is_subset
        True if it is a subset, False otherwise.

    """
    if len(set(query_genes)) != len(query_genes):
        print("Warning: Your genes are not unique.")

    if set(ref_genes).issubset(set(query_genes)):
        print("All ref genes are in query dataset. Can use pretrained models.")
        is_subset = True
    else:
        print("Not all reference genes are in query dataset. Retraining models.")
        is_subset = False
    return is_subset


def make_batch_covariate(adata, batch_keys):
    """
    Combines all the batches in batch_keys into a single batch.
    Saves results into adata.obs['_batch']

    Parameters
    ----------
    adata
        Anndata object
    batch_keys
        List of keys in adat.obs corresponding to batches
    """
    adata.obs["_batch"] = ""
    for key in batch_keys:
        v1 = adata.obs["_batch"].values
        v2 = adata.obs[key].values
        adata.obs["_batch"] = [a + b for a, b in zip(v1, v2)]


def process_query(
    query_adata,
    tissue,
    save_folder,
    ref_adata_path,
    ref_labels_key="_labels_annotation",
    ref_batch_keys=["Method", "Donor"],
    ref_layers_key=None,
    ref_cell_ontology_key="final_annotation_cell_ontology_id",
    query_labels_key=None,
    query_batch_key=None,
    query_layers_key=None,
    pretrained_scvi_path=None,
    unknown_celltype_label="unknown",
    training_mode="online",
):
    """
    Processes the query dataset in preperation for the annotation pipeline.
    

    Parameters
    ----------
    query_adata
        AnnData of query cells
    tissue
        Tissue to do the following
    save_folder
        Folder to save data to
    ref_adata_path
        Path to reference AnnData
    ref_labels_key
        Key in obs field of reference AnnData for labels
    ref_batch_keys
        List of Keys (or None) in obs field of reference AnnData to
        use for labels
    ref_layers_key:
        If not None, will use data from ref_adata.layers[ref_layers_key]
    ref_cell_ontology_key
        Key in obs field of reference AnnData for ontology ids
    query_batch_key
        Key in obs field of query adata for batch information.
    query_layers_key
        If not None, will use data from query_adata.layers[query_layers_key]. 
    query_labels_key
        Key in obs field of query adata for label information. 
        This is only used for training scANVI. 
        Make sure to set unknown_celltype_label to mark unlabelled cells.
    unknown_celltype_label
        If query_labels_key is not None, cells with label unknown_celltype_label
        will be treated as unknown and will be predicted by the model.
    pretrained_scvi_path
        Path to pretrained scvi model
    training_mode
        If training_mode=='offline', will train scVI and scANVI from scratch. Else if
        training_mode=='online' and all the genes in the pretrained models are present
        in query adata, will train the scARCHES version of scVI and scANVI, resulting in
        faster training times.
    
    Returns
    -------
    adata 
        AnnData object that is setup for use with the annotation pipeline

    """
    if query_adata.n_obs == 0:
        raise ValueError("Input query anndata has no cells.")

    if training_mode == "online":
        # check if you can use pretrained model
        pretrained_genes = get_pretrained_model_genes(pretrained_scvi_path)
        query_genes = query_adata.var_names.to_numpy().astype("str")
        is_subset = check_genes_is_subset(pretrained_genes, query_genes)
    elif training_mode == "offline":
        is_subset = False

    if is_subset is True:
        print("Training online scvi and scanvi.")
        ref_adata_path = os.path.join(pretrained_scvi_path, "adata.h5ad")
    elif is_subset is False:
        print("Is not subset, training offline.")
        training_mode = "offline"

    ref_adata = anndata.read(ref_adata_path)
    
    if ref_batch_keys is not None:
        make_batch_covariate(ref_adata, ref_batch_keys)
        ref_batch_key = "_batch"

    ref_adata.obs["_labels_annotation"] = ref_adata.obs[ref_labels_key]
    ref_adata.obs["_batch_annotation"] = ref_adata.obs[ref_batch_key]
    ref_adata.obs["_dataset"] = "ref"

    if ref_layers_key is not None:
        ref_adata.layers["scvi_counts"] = ref_adata.layers[ref_layers_key]
    else:
        ref_adata.layers["scvi_counts"] = ref_adata.X

    # subsample the reference cells used for training certain models
    ref_adata.obs["_ref_subsample"] = False
    ref_subsample_idx = subsample_dataset(
        ref_adata, ref_labels_key, ignore_label=[unknown_celltype_label]
    )
    ref_adata.obs["_ref_subsample"][ref_subsample_idx] = True

    ref_adata_not_manual = ref_adata[ref_adata.obs["Manually Annotated"] == "False"].copy()
    ref_adata = ref_adata[ref_adata.obs["Manually Annotated"] == "True"].copy()
    
    if query_batch_key is not None:
        query_adata.obs["_batch_annotation"] = (
            query_adata.obs[query_batch_key].astype("str") + "_query"
        )
    else:
        query_adata.obs["_batch_annotation"] = "query"

    query_adata.obs["_dataset"] = "query"
    query_adata.obs["_ref_subsample"] = False
    query_adata.obs[ref_cell_ontology_key] = unknown_celltype_label

    if query_layers_key is not None:
        query_adata.layers["scvi_counts"] = query_adata.layers[query_layers_key]
    else:
        query_adata.layers["scvi_counts"] = query_adata.X

    if query_labels_key is not None:
        query_labels = query_adata.obs[query_labels_key]
        known_cell_idx = np.where(query_labels != unknown_celltype_label)[0]
        if len(known_cell_idx) != 0:
            query_adata.obs["_labels_annotation"][known_cell_idx] = query_labels[
                known_cell_idx
            ]
        else:
            query_adata.obs["_labels_annotation"] = unknown_celltype_label
    else:
        query_adata.obs["_labels_annotation"] = unknown_celltype_label
    
    if training_mode == "online":
        query_adata = query_adata[:, ref_adata.var_names].copy()
        adata = anndata.concat((ref_adata, query_adata))
    elif training_mode == "offline":
        adata = anndata.concat((ref_adata, query_adata))
    
    adata.uns["_training_mode"] = training_mode

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10, zero_center=False)

    sc.pp.highly_variable_genes(
        adata, n_top_genes=4000, subset=True, layer="scvi_counts", flavor="seurat_v3"
    )
    sc.tl.pca(adata)
    scvi.data.setup_anndata(
        adata,
        batch_key="_batch_annotation",
        labels_key="_labels_annotation",
        layer="scvi_counts",
    )
    
    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
        
    ref_query_results_fn = os.path.join(save_folder, 'annotated_query_plus_ref.h5ad')
    anndata.concat((query_adata, ref_adata), join='outer').write(ref_query_results_fn)
    
    query_results_fn = os.path.join(save_folder, 'annotated_query.h5ad')
    query_adata.write(query_results_fn)
    return adata


def compute_consensus(adata, prediction_keys):
    """
    Computes consensus prediction and statistics between all methods. 
    
    Parameters
    ----------
    adata
        AnnData object
    prediction_keys
        Keys in adata.obs for for predicted values
    
    Returns
    -------
    Saves the consensus prediction in adata.obs['consensus_prediction']
    Saves the consensus percentage between methods in adata.obs['consensus_percentage']
    """
    consensus_prediction = adata.obs[prediction_keys].apply(majority_vote, axis=1)
    adata.obs["consensus_prediction"] = consensus_prediction

    agreement = adata.obs[prediction_keys].apply(majority_count, axis=1)
    agreement *= 100
    adata.obs["consensus_percentage"] = agreement.values.round(2).astype(str)


def majority_vote(x):
    a, b = np.unique(x, return_counts=True)
    return a[np.argmax(b)]


def majority_count(x):
    a, b = np.unique(x, return_counts=True)
    return np.max(b) / np.sum(b)


def annotate_data(
    adata, methods, save_path, pretrained_scvi_path=None, pretrained_scanvi_path=None
):
    if not os.path.exists(save_path):
        os.mkdir(save_path)
        
    ref_query_results_fn = os.path.join(save_path, "annotated_query_plus_ref.h5ad")
    query_results_fn = os.path.join(save_path, "annotated_query.h5ad")

    if "bbknn" in methods:
        run_bbknn(adata, batch_key="_batch_annotation")
        run_knn_on_bbknn(
            adata, labels_key="_labels_annotation", result_key="knn_on_bbknn_pred"
        )
        save_results(
            adata,
            ref_query_results_fn,
            obs_keys=["knn_on_bbknn_pred"],
            obsm_keys=["bbknn_umap"],
        )
        save_results(
            adata,
            query_results_fn,
            obs_keys=["knn_on_bbknn_pred"],
            obsm_keys=["bbknn_umap"],
        )

    if "scvi" in methods:
        training_mode = adata.uns["_training_mode"]
        scvi_obsm_latent_key = "X_scvi_" + training_mode

        run_scvi(
            adata,
            max_epochs=None,
            n_latent=50,
            dropout_rate=0.1,
            dispersion="gene-batch",
            obsm_latent_key=scvi_obsm_latent_key,
            pretrained_scvi_path=pretrained_scvi_path,
        )
        knn_pred_key = "knn_on_scvi_{}_pred".format(training_mode)
        run_knn_on_scvi(adata, obsm_key=scvi_obsm_latent_key, result_key=knn_pred_key)
        save_results(
            adata,
            ref_query_results_fn,
            obs_keys=[knn_pred_key],
            obsm_keys=[scvi_obsm_latent_key, scvi_obsm_latent_key + "_umap"],
        )
        save_results(
            adata,
            query_results_fn,
            obs_keys=[knn_pred_key],
            obsm_keys=[scvi_obsm_latent_key, scvi_obsm_latent_key + "_umap"],
        )

    if "scanvi" in methods:
        training_mode = adata.uns["_training_mode"]
        obsm_latent_key = "X_scanvi_{}".format(training_mode)
        predictions_key = "scanvi_{}_pred".format(training_mode)
        run_scanvi(
            adata,
            max_epochs=None,
            n_latent=100,
            dropout_rate=0.1,
            obsm_latent_key=obsm_latent_key,
            obs_pred_key=predictions_key,
            pretrained_scanvi_path=pretrained_scanvi_path,
        )

        save_results(
            adata,
            ref_query_results_fn,
            obs_keys=[predictions_key],
            obsm_keys=[obsm_latent_key],
        )
        save_results(
            adata,
            query_results_fn,
            obs_keys=[predictions_key],
            obsm_keys=[obsm_latent_key],
        )
    if "svm" in methods:
        run_svm_on_hvg(adata)
        save_results(adata, ref_query_results_fn, obs_keys=["svm_pred"])
        save_results(adata, query_results_fn, obs_keys=["svm_pred"])

    if "rf" in methods:
        run_rf_on_hvg(adata)
        save_results(adata, ref_query_results_fn, obs_keys=["rf_pred"])
        save_results(adata, query_results_fn, obs_keys=["rf_pred"])

    if "onclass" in methods:
        run_onclass(
            adata=adata,
            layer="scvi_counts",
            max_iter=20,
            cl_obo_file="cl.obo",
            cl_ontology_file="cl.ontology",
            nlp_emb_file="cl.ontology.nlp.emb",
        )
        save_results(adata, ref_query_results_fn, obs_keys=["onclass_pred"])
        save_results(adata, query_results_fn, obs_keys=["onclass_pred"])

    if "scanorama" in methods:
        run_scanorama(adata, batch_key="_batch_annotation")
        run_knn_on_scanorama(adata)
        save_results(
            adata,
            ref_query_results_fn,
            obs_keys=["knn_on_scanorama_pred"],
            obsm_keys=["scanorama_umap"],
        )
        save_results(
            adata,
            query_results_fn,
            obs_keys=["knn_on_scanorama_pred"],
            obsm_keys=["scanorama_umap"],
        )

    # Here we compute the consensus statistics
    all_prediction_keys = [
        "knn_on_bbknn_pred",
        "knn_on_scvi_online_pred",
        "knn_on_scvi_offline_pred",
        "scanvi_online_pred",
        "scanvi_offline_pred",
        "svm_pred",
        "rf_pred",
        "onclass_pred",
        "knn_on_scanorama_pred",
    ]

    obs_keys = adata.obs.keys()
    pred_keys = [key for key in obs_keys if key in all_prediction_keys]
    compute_consensus(adata, pred_keys)
    
    save_results(
        adata,
        ref_query_results_fn,
        obs_keys=["consensus_prediction", "consensus_percentage"],
    )
    save_results(
        adata,
        query_results_fn,
        obs_keys=["consensus_prediction", "consensus_percentage"],
    )
    make_agreement_plots(adata, methods=pred_keys, save_folder=save_path)

    
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
            results.obs[key] = adata[results.obs_names].obs[key]
        for key in obsm_keys:
            results.obsm[key] = adata[results.obs_names].obsm[key]
        results.write(results_adata_path, compression)
    else:
        adata.write(results_adata_path, compression)


@try_method("Integrating data with bbknn")
def run_bbknn(adata, batch_key="_batch"):
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
    adata.obsm["bbknn_umap"] = adata.obsm["X_umap"]
    del adata.obsm["X_umap"]
    return adata


@try_method("Classifying with knn on bbknn distances")
def run_knn_on_bbknn(
    adata,
    labels_key="_labels_annotation",
    result_key="knn_on_bbknn_pred",
):
    distances = adata.obsp["distances"]

    ref_idx = adata.obs["_dataset"] == "ref"
    query_idx = adata.obs["_dataset"] == "query"

    ref_dist_idx = np.where(ref_idx == True)[0]
    query_dist_idx = np.where(query_idx == True)[0]

    train_y = adata[ref_idx].obs[labels_key].to_numpy()
    train_distances = distances[ref_dist_idx, :][:, ref_dist_idx]

    knn = KNeighborsClassifier(n_neighbors=2, metric="precomputed")
    knn.fit(train_distances, y=train_y)

    test_distances = distances[query_dist_idx, :][:, ref_dist_idx]
    knn_pred = knn.predict(test_distances)

    # save_results. ref cells get ref annotations, query cells get predicted
    adata.obs[result_key] = adata.obs[labels_key]
    adata.obs[result_key][query_idx] = knn_pred
    print('Saved knn on bbknn results to adata.obs["{}"]'.format(result_key))


@try_method("Classifying with random forest")
def run_rf_on_hvg(
    adata,
    labels_key="_labels_annotation",
    save_key="rf_pred",
    layers_key="scvi_counts",
):
    train_idx = adata.obs["_ref_subsample"]
    test_idx = adata.obs["_dataset"] == "query"

    train_x = adata[train_idx].layers[layers_key]
    train_y = adata[train_idx].obs[labels_key].to_numpy()
    test_x = adata[test_idx].layers[layers_key]

    print("Training random forest classifier with {} cells".format(len(train_y)))
    rf = RandomForestClassifier()
    rf.fit(train_x, train_y)
    rf_pred = rf.predict(test_x)

    adata.obs[save_key] = adata.obs[labels_key]
    adata.obs[save_key][test_idx] = rf_pred


@try_method("Running OnClass")
def run_onclass(
    adata,
    cl_obo_file,
    cl_ontology_file,
    nlp_emb_file,
    labels_key="_labels_annotation",
    layer=None,
    save_key="onclass_pred",
    n_hidden=500,
    max_iter=20,
    save_model="onclass_model",
    shard_size=50000,
):
    celltype_dict, clid_2_name = make_celltype_to_cell_ontology_id_dict(cl_obo_file)
    cell_ontology_obs_key = make_cell_ontology_id(adata, labels_key, celltype_dict)

    train_model = OnClassModel(
        cell_type_nlp_emb_file=nlp_emb_file, cell_type_network_file=cl_ontology_file
    )

    train_idx = adata.obs["_dataset"] == "ref"
    test_idx = adata.obs["_dataset"] == "query"

    if layer is None:
        train_X = adata[train_idx].X.todense()
        test_X = adata[test_idx].X.todense()
    else:
        train_X = adata[train_idx].layers[layer].todense()
        test_X = adata[test_idx].layers[layer].todense()

    train_genes = adata[train_idx].var_names
    train_Y = adata[train_idx].obs[cell_ontology_obs_key]

    test_adata = adata[test_idx]

    _ = train_model.EmbedCellTypes(train_Y)

    model_path = "OnClass"

    train_model.BuildModel(ngene=adata.n_vars)
    train_model.Train(
        train_X, train_Y, save_model=model_path, genes=train_genes, max_iter=max_iter
    )

    test_adata.obs[save_key] = "na"

    if test_adata.n_obs > shard_size:
        for i in range(0, test_adata.n_obs, shard_size):
            tmp_X = test_X[i : i + shard_size]
            onclass_pred = train_model.Predict(tmp_X)
            pred_label_str = [train_model.i2co[l] for l in onclass_pred[2]]
            pred_label_str = [clid_2_name[i] for i in pred_label_str]
            test_adata.obs[save_key][i : i + shard_size] = pred_label_str
    else:
        onclass_pred = train_model.Predict(test_X)
        pred_label_str = [train_model.i2co[l] for l in onclass_pred[2]]
        pred_label_str = [clid_2_name[i] for i in pred_label_str]
        test_adata.obs[save_key] = pred_label_str

    adata.obs[save_key] = adata.obs[labels_key].astype(str)
    adata.obs[save_key][test_adata.obs_names] = test_adata.obs[save_key]

    return adata


@try_method("Classifying with SVM")
def run_svm_on_hvg(
    adata,
    labels_key="_labels_annotation",
    save_key="svm_pred",
    layers_key="scvi_counts",
):
    train_idx = adata.obs["_ref_subsample"]
    test_idx = adata.obs["_dataset"] == "query"

    train_x = adata[train_idx].layers[layers_key]
    train_y = adata[train_idx].obs[labels_key].to_numpy()
    test_x = adata[test_idx].layers[layers_key]

    clf = svm.LinearSVC(max_iter=1000)
    clf.fit(train_x, train_y)
    svm_pred = clf.predict(test_x)

    # save_results
    adata.obs[save_key] = adata.obs[labels_key]
    adata.obs[save_key][test_idx] = svm_pred


@try_method("Running scVI")
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

    # temporary scvi hack
    tmp_mappings = adata.uns["_scvi"]["categorical_mappings"]["_scvi_labels"]
    model.scvi_setup_dict_["categorical_mappings"]["_scvi_labels"] = tmp_mappings

    adata.obsm[obsm_latent_key] = model.get_latent_representation(adata)

    sc.pp.neighbors(adata, use_rep=obsm_latent_key)
    sc.tl.umap(adata)
    adata.obsm[obsm_latent_key + "_umap"] = adata.obsm["X_umap"]
    del adata.obsm["X_umap"]

    if save_folder is not None:
        model.save(save_folder, overwrite=overwrite, save_anndata=save_anndata)


@try_method("Classifying with knn on scVI latent space")
def run_knn_on_scvi(
    adata,
    labels_key="_labels_annotation",
    obsm_key="X_scvi",
    result_key="knn_on_scvi_pred",
):
    if obsm_key not in adata.obsm.keys():
        raise ValueError("Please train scVI first or pass in a valid obsm_key.")

    print(
        "Training knn on scvi latent space. "
        + 'Using latent space in adata.obsm["{}"]'.format(obsm_key)
    )

    ref_idx = adata.obs["_dataset"] == "ref"
    query_idx = adata.obs["_dataset"] == "query"

    train_X = adata[ref_idx].obsm[obsm_key]
    train_Y = adata[ref_idx].obs[labels_key].to_numpy()

    test_X = adata[query_idx].obsm[obsm_key]

    knn = KNeighborsClassifier(n_neighbors=15, weights="uniform")
    knn.fit(train_X, train_Y)
    knn_pred = knn.predict(test_X)

    # save_results
    adata.obs[result_key] = adata.obs[labels_key]
    adata.obs[result_key][query_idx] = knn_pred


@try_method("Classifying with knn on scanorama latent space")
def run_knn_on_scanorama(
    adata,
    labels_key="_labels_annotation",
    obsm_key="X_scanorama",
    result_key="knn_on_scanorama_pred",
):
    print("Running knn on scanorama")
    if obsm_key not in adata.obsm.keys():
        print("Please run scanorama first or pass in a valid obsm_key.")

    ref_idx = adata.obs["_dataset"] == "ref"
    query_idx = adata.obs["_dataset"] == "query"

    train_X = adata[ref_idx].obsm[obsm_key]
    train_Y = adata[ref_idx].obs["_labels_annotation"].to_numpy()

    test_X = adata[query_idx].obsm[obsm_key]

    knn = KNeighborsClassifier(n_neighbors=15, weights="uniform")
    knn.fit(train_X, train_Y)
    knn_pred = knn.predict(test_X)

    # save_results
    adata.obs[result_key] = adata.obs[labels_key]
    adata.obs[result_key][query_idx] = knn_pred


@try_method("Running scanorama")
def run_scanorama(adata, batch_key="_batch"):
    # TODO add check if in colab and n_genes > 120000
    # throw warning
    adatas = [adata[adata.obs[batch_key] == i] for i in np.unique(adata.obs[batch_key])]
    scanorama.integrate_scanpy(adatas, dimred=50)
    tmp_adata = anndata.concat(adatas)
    adata.obsm["X_scanorama"] = tmp_adata[adata.obs_names].obsm["X_scanorama"]

    print("Computing umap on scanorama")
    sc.pp.neighbors(adata, use_rep="X_scanorama")
    sc.tl.umap(adata)
    adata.obsm["scanorama_umap"] = adata.obsm["X_umap"]
    del adata.obsm["X_umap"]


@try_method("Running scANVI")
def run_scanvi(
    adata,
    unlabeled_category="unknown",
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
            unlabeled_category=unlabeled_category,
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
