# this is hot to set up the anndata object to run scVI 
scvi.data.setup_anndata(
    adata, batch_key="_batch", labels_key="_labels", layer="counts"
)

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