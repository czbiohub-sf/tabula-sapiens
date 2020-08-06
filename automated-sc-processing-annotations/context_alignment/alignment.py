import os

import argparse
import scanpy as sc
import numpy as np
import IPython.display as IP
# import igraph
# import louvain
# import leidenalg
from scvi.inference import UnsupervisedTrainer
from scvi.dataset import AnnDatasetFromAnnData
from scvi.models import VAE

import torch
import scanorama

# default scVI training parameters
use_cuda = True
retrain = False


def main(
        input_file_path,
        output_file_path,
        batch='batch',
        batch_correction='bbknn',
        tissue='all',
        # min_genes=200,
        # min_counts=2500,
        n_neighbors=25,
        n_pcs=20,
        cluster_resolution=5
):
    print("subsetting sapiens data")
    sapiens = sc.read_h5ad("../data/OnClass_data/data_used_for_training/tabula-muris-senis-facs_cell_ontology.h5ad")
    # sapiens[sapiens.obs['tissue'] == tissue].copy()

    print('reading data in')
    adata = sc.read_h5ad(input_file_path)
    IP.display(adata)

    # adata.X = adata.raw.X

    adata = adata.concatenate(sapiens)
    IP.display(adata)
    if tissue != 'all':
        sapiens = sapiens[sapiens.obs['tissue']==tissue].copy()
        if sapiens.n_obs == 0:
            raise ValueError('{} is not a tissue in the atlas.'.format(tissue))

    if batch_correction == 'bbknn':
        print('normalization & scaling')
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        adata = sc.pp.filter_genes_dispersion(adata, subset=False, min_disp=.5, max_disp=None,
                              min_mean=.0125, max_mean=10, n_bins=20, n_top_genes=None,
                              log=True, copy=True)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10, zero_center=False)
        # adata.uns['tissue_colors'] = list(tissue_color_dict.values())

        print('pca')
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pl.pca_overview(adata)

        print('neighs')
        # sc.pp.neighbors(adata)
        sc.external.pp.bbknn(adata, metric='euclidean',
                     batch_key=batch,
                     approx=True,
                     n_pcs=n_pcs, trim=None, n_trees=10,
                     use_faiss=True, set_op_mix_ratio=1.0, local_connectivity=1)

        print('umap computing')
        sc.tl.umap(adata, n_components=2)
        print('test')

        print('clustering')
        sc.tl.louvain(adata, resolution=cluster_resolution)
        sc.tl.leiden(adata, resolution=cluster_resolution)

        print('save h5ad and launch cellxgene')
        adata.write(output_file_path + "bbknn.h5ad")

    if batch_correction == 'scvi':
        model_file = output_file_path + 'scVI.model.pkl'
        print('Define scVI Model')
        data = AnnDatasetFromAnnData(adata, batch_label='batch')
        vae = VAE(data.nb_genes, n_batch=data.n_batches,
              n_layers=3, n_latent=50, dispersion='gene-batch',
              reconstruction_loss='zinb')
        trainer = UnsupervisedTrainer(
            vae,
            data,
            train_size=0.99,
            use_cuda=use_cuda,
            frequency=5,
            data_loader_kwargs={"pin_memory": False},
            n_epochs_kl_warmup=10
        )
        print('train scVI')
        if retrain:
            trainer.train(n_epochs=1)
            torch.save(trainer.model.state_dict(), model_file)
        else:
            if os.path.isfile(model_file):
                trainer.model.load_state_dict(torch.load(model_file))
                trainer.model.eval()
            else:
                trainer.train(n_epochs=1)
                torch.save(trainer.model.state_dict(), model_file)

        posterior = trainer.create_posterior(
            trainer.model, data, indices=np.arange(len(data))
        ).sequential()
        latent = posterior.get_latent()
        adata.obsm['X_scvi'] = latent[0]
        print('umap computing')
        sc.pp.neighbors(adata, n_neighbors=100, n_pcs=30, use_rep="X_scvi")
        sc.tl.umap(adata)

        print('clustering')
        sc.tl.louvain(adata, resolution=cluster_resolution)
        sc.tl.leiden(adata, resolution=cluster_resolution)

        print('save h5ad and launch cellxgene')
        adata.write(output_file_path + "scVI.h5ad")

    if batch_correction == 'scanorama':
        print('normalization & scaling & gene selection')
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key=batch)
        var_select = adata.var.highly_variable_nbatches > 1
        var_genes = var_select.index[var_select]

        adatas = [adata[adata.obs[batch] == i, var_genes] for i in np.unique(adata.obs[batch])]
        integrated = scanorama.integrate_scanpy(adatas, dimred=50)
        integrated = np.concatenate(integrated)
        adata.obsm['X_scanorama'] = integrated

        print('umap computing')
        sc.pp.neighbors(adata, n_neighbors=100, n_pcs=50, use_rep="X_scanorama")
        sc.tl.umap(adata)

        print('clustering')
        sc.tl.louvain(adata, resolution=cluster_resolution)
        sc.tl.leiden(adata, resolution=cluster_resolution)

        print('save h5ad and launch cellxgene')
        adata.write(output_file_path + "scanorama.h5ad")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Processes raw data using cell size factors and log1p normalization.
        ''')
    parser.add_argument(
        '--input', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--output', required=True,
        help='directory where batched corrected file will be written')
    parser.add_argument(
        '--batch_correction', required=False,
        help='batch correction method. One of: bbknn, scvi, scanorama')
    args = parser.parse_args()
    main(args.input, args.output, batch_correction = args.batch_correction)
