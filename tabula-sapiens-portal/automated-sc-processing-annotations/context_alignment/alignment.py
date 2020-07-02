import argparse
import scanpy as sc
# import numpy
import IPython.display as IP
# import igraph
# import louvain
# import leidenalg


def main(
    input_file_path,
    output_file_path,
    batch = 'batch',
    batch_correction = 'bbknn',
    tissue = 'all',
    # min_genes=200,
    # min_counts=2500,
    n_neighbors=25,
    n_pcs=20,
    cluster_resolution=5
):
    print("subsetting sapiens data")
    sapiens = sc.read_h5ad("../data/tabula-sapiens.h5ad")
    if tissue != 'all':
        try:
            sapiens = sapiens[sapiens.obs['tissue']==tissue].copy()

            print('reading data in')
            adata = sc.read_h5ad(input_file_path)
            IP.display(adata)

            adata.X = adata.raw.X

            adata = adata.concatenate(sapiens)
            IP.display(adata)

            if batch_correction == 'bbknn':
                print('normalization & scaling')
                sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
                adata = sc.pp.filter_genes_dispersion(adata, subset = False, min_disp=.5, max_disp=None, 
                                        min_mean=.0125, max_mean=10, n_bins=20, n_top_genes=None, 
                                        log=True, copy=True)
                sc.pp.log1p(adata)
                sc.pp.scale(adata, max_value=10, zero_center=False)
                # adata.uns['tissue_colors'] = list(tissue_color_dict.values())

                print('pca')
                sc.tl.pca(adata,svd_solver='arpack')
                sc.pl.pca_overview(adata)

                print('neighs')
                # sc.pp.neighbors(adata)
                sc.external.pp.bbknn(adata, 
                n_neighbors=n_neighbors,
                                    batch_key=batch, 
                                    approx=True, metric='angular',
                                    n_pcs=n_pcs, trim=None, n_trees=10, 
                                    use_faiss=True, set_op_mix_ratio=1.0, local_connectivity=1)

                print('umap computing')
                sc.tl.umap(adata,n_components=2)

                print('clustering')
                sc.tl.louvain(adata, resolution = cluster_resolution)
                sc.tl.leiden(adata, resolution = cluster_resolution)

                print('save h5ad and launch cellxgene')
                adata.write(output_file_path)

                
        except:
            print("Tissue selected not in the atlas")

    
    


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
    args = parser.parse_args()

    main(args.input, args.output)
