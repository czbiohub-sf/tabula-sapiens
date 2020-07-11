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
    min_genes=200,
    min_counts=2500,
    n_neighbors=25,
    n_pcs=20,
    cluster_resolution=5
):

    print('reading data in')
    adata = sc.read_h5ad(input_file_path)
    IP.display(adata)

    print('pre-processing')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)

    IP.display(adata)

    print('normalization & scaling')
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    adata_nonans = sc.pp.filter_genes_dispersion(
        adata, subset=False, min_disp=.5, max_disp=None,
        min_mean=.0125, max_mean=10, n_bins=20, n_top_genes=None,
        log=True, copy=True)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10, zero_center=False)

    # print('pre-processing visualizations')  # optional step
    # sc.pl.violin(
    #     adata_nonans, [
    #         'n_genes', 'n_counts'], jitter=0.4, multi_panel=True)
    # sc.pl.scatter(adata_nonans, x='n_counts', y='n_genes')

    # sc.pl.violin(adata_nonans, ['n_counts', 'n_genes'],
    #              groupby='tissue', size=2, log=True)
    adata.raw = adata

    print('pca')
    sc.tl.pca(adata)
    # # optional
    # sc.pl.pca_overview(adata)

    print('neighs')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    print('umap computing')
    sc.tl.umap(adata, n_components=2)

    print('clustering')
    sc.tl.louvain(adata, resolution=cluster_resolution)
    sc.tl.leiden(adata, resolution=cluster_resolution)

    print('save h5ad and launch cellxgene')
    adata.write(output_file_path)


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
        help='directory where arrow files should be written')
    args = parser.parse_args()

    main(args.input, args.output)
