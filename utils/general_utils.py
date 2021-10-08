import scanpy as sc
import numpy as np
# create a list of letters
def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in range(ord(c1), ord(c2) + 1):
        yield chr(c)

# create list of wells in range
def well_range(c1,n1,c2,n2):
    well_order_by_columns = [
            f"{w}{n:1}" for n in range(n1, n2) for w in char_range(c1, c2)
        ]
    return well_order_by_columns

# transform all categorical columns into strings
def remove_cats(adata):

    cat_columns = adata.obs.select_dtypes(['category']).columns
    adata.obs[cat_columns] = adata.obs[cat_columns].astype(str)
    
    return adata

# balanced downsampling 
import logging as logg
def downsample_to_smallest_category(
        dataMatrix, 
        classOfInterest = 'classification_group',
        random_state = None,
        min_cells = 15,
        keep_small_categories = True
) -> sc.AnnData:
    """
    returns an annData object in which all categories in 'classOfInterest' have
    the same size
    classOfInterest
        column with the categories to downsample
    min_cells
        Minimum number of cells to downsample.
        Categories having less than `min_cells` are discarded unless
        keep_small_categories is True
    keep_small_categories
        Be default categories with less than min_cells are discarded.
        Set to true to keep them
    """

    counts = dataMatrix.obs[classOfInterest].value_counts(sort=False)
    if len(counts[counts < min_cells]) > 0 and keep_small_categories is False:
        logg.warning(
            "The following categories have less than {} cells and will be "
            "ignored: {}".format(min_cells, dict(counts[counts < min_cells]))
        )
    min_size = min(counts[counts >= min_cells])
    sample_selection = None
    for sample, num_cells in counts.items():
        if num_cells <= min_cells:
            if keep_small_categories:
                sel = dataMatrix.obs.index.isin(
                    dataMatrix.obs[dataMatrix.obs[classOfInterest] == sample].index)
            else:
                continue
        else:
            sel = dataMatrix.obs.index.isin(
                dataMatrix.obs[dataMatrix.obs[classOfInterest] == sample]
                .sample(min_size, random_state=random_state)
                .index
            )
        if sample_selection is None:
            sample_selection = sel
        else:
            sample_selection |= sel
    logg.info(
        "The cells in category {!r} had been down-sampled to have each {} cells. "
        "The original counts where {}".format(classOfInterest, min_size, dict(counts))
    )
    return dataMatrix[sample_selection].copy()

# remove duplicated barcodes from 10X
def FindUniqueCells(tenx): 
    barcodes = tenx.obs.index
    barcodes = [x.split('-')[0] for x in barcodes]

    unique_barcodes, count = np.unique(barcodes, return_counts=True)
    is_unique = [x in unique_barcodes[count==1] for x in barcodes]
    return(is_unique)
