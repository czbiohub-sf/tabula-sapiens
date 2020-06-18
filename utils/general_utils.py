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

# convert ensembl to gene symbol
def convert_ensembl_symbol(adata,species='human',idcol = 'ensembls'):
    import mygene
    mg = mygene.MyGeneInfo()

    mygene_converter = mg.querymany(list(adata.var[idcol]),scopes='all', species=species, as_dataframe=True)
    mygene_converter.loc[mygene_converter['notfound']==True,'symbol'] = mygene_converter.loc[mygene_converter['notfound']==True].index

    adata.var = adata.var.merge(
        mygene_converter.reset_index(),left_on='ensembls',right_on='query').sort_values(
        by='_score',ascending=False).drop_duplicates(
        'ensembl_id').set_index('symbol')
    
    return adata

# transform all categorical columns into strings
def remove_cats(adata):

    cat_columns = adata.obs.select_dtypes(['category']).columns
    adata.obs[cat_columns] = adata.obs[cat_columns].astype(str)
    
    return adata

# balanced downsampling 
import logging as logg
def downsample_to_smallest_category(
        adata,
        column="sample_short",
        random_state=None,
        min_cells=15,
        keep_small_categories=False
) -> sc.AnnData:
    """
    returns an annData object in which all categories in 'column' have
    the same size

    column
        column with the categories to downsample
    min_cells
        Minimum number of cells to downsample.
        Categories having less than `min_cells` are discarded unless
        keep_small_categories is True
    keep_small_categories
        Be default categories with less than min_cells are discarded.
        Set to true to keep them
    """
    counts = adata.obs[column].value_counts(sort=False)
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
                sel = adata.obs.index.isin(
                    adata.obs[adata.obs[column] == sample].index)
            else:
                continue
        else:
            sel = adata.obs.index.isin(
                adata.obs[adata.obs[column] == sample]
                .sample(min_size, random_state=random_state)
                .index
            )
        if sample_selection is None:
            sample_selection = sel
        else:
            sample_selection |= sel
    logg.info(
        "The cells in category {!r} had been down-sampled to have each {} cells. "
        "The original counts where {}".format(column, min_size, dict(counts))
    )
    return adata[sample_selection].copy()
