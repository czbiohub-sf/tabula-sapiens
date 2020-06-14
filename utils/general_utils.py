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