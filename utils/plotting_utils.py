#WIP
def tiss_cell_fractions(adata,
                        technology_col='method',
                        groupby='tissue',
                        category='cell_ontology_class',
                        dataset="Pilot2"):
    
    if technology_col:
        adata = adata[adata.obs[technology_col]=='10X'].copy()
    
    for t in list(set(adata.obs[groupby])):
        print(t)
        tiss = adata[adata.obs[groupby] == t].copy()
        
        aux = tiss.obs.groupby([category]).count()
        aux = pd.DataFrame(aux).reset_index()
        aux['sample'] = aux[aux.columns[1]]/aux[aux.columns[1]].sum()

        f, ax = plt.subplots(figsize=(15,10)) 
        g = sns.barplot(data = aux, y = 'sample',x = category, ax = ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90);
        ax.set(xlabel= dataset +" "+ t +' cell types', ylabel='Relative abundance in 10X data');
        plt.tight_layout()
        plt.savefig("./cell_fractions/"+save+"_"+t+'.pdf')