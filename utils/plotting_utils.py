import pandas as pd

# automcatically regenerate tissue cell fractions plots
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


# plot a sankey diagram using a pandas dataframe
def genSankey(df,cat_cols=[],cat_cols_color = [],value_cols='',title='Sankey Diagram'):
    
    """
    returns a fig object to use with plotly

    cat_cols
        columns to be used as nodes
    cat_cols_color
        List of lists with nodes colors.
        Provide a list with same length as set(category) or a single color
    value_cols
        size of the links
    title
        plot title
    """
    
    
    labelList = []
    colorNumList = []
    for catCol in cat_cols:
        labelListTemp =  list(set(df[catCol].values))
        colorNumList.append(len(labelListTemp))
        labelList = labelList + labelListTemp

    # remove duplicates from labelList
    labelList = list(dict.fromkeys(labelList))

    # define colors based on number of levels
    colorList = []
    if cat_cols_color == []:
        # maximum of 6 value cols -> 6 colors
        colorPalette = ['#4B8BBE','#306998','#FFE873','#FFD43B','#646464']
        for idx, colorNum in enumerate(colorNumList):
            colorList = colorList + [colorPalette[idx]]*colorNum
    else:
        position = 0
        for idx, colorNum in enumerate(colorNumList):
            if type(cat_cols_color[idx]) is dict:#len(cat_cols_color[idx])==colorNum:
                idx_colors = cat_cols_color[idx]
                for i in range(position,position+colorNum):
    #                 colorList = colorList + idx_colors[labelList[i]]
                    colorList.append(idx_colors[labelList[i]])
            else:
                colorList = colorList + [cat_cols_color[idx]]*colorNum
        #         colorList = colorList + [colorPalette[idx]]*colorNum

            position = position + colorNum

    # transform df into a source-target pair
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            tempDf.columns = ['source','target','count']
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
        sourceTargetDf = sourceTargetDf.groupby(['source','target']).agg({'count':'sum'}).reset_index()

    # add index for source-target pair
    sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
    sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))

    # creating the sankey diagram
    data = dict(
        type='sankey',
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(
            color = "black",
            width = 0.5
          ),
          label = labelList,
          color = colorList
        ),
        link = dict(
          source = sourceTargetDf['sourceID'],
          target = sourceTargetDf['targetID'],
          value = sourceTargetDf['count']
        )
      )

    layout =  dict(
        title = title,
        font = dict(
          size = 10
        )
    )

    fig = dict(data=[data], layout=layout)
    return fig