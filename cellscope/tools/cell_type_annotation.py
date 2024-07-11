from io import StringIO
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from typing import Optional
import matplotlib.pyplot as plt
import bioquest as bq


def add_label(
    adata: anndata.AnnData, 
    annotation: pd.DataFrame,
    reference_key: str, 
    cell_type_key: str = 'CellType',
    ) -> Optional[anndata.AnnData]:

    '''
    labeled(adata,cluster_names=new_cluster_names,reference_key='leiden',cell_type_key='CellType')
    '''
    _annot_df = pd.read_csv(annotation,sep='\t', dtype='object')
    _annot_df['Cluster']=_annot_df.Cluster.str.split(',')
    _annot_df=_annot_df.explode('Cluster')

    _ref_df = adata.obs.loc[:, [reference_key]]
    adata.obs[cell_type_key] = pd.merge(_ref_df, _annot_df, on=reference_key, how='left')['CellType'].values

    return adata[~adata.obs['CellType'].isna(),]



