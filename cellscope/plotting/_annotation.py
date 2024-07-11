import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, Union, Optional,Sequence
import pandas as pd
import scanpy as sc
import bioquest
import pathlib
import functools
from ._emmbeding import dimplot

def plot_marker(
    adata: sc.AnnData, 
    annotation: Union[pathlib.PosixPath, str],
    outdir:Union[pathlib.PosixPath,str] = pathlib.Path().absolute(),
    cell_type_key: str = 'CellType',
    marker_key:str = 'Marker',
    formats:tuple = ('pdf','png'),
    palette:str='Reds'
    ):
    bioquest.tl.mkdir(outdir)
    _saveimg = bioquest.tl._saveimg(formats=formats,outdir=outdir,dpi=300)

    _annot_df = pd.read_csv(annotation,sep='\t', dtype='object')

    marker_dict = {}
    for index,row in _annot_df.iterrows():
        marker_dict[row[cell_type_key]] = row[marker_key].split(',')

    sc.pl.dotplot(adata,var_names=marker_dict,groupby=cell_type_key,show=False,cmap=palette)
    _saveimg('MakerDotPlot')
    sc.pl.stacked_violin(adata,var_names=marker_dict,groupby=cell_type_key,show=False,cmap=palette)
    _saveimg('MakerStackedViolin')
    sc.pl.matrixplot(adata,var_names=marker_dict,groupby=cell_type_key,show=False,cmap=palette)
    _saveimg('MakerMatrixplot')
    sc.pl.heatmap(adata,var_names=marker_dict,groupby=cell_type_key,show=False,cmap=palette)
    _saveimg('MakerHeatmap')
    for k,v in marker_dict:
        dimplot(adata,reduction='X_umap',outdir=outdir,filename=f'MakerDimPlot_{k}',color=v)

def auc_heatmap(adata,marker,out_prefix,ref_key="Cluster",figsize=(12,6),use_raw=True):
    import decoupler
    net=marker.melt(var_name="source",value_name="target").dropna()
    decoupler.run_aucell(adata,net,source="source",target="target",min_n=1,seed=1314,use_raw=use_raw)
    dt2=adata.obsm["aucell_estimate"].groupby(by=adata.obs.loc[:,ref_key]).agg(np.mean)
    import seaborn
    seaborn.clustermap(dt2.T,method='complete',z_score=0,cmap="viridis",figsize=figsize);
    plt.savefig(f"{out_prefix}.pdf",bbox_inches='tight')
    dt2.index.name="CellType"
    dt2.to_csv(f"{out_prefix}_score.csv.gz")

def score_heatmap(adata,marker_df,reference_key="Cluster",figsize=(9,6),return_score=False,save_fig=False):
    obs = adata.obs
    markers_dict = {x:np.intersect1d(marker_df.loc[:,x].dropna(),adata.raw.var_names) for x in  marker_df.columns}
    for x in markers_dict.keys():
        sc.tl.score_genes(adata,gene_list=markers_dict[x],score_name=f"{x}_Marker_Score")
    dt = bq.tl.select(adata.obs,columns=[reference_key],pattern="_Marker_Score$")
    adata.obs = obs
    a=dt.groupby(by=reference_key).apply(np.mean,axis=0)
    a.columns = bq.st.removes(string=a.columns,pattern=r"_Marker_Score$")
    import seaborn as sns
    sns.clustermap(a.T,method='complete',standard_scale=0,cmap="viridis",figsize=figsize);
    if return_score:
        return dt
    if save_fig:
        plt.savefig(f"{save_fig}/anno_heatmap.pdf",bbox_inches='tight')
    a.to_csv("score.csv.gz")