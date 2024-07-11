import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, Union, Optional,Sequence
import pandas as pd
import scanpy as sc
import bioquest
import pathlib
import functools

def dimplot(
    adata: sc.AnnData,
    reduction:str,
    filename:str,
    dim1label:str='UMAP1',
    dim2label:str='UMAP2',
    color:Union[str,list[str],Tuple[str],Sequence[str]]="CellType",
    formats:Union[Tuple[str]]=("pdf", "png"),
    frameon:bool=False,
    outdir: Union[pathlib.PosixPath, str] = pathlib.Path().absolute(),
    width=7,
    height=6,
    dpi:int=300,
    **kwds
):
    _saveimg = bioquest.tl.saveimg(formats=formats, outdir=outdir, dpi=dpi)
    plt.rcParams["figure.figsize"] = (width, height)
    ax = sc.pl.embedding(
        adata,
        basis=reduction,
        color=color,
        show=False,
        frameon=frameon,
        **kwds
    )
    axes = ax if isinstance(ax,list) else [ax]
    for x in axes:
        x.arrow(
            -7,
            -12,
            5 * height / width,
            0,
            head_width=0.5,
            head_length=0.5,
            width=0.1,
            color="black",
        )
        x.arrow(
            -7,
            -12,
            0,
            5,
            head_width=0.5,
            head_length=0.5,
            width=0.1 * height / width,
            color="black",
        )
        x.text(-6.5, -13.5, dim1label, fontdict=dict(weight="bold", color="black"))
        x.text(-8, -11.3, dim2label, fontdict={"fontweight": "bold", "rotation": "vertical"})
    _saveimg(filename=filename, figsize=(width, height))

def plot_batch_effect(
    adata: sc.AnnData,
    *,
    cluster_key: str = 'Sample',
    batch_key: str = "Sample",
    n_jobs:int=8,
    resolution: float = 1.0,
    outdir: Union[pathlib.PosixPath, str] = pathlib.Path().absolute(),
    legend_loc: str = "right margin",  # right margin, on data,
    legend_fontsize: str = "small",  # [‘xx-small’, ‘x-small’, ‘small’, ‘medium’, ‘large’, ‘x-large’, ‘xx-large’]
    n_pcs: int = 20,
    pca_use_hvg: bool = True,
    n_neighbors: int = 15,
    dpi: int = 300,
    formats: Tuple[str] = ("pdf", "png"),
) -> Optional[sc.AnnData]:
    bioquest.tl.mkdir(outdir)
    _saveimg = bioquest.tl.saveimg(formats=formats, outdir=outdir, dpi=dpi)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50, use_highly_variable=pca_use_hvg)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
    _saveimg(f"PCA_variance_ratio")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.tsne(adata,n_jobs=n_jobs)
    # batch_key='Sample';legend_loc='right margin';reduction='umap';outdir=f"{OUTDIR}/batch_effect_before_integratation"
    # legend_fontsize='small';cluster_key='Cluster_before_integratation'
    sc.tl.leiden(adata, key_added=cluster_key, resolution=resolution)

    _dimplot=functools.partial(dimplot,adata=adata,legend_fontsize=legend_fontsize,outdir=outdir,legend_loc=legend_loc)
    _dimplot(reduction='umap',color=batch_key,filename=f"{batch_key}_UMAP")
    _dimplot(reduction='tsne',color=batch_key,filename=f"{batch_key}_TSNE",dim1label='TSNE1',dim2label='TSNE2')
    _dimplot(reduction='umap',color=cluster_key,filename=f"{cluster_key}_UMAP")
    _dimplot(reduction='tsne',color=cluster_key,filename=f"{cluster_key}_TSNE",dim1label='TSNE1',dim2label='TSNE2')
