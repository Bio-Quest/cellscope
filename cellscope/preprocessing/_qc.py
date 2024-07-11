import itertools
import re
import pathlib
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import functools
import bioquest
from typing import Tuple, Union, Optional, List

def flag_gene_family(
    adata: sc.AnnData,
    *,
    gene_family_name: str,
    gene_family_pattern: str = None,
    gene_list: list = None,
) -> None:
    """
    Flags a gene or gene_family in .var with boolean. (e.g all mitochondrial genes).
    Please only choose gene_family prefix or gene_list

    Parameters
    ----------
        adata
            AnnData object
        gene_family_name
            name of columns in .var where you want to store informationa as a boolean
        gene_family_pattern
            pattern of the gene family (eg. mt- for all mitochondrial genes in mice)
        gene_list
            list of genes to flag in `.var`
    Returns
    -------
        adds the boolean column in `.var`

    """
    if gene_family_pattern:
        adata.var[gene_family_name] = adata.var.index.str.contains(pat=gene_family_pattern,flags=re.IGNORECASE,regex=True)
    if gene_list:
        adata.var[gene_family_name] = adata.var.index.isin(gene_list)


def fastqc(
    adata: sc.AnnData,
    *,
    qc_vars: list,
    sample: str = "Sample",
    outdir: Union[pathlib.PosixPath,str] = pathlib.Path().absolute(),
    min_genes: int = 200,
    min_cells: int = 3,
    percent_top: str = (50,),
    log1p: bool = True,
    dpi: int = 300,
    inplace: bool = True,
    formats: Union[str, Tuple[str, ...]] = ("pdf", "png"),
) -> Optional[sc.AnnData]:
    """
    fastqc
    sk.pp.fastqc(adata_spatial,sample="Sample",outdir=outdir,mitochondrion=True)
    """
    bioquest.tl.mkdir(outdir)
    _saveimg = bioquest.tl.saveimg(
        formats=formats, outdir=outdir, dpi=dpi
    )
    n_samples = len(adata.obs[sample])
    _adata = adata if inplace else adata.copy()
    # Cell Cycle Phase Classification
    # cell_phase(_adata)
    if min_cells:
        sc.pp.filter_cells(_adata, min_genes=min_genes)
    if min_cells:
        sc.pp.filter_genes(_adata, min_cells=min_cells)

    sc.pp.calculate_qc_metrics(
        _adata, qc_vars=qc_vars, percent_top=percent_top, log1p=log1p, inplace=True
    )

    # /* scatter plot */
    if qc_vars:
        ks = ["total_counts", "n_genes_by_counts"] + [
            f"pct_counts_{x}" for x in qc_vars
        ]
        iters = list(itertools.product(ks[0:2], ks[2::])) + [
            ("total_counts", "n_genes_by_counts")
        ]
        _n = len(iters)
        _, axes = plt.subplots(1, _n, figsize=(5 * _n, 3))
        for i, (x, y) in enumerate(iters):
            sc.pl.scatter(
                _adata,
                x=x,
                y=y,
                color=sample,
                ax=axes[i],
                legend_loc="none",
                show=False,
            )
        plt.subplots_adjust(wspace=0.5)
        _saveimg("QC_Scatter")
    else:
        sc.pl.scatter(
            _adata, x=ks[0], y=ks[1], color=sample, legend_loc="none", show=False
        )
        _saveimg("QC_Scatter")

    # /* violin plot */
    _, axes = plt.subplots(1, _n, figsize=(6 * _n, 6))
    for a, k in enumerate(ks):
        sc.pl.violin(
            _adata,
            keys=k,
            groupby=sample,
            rotation=90,
            jitter=False,
            show=False,
            stripplot=False,
            ax=axes[a],
        )
    plt.subplots_adjust(wspace=0.5)
    _saveimg("QC_Violin")

    return None if inplace else _adata


def mad_filter(adata: sc.AnnData, *metric_nmad: Tuple[str, int], **kwds) -> sc.AnnData:
    """
    use median ± n * mad as indicators to filter out outliers
    """
    batch_key = kwds.get("batch_key", None)

    def is_outlier(
        adata, metric: str, nmads: int, batch_key: Optional[str] = batch_key
    ):
        def helper(m: pd.Series, nmads: int = nmads):
            from scipy.stats import median_abs_deviation

            median_ = m.median()
            n_mad = median_abs_deviation(m) * nmads
            return np.logical_or(m < median_ - n_mad, m > median_ + n_mad)

        if batch_key:
            batches = adata.obs.loc[:, batch_key]
            return adata.obs.loc[:, metric].groupby(batches).apply(helper).droplevel(0)
        return helper(adata.obs.loc[:, metric])

    outliers = functools.reduce(
        np.logical_or, [[is_outlier(adata, x, y)] for x, y in metric_nmad]
    )
    return adata[~outliers]
