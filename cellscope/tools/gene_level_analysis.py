# Gene-level analysis (DEGs, Enrichment, GRN)
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def aucell(adata: sc.AnnData, gene_list: list, use_raw: bool = True):
    import decoupler
    net=pd.DataFrame(dict(source=np.array("AAA"),target=gene_list))
    decoupler.run_aucell(adata,net,source="source",target="target",min_n=1,seed=0,use_raw=use_raw)
    return adata.obsm["aucell_estimate"].loc[:,"AAA"].to_list()


def get_rank_array(adata, key, rank_name="rank_genes_groups"):
    deg_df = adata.uns[rank_name]
    return np.array(deg_df[key].tolist()).flatten()

def find_all_markers(
    adata: sc.AnnData, groupby:str='Cluster', use_raw:bool=True
):
    sc.tl.rank_genes_groups(adata, groupby=groupby,reference='rest',use_raw=use_raw, pts=True, method='wilcoxon',key_added="AllMakers")
    identy = list(adata.uns['AllMakers']["names"].dtype.names)
    degs = pd.DataFrame(
        {
            "Score": get_rank_array(adata, "scores", rank_name='AllMakers'),
            "Pvalue": get_rank_array(adata, "pvals", rank_name='AllMakers'),
            "Padj": get_rank_array(adata, "pvals_adj", rank_name='AllMakers'),
            "LogFC": get_rank_array(adata, "logfoldchanges", rank_name='AllMakers'),
            "Feature": get_rank_array(adata, "names", rank_name='AllMakers'),
        }
    )
    if use_raw:
        degs.insert(loc=0, value=identy * adata.raw.n_vars, column="Identy")
    else:
        degs.insert(loc=0, value=identy * adata.n_vars, column="Identy")

    pts = adata.uns['AllMakers']["pts"].melt(ignore_index=False, value_name="PTS")
    pts_rest = adata.uns['AllMakers']["pts_rest"].melt(
        ignore_index=False, value_name="PTS_Rest"
    )
    pts_rest.drop(columns="variable", inplace=True)
    pts = pd.concat([pts, pts_rest], axis=1)
    pts.insert(
        loc=0, column="TEMP", value=pts.index.values + "XXX" + pts.variable.values
    )
    pts.reset_index(drop=True, inplace=True)

    degs.insert(
        loc=0, column="TEMP", value=degs.Feature.values + "XXX" + degs.Identy.values
    )
    degs = pd.merge(pts, degs, on="TEMP")
    degs.drop(columns=["variable", "TEMP"], inplace=True)

    degs.set_index("Feature", drop=True, inplace=True)
    degs.sort_values(by="Identy", inplace=True)
    return degs

def deseq(
    adata: sc.AnnData,
    cluster_key: str = "CellType",
    sample_key: str = "Sample",
    condition_key: str = "Condition",
    ref_level: str = "Control",
    n_jobs: int = 20,
) -> dict:
    """
    PseudoBulk DESeq2 for scRNA-seq
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import functools

    # remove outliers
    adata = adata[:, adata.X.sum(axis=0) > 10]

    clinical = (
        adata.obs.loc[:, [sample_key, condition_key]]
        .drop_duplicates(subset=sample_key)
        .set_index(sample_key)
        .rename(columns={condition_key: "condition"})
    )
    clinical.index.name = None
    case_level = np.setdiff1d(clinical.condition, ref_level)[0]

    count_df = adata.to_df()
    clusters = adata.obs.loc[:, cluster_key]
    samples = adata.obs.loc[:, sample_key]

    agg_sum = count_df.groupby([clusters, samples]).agg("sum")

    def run_deseq(i, agg_df=agg_sum, clinical=clinical, case=case_level, ref=ref_level,n_jobs=n_jobs):
        counts = agg_df.loc[i, :]
        counts.index.name = None
        _clinical = clinical.loc[counts.index, :]

        # 构建DeseqDataSet 对象
        dds = DeseqDataSet(
            counts=counts,
            clinical=_clinical,
            ref_level=["condition", ref],
            design_factors="condition",
            refit_cooks=True,
            n_cpus=n_jobs,
        )
        # 离散度和log fold-change评估.
        dds.deseq2()
        # 差异表达统计检验分析
        stat_res = DeseqStats(
            dds,
            alpha=0.05,
            cooks_filter=True,
            contrast=["condition", case, ref],
            independent_filter=True,
            n_cpus=n_jobs,
            joblib_verbosity=0,
        )
        stat_res.summary()
        return stat_res.results_df.reset_index().rename(
            columns=dict(
                index="GeneID", log2FoldChange="LFC", padj="FDR", pvalue="Pvalue"
            )
        )

    return {x: run_deseq(x) for x in clusters.unique()}
