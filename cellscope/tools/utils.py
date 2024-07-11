import pathlib
import scanpy as sc
import anndata
from typing import Optional, Union

def subset(
        adata: anndata.AnnData,
        subsets: dict,
        inplace: bool = False
) -> Optional[anndata.AnnData]:
    """
    filter/subset a AnnData according to subsets conditions
    """
    _a = adata if inplace else adata.copy()
    for k in subsets:
        v = subsets.get(k)
        if isinstance(v, list):
            _lg = _a.obs[k].isin(v)
            _a = _a[_lg, :]
        else:
            _lg = _a.obs[k].apply(lambda x: eval(v))
            _a = _a[_lg, :]
    return None if inplace else _a

class Chrom_size():
    hg38 = {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
        "chr4": 190214555,
        "chr5": 181538259,
        "chr6": 170805979,
        "chr7": 159345973,
        "chr8": 145138636,
        "chr9": 138394717,
        "chr10": 133797422,
        "chr11": 135086622,
        "chr12": 133275309,
        "chr13": 114364328,
        "chr14": 107043718,
        "chr15": 101991189,
        "chr16": 90338345,
        "chr17": 83257441,
        "chr18": 80373285,
        "chr19": 58617616,
        "chr20": 64444167,
        "chr21": 46709983,
        "chr22": 50818468,
        "chrX": 156040895,
        "chrY": 57227415
    }
    mm10 = {
        "chr1": 195471971,
        "chr2": 182113224,
        "chr3": 160039680,
        "chr4": 156508116,
        "chr5": 151834684,
        "chr6": 149736546,
        "chr7": 145441459,
        "chr8": 129401213,
        "chr9": 124595110,
        "chr10": 130694993,
        "chr11": 122082543,
        "chr12": 120129022,
        "chr13": 120421639,
        "chr14": 124902244,
        "chr15": 104043685,
        "chr16": 98207768,
        "chr17": 94987271,
        "chr18": 90702639,
        "chr19": 61431566,
        "chrX": 171031299,
        "chrY": 91744698,
    }
