import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, Union, Optional
import pandas as pd
import scanpy as sc

def cell_ratio(
    adata, 
    x:str, 
    y:str,
    *,
    palette = None,
    normalize = True,
    od=None,
    legend=True,
    figsize=(6, 3)
    ):

    df = adata.obs.loc[:,[x,y]]
    x_items = sorted(df[x].unique().tolist())
    y_items = sorted(df[y].unique().tolist())

    if palette is None:
        palette = dict(zip(y_items, adata.uns[f'{y}_colors']))

    heights = []
    for x_item in x_items:
        tmp_result = []
        x_item_counter = df[df[x]==x_item][y].value_counts().to_dict()
        for y_item in y_items:
            tmp_result.append(x_item_counter.get(y_item, 0))
        heights.append(tmp_result)
    heights = np.asarray(heights)

    if normalize:
        heights = heights/np.sum(heights, axis=0)
    heights = (heights.T/np.sum(heights, axis=1)).T

    plt.figure(figsize=figsize)
    _last = np.matrix([0.]* heights.shape[0])
    for i, y_item in enumerate(y_items):
        p = plt.bar(range(0, heights.shape[0]), heights[:, i],
                    bottom=np.asarray(_last)[0],
                    color=palette.get(y_item, 'b'),
                    label=y_item
                    )
        _last = _last + np.matrix(heights[:, i])
    plt.xticks(range(0, len(x_items)),labels=x_items,rotation=90)
    plt.ylim((0, 1))
    if legend:
        plt.legend()
        ax = plt.gca()
        ax.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    if od is not None:
        plt.savefig(od + f'/{x}_{y}_cell_ratio.pdf', dpi=300,bbox_inches="tight")
    return pd.DataFrame(heights,columns=y_items,index=x_items)