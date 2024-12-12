import anndata as ad
import bbknn
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import sys


if __name__ == "__main__":
    dataset = sys.argv[1]

    if dataset == "TNK":
        fname = "/WUH721818ALE6L4/projects/TBK1/02_TNK_subcluster/TNK_bbknn.h5ad"
        fname_out = "/WUH721818ALE6L4/projects/TBK1/02_TNK_subcluster/TNK_bbknn_cellinfo_res12sub04.csv"
        adata = sc.read_h5ad(fname)
        sc.tl.leiden(
            adata, resolution=0.3,
            neighbors_key="bkn_aft",
            key_added=f"bkn_leiden_1.2_sub4",
            restrict_to=("bkn_leiden_1.2", ["4"])
        )
        adata.obs.loc[:, "bkn_leiden_1.2_sub4"].to_csv(fname_out)
