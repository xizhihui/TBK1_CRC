import anndata as ad
import bbknn
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import sys


if __name__ == "__main__":
    # check inputs and outputs ----
    cellinfo_file = sys.argv[1]
    outprefix = sys.argv[2]
    outfile = [
        outprefix + "_bbknn_cellinfo.csv",
        outprefix + "_bbknn.h5ad"
    ]

    if not os.path.exists(cellinfo_file):
        raise FileNotFoundError(f"{cellinfo_file} not found.")

    if os.path.exists(outfile[0]) or os.path.exists(outfile[1]):
        outfile = ", ".join(outfile)
        raise FileExistsError(f"{outfile} may exist, please check!")


    # load the QCed matrix and cellinfo ----
    mtx_qc = "/WUH721818ALE6L4/projects/TBK1/01_main_cluster/matrix_after_QC"
    adata = sc.read_10x_mtx(mtx_qc, var_names="gene_symbols")
    cellinfo = pd.read_csv(cellinfo_file, sep="\t", index_col=0, header=0)
    adata.obs.loc[cellinfo.index, cellinfo.columns] = cellinfo


    # preprocess ----
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    sc.tl.pca(adata)

    # bbknn with redge regression ----
    bbknn.bbknn(adata, batch_key='sample', key_added="bkn_pre")
    sc.tl.leiden(adata, resolution=1, neighbors_key="bkn_pre", key_added="bkn_pre_leiden_1")
    bbknn.ridge_regression(adata, batch_key=['sample'], confounder_key=['bkn_pre_leiden_1'])
    sc.pp.pca(adata)
    bbknn.bbknn(adata, batch_key='sample', key_added="bkn_aft")
    sc.tl.umap(adata, neighbors_key="bkn_aft")

    # get clusters ----
    clusters = []
    for resolution in np.array(range(4, 21, 2)) / 10:
        sc.tl.leiden(
            adata, resolution=resolution, neighbors_key="bkn_aft", 
            key_added=f"bkn_leiden_{resolution}"
        )
        clusters.append(f"bkn_leiden_{resolution}")
    clustersdf = adata.obs.loc[:, clusters]
    clustersdf.loc[:, ["bkn50UMAP_1", "bkn50UMAP_2"]] = adata.obsm["X_umap"]

    # save object and results ----
    clustersdf.to_csv(outfile[0])
    adata.write_h5ad(outfile[1])
