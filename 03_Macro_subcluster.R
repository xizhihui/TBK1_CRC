dir.create("/WUH721818ALE6L4/projects/TBK1/03_Macro_subcluster", showWarnings = FALSE)
setwd("/WUH721818ALE6L4/projects/TBK1/03_Macro_subcluster")


suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(ggsci)
  library(magrittr)
})

sces <- readRDS("/WUH721818ALE6L4/projects/TBK1/01_main_cluster/main_harmony.rds")
sces_cellinfo <- readRDS("/WUH721818ALE6L4/projects/TBK1/01_main_cluster/main_final_cellinfo.rds")
sces@meta.data[, colnames(sces_cellinfo$metas)] <- sces_cellinfo$metas

macro <- subset(sces, subset = celltypes == "macro_Macro")
meta2keep <- colnames(sces@meta.data)[1:6]
meta2keep <- c(meta2keep, "sample2", "group")
macro <- rusinglecell::clean_seurat(macro, metas = meta2keep)
write.table(macro@meta.data[, 1:6], file = "cellinfo_Macro_raw.tsv", sep = "\t", quote = FALSE)

# main cluster ----
basic_clustering <- function (sce, resolutions = c(1.2, 1, 0.8, 0.6), reduction = "pca", dims = 30) {
  if (length(dims) == 1) {
    dims <- 1:dims
  }

  sce %>% 
    FindNeighbors(dims = dims, reduction = reduction) %>% 
    FindClusters(resolution = resolutions) %>%
    RunUMAP(dims = dims, reduction = reduction) %>% 
    RunTSNE(dims = dims, reduction = reduction, check_duplicates = FALSE)
}


## global params ----
resolutions <- seq(0.4, 2, 0.2)
vars2regress <- c("nFeature_RNA", "nCount_RNA", "pct.mito", "pct.ribo")


## batch correction: harmony ----
schm <- SCTransform(macro, vars.to.regress = vars2regress) %>%
    rusinglecell::check_vargenes(remove = TRUE) %>%
    RunPCA() %>%
    harmony::RunHarmony(
        group.by.vars = "sample",
        assay.use = "SCT",
        max.iter.harmony = 30
    ) %>%
    rusinglecell::get_elbow_pc(name = "harmony") %>%
    basic_clustering(
        resolution = resolutions,
        reduction = "harmony",
        dims = .@misc$harmony.pcs
    )
ggsave("elbow_harmony.png", width = 9, height = 6, plot = schm@misc$harmony.elbow)
saveRDS(schm, file = "macro_harmony.rds")


## batch correction: multi-cca ----
sccc <- SplitObject(macro, "sample") %>%
    lapply(SCTransform, vars.to.regress = vars2regress)
sccc_feat <- SelectIntegrationFeatures(sccc)
sccc <- PrepSCTIntegration(sccc, anchor.features = sccc_feat)
sccc_anchor <- FindIntegrationAnchors(sccc, anchor.features = sccc_feat, normalization.method = "SCT")
sccc <- IntegrateData(sccc_anchor, normalization.method = "SCT") %>%
    RunPCA() %>%
    rusinglecell::get_elbow_pc(name = "pca") %>%
    basic_clustering(
        resolution = resolutions,
        reduction = "pca",
        dims = .@misc$pca.pcs
    )
ggsave("elbow_SeuratCCA.png", width = 9, height = 6, plot = sccc@misc$pca.elbow)
saveRDS(sccc, file = "macro_SeuratCCA.rds")


## batch correction: bbknn ----
bbknn_args = c(
    "python3",
    "/WUH721818ALE6L4/projects/TBK1/run_bbknn.py",
    "/WUH721818ALE6L4/projects/TBK1/03_Macro_subcluster/cellinfo_Macro_raw.tsv",
    "/WUH721818ALE6L4/projects/TBK1/03_Macro_subcluster/macro"
)
system(paste0(bbknn_args, collapse = " "))


# check clusters ----
GrayRed <- c("lightgray", "red")
DefaultAssay(schm) <- "RNA"
schm <- NormalizeData(schm)
DefaultAssay(sccc) <- "RNA"
sccc <- NormalizeData(sccc)
bkn <- read.csv("/WUH721818ALE6L4/projects/TBK1/03_Macro_subcluster/macro_bbknn_cellinfo.csv")
schm[["bkn50umap"]] <- bkn[, c("bkn50UMAP_1", "bkn50UMAP_2")] %>% 
    as.matrix() %>% set_rownames(bkn$X) %>%
    CreateDimReducObject()
schm@meta.data[, colnames(bkn)[2:10] <- bkn[, colnames(bkn)[2:10]]


DimPlot(schm, group.by = "sample", reduction = "umap") + scale_color_igv()
ggsave("umap_harmony_sample.png", width = 7, height = 6)
clustersA <- grep("SCT_snn_res", colnames(schm@meta.data), value = T)
DimPlot(schm, group.by = clustersA, label = T) & scale_color_igv() & NoLegend()
ggsave("umap_harmony_clusters.png", width = 12, height = 12)
mks_harmony <- lapply(setNames(nm = clustersA), rusinglecell::get_markers, sceobj = schm)

DimPlot(schm, group.by = "sample", reduction = "bkn50umap") + scale_color_igv()
ggsave("umap_bbknn_sample.png", width = 7, height = 6)
clustersB <- grep("bkn_leiden", colnames(schm@meta.data), value = T)
DimPlot(schm, group.by = clustersB, label = T) & scale_color_igv() & NoLegend()
ggsave("umap_bbknn_clusters.png", width = 12, height = 12)
mks_bbknn <- lapply(setNames(nm = clustersB), rusinglecell::get_markers, sceobj = schm)

DimPlot(sccc, group.by = "sample", reduction = "umap") + scale_color_igv()
ggsave("umap_SeuratCCA_sample.png", width = 7, height = 6)
clustersC <- grep("integrated_snn_res", colnames(sccc@meta.data), value = T)
DimPlot(sccc, group.by = clustersC, label = T) & scale_color_igv() & NoLegend()
ggsave("umap_SeuratCCA_clusters.png", width = 12, height = 12)
mks_seuratcca <- lapply(setNames(nm = clustersC), rusinglecell::get_markers, sceobj = sccc)


## check the features, cluster markers in each clusters mannually ----
View(mks_harmony[[1]])
View(mks_bbknn[[1]])
View(mks_seuratcca[[1]])

# c("Cd3d", "Cd3e", "Nkg7") %>%
# c("APC", "Malat1") %>%
# c("Arg1", "Ccl6", "Cxcl3", "Ccl9", "Spp1", "Pf4", "Cxcl1", "Hilpda", "Malat1") %>%
# c("C1qa", "C1qb", "C1qc", "Apoe", "Ctsl", "Ctsb", "Lgmn", "Ccl7") %>%
# c("Il1b", "Plac8", "Isg15","Ifitm3", "Ifitm2", "Ifi27l2a", "Ifitm6", "Tnf", "Hp") %>%
c("H2-Ab1", "H2-Aa", "Cd74", "H2-Eb1", "H2-DMb1", "Il1b") %>%
  FeaturePlot(schm, ., reduction = "bkn50umap", cols = GrayRed)
  # FeaturePlot(schm, ., reduction = "hm18umap", cols = GrayRed)
  # FeaturePlot(sccc, ., reduction = "cca13umap", cols = GrayRed)


## final cluster with BBKNN ----
schm$finalcluster <- as.numeric(as.character(schm$bkn_0.6))
# CD3+ cluster appeared starting from the resolution 0.8 (cluster 9).
schm$finalcluster[schm$bkn_0.8 == "9"] <- max(schm$finalcluster) + 1
schm$finalcluster <- rusinglecell::resortcluster(as.character(schm$finalcluster))

subtypes <- setNames(nm = as.character(0:7))
subtypes["4"] <- "Macro_C1_C1qa"
subtypes["0"] <- "Macro_C2_Cxcl3"
subtypes["1"] <- "Macro_C3_H2-Eb1"
subtypes["5"] <- "Macro_C4_Isg15"
subtypes["6"] <- "Macro_C5_Ifitm6"
subtypes["2"] <- "Macro_C6_Malat1"
subtypes["3"] <- "Macro_C6_Malat1"
subtypes["7"] <- "doublets"

schm$subtypes <- subtypes[as.character(schm$finalcluster)]
schm$subtypes <- factor(schm$subtypes, sort(unique(subtypes)))
subtypes_cols <- setNames(
    c("lightgray", as.character(BuenColors::jdb_palette("corona", 6))), 
    sort(unique(subtypes))
)
schm[["umap"]] <- Embeddings(schm, "bkn50umap") %>%
  magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>%
  CreateDimReducObject()


# save the cellinfo ----
obj2save <- list()
obj2save$metas <- schm@meta.data[, c("finalcluster", "subtypes", "sample2", "group")]
obj2save$celltypes_cols <- subtypes_cols
obj2save$umap <- schm@reductions$umap
saveRDS(obj2save, file = "macro_final_cellinfo.rds")


# Figures ----
schm2 <- subset(schm, subset = subtypes != "doublets")

## clusters and features ----
UMAPPlot(schm2, group.by = "subtypes", cols = subtypes_cols) + ggtitle("Mono/Macrophage cells")
ggsave("final_umap_celltypes.pdf", width = 6.5, height = 5)
UMAPPlot(schm2, group.by = "subtypes", cols = subtypes_cols, split.by = "group") + ggtitle(NULL)
ggsave("final_umap_celltypes_splitBy_group.pdf", width = 13, height = 4)
UMAPPlot(schm2, label = T) + scale_color_igv()
ggsave("final_umap_cluster.pdf", width = 5.5, height = 5)

aimgenes <- c(
  "C1qa", "C1qb", 
  "Trem2", "Apoe", "Pf4", "Spp1", "Arg1", 
  "Cxcl3", "Cxcl1",
  "H2-Eb1", "H2-Aa", "Cd74", "Il1b", "Tnf",
  "Isg15", "Cxcl9", "Cxcl10", "Ifitm2", "Ifitm3", "Ifitm6",
  "Malat1", "Neat1", "mt-Co1", "mt-Co3"
)
VlnPlot(
    schm2, aimgenes, stack = T, flip = T, 
    group.by = "subtypes", cols = subtypes_cols, fill.by = "ident"
) + NoLegend()
ggsave("final_violin_features.pdf", width = 5, height = 10)

DotPlot(schm2, features = aimgenes, group.by = "subtypes") + 
  labs(x = NULL, y = NULL) + RotatedAxis() +
  scale_color_gradient2(low = "steelblue", high = "red")
ggsave("final_dotplot_markers.pdf", width = 10, height = 3.5)

## cell proportion or preference ----
cc <- rubasic::count_ratio2(schm2@meta.data, "sample2", "subtypes") %>%
  dplyr::rename(sample = sample2) %>%
  dplyr::mutate(group = factor(gsub("_[12]$", "", sample), levels(schm$group)))
ccplt <- rusinglecell::plot_cellcount(cc, "subtypes", cols = subtypes_cols, combine = FALSE, addlabel = T)
ccplt$ratio + theme_bw() + labs(y = "Cell proportion in total macrophages")
ggsave("final_cell_proportion_in_sample.pdf", width = 8, height = 4)

roe <- rusinglecell::odds_ratio(droplevels(schm2$subtypes), schm2$group)
rusinglecell::plot_roe(roe, "group") +
  theme(
    axis.title.y = element_text(family = "Arial", size = 14, face = "bold"),
    axis.text = element_text(family = "Arial", size = 12),
    legend.text = element_text(family = "Arial", size = 12),
    legend.title = element_text(family = "Arial", size = 14),
    plot.title = element_text(family = "Arial", size = 16, face = "bold", hjust = 0.5)
  )
ggsave("final_cell_proportion_in_group_preference.png", width = 6, height = 4)
