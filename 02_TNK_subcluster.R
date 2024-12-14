dir.create("/WUH721818ALE6L4/projects/TBK1/02_TNK_subcluster", showWarnings = FALSE)
setwd("/WUH721818ALE6L4/projects/TBK1/02_TNK_subcluster")


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

tnk <- subset(sces, subset = celltypes %in% c("T cell", "NK cell"))
meta2keep <- colnames(sces@meta.data)[1:6]
meta2keep <- c(meta2keep, "sample2", "group")
tnk <- rusinglecell::clean_seurat(tnk, metas = meta2keep)
write.table(tnk@meta.data[, 1:6], file = "cellinfo_TNK_raw.tsv", sep = "\t", quote = FALSE)

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
schm <- SCTransform(tnk, vars.to.regress = vars2regress) %>%
    rusinglecell::check_vargenes(
        dtypes = c("mito", "ribo", "igv", "igc", "trabv", "heat"),
        remove = TRUE
    ) %>%
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
saveRDS(schm, file = "TNK_harmony.rds")


## batch correction: multi-cca ----
sccc <- SplitObject(tnk, "sample") %>%
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
saveRDS(sccc, file = "TNK_SeuratCCA.rds")


## batch correction: bbknn ----
bbknn_args = c(
    "python3",
    "/WUH721818ALE6L4/projects/TBK1/run_bbknn.py",
    "/WUH721818ALE6L4/projects/TBK1/02_TNK_subcluster/cellinfo_TNK_raw.tsv",
    "/WUH721818ALE6L4/projects/TBK1/02_TNK_subcluster/TNK"
)
system(paste0(bbknn_args, collapse = " "))


# check clusters ----
GrayRed <- c("lightgray", "red")
DefaultAssay(schm) <- "RNA"
schm <- NormalizeData(schm)
DefaultAssay(sccc) <- "RNA"
sccc <- NormalizeData(sccc)
bkn <- read.csv("/WUH721818ALE6L4/projects/TBK1/02_TNK_subcluster/TNK_bbknn_cellinfo.csv")
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

# c("Cd3d", "Ccr7", "Cd8a", "Cd4", "Foxp3", "Mki67", "Ncr1", "Klrb1c", "mt-Co1") %>%
# c("Gzma", "Ccl5", "Gzmb", "Gzmk") %>%
# c("Ifitm1", "Ifitm2", "Lgals3", "S100a4") %>%
# c("Ifitm3", "Ifi27l2a", "Il7r", "Ccr7") %>%
# paste0("Gzm", letters[1:10]) %>%
# c("Nkg7", "Prf1") %>%
# c("Pdcd1", "Ctla4", "Tigit", "Havcr2", "Lag3") %>%
# c("Isg15", "mt-Co1", "mt-Co2", "Malat1") %>%
# c("Ifng", "Ifitm1", "Ifitm2", "Ifitm3", "Isg15", "Ifi27l2a", "Cxcl9", "Cxcl10") %>%
c("Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Xcl1") %>%
  FeaturePlot(schm, ., reduction = "bkn50umap", cols = GrayRed)
  # FeaturePlot(schm, ., reduction = "hm19umap", cols = GrayRed)
  # FeaturePlot(sccc, ., reduction = "cca12umap", cols = GrayRed)


## final cluster with BBKNN ----
schm$finalcluster <- as.numeric(as.character(schm$bkn_leiden_1.2))
# split the NKT from T cells
bkn_subcluster_arg <- c(
    "python3",
    "/WUH721818ALE6L4/projects/TBK1/run_bbknn_subcluster.py",
    "TNK"
)
system(paste0(bkn_subcluster_arg, collapse = " "))
bknRes12sub4 <- read.csv("./TNK_bbknn_cellinfo_res12sub04.csv")
schm$bkn_leiden_1.2_sub4 <- bknRes12sub4$bkn_leiden_1.2_sub4
schm$finalcluster[schm$bkn_leiden_1.2_sub4 == "4,1"] <- max(schm$finalcluster) + 1
# merge cluster 8 with cluster 0 due to multiple common cluster markers
schm$finalcluster[schm$bkn_leiden_1.2 == "8"] <- 0
schm$finalcluster <- rusinglecell::resortcluster(as.character(schm$finalcluster))


subtypes <- setNames(nm = as.character(0:14))
subtypes["1"] <- "CD8T_C1_Ccl5"
subtypes["0"] <- "CD8T_C2_Gzma"
subtypes["2"] <- "CD8T_C3_Ifitm1"
subtypes["7"] <- "CD8T_C4_Gzmf"
subtypes["8"] <- "CD8T_C5_Ifng"
subtypes["14"] <- "CD8T_C6_Isg15"
subtypes["10"] <- "CD8T_C7_Mki67"
subtypes["9"] <- "CD8T_C8_Ccr7"
subtypes["4"] <- "CD8T_C9_Mito"
subtypes["11"] <- "CD8T_C9_Mito"
subtypes["13"] <- "CD4T_C1_Treg"
subtypes["3"] <- "CD4T_C2_Ifitm2"
subtypes["5"] <- "NK_C1_Ccl5"
subtypes["6"] <- "NK_C2_Gzmc"
subtypes["12"] <- "NKT"

tnk_cols <- setNames(
  c(
    "steelblue", "#474FB7",
    "#0F3919", "#5B7739", "#4DB02F", "#B6BD59", 
    "#BFC699", "#E7E7B0", "#88A088", "#B1BFB1", "#C3CD9E",
    "#A5311F", "#CD7F2E", "#D7B8C0"
  ),
  nm = unique(subtypes) %>% sort()
)

schm[["umap"]] <- Embeddings(schm, "bkn50umap") %>%
  magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>%
  CreateDimReducObject()
schm$subtypes <- subtypes[as.character(schm$finalcluster)]
schm$finalclusterOrd <- factor(schm$finalcluster, as.character(sort(subtypes) %>% names()))


# save the cellinfo ----
obj2save <- list()
obj2save$metas <- schm@meta.data[, c("finalcluster", "finalclusterOrd", "subtypes", "sample2", "group")]
obj2save$celltypes_cols <- tnk_cols
obj2save$celtypes_ord <- levels(schm$finalclusterOrd)
obj2save$umap <- schm@reductions$umap
saveRDS(obj2save, file = "TNK_final_cellinfo.rds")


# Figures ----
## clusters and features ----
UMAPPlot(schm, group.by = "subtypes", cols = tnk_cols) + ggtitle("T/NK cells")
ggsave("final_umap_celltypes.pdf", width = 6.5, height = 5)
UMAPPlot(schm, group.by = "subtypes", cols = tnk_cols, split.by = "group") + ggtitle(NULL)
ggsave("final_umap_celltypes_splitBy_group.pdf", width = 13, height = 4)
UMAPPlot(schm, label = T) + scale_color_igv()
ggsave("final_umap_cluster.pdf", width = 5.5, height = 5)

# c("Cd3d", "Cd3g", "Cd8a", "Cd8b1", "Cd4", "Foxp3", "Ncr1", "Klrb1c") %>%
  c("Nkg7", "Ccl5", "Ifitm2", "Ifitm1", "Ifng", "Isg15", "Ifit3") %>%
  c("Gzma", "Gzmb", "Gzmc", "Gzmd", "Gzme", "Gzmf") %>%
  c("Mki67", "Top2a", "mt-Co1", "mt-Co3", "Ccl3", "Ccl4") %>%
  VlnPlot(
    schm, features = ., group.by = "subtypes", pt.size = 0,
    flip = T, stack = T, cols = tnk_cols, fill.by = "ident"
  ) + NoLegend()
ggsave("final_violin_features.pdf", width = 8, height = 12)

c("Nkg7", "Ccl5", "Ifitm2", "Ifitm1", "Ifng", "Isg15", "Ifit3") %>%
  c("Gzma", "Gzmb", "Gzmc", "Gzmd", "Gzme", "Gzmf") %>%
  c("Mki67", "Top2a", "mt-Co1", "mt-Co3", "Ccl3", "Ccl4", "Xcl1") %>%
  DotPlot(schm, features = ., group.by = "subtypes", cols = GrayRed) + 
  ThemeDot() + scale_color_gradient2(low = "navy", high = "red")
ggsave("final_dotplot_subtype_marker.pdf", width = 9, height = 6)


## cell proportion or preference ----
cc <- rubasic::count_ratio2(schm@meta.data, "sample2", "subtypes") %>%
  dplyr::rename(sample = sample2) %>%
  dplyr::mutate(
    group = factor(gsub("_[12]$", "", sample), levels(schm$group)),
    sample = factor(sample, levels(schm$sample2))
  )
write.csv(cc, "final_cell_proportion.csv")
ggplot(cc, aes(x = subtypes, y = ratio, color = group)) +
  geom_boxplot() + 
  geom_point(position=position_dodge(width = 0.75),aes(group = group)) + 
  theme_bw() + RotatedAxis() +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black")) +
  scale_color_aaas() + scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "Cell proportion in total T or NK cells", color = NULL) +
  geom_vline(xintercept = 0.5 + 1:13, color = "lightgray", linetype = "dashed")
ggsave("final_cell_proport_in_TNK.pdf", width = 7, height = 3.5)

roe_group <- rusinglecell::odds_ratio(schm$subtypes, schm$group)
roe_group$roe <- roe_group$roe[, levels(schm$group)]
rusinglecell::plot_roe(roe_group, "group") +
  theme(
    axis.title.y = element_text(family = "Arial", size = 14, face = "bold"),
    axis.text = element_text(family = "Arial", size = 12),
    legend.text = element_text(family = "Arial", size = 12),
    legend.title = element_text(family = "Arial", size = 14),
    plot.title = element_text(family = "Arial", size = 16, face = "bold", hjust = 0.5)
  )
ggsave("final_cell_proport_in_TNK_preference.pdf", width = 6, height = 6)


## exhaustion signature ----
extgenes <- c("Pdcd1", "Tigit", "Lag3", "Havcr2", "Ctla4")
schm <- AddModuleScore(schm, list(ext = extgenes))
schm$Exhaust <- schm$Cluster1
aimcls <- c(
    "CD8T_C1_Ccl5", "CD8T_C2_Gzma", "CD8T_C4_Gzmf", 
    "CD8T_C5_Ifng", "NK_C1_Ccl5", "NK_C2_Gzmc", "NKT"
)
tmpdf <- FetchData(tnk, c("subtypes", "Exhaust", "group", extgenes)) %>%
  tidyr::pivot_longer(cols = c(extgenes, "Exhaust"), values_to = "score", names_to = "gene") %>%
  dplyr::mutate(
    gene = ifelse(gene == "Exhaust", "Exhaustion score", gene),
    gene = factor(gene, c(extgenes, "Exhaustion score"))
  ) %>%
  dplyr::filter(subtypes %in% aimcls)

tmpdf %>%
  ggplot(aes(x = subtypes, y = score, fill = subtypes)) +
  geom_violin(scale = "width", draw_quantiles = c(0.25,0.5,0.75)) +
  theme_bw() + RotatedAxis() + 
  scale_fill_manual(values = tnk_cols) + NoLegend() +
  facet_wrap(~ gene, ncol = 3, scales = "free_y") + theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = NA, color = NA)
  ) + labs(x = NULL, y = "Expression level or Module score")
ggsave("final_exhaust_score_in_aim_cluster.pdf", width = 9, height = 5)


## Comb vs others in total NK ----
mk_combVsOther_nk <- FindMarkers(
  subset(schm, subset = subtypes %in% c("NK_C1_Ccl5", "NK_C2_Gzmc")),
  ident.1 = "Comb",
  logfc.threshold = 0, min.pct = 0, group.by = "group"
)
mk_combVsOther_nk$gene <- rownames(mk_combVsOther_nk)
mk_combVsOther_nk$subtype <- "NK"
write.csv(mk_combVsOther_nk, file = "final_DEGs_Comb_vs_Other_in_NK.csv")
rubasic::volcanoplot(
  mk_combVsOther_nk, "avg_log2FC", 0.5, "p_val_adj",
  label.by = "fc", p.inflimit = TRUE
) + labs(x = "avg_log2FC", title = "Comb vs Others in NK")
ggsave("final_DEGs_Comb_vs_Other_in_NK-volcano.pdf", width = 6, height = 6)
