dir.create("/WUH721818ALE6L4/projects/TBK1/01_main_cluster", showWarnings = FALSE)
setwd("/WUH721818ALE6L4/projects/TBK1/01_main_cluster")


suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(ggsci)
  library(magrittr)
})


samples <- unlist(strsplit("vsv51_2 vsv51_1 GSK8612_2 GSK8612_1 Control_2 Control_1 Comb_2 Comb_1", " "))
samples <- setNames(nm = samples)


# quality control ----
load_data <- function(sname, fpath) {
  cnt <- Read10X(fpath, gene.column = 1)
  obj <- CreateSeuratObject(cnt, sname)
  obj$sample <- sname
  obj <- PercentageFeatureSet(obj, pattern = "^[Mm][Tt]-", col.name = "pct.mito")
  obj <- PercentageFeatureSet(obj, pattern = "^[Rr][Pp][LlSs]", col.name = "pct.ribo")
  obj
}

basic_flow <- function(obj) {
  NormalizeData(obj) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 0.8) %>%
    RunUMAP(dims = 1:30)
}

save_matrix <- function(sces, outpath) {
    dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
    fpaths <- file.path(outpath, c("matrix.mtx", "barcodes.tsv", "features.tsv"))
    Matrix::writeMM(sces@assays$RNA@count, file = fpaths[1])
    writeLines(colnames(sces), con = fpaths[2])
    data.frame(
        gene = rownames(sces),
        symbol = rownames(sces),
        biotype = "Gene Expression"
    ) %>%
        write.table(
            file = fpaths[3], sep = "\t", row.names = FALSE, 
            col.names = FALSE, quote = FALSE
        )
    system(sprintf("gzip %s/*", outpath))
}

scelist <- lapply(samples, function(sname) {
    load_data(sname, sprintf("../00_DNBC4tools/%s/02.count/filter_matrix", sname)) %>%
        basic_flow() %>%
        rusinglecell::run_doubletfinder()
})

sces <- lapply(
  scelist,
  subset,
  subset = nFeature_RNA >= 500 &
    nFeature_RNA <= 6000 &
    nCount_RNA <= 30000 &
    doublets == "Singlet" &
    pct.mito <= 10
)

sces <- lapply(
    sces,
    rusinglecell::clean_seurat,
    metas = colnames(sces[[1]]@meta.data)[1:6]
)
sces <- merge(sces[[1]], sces[-1])
save_matrix(sces, "matrix_after_QC")
write.table(sces@meta.data[, 1:6], file = "cellinfo_main_raw.tsv", sep = "\t", quote = FALSE)


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
schm <- SCTransform(sces, vars.to.regress = vars2regress) %>%
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
saveRDS(schm, file = "main_harmony.rds")


## batch correction: multi-cca ----
sccc <- SplitObject(sces, "sample") %>%
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
saveRDS(sccc, file = "main_SeuratCCA.rds")


## batch correction: bbknn ----
bbknn_args = c(
    "python3",
    "/WUH721818ALE6L4/projects/TBK1/run_bbknn.py",
    "/WUH721818ALE6L4/projects/TBK1/01_main_cluster/cellinfo_main_raw.tsv",
    "/WUH721818ALE6L4/projects/TBK1/01_main_cluster/main"
)
system(paste0(bbknn_args, collapse = " "))


# check clusters ----
GrayRed <- c("lightgray", "red")
DefaultAssay(schm) <- "RNA"
schm <- NormalizeData(schm)
DefaultAssay(sccc) <- "RNA"
sccc <- NormalizeData(sccc)
bkn <- read.csv("/WUH721818ALE6L4/projects/TBK1/01_main_cluster/main_bbknn_cellinfo.csv")
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


ctmain <- list(
  "T cell" = c("Cd3d", "Cd3e", "Cd8a", "Cd8b1", "Cd4", "Foxp3"),
  "NK cell" = c("Ncr1", "Prf1", "Klrb1c"),
  "B cell" = c("Ms4a1", "Cd79a", "Cd79b"),
  "pDC" = c("Ighm", "Iglc3", "Siglech", "Irf8"),
  "DC" = c("Fscn1", "Ccr7", "H2-Aa", "H2-Eb1", "Xcr1", "Clec9a", "Clec10a"),
  "Monocyte/Macrophage" = c("Aif1", "Lyz2", "Csf1r", "Arg1", "Cd68", "Cd86", "C1qa"),
  "Neutrophil" = c("G0s2", "Csf3r", "Fpr1", "S100a9", "S100a8"),
  "Basophil" = c("Cd200r3", "Fcer1a", "Cd69"),
  "Erythrocytes" = c("Hba-a1", "Hba-a2", "Hbb-bs"),
  "Fibroblast" = c("Sparc", "Cald1", "Col3a1", "Fbln2")
)

for (gene in unlist(ctmain)) {
    FeaturePlot(schm, gene, cols = GrayRed, reduction = "umap")
    ggsave(sprintf("feat_harmony_%s.png", gene), width = 4, height = 4)

    FeaturePlot(schm, gene, cols = GrayRed, reduction = "bkn50umap")
    ggsave(sprintf("feat_bbknn_%s.png", gene), width = 4, height = 4)

    FeaturePlot(sccc, gene, cols = GrayRed, reduction = "umap")
    ggsave(sprintf("feat_SeuratCCA_%s.png", gene), width = 4, height = 4)
}

## check the features, cluster markers in each clusters mannually ----
View(mks_harmony[[1]])
View(mks_bbknn[[1]])
View(mks_seuratcca[[1]])

## final cluster with BBKNN ----
# select the best cluster with minor resolution to distinguish the cell types.
schm$finalcluster <- as.numeric(as.character(schm$bkn_leiden_0.6))
# split the migratory DC from the main DC
schm$finalcluster[schm$bkn_leiden_2 %in% "20"] <- max(schm$finalcluster) + 1
# split the Treg from the main CD4T
schm$finalcluster[schm$bkn_leiden_2 %in% "23"] <- max(schm$finalcluster) + 1
# split the CD3+CSF3R+ from the CSF3R+ cells
schm$finalcluster[schm$bkn_leiden_2 %in% "25"] <- max(schm$finalcluster) + 1
# resort the final clusters and force a new order with cell type
schm$finalcluster <- rusinglecell::resortcluster(as.character(schm$finalcluster))
schm$finalclusterOrd <- factor(
  schm$finalcluster,
  as.character(c(0,4,7,13,10,14,6,16,15,12,11,1,2,3,9,8,17,18,5))
)

schm[["umap"]] <- Embeddings(schm, "bkn50umap") %>%
  magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>%
  CreateDimReducObject()

celltypes <- setNames(nm = as.character(0:18))
celltypes[as.character(c(0,4,7,13,10))] <- "T cell"
celltypes[as.character(c(6))] <- "NK cell"
celltypes[as.character(c(16))] <- "B cell"
celltypes[as.character(c(15))] <- "pDC"
celltypes[as.character(c(11,12))] <- "DC"
celltypes[as.character(c(1,2,3,9))] <- "Mono_Macro"
celltypes[as.character(c(8))] <- "Neutrophil"
celltypes[as.character(c(17))] <- "Basophil"
celltypes[as.character(c(18))] <- "Erythrocytes"
celltypes[as.character(c(5))] <- "Fibroblast"
celltypes[as.character(c(14))] <- "Doublets"
celltypes_ords <- c(
  "T cell", "NK cell", "B cell", "pDC", "DC", "Mono_Macro",
  "Neutrophil", "Basophil", "Erythrocytes", "Fibroblast",
  "Doublets"
)
celltypes_cols <- setNames(
  BuenColors::jdb_palette("corona", 11) %>% as.character() %>%
    setdiff("#7f7f7f") %>% c("lightgray"),
  nm = celltypes_ords
)
schm$celltypes <- celltypes[as.character(schm$finalcluster)]
schm$celltypes <- factor(schm$celltypes, celltypes_ords)

schm$sample2 <- plyr::mapvalues(
  schm$sample,
  c("GSK8612_1", "GSK8612_2", "vsv51_1", "vsv51_2"),
  c("GSK_1", "GSK_2", "VSV51_1", "VSV51_2")
) %>% factor(c("Control_1", "Control_2", "GSK_1", "GSK_2", "VSV51_1", "VSV51_2", "Comb_1", "Comb_2"))
schm$group <- factor(gsub("_[12]$", "", schm$sample2), c("Control", "GSK", "VSV51", "Comb"))

# save the cellinfo ----
obj2save <- list()
obj2save$metas <- sces@meta.data[, c("finalcluster", "finalclusterOrd", "celltypes", "sample2", "group")]
obj2save$celltypes_cols <- celltypes_cols
obj2save$celtypes_ord <- celltypes_ords
obj2save$umap <- sces@reductions$umap
saveRDS(obj2save, file = "main_final_cellinfo.rds")


# Figures ----
## clusters and features ----
UMAPPlot(schm, group.by = "finalcluster", label = T, repel = T) + 
  scale_color_igv() + ggtitle(NULL)
ggsave("umap_cluster.pdf", width = 6, height = 5)

DotPlot(schm, features = ctmain, group.by = "finalclusterOrd", cols = c("navy", "red")) +
  RotatedAxis() + scale_color_gradient2(low = "navy", high = "red") +
  theme(
    plot.background = element_rect(fill = "white"), 
    strip.text = element_text(size = 8, face = "bold")
  )
ggsave(
  "final_dotplot_cellmarker_in_clusterOrdered.pdf", 
  width = 15, height = 5
)

DotPlot(schm, features = ctmain, group.by = "finalcluster", cols = c("navy", "red")) +
  RotatedAxis() + scale_color_gradient2(low = "navy", high = "red") +
  theme(
    plot.background = element_rect(fill = "white"), 
    strip.text = element_text(size = 8, face = "bold")
  )
ggsave(
  "final_dotplot_cellmarker_in_cluster.pdf", 
  width = 15, height = 5
)

UMAPPlot(schm, group.by = "celltypes", cols = celltypes_cols) + ggtitle(NULL)
ggsave("final_umap_celltype.pdf", width = 6.5, height = 5)

DotPlot(schm, group.by = "celltypes", features = ctmain2) +
  RotatedAxis() + scale_color_gradient2(low = "navy", high = "red") +
  theme(
    plot.background = element_rect(fill = "white"), 
    strip.text = element_text(size = 8, face = "bold")
  )
ggsave(
  "final_dotplot_cellmarker_in_celltype.pdf",
  width = 15, height = 4
)

UMAPPlot(schm, group.by = "celltypes", cols = celltypes_cols, split.by = "group") + ggtitle(NULL)
ggsave("final_umap_celltype_splitBy_group.pdf", width = 13, height = 3)

# cell proportion ----
cc <- schm@meta.data %>%
  dplyr::filter(!celltypes %in% c("Fibroblast", "Doublets")) %>%
  rubasic::count_ratio2("sample", "celltypes")
write.csv(cc, file = "final_cell_proportion.csv")
ccplt <- rusinglecell::plot_cellcount(
  cc, "celltypes", cols = celltypes_cols,
  addlabel = TRUE,
  tolabel = c("T cell", "NK cell", "DC", "Mono_Macro", "Neutrophil"),
  combine = FALSE
)
ccplt$ratio + theme_bw() + labs(y = "Cell proportion in immune cells") +
  RotatedAxis() + theme(axis.text = element_text(color = "black"))
ggsave("final_cell_proportion.pdf", width = 8, height = 6)
