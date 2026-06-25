library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


obj_name <- load("~/Desktop/Loss_of_Y_analysis/DCM_public/GSE183852_DCM_Cells.Robj")
obj_name

table(HDCM$orig.ident)
table(HDCM$Condition)
table(HDCM$Sex)
table(HDCM$Names)

HDCM <- subset(HDCM, subset = Condition == c("DCM"))
HDCM <- subset(HDCM, subset = Sex == c("Male"))

HDCM1 <- subset(HDCM, subset = orig.ident == c("HDCM1"))
HDCM4 <- subset(HDCM, subset = orig.ident == c("HDCM4"))
HDCM6 <- subset(HDCM, subset = orig.ident == c("HDCM6"))

HDCM1 <- AddMetaData(HDCM1, "HDCM1", col.name = "sample")
HDCM4 <- AddMetaData(HDCM4, "HDCM4", col.name = "sample")
HDCM6 <- AddMetaData(HDCM6, "HDCM6", col.name = "sample")


HDCM1[["percent.mt"]] <- PercentageFeatureSet(HDCM1, pattern = "^MT-")
HDCM4[["percent.mt"]] <- PercentageFeatureSet(HDCM4, pattern = "^MT-")
HDCM6[["percent.mt"]] <- PercentageFeatureSet(HDCM6, pattern = "^MT-")

HDCM1 <- subset(HDCM1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
HDCM4 <- subset(HDCM4, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
HDCM6 <- subset(HDCM6, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

HDCM1 <- RenameCells(HDCM1, new.names = paste0(colnames(HDCM1), "_", "HDCM1"))
HDCM4 <- RenameCells(HDCM4, new.names = paste0(colnames(HDCM4), "_", "HDCM4"))
HDCM6 <- RenameCells(HDCM6, new.names = paste0(colnames(HDCM6), "_", "HDCM6"))

data.DCM1_LVN <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/N1/LVN")
data.DCM1_LVP <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/N1/LVP")
data.DCM1_RVN <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/N1/RVN")
data.DCM1_RVP <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/N1/RVP")
data.DCM2_LVN <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM2/LVN")
data.DCM2_LVP <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM2/LVP")
data.DCM2_RVN <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM2/RVN")
data.DCM2_RVP <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM2/RVP")
data.DCM3_LVN <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM3/LVN")
data.DCM3_LVP <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM3/LVP")
data.DCM3_RVN <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM3/RVN")
data.DCM3_RVP <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/DCM_public/GSE145154_RAW/DCM3/RVP")


CM.list <- list(
  DCM1_LVN  = data.DCM1_LVN,
  DCM1_LVP  = data.DCM1_LVP,
  DCM1_RVN  = data.DCM1_RVN,
  DCM1_RVP  = data.DCM1_RVP,
  DCM2_LVN  = data.DCM2_LVN,
  DCM2_LVP  = data.DCM2_LVP,
  DCM2_RVN  = data.DCM2_RVN,
  DCM2_RVP  = data.DCM2_RVP,
  DCM3_LVN  = data.DCM3_LVN,
  DCM3_LVP  = data.DCM3_LVP,
  DCM3_RVN  = data.DCM3_RVN,
  DCM3_RVP  = data.DCM3_RVP
)


process_sample <- function(counts, sample_name) {
  
  obj <- CreateSeuratObject(
    counts = counts, 
    project = "CM_project",
    min.cells = 3,
    min.features = 200
  )
  
  obj <- AddMetaData(obj, sample_name, col.name = "sample")
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
  
  obj <- RenameCells(obj, new.names = paste0(colnames(obj), "_", sample_name))
  
  return(obj)
}


for (s in names(CM.list)) {
  message("Processing: ", s)
  CM.list[[s]] <- process_sample(CM.list[[s]], s)
}

DCM.list <- list(
  HDCM1     = HDCM1,
  HDCM4     = HDCM4,
  HDCM6     = HDCM6
)

CM.list <- c(CM.list, DCM.list)


length(CM.list) 

gene_list <- read.csv("~/Desktop/yoshida_test/Loss_of_Y/MSY_gene.txt")
class(gene_list)
gene_list <- gene_list$Gene.name

gene_list <- c("RPS4Y1","ZFY","USP9Y","DDX3Y","KDM5D","EIF1AY")

fetch_counts_with_zeros <- function(seurat_obj, gene_list, slot = "counts") {
  
  mat <- GetAssayData(seurat_obj, slot = slot)
  
  out <- matrix(0, nrow = ncol(mat), ncol = length(gene_list))
  colnames(out) <- gene_list
  rownames(out) <- colnames(mat)
  
  present <- intersect(gene_list, rownames(mat))
  
  if (length(present) > 0) {
    out[, present] <- t(as.matrix(mat[present, , drop = FALSE]))
  } else {
    warning("指定した gene_list がどれもこのオブジェクトに存在しません。すべて 0 とみなされます。")
  }
  
  return(as.data.frame(out))
}

judge_LoY <- function(obj_list, male_samples, gene_list, min_features = 500) {
  
  expr_list <- lapply(male_samples, function(s) {
    obj <- obj_list[[s]]
    
    df  <- fetch_counts_with_zeros(obj, gene_list, slot = "counts")

    meta <- obj@meta.data[colnames(obj), , drop = FALSE]
    df$nFeature_RNA <- meta$nFeature_RNA
    df$sample <- s
    
    return(df)
  })
  
  # cell × (gene_list + nFeature_RNA + sample)
  expr_all <- do.call(rbind, expr_list)
  
  expr_filt <- expr_all[expr_all$nFeature_RNA >= min_features, , drop = FALSE]
  
  if (nrow(expr_filt) == 0) {
    warning("QC (nFeature_RNA >= min_features) を満たす男性細胞が 1 つもありません。")
    return(character(0))
  }
  
  is_LoY <- rowSums(expr_filt[, gene_list, drop = FALSE] > 0) == 0

  LoY_cells <- rownames(expr_filt)[is_LoY]
  
  return(LoY_cells)
}

sample_names <- names(CM.list)

CM_LoY <- judge_LoY(CM.list, names(CM.list), gene_list,min_features = 500)

LoY_cells <- c(CM_LoY)

add_LoY_metadata <- function(AF.list, LoY_cells) {
  
  for (s in names(AF.list)) {
    obj <- AF.list[[s]]
    
    obj[["LoY"]] <- "no_LoY"
    
    common <- intersect(colnames(obj), LoY_cells)
    if (length(common) > 0) {
      obj$LoY[common] <- "LoY"
    }
    
    AF.list[[s]] <- obj
  }
  
  return(AF.list)
}

CM.list <- add_LoY_metadata(CM.list, LoY_cells)

combined <- Reduce(function(x, y) merge(x, y), CM.list)

combined <- NormalizeData(combined, verbose = FALSE)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)

combined <- ScaleData(combined, verbose = FALSE)

combined <- SCTransform(combined, vars.to.regress = c("percent.mt"))

combined <- RunPCA(combined, assay = "SCT", verbose = FALSE)

library(Rcpp)
library(harmony)

combined <- RunHarmony(combined, group.by.vars = "sample", assay.use = "SCT",
                       reduction.use = "pca",plot_convergence = TRUE)


ElbowPlot(combined, ndims = 50)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:50)
combined <- FindClusters(combined, resolution = 0.05)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:50)

DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined, reduction = "umap", split.by = "LoY")
#DimPlot(combined, reduction = "umap", split.by = "sample")

table(combined$orig.ident)

saveRDS(combined, "~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_sample_all_cell_DCM.rds")
combined <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_sample_all_cell.rds")

##Heatmap
combined_subset <- subset(combined, cells = Cells(combined[["RNA"]]), downsample = 2500)
markers <- FindAllMarkers(combined_subset, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)

combined_subset <- ScaleData(combined_subset, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_subset, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=88)


# Myeloid----
VlnPlot(combined, assay = "RNA",features = c("CD68","CD14","CSF3R","CSF1R","IL1B","S100A8", "S100A9","MKI67","FCGR3A") ,pt.size = 0)
VlnPlot(combined, assay = "RNA",features = c("CD163","CD14","PTPRC") ,pt.size = 0)

# Bcells-----
VlnPlot(combined, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") ,pt.size = 0)

VlnPlot(combined, features = c("CDH5","ACTA2","VWF","MKI67","MYH11","CD34","IGHG4") ,pt.size = 0)

#SMC
VlnPlot(combined, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DRB5","HLA-DPB1") ,pt.size = 0)
VlnPlot(combined, assay = "RNA", features = c("ACTA2","HAND2","ERBB4","SERTAD4","DLX5","SPARCL1","TAGLN","MYOCD") ,pt.size = 0)


# Tcells-----
VlnPlot(combined, features = c("CD3E","CD4","CD8A","IL7R","LEF1","GZMK","NKG7") ,pt.size = 0)

VlnPlot(combined, features = c("GLS"), split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))


#cDC1
VlnPlot(combined, pt.size = 0 , features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
#cDC2
VlnPlot(combined, pt.size = 0 , features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
#pDC
VlnPlot(combined, pt.size = 0 , features = c("IL3RA","CLEC4C") )

###Myeloid

combined_Mye <- subset(combined, idents = c("1"))
DefaultAssay(combined_Mye) <- "RNA"
combined_Mye <- NormalizeData(combined_Mye, verbose = FALSE)
combined_Mye <- FindVariableFeatures(combined_Mye, selection.method = "vst", nfeatures = 3000)
combined_Mye <- ScaleData(combined_Mye, verbose = FALSE)
#combined_Mye <- SCTransform(combined_Mye, vars.to.regress = c("percent.mt"))
#combined_Mye <- RunPCA(combined_Mye, assay = "SCT", verbose = FALSE)
combined_Mye <- RunPCA(combined_Mye, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_Mye <- RunHarmony(combined_Mye, group.by.vars = "sample", assay.use = "RNA",
                           reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_Mye, ndims = 50, reduction = "harmony")
combined_Mye <- FindNeighbors(combined_Mye, reduction = "harmony", dims = 1:20)
combined_Mye <- FindClusters(combined_Mye, resolution = 0.05)
combined_Mye <- RunUMAP(combined_Mye, reduction = "harmony", dims = 1:20)

DimPlot(combined_Mye, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_Mye, reduction = "umap", split.by = "LoY")

table(combined_Mye$LoY)

###Heatmap
markers <- FindAllMarkers(combined_Mye, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

combined_Mye_heat <- ScaleData(combined_Mye, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_Mye_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=80)


# Myeloid----
VlnPlot(combined_Mye, assay = "RNA",features = c("CD68","CD14","CSF3R","CSF1R","IL1B","S100A8", "S100A9","MKI67","FCGR3A") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA",features = c("CD163","CD14","PTPRC") ,pt.size = 0)

# Bcells-----
VlnPlot(combined_Mye, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") ,pt.size = 0)

VlnPlot(combined, features = c("CDH5","ACTA2","VWF","MKI67","MYH11","CD34","IGHG4") ,pt.size = 0)

#SMC
VlnPlot(combined_Mye, assay = "RNA", features = c("ACTA2","HAND2","ERBB4","SERTAD4","DLX5","SPARCL1","TAGLN","MYOCD") ,pt.size = 0)
# Tcells-----
VlnPlot(combined_Mye, features = c("CD3E","CD4","CD8A","IL7R","LEF1","GZMK","NKG7") ,pt.size = 0)

#cDC1
VlnPlot(combined_Mye, pt.size = 0 , features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
#cDC2
VlnPlot(combined_Mye, pt.size = 0 , features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
#pDC
VlnPlot(combined, pt.size = 0 , features = c("IL3RA","CLEC4C") )

##2:NK

combined_Mye2 <- subset(combined_Mye, idents = c("2"),invert=TRUE)
DefaultAssay(combined_Mye2) <- "RNA"
combined_Mye2 <- NormalizeData(combined_Mye2, verbose = FALSE)
combined_Mye2 <- FindVariableFeatures(combined_Mye2, selection.method = "vst", nfeatures = 3000)
combined_Mye2 <- ScaleData(combined_Mye2, verbose = FALSE)
#combined_Mye2 <- SCTransform(combined_Mye2, vars.to.regress = c("percent.mt"))
#combined_Mye2 <- RunPCA(combined_Mye2, assay = "SCT", verbose = FALSE)
combined_Mye2 <- RunPCA(combined_Mye2, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_Mye2 <- RunHarmony(combined_Mye2, group.by.vars = "sample", assay.use = "RNA",
                            reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_Mye2, ndims = 50, reduction = "harmony")
combined_Mye2 <- FindNeighbors(combined_Mye2, reduction = "harmony", dims = 1:30)
combined_Mye2 <- FindClusters(combined_Mye2, resolution = 0.1)
combined_Mye2 <- RunUMAP(combined_Mye2, reduction = "harmony", dims = 1:30)

DimPlot(combined_Mye2, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_Mye2, reduction = "umap", split.by = "LoY")

table(combined_Mye2$LoY)

###Heatmap
markers <- FindAllMarkers(combined_Mye2, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

combined_Mye2_heat <- ScaleData(combined_Mye2, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_Mye2_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=70)

##3:myocardial,4:vascular wall

combined_Mye3 <- subset(combined_Mye2, idents = c("3","4"),invert=TRUE)
DefaultAssay(combined_Mye3) <- "RNA"
combined_Mye3 <- NormalizeData(combined_Mye3, verbose = FALSE)
combined_Mye3 <- FindVariableFeatures(combined_Mye3, selection.method = "vst", nfeatures = 3000)
combined_Mye3 <- ScaleData(combined_Mye3, verbose = FALSE)
#combined_Mye3 <- SCTransform(combined_Mye3, vars.to.regress = c("percent.mt"))
#combined_Mye3 <- RunPCA(combined_Mye3, assay = "SCT", verbose = FALSE)
combined_Mye3 <- RunPCA(combined_Mye3, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_Mye3 <- RunHarmony(combined_Mye3, group.by.vars = "sample", assay.use = "RNA",
                            reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_Mye3, ndims = 50, reduction = "harmony")
combined_Mye3 <- FindNeighbors(combined_Mye3, reduction = "harmony", dims = 1:40)
combined_Mye3 <- FindClusters(combined_Mye3, resolution = 0.3)
combined_Mye3 <- RunUMAP(combined_Mye3, reduction = "harmony", dims = 1:40)

combined_Mye3 <- subset(combined_Mye3, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)

DimPlot(combined_Mye3, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_Mye3, reduction = "umap", split.by = "LoY")

table(combined_Mye3$LoY)
table(combined_Mye3$sample)

###Heatmap
markers <- FindAllMarkers(combined_Mye3, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

combined_Mye3_heat <- ScaleData(combined_Mye3, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_Mye3_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=70)

combined_Mye3 <- AddMetaData(combined_Mye3, Idents(combined_Mye3),col.name = "cell_annotation_final")

saveRDS(combined_Mye3, "~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_DCM_sample_Myeloid_cell_gene3000_260414.rds")
combined_Mye3 <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_DCM_sample_Myeloid_cell_gene3000_260414.rds")

VlnPlot(merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

merge <- merge(combined_Mye,combined_Mye3)

Idents(combined_Mye) <- "LoY"

#Dotplot
combined_Mye3_trans <- SCTransform (combined_Mye3)

cd_genes <- c("CD14","CD68","CSF3R","FCGR3B","CSF1R","FCGR3A","HLA-DQA1","HLA-DQB1","CLEC9A","CLEC10A","S100A8","S100A12","VCAN","LYZ","IL1B","AREG","EREG","HBEGF","LYVE1","FOLR2","MRC1","TNF","CCL3","CCL4","TREM2","SPP1","APOE","MKI67","TOP2A")
DotPlot(combined_Mye3_trans,features = cd_genes)+RotatedAxis()+coord_flip()



##annotation
#0:TNF+ inflammtory mac
#1:LYVE1+ mac
#2:IL1B+ inflammatory mono 
#3:cDC2
#4:TREM2+ mac
#5:noncalssical mono
#6:proliferative mac
#7:cDC1

##IL1B
VlnPlot(combined_Mye3, assay = "RNA",features = c("IL1B","EREG","AREG","HBEGF","CXCL2","VEGFA","CD44","IL1RN","CXCL3") ,pt.size = 0)
VlnPlot(combined_Mye3, assay = "RNA",features = c("TNF","CXCL2") ,pt.size = 0)

VlnPlot(combined_Mye3, assay = "RNA",features = c("FN1","SPP1") ,pt.size = 0,split.by = "LoY")

#TREM2+
VlnPlot(combined_Mye3, assay = "RNA",features = c("TREM2","SPP1","LPL","APOE","ABCG1") ,pt.size = 0)

#LYVE1+
VlnPlot(combined_Mye3, assay = "RNA",features = c("LYVE1","FOLR2","IGF1","C1QA") ,pt.size = 0)
##proliferative
VlnPlot(combined_Mye3, assay = "RNA", features = c("SPC25","MKI67","TOP2A","BIRC5","FABP5") ,pt.size = 0)


FeaturePlot(combined_Mye3, features = "TGFB1")
VlnPlot(combined_Mye3, assay = "RNA",features = c("TGFB1","FN1") ,pt.size = 0)

combined_Mye3 <- RenameIdents(
  combined_Mye3,
  `0` = "TNF+ inflammatory macs",
  `1` = "LYVE1+ macs",
  `2` = "IL1B+ inflammatory monos / macs",
  `3` = "cDC2",
  `4` = "TREM2+ macs",
  `5` = "Non classical monos",
  `6` = "proliferative macs",
  `7` = "cDC1"
)

Idents(combined_Mye3) <- factor(Idents(combined_Mye3), levels = c("Non classical monos", "cDC1", "cDC2", "IL1B+ inflammatory monos / macs","LYVE1+ macs","TNF+ inflammatory macs","TREM2+ macs","proliferative macs"))

#Idents(combined_Mye3) <- "cell_annotation_final"

combined <- subset(combined, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)
DimPlot(combined)


tab <- table(
  cluster = Idents(combined_Mye3),   # seu$seurat_clusters
  state   = combined_Mye3$LoY      # "A" / "B"
)

tab <- table(
  cluster = combined_Mye3$sample,   # seu$seurat_clusters
  state   = combined_Mye3$LoY      # "A" / "B"
)

tab

prop_tab <- prop.table(tab, margin = 1)

prop_A <- prop_tab[, "LoY"]

prop_A

prop_df <- data.frame(
  cluster   = rownames(prop_tab),
  frac_A    = prop_tab[, "LoY"],
  n_total   = rowSums(tab),
  n_A       = tab[, "LoY"],
  n_B       = tab[, "no_LoY"],
  percent   = tab[, "LoY"] / rowSums(tab) * 100
)

prop_df

#cluster     frac_A n_total  n_A  n_B   percent
#0       0 0.14328771    7963 1141 6822 14.328771
#1       1 0.28229568    4931 1392 3539 28.229568
#2       2 0.19228296    3110  598 2512 19.228296
#3       3 0.11386139    2828  322 2506 11.386139
#4       4 0.13750000    2400  330 2070 13.750000
#5       5 0.14467515    1493  216 1277 14.467515
#6       6 0.09236948     249   23  226  9.236948
#7       7 0.13970588     136   19  117 13.970588


#####Myeloid volcanoplot 
test3 <- subset(combined_Mye3, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)

DimPlot(test3, reduction = "umap")
table(combined_Mye3$sample)
VlnPlot(test3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


test <- subset(test3, idents = "IL1B+ inflammatory monos / macs")
test <- subset(combined_Mye3, idents = "TREM2+ macs")
test <- subset(combined_Mye3, idents = "Neutrophils")
test <- subset(combined_Mye3, idents = "LYVE1+ macs")
test <- test3
test <- subset(combined_Mye3, idents = c("Neutrophils","cDC1","cDC2","Non classical monos","classical monos"),invert=TRUE)
test <- subset(combined_Mye3, idents = c("IL1B+ inflammatory mono / mac","classical mono","TNF+ inflammatory mac"))
test <- subset(combined_Mye3, idents = c("Non classical monos"))
test <- subset(combined_Mye3, idents = c("LYVE1+ macs","TREM2+ macs","IL1B+ inflammatory monos / macs","TNF+ inflammatory macs"))
test <- subset(combined_Mye3, idents = c("LYVE1+ macs","TREM2+ macs"))
test <- subset(combined_Mye3, idents = c("TNF+ inflammatory macs"))

test <- subset(combined_Mye3, subset = sample %in% c("Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi"))

DimPlot(test, reduction = "umap")

DefaultAssay(test) <- "RNA"
table(test$LoY)

test_LoY <- subset(test, subset = LoY == "LoY")
test_no_LoY <- subset(test, subset = LoY == "no_LoY")

VlnPlot(test_LoY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(test_no_LoY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(combined_Mye, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


table(test$LoY)

genes_to_keep <- setdiff(rownames(test), c(gene_list))

#subset
test[["RNA"]]@counts     <- test[["RNA"]]@counts[genes_to_keep, ]
test[["RNA"]]@data       <- test[["RNA"]]@data[genes_to_keep, ]
#test[["RNA"]]@scale.data <- test[["RNA"]]@scale.data[genes_to_keep, ]


Idents(test) <- "LoY"

Aur.table_chipvsnon <- FindMarkers(test, ident.1 = c("LoY"), ident.2 =c("no_LoY"), verbose = FALSE, logfc.threshold = 0)

Aur.table_chipvsnon$logp <- -log10(Aur.table_chipvsnon$p_val)

Aur.table_chipvsnon_filtered_left = subset(Aur.table_chipvsnon, (logp>=1 & avg_log2FC <= -0.34) | (logp>=5 & avg_log2FC <= -0.25))
Aur.table_chipvsnon_filtered_right = subset(Aur.table_chipvsnon, (logp>=4 & avg_log2FC >= 0.15) | (logp>=2 & avg_log2FC >= 0.35))

Aur.table_chipvsnon_filtered_left = subset(Aur.table_chipvsnon, (logp>=1 & avg_log2FC <= -0.7) | (logp>=2.6 & avg_log2FC <= -0.1))
Aur.table_chipvsnon_filtered_right = subset(Aur.table_chipvsnon, (logp>=1.0 & avg_log2FC >= 0.15) | (logp>=1 & avg_log2FC >= 0.15))

genes.to.label.left <- rownames(Aur.table_chipvsnon_filtered_left)
genes.to.label.right <- rownames(Aur.table_chipvsnon_filtered_right)
genes.to.label.right <- c(genes.to.label.right,"CXCR2","S100A12","IL17RA","TGFBR2")
genes.to.label.right <- c(genes.to.label.right,"TGFB1","FN1")

#genes.to.label.left <- "IL1B"
genes.to.label.right <- c(genes.to.label.right,"S100A8","S100A9","S100A6","S100A12","CXCL3","CXCL1")
genes.to.label.left <- c(genes.to.label.left,"CXCL3","IL1B","CXCL8")
genes.to.label.left <- c(genes.to.label.left,"CXCL3")

genes.to.label.right <- c(genes.to.label.right,"RETN","MS4A7","CXCL3","CXCL1")
genes.to.label.left <- c(genes.to.label.left,"CXCL8")


p1 <- ggplot(Aur.table_chipvsnon, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

Aur.table_chipvsnon$gene <- rownames(Aur.table_chipvsnon)

library(ggplot2)
library(ggrepel)

p1 <- ggplot(Aur.table_chipvsnon,
             aes(x = avg_log2FC, y = logp)) +
  geom_point()

p1 <- p1 +
  geom_label_repel(
    data = subset(Aur.table_chipvsnon, gene %in% genes.to.label.right),
    aes(label = gene),
    color = "red",
    size = 5.5,
    fill  = "white",  
    box.padding   = 0.3,
    point.padding = 0.2,
    label.size    = 0
  ) +
  geom_label_repel(
    data = subset(Aur.table_chipvsnon, gene %in% genes.to.label.left),
    aes(label = gene),
    color = "blue",
    size = 5.5,
    fill  = "white",
    box.padding   = 0.3,
    point.padding = 0.2,
    label.size    = 0
  )
p1



##GO
type0.markers <- Aur.table_chipvsnon[Aur.table_chipvsnon$p_val < 0.05&Aur.table_chipvsnon$avg_log2FC > 0.2,]
type1.markers <- Aur.table_chipvsnon[Aur.table_chipvsnon$p_val < 0.05&Aur.table_chipvsnon$avg_log2FC < -0.2,]
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

library(clusterProfiler)
type0.gene_SYMBOLs <- rownames(type0.markers)
type0.gene_IDs <- AnnotationDbi::select(hs, keys=type0.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID

type1.gene_SYMBOLs <- rownames(type1.markers)
type1.gene_IDs <- AnnotationDbi::select(hs, keys=type1.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
# compare cluster
type0.gene_SYMBOLs <- rownames(type0.markers)
type0.gene_IDs <- AnnotationDbi::select(hs, keys=type0.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
type1.gene_SYMBOLs <- rownames(type1.markers)
type1.gene_IDs <- AnnotationDbi::select(hs, keys=type1.gene_SYMBOLs, columns = c("ENTREZID", "SYMBOL"), keytype="SYMBOL")$ENTREZID
genelist <- list(type0.gene_IDs, type1.gene_IDs)
names(genelist) <- c("type0", "type1")

# GO:BP
cgBP <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="BP",OrgDb='org.Hs.eg.db',pvalueCutoff = 1)
dotplot(cgBP,showCategory = 10) + 
  theme(axis.text.y = element_text(size = 9, lineheight = 0.7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

dotplot(cgBP, showCategory = 6) +
  scale_color_continuous(
    low = "red",
    high = "blue", 
    limits = c(0.004, 0.01),
    oob = squish   
  ) +
  theme(axis.text.y = element_text(size = 9, lineheight = 0.7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

df <- cgBP@compareClusterResult

df_keep <- bind_rows(
  df %>% group_by(Cluster) %>% slice_min(order_by = p.adjust, n = 9, with_ties = FALSE) %>% ungroup(),
  df %>% filter(Description == "neutrophil chemotaxis")
) %>% distinct()

cgBP_sub <- cgBP
cgBP_sub@compareClusterResult <- df_keep

dotplot(cgBP_sub, showCategory = 10) +
  theme(axis.text.y = element_text(size = 9, lineheight = 0.7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))



combined_Mye4 <- subset(combined_Mye3, idents = c("0","1","3","4"))
DefaultAssay(combined_Mye4) <- "RNA"
combined_Mye4 <- NormalizeData(combined_Mye4, verbose = FALSE)
combined_Mye4 <- FindVariableFeatures(combined_Mye4, selection.method = "vst", nfeatures = 3000)
combined_Mye4 <- ScaleData(combined_Mye4, verbose = FALSE)
combined_Mye4 <- SCTransform(combined_Mye4, vars.to.regress = c("percent.mt"))
combined_Mye4 <- RunPCA(combined_Mye4, assay = "SCT", verbose = FALSE)
#combined_Mye4 <- RunPCA(combined_Mye4, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_Mye4 <- RunHarmony(combined_Mye4, group.by.vars = "sample", assay.use = "RNA",
                            reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_Mye4, ndims = 50, reduction = "harmony")
combined_Mye4 <- FindNeighbors(combined_Mye4, reduction = "harmony", dims = 1:20)
combined_Mye4 <- FindClusters(combined_Mye4, resolution = 0.3)
combined_Mye4 <- RunUMAP(combined_Mye4, reduction = "harmony", dims = 1:20)

DimPlot(combined_Mye4, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_Mye4, reduction = "umap", split.by = "LoY")

table(combined_Mye4$LoY)

###Heatmap
markers <- FindAllMarkers(combined_Mye4, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

combined_Mye4_heat <- ScaleData(combined_Mye4, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_Mye4_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p



###Phate
# r-reti2 pip install phate
library(reticulate)

use_condaenv("r-reti2")

if (!suppressWarnings(require(devtools))) install.packages("devtools")
reticulate::py_install("phate", pip=TRUE)
devtools::install_github("KrishnaswamyLab/phateR")

library(phateR)

combined_Mye3 <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_DCM_sample_Myeloid_cell_gene3000_260414.rds")
seurat_obj <- subset(combined_Mye3, idents = c("IL1B+ inflammatory monos / macs","LYVE1+ macs","TREM2+ macs","TNF+ inflammatory macs"))
seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)

Idents(seurat_obj) <- "cell_annotation_final"
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)


seurat_obj <- subset(seurat_obj,idents = c("classical mono","IL1B+ inflammatory mac"), invert=TRUE)

DefaultAssay(seurat_obj) <- "RNA"
set.seed(123)
pca_coords <- seurat_obj@reductions$harmony@cell.embeddings
pca_coords <- as.matrix(pca_coords)

#pca_coords <- t(pca_coords)

dim(pca_coords)
phate_emb <- phate(pca_coords, knn=100,t=55)

mat <- phate_emb$embedding
dim(mat)
rownames(mat) <- Cells(seurat_obj)
colnames(mat) <- c("PHATE_1", "PHATE_2")

seurat_obj[["phate"]] <- CreateDimReducObject(
  embeddings = mat,
  key = "PHATE_",
  assay = DefaultAssay(seurat_obj)
)

DimPlot(seurat_obj,reduction = "phate",label = F, cols = c("#00BFC4","#00A9FF","#C77CFF","#FF61CC"))
DimPlot(seurat_obj,reduction = "phate",label = F)

FeaturePlot(seurat_obj,reduction = "phate", features = 'IL1B')

seurat_obj <- AddMetaData(seurat_obj,Idents(seurat_obj),col.name = "cell_annotation_final")
Idents(seurat_obj) <- "LoY"
Idents(seurat_obj) <- "cell_annotation_final"

levels(seurat_obj)


DimPlot(seurat_obj, reduction = "phate",cells.highlight = WhichCells(seurat_obj, idents = c("LoY")), pt.size = 0.01)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot")

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(phateR)
library(ggplot2)

sce <- as.SingleCellExperiment(seurat_obj)

X_harmony <- Embeddings(seurat_obj, "harmony")[, 1:50]

X_harmony <- X_harmony[colnames(sce), , drop = FALSE]

reducedDim(sce, "Harmony") <- X_harmony

colData(sce)$cluster <- Idents(seurat_obj)[colnames(sce)]

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim = "Harmony"
)

pt <- slingPseudotime(sce)[, 2]
seurat_obj$SLING <- pt

#phate_emb <- phate(pca_coords)
set.seed(123)

ph <- phate(pca_coords, knn=100,t=55)
ph_df <- as.data.frame(ph$embedding)
colnames(ph_df) <- c("PHATE1","PHATE2")
ph_df$SLING <- pt

ggplot(ph_df, aes(PHATE1, PHATE2, color = SLING)) +
  geom_point(size = 0.4, alpha = 0.8) +
  theme_classic() +
  labs(title="PHATE colored by Slingshot pseudotime", color="pseudotime")

mat <- ph$embedding
dim(mat)
rownames(mat) <- Cells(seurat_obj)
colnames(mat) <- c("PHATE_1", "PHATE_2")

seurat_obj[["phate"]] <- CreateDimReducObject(
  embeddings = mat,
  key = "PHATE_",
  assay = DefaultAssay(seurat_obj)
)
Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = c("IL1B+ inflammatory monos / macs","LYVE1+ macs","TNF+ inflammatory macs","TREM2+ macs","proliferative macs"))

DimPlot(seurat_obj,reduction = "phate",label = F, cols = c("#00BE67","#00BFC4","#00A9FF","#C77CFF"))


umap_df <- as.data.frame(Embeddings(seurat_obj, "phate"))
colnames(umap_df)[1:2] <- c("PHATE_1", "PHATE_2")

umap_df$SLING <- seurat_obj$SLING

ggplot(umap_df, aes(PHATE_1, PHATE_2, color = SLING)) +
  geom_point(size = 0.4, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "UMAP colored by Slingshot pseudotime",
    color = "pseudotime"
  )

seurat_obj <- AddMetaData(seurat_obj,Idents(seurat_obj),col.name = "cell_annotation_final")
Idents(seurat_obj) <- "LoY"
Idents(seurat_obj) <- "cell_annotation_final"

DimPlot(seurat_obj, reduction = "phate",cells.highlight = WhichCells(seurat_obj, idents = c("LoY")))

saveRDS(seurat_obj, "~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_DCM_sample_Myeloid_cell_gene3000_260420_PHATE_analysis.rds")

seurat_obj <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_DCM_sample_Myeloid_cell_gene3000_260420_PHATE_analysis.rds")

Idents(combined_Mye) <- "cell_annotation_final"

combined_Mye$cell_annotation_final <- recode(
  combined_Mye$cell_annotation_final,
  `LYVE1+ resident mac` = "LYVE1+ mac"
)

table(combined_Mye$cell_annotation_final)

combined_Mye <- AddMetaData(combined_Mye, paste(Idents(combined_Mye), combined_Mye$LoY),col.name = "cluster_LoY")
table(combined_Mye$cluster_LoY)



test3 <- subset(combined_Mye3, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)
test3 <- combined

id <- 1:ncol(test3)
name <- Idents(test3)
LoY <- test3$LoY

data = data.frame(id, name, LoY, stringsAsFactors = F)

library(dplyr)
library(ggplot2)
library(scales)

df_sum <- data %>%
  dplyr::count(name, LoY) %>%
  group_by(name) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

df_sum <- df_sum %>%
  mutate(LoY = factor(LoY, levels = c("no_LoY", "LoY")))

df_sum$name2 <- factor(df_sum$name, levels = c("proliferative macs","TREM2+ macs","TNF+ inflammatory macs","LYVE1+ macs","IL1B+ inflammatory monos / macs","cDC2","cDC1","Non classical monos"))  # ←順序を逆に


ggplot(df_sum, aes(x = name, y = p, fill = LoY)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = percent(p, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("LoY" = "#F8766D", "no_LoY" = "#00BFC4")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  coord_flip() +
  theme_minimal()




get_MSY_sum_per_cell <- function(seurat_obj, gene_list, slot = "counts") {
  
  mat <- GetAssayData(seurat_obj, slot = slot)

  present <- intersect(gene_list, rownames(mat))
  if (length(present) == 0) {
    stop("指定した gene_list の遺伝子が、このオブジェクトには1つも存在しません。")
  }
  
  msy_mat <- mat[present, , drop = FALSE]  # genes × cells
  

  MSY_sum   <- Matrix::colSums(msy_mat)
  MSY_mean  <- MSY_sum / length(present)
  MSY_n_det <- Matrix::colSums(msy_mat > 0)
  

  meta <- seurat_obj@meta.data
  meta$cell <- rownames(meta)
  
  df <- data.frame(
    cell        = colnames(seurat_obj),
    MSY_sum     = as.numeric(MSY_sum[colnames(seurat_obj)]),
    MSY_mean    = as.numeric(MSY_mean[colnames(seurat_obj)]),
    MSY_n_gene  = as.numeric(MSY_n_det[colnames(seurat_obj)]),
    meta[colnames(seurat_obj), , drop = FALSE],
    row.names   = colnames(seurat_obj)
  )
  
  return(df)
}

msy_df <- get_MSY_sum_per_cell(test, gene_list, slot = "counts")

library(ggplot2)

ggplot(msy_df, aes(x = MSY_sum)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    x = "MSY_sum（MSY遺伝子の合計カウント）",
    y = "細胞数",
    title = "MSY_sum の分布（全細胞）"
  )


test2 <- subset(combined_Mye3, subset = LoY == "no_LoY")
test3 <- subset(test2, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)
VlnPlot(test3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

df <- fetch_counts_with_zeros(test2, gene_list, slot = "counts")

df["sum"] <- rowSums(df)

freq <- table(df["sum"])      
pct <- prop.table(table(df$sum)) * 100

barplot(
  pct,
  xlab = "LoY gene counts",
  ylab = "Percent (%)",
  ylim = c(0, max(pct) * 1.1)
)

mean_loy <- mean(df$sum, na.rm = TRUE)
median_loy <- median(df$sum, na.rm = TRUE)

cat("Mean of LoY gene counts:", mean_loy, "\n")
cat("Median of LoY gene counts:", median_loy, "\n")

####DCM
#Mean of LoY gene counts: 3.591966 
#Median of LoY gene counts: 3 

fraction_percent <- mean(df$sum > 0) * 100

fraction_percent


test3 <- subset(combined_Mye3, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)
test3 <- subset(test3, idents = c("LYVE1+ macs"))
table(test3$LoY)



combined_Mye <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/With_public_data/All_male_sample_Myeloid_260305_non_doublets_final4.rds")
combined_Mye3 <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_DCM_sample_Myeloid_cell_gene3000_260414.rds")

combined_Mye <- AddMetaData(combined_Mye, paste(combined_Mye$sample, combined_Mye$LoY),col.name = "sample_LoY")
combined_Mye3 <- AddMetaData(combined_Mye3, paste(combined_Mye3$sample, combined_Mye3$LoY),col.name = "sample_LoY")

meta <- combined_Mye3@meta.data %>%
  tibble::rownames_to_column("cell")

library(lme4)

meta2 <- meta %>%
  mutate(
    LoY_bin = ifelse(LoY == "LoY", 1, 0),
    seurat_clusters = factor(cell_annotation_final),
    sample = factor(sample)
  )

fit <- glmer(
  LoY_bin ~ seurat_clusters + (1 | sample),
  data = meta2,
  family = binomial
)

summary(fit)

target_cluster <- "LYVE1+ macs"

meta3 <- meta2 %>%
  mutate(target_vs_other = factor(ifelse(seurat_clusters == target_cluster, "target", "other"),
                                  levels = c("other", "target")))

fit_target <- glmer(
  LoY_bin ~ target_vs_other + (1 | sample),
  data = meta3,
  family = binomial
)

summary(fit_target)

exp(fixef(fit_target))



######all cell UMAP

#####annotation
#0  NK/cytotoxic lymphocytes
#1  Inflammatory macrophages / monocyte-derived macrophages
#2  Blood endothelial cells
#3  Fibroblasts
#4  Pericytes / mural cells
#5  Ventricular cardiomyocytes
#6  Plasma cells
#7  Cycling cells proliferative cell
#8  Lymphatic endothelial cells
#9  Schwann cells / peripheral glial cells
#10 Mast cells

combined <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_sample_all_cell_DCM.rds")

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)

combined <- subset(combined, subset = nCount_RNA > 3000)


combined <- AddMetaData(combined, Idents(combined), col.name = "first_annotation")

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:50)
combined <- FindClusters(combined, resolution = 0.05)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:50)


DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined, reduction = "umap", split.by = "LoY")

Idents(combined) <- "first_annotation"

table(combined$LoY)
table(combined$sample)

combined_subset <- subset(combined, cells = Cells(combined[["RNA"]]), downsample = 2500)
markers <- FindAllMarkers(combined, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)

combined_subset <- ScaleData(combined_subset, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_subset, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p


saveRDS(combined, "~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_sample_all_cell_DCM_final.rds")
combined <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_sample_all_cell_final.rds")

##### immune cell

combined <- AddMetaData(combined, Idents(combined), col.name = "first_annotation")

combined_immune <- subset(combined, idents = c("0","1","6"))
DefaultAssay(combined_immune) <- "RNA"
combined_immune <- NormalizeData(combined_immune, verbose = FALSE)
combined_immune <- FindVariableFeatures(combined_immune, selection.method = "vst", nfeatures = 3000)
combined_immune <- ScaleData(combined_immune, verbose = FALSE)
combined_immune <- SCTransform(combined_immune, vars.to.regress = c("percent.mt"))
combined_immune <- RunPCA(combined_immune, assay = "SCT", verbose = FALSE)
#combined_immune <- RunPCA(combined_immune, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_immune <- RunHarmony(combined_immune, group.by.vars = "sample", assay.use = "SCT",
                              reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_immune, ndims = 50, reduction = "harmony")
combined_immune <- FindNeighbors(combined_immune, reduction = "harmony", dims = 1:15)
combined_immune <- FindClusters(combined_immune, resolution = 0.09)
combined_immune <- RunUMAP(combined_immune, reduction = "harmony", dims = 1:15)

DimPlot(combined_immune, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_immune, reduction = "umap", split.by = "LoY")

table(combined_immune$LoY)

combined_subset <- subset(combined_immune, cells = Cells(combined_immune[["RNA"]]), downsample = 6000)
markers <- FindAllMarkers(combined_immune, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)

combined_subset <- ScaleData(combined_subset, assay = "RNA", features = top10$gene)
p1 <- DoHeatmap(combined_subset, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p1

print(top10, n=49)

VlnPlot(combined_immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)

##### immune cell

combined_immune2 <- subset(combined_immune, idents = c("4"), invert=TRUE)
DefaultAssay(combined_immune2) <- "RNA"
combined_immune2 <- NormalizeData(combined_immune2, verbose = FALSE)
combined_immune2 <- FindVariableFeatures(combined_immune2, selection.method = "vst", nfeatures = 3000)
combined_immune2 <- ScaleData(combined_immune2, verbose = FALSE)
combined_immune2 <- SCTransform(combined_immune2, vars.to.regress = c("percent.mt"))
combined_immune2 <- RunPCA(combined_immune2, assay = "SCT", verbose = FALSE)
#combined_immune2 <- RunPCA(combined_immune2, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_immune2 <- RunHarmony(combined_immune2, group.by.vars = "sample", assay.use = "SCT",
                               reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_immune2, ndims = 50, reduction = "harmony")
combined_immune2 <- FindNeighbors(combined_immune2, reduction = "harmony", dims = 1:15)
combined_immune2 <- FindClusters(combined_immune2, resolution = 0.11)
combined_immune2 <- RunUMAP(combined_immune2, reduction = "harmony", dims = 1:15)

DimPlot(combined_immune2, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_immune2, reduction = "umap", split.by = "LoY")

table(combined_immune2$LoY)

combined_subset <- subset(combined_immune2, cells = Cells(combined_immune2[["RNA"]]), downsample = 6000)
markers <- FindAllMarkers(combined_immune2, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)

combined_subset <- ScaleData(combined_subset, assay = "RNA", features = top10$gene)
p1 <- DoHeatmap(combined_subset, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p1

print(top10, n=49)

VlnPlot(combined_immune2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)

saveRDS(combined_immune2, "~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_sample_immune_cell_DCM_final2.rds")
combined_immune2 <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_sample_immune_cell_final2.rds")


#####annotation
#0  macrophage
#1  NK/T
#2  non classical mono
#3  classical mono
#4  B
#5  plasma


##LoY graph
id <- 1:ncol(combined_immune2)
name <- Idents(combined_immune2)
LoY <- combined_immune2$LoY

data = data.frame(id, name, LoY, stringsAsFactors = F)

library(dplyr)
library(ggplot2)
library(scales)

df_sum <- data %>%
  dplyr::count(name, LoY) %>%
  group_by(name) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

df_sum <- df_sum %>%
  mutate(LoY = factor(LoY, levels = c("no_LoY", "LoY")))

df_sum$name2 <- factor(df_sum$name, levels = c(5,4,3,2,1,0))  

ggplot(df_sum, aes(x = name2, y = p, fill = LoY)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = percent(p, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("LoY" = "#F8766D", "no_LoY" = "#00BFC4")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  coord_flip() +
  theme_minimal()




##########  vlnplot like
DimPlot(combined_Mye3)
combined_Mac <- subset(combined_Mye3, idents= c("TNF+ inflammatory macs","LYVE1+ macs","TREM2+ macs"))

combined_Mac$cell_type_final <- recode(
  combined_Mac$cell_annotation_final,
  `LYVE1+ macs` = "Remodeling-associated mac",
  `TNF+ inflammatory macs` = "Inflammatory mac",
  `TREM2+ macs` = "Remodeling-associated mac"
)###Heatmap

Idents(combined_Mac) <- "cell_type_final"
DimPlot(combined_Mac)

df <- data.frame(
  cell    = colnames(combined_Mac),
  name    = as.character(Idents(combined_Mac)),
  LoY     = combined_Mac$LoY,
  sample  = combined_Mac$sample,
  stringsAsFactors = FALSE
)

df <- df %>%
  mutate(
    LoY = factor(LoY, levels = c("no_LoY", "LoY")),
    name = factor(
      name,
      levels = c(
        "Remodeling-associated mac",
        "Inflammatory mac"
      )
    )
  )


df_prop <- df %>%
  group_by(sample, name) %>%
  summarise(
    n_total = n(),
    n_LoY = sum(LoY == "LoY", na.rm = TRUE),
    prop_LoY = n_LoY / n_total,
    .groups = "drop"
  )

head(df_prop)

boxplot(prop_LoY~name,data=df_prop, outline=F, width=c(0.6,0.6), col=c("#00B0F0", "#FFC000"))
stripchart(prop_LoY~name,              # Data
           data=df_prop,
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = "black",           # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE) 

ggplot(df_prop, aes(x = name, y = prop_LoY)) +
  geom_violin(
    aes(fill = name),
    trim = FALSE,
    scale = "width",
    color = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.9,
    color = "black"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "LoY fraction per sample"
  ) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14)
  )


result <- wilcox.test(prop_LoY~name,data=df_prop, paired = TRUE)
print(result)





####age,LoY proportion
library(ggplot2)
library(ggpubr)
library(scales)

age <- c(79,75,72,70,64,62,58,60,62,55,74,60)
LoY_fraction <- c(27.61,14.58,3.75,14.00,18.58,8.13,9.08,5.74,8.10,8.49,5.20,7.55)

df <- data.frame(
  age = age,
  LoY_fraction = LoY_fraction
)

library(ggplot2)

ggplot(df, aes(x = age, y = LoY_fraction)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 14) +
  labs(
    x = "Age",
    y = "LoY fraction",
    title = "Correlation between age and LoY fraction"
)

library(ggplot2)
library(ggpubr)

ggplot(df, aes(x = age, y = LoY_fraction)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "pearson") +
  theme_classic(base_size = 14) +
  labs(
    x = "Age",
    y = "LoY fraction",
    title = "Correlation between age and LoY fraction"
  )

## R=0.47





