library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

AF.list<-list()

##Our data 9 
data.Aur2 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur2_without_CellPlex/outs/filtered_feature_bc_matrix")
data.Aur3Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur3/outs/per_sample_outs/Aur3_Endo/count/sample_feature_bc_matrix")
data.Aur3Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur3/outs/per_sample_outs/Aur3_Epi/count/sample_feature_bc_matrix")
data.Aur4Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4/outs/per_sample_outs/Aur4_Endo/count/sample_feature_bc_matrix")
data.Aur4Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4/outs/per_sample_outs/Aur4_Epi/count/sample_feature_bc_matrix")
data.Aur12 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur12/outs/filtered_feature_bc_matrix")
data.Aur13 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur13/outs/filtered_feature_bc_matrix")
data.Aur10 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur10/outs/filtered_feature_bc_matrix")
data.Aur16 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur16/outs/filtered_feature_bc_matrix")
data.Aur19 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur19/outs/filtered_feature_bc_matrix")

##science
data.sci1 <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE197518_RAW/GSM5919345_1_MR_AF_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
data.sci2 <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE197518_RAW/GSM5919346_2_MR_AF_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
data.sci3 <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE197518_RAW/GSM5919347_3_MR_AF_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
data.sci4 <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE197518_RAW/GSM5919348_4_MR_AF_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
data.sci7 <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE197518_RAW/GSM5919349_7_MR_AF_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
data.sci13 <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE197518_RAW/GSM5919350_13_MR_AF_filtered_feature_bc_matrix/filtered_feature_bc_matrix")
data.sci14 <- Read10X(data.dir = "~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE197518_RAW/GSM5919351_14_MR_AF_filtered_feature_bc_matrix/filtered_feature_bc_matrix")


##nature
data.nat1 <- Read10X_h5("~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE263154_RAW/GSM8186782_filtered_feature_bc_matrix_RAA1.h5")
data.nat2 <- Read10X_h5("~/Desktop/Loss_of_Y_analysis/public_AF_dataset/GSE263154_RAW/GSM8186788_filtered_feature_bc_matrix_RAA2.h5")

Aur.list <- list(
  Aur2      = data.Aur2,
  Aur3Endo  = data.Aur3Endo$`Gene Expression`,
  Aur3Epi   = data.Aur3Epi$`Gene Expression`,
  Aur4Endo  = data.Aur4Endo$`Gene Expression`,
  Aur4Epi   = data.Aur4Epi$`Gene Expression`,
  Aur12     = data.Aur12,
  Aur13     = data.Aur13,
  Aur10     = data.Aur10,
  Aur16     = data.Aur16,
  Aur19     = data.Aur19
)

Sci.list <- list(
  # science (7 datasets)
  sci1 = data.sci1,
  sci2 = data.sci2,
  sci3 = data.sci3,
  sci4 = data.sci4,
  sci7 = data.sci7,
  sci13 = data.sci13,
  sci14 = data.sci14
)


Nat.list <- list(
  # nature (2 datasets)
  nat1 = data.nat1$`Gene Expression`,
  nat2 = data.nat2$`Gene Expression`
)

raw.list <- list(
  Aur2      = data.Aur2,
  Aur3Endo  = data.Aur3Endo$`Gene Expression`,
  Aur3Epi   = data.Aur3Epi$`Gene Expression`,
  Aur4Endo  = data.Aur4Endo$`Gene Expression`,
  Aur4Epi   = data.Aur4Epi$`Gene Expression`,
  Aur12     = data.Aur12,
  Aur13     = data.Aur13,
  Aur10     = data.Aur10,
  Aur16     = data.Aur16,
  Aur19     = data.Aur19,
  
  # science (7 datasets)
  sci1 = data.sci1,
  sci2 = data.sci2,
  sci3 = data.sci3,
  sci4 = data.sci4,
  sci7 = data.sci7,
  sci13 = data.sci13,
  sci14 = data.sci14,
  
  # nature (2 datasets)
  nat1 = data.nat1$`Gene Expression`,
  nat2 = data.nat2$`Gene Expression`
)



process_sample <- function(counts, sample_name) {
  
  # Create a Seurat object
  obj <- CreateSeuratObject(
    counts = counts, 
    project = "AF_project",
    min.cells = 3,
    min.features = 200
  )
  
  # Sample name
  obj <- AddMetaData(obj, sample_name, col.name = "sample")
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
  
  obj <- RenameCells(obj, new.names = paste0(colnames(obj), "_", sample_name))
  
  return(obj)
}



for (s in names(Aur.list)) {
  message("Processing: ", s)
  Aur.list[[s]] <- process_sample(Aur.list[[s]], s)
}

for (s in names(Sci.list)) {
  message("Processing: ", s)
  Sci.list[[s]] <- process_sample(Sci.list[[s]], s)
}

for (s in names(raw.list)) {
  message("Processing: ", s)
  AF.list[[s]] <- process_sample(raw.list[[s]], s)
}


length(AF.list)  

gene_list <- read.csv("~/Desktop/yoshida_test/Loss_of_Y/MSY_gene.txt")
class(gene_list)
gene_list <- gene_list$Gene.name


## Function to retrieve counts of MSY genes, assigning 0 to genes that are absent -----------------
fetch_counts_with_zeros <- function(seurat_obj, gene_list, slot = "counts") {
  
  # Extract raw counts
  mat <- GetAssayData(seurat_obj, slot = slot)
  
  # Output matrix (cell × gene_list)
  out <- matrix(0, nrow = ncol(mat), ncol = length(gene_list))
  colnames(out) <- gene_list
  rownames(out) <- colnames(mat)
  
  # Replace only genes that are present in gene_list
  present <- intersect(gene_list, rownames(mat))
  
  if (length(present) > 0) {
    out[, present] <- t(as.matrix(mat[present, , drop = FALSE]))
  } else {
    warning("None of the specified genes in gene_list are present in this object. All values will be treated as 0.")
  }
  
  return(as.data.frame(out))
}

## LoY classification: define cells as LoY if all MSY gene counts are 0 -----------------------
## obj_list    : named list of sample IDs -> Seurat objects
## male_samples: vector of male sample IDs to be evaluated for LoY
## gene_list   : MSY gene panel (e.g., c("RPS4Y1","ZFY","USP9Y","DDX3Y","KDM5D","EIF1AY"))
## min_features: QC threshold (cells with nFeature_RNA < min_features are excluded from classification)
judge_LoY <- function(obj_list, male_samples, gene_list, min_features = 500) {
  
  # ① Collect cell × MSY gene counts from male samples ---------------------------
  expr_list <- lapply(male_samples, function(s) {
    obj <- obj_list[[s]]
    
    # Extract raw counts of MSY genes (missing genes are set to 0)
    df  <- fetch_counts_with_zeros(obj, gene_list, slot = "counts")
    
    # Metadata (add nFeature_RNA)
    meta <- obj@meta.data[colnames(obj), , drop = FALSE]
    df$nFeature_RNA <- meta$nFeature_RNA
    df$sample <- s
    
    return(df)
  })
  
  # cell × (gene_list + nFeature_RNA + sample)
  expr_all <- do.call(rbind, expr_list)
  
  # ② QC: remove cells with too few nFeature_RNA ------------------------------------
  expr_filt <- expr_all[expr_all$nFeature_RNA >= min_features, , drop = FALSE]
  
  if (nrow(expr_filt) == 0) {
    warning("No male cells passed QC (nFeature_RNA >= min_features).")
    return(character(0))
  }
  
  # ③ LoY classification: all genes in gene_list are 0 (the number of genes with counts > 0 is 0) ----------
  is_LoY <- rowSums(expr_filt[, gene_list, drop = FALSE] > 0) == 0
  
  # Return cell names (cell barcodes) classified as LoY
  LoY_cells <- rownames(expr_filt)[is_LoY]
  
  return(LoY_cells)
}

male_samples_Aur <- c("Aur10","Aur12","Aur13","Aur19","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi")  # Aur male

male_samples_sci <- c("sci1", "sci2","sci7", "sci13","sci14")  # Sci male

male_samples_Nat <- c("nat1", "nat2")  # Nat


Aur_LoY <- judge_LoY(AF.list, male_samples_Aur, gene_list,min_features = 500)
Sci_LoY <- judge_LoY(AF.list, male_samples_sci, gene_list,min_features = 500)
Nat_LoY <- judge_LoY(AF.list, male_samples_Nat, gene_list,min_features = 500)

LoY_cells <- c(Aur_LoY, Sci_LoY, Cli_LoY, Nat_LoY)

## Function: add LoY metadata to AF.list all at once
add_LoY_metadata <- function(AF.list, LoY_cells) {
  
  for (s in names(AF.list)) {
    obj <- AF.list[[s]]
    
    # Initialization (all FALSE)
    obj[["LoY"]] <- "no_LoY"
    
    # Set only overlapping cells to TRUE
    common <- intersect(colnames(obj), LoY_cells)
    if (length(common) > 0) {
      obj$LoY[common] <- "LoY"
    }
    
    AF.list[[s]] <- obj
  }
  
  return(AF.list)
}

AF.list <- add_LoY_metadata(AF.list, LoY_cells)

AF.list[["Aur2"]] <- NULL
AF.list[["Aur16"]] <- NULL
AF.list[["sci3"]] <- NULL
AF.list[["sci4"]] <- NULL

combined <- Reduce(function(x, y) merge(x, y), AF.list)

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
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.05)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)

DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined, reduction = "umap", split.by = "LoY")
#DimPlot(combined, reduction = "umap", split.by = "sample")

table(combined$sample)

combined_subset <- subset(combined, cells = Cells(combined[["RNA"]]), downsample = 2500)
markers <- FindAllMarkers(combined_subset, assay = "RNA", only.pos = TRUE)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
combined_subset <- ScaleData(combined_subset, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_subset, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

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
VlnPlot(combined, features = c("CD3E","CD4","CD8A","IL7R","LEF1","GZMK") ,pt.size = 0)

VlnPlot(combined, features = c("GLS"), split.by = "stim" , cols = c("sa" ="deepskyblue", "sym"="red"))


#cDC1
VlnPlot(combined, pt.size = 0 , features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
#cDC2
VlnPlot(combined, pt.size = 0 , features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
#pDC
VlnPlot(combined, pt.size = 0 , features = c("IL3RA","CLEC4C") )


### Extract myeloid cells

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
combined_Mye <- FindNeighbors(combined_Mye, reduction = "harmony", dims = 1:45)
combined_Mye <- FindClusters(combined_Mye, resolution = 0.05)
combined_Mye <- RunUMAP(combined_Mye, reduction = "harmony", dims = 1:45)

DimPlot(combined_Mye, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_Mye, reduction = "umap", split.by = "LoY")
DimPlot(combined_Mye, reduction = "umap", split.by = "annotation")

Idents(combined_Mye) <- "cluster_global"

table(combined_Mye$LoY)
table(combined_Mye$sample)

VlnPlot(combined_Mye, assay = "RNA",features = c("CD68","CD14","FCGR3A") ,pt.size = 0)

###Heatmap
markers <- FindAllMarkers(combined_Mye, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

combined_Mye_heat <- ScaleData(combined_Mye, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_Mye_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=80)


### Extract macrophages, monocytes, and DCs
combined_macmod <- subset(combined_Mye, idents = c("0","1","2")) 
DefaultAssay(combined_macmod) <- "RNA"
combined_macmod <- NormalizeData(combined_macmod, verbose = FALSE)
combined_macmod <- FindVariableFeatures(combined_macmod, selection.method = "vst", nfeatures = 3000)
combined_macmod <- ScaleData(combined_macmod, verbose = FALSE)
#combined_macmod <- SCTransform(combined_macmod, vars.to.regress = c("percent.mt"))
#combined_macmod <- RunPCA(combined_macmod, assay = "SCT", verbose = FALSE)
combined_macmod <- RunPCA(combined_macmod, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_macmod <- RunHarmony(combined_macmod, group.by.vars = "sample", assay.use = "RNA",
                              reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_macmod, ndims = 50, reduction = "harmony")
combined_macmod <- FindNeighbors(combined_macmod, reduction = "harmony", dims = 1:40)
combined_macmod <- FindClusters(combined_macmod, resolution = 0.30)
combined_macmod <- RunUMAP(combined_macmod, reduction = "harmony", dims = 1:40)

DimPlot(combined_macmod, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_macmod, reduction = "umap", split.by = "LoY")
DimPlot(combined_macmod, reduction = "umap", split.by = "annotation")

VlnPlot(combined_macmod, assay = "RNA",features = c("CD68","CD14","FCGR3A") ,pt.size = 0)

table(combined_macmod$LoY)
table(combined_macmod$sample)

###Heatmap
markers <- FindAllMarkers(combined_macmod, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

combined_macmod_heat <- ScaleData(combined_macmod, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_macmod_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=70)

FeaturePlot(combined_macmod, features=c('CD14'), min.cutoff=0, max.cutoff='q90')


### Extract macrophages and monocytes
combined_macmod2 <- subset(combined_macmod, idents = c("0","1","2","6")) 
DefaultAssay(combined_macmod2) <- "RNA"
combined_macmod2 <- NormalizeData(combined_macmod2, verbose = FALSE)
combined_macmod2 <- FindVariableFeatures(combined_macmod2, selection.method = "vst", nfeatures = 3000)
combined_macmod2 <- ScaleData(combined_macmod2, verbose = FALSE)
#combined_macmod2 <- SCTransform(combined_macmod2, vars.to.regress = c("percent.mt"))
#combined_macmod2 <- RunPCA(combined_macmod2, assay = "SCT", verbose = FALSE)
combined_macmod2 <- RunPCA(combined_macmod2, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_macmod2 <- RunHarmony(combined_macmod2, group.by.vars = "sample", assay.use = "RNA",
                               reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_macmod2, ndims = 50, reduction = "harmony")
combined_macmod2 <- FindNeighbors(combined_macmod2, reduction = "harmony", dims = 1:50)
combined_macmod2 <- FindClusters(combined_macmod2, resolution = 0.5)
combined_macmod2 <- RunUMAP(combined_macmod2, reduction = "harmony", dims = 1:50)

DimPlot(combined_macmod2, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_macmod2, reduction = "umap", split.by = "LoY")
DimPlot(combined_macmod2, reduction = "umap", split.by = "annotation")

table(combined_macmod2$LoY)
table(combined_macmod2$sample)

###Heatmap
markers <- FindAllMarkers(combined_macmod2, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

combined_macmod2_heat <- ScaleData(combined_macmod2, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_macmod2_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=120)

VlnPlot(combined_macmod2, assay = "RNA",features = c("CD52","CD68","CD14") ,pt.size = 0)

#IL1B+
VlnPlot(combined_macmod2, assay = "RNA",features = c("IL1B","EREG","AREG","HBEGF","CXCL2","VEGFA","CD44","IL1RN","CXCL3") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("IL1B","EREG","AREG","HBEGF","CXCL2","CXCL3") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("TNF","CXCL2") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("IL1B","AREG") ,pt.size = 0)

#TREM2+
VlnPlot(combined_macmod2, assay = "RNA",features = c("TREM2","SPP1","LPL","APOE","ABCG1") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("TREM2","SPP1") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("APOE","C1QC") ,pt.size = 0)

#LYVE1+
VlnPlot(combined_macmod2, assay = "RNA",features = c("LYVE1","FOLR2","IGF1","C1QA") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("LYVE1","FOLR2") ,pt.size = 0)

#No mast cells
VlnPlot(combined_macmod2, assay = "RNA",features = c("HDC","KIT") ,pt.size = 0)

#No Mreg DC
VlnPlot(combined_macmod2, assay = "RNA",features = c("FSCN1","CCR7","HLA-DRA","HLA-DPA1") ,pt.size = 0)

#Classical monocytes
VlnPlot(combined_macmod2, assay = "RNA",features = c("LYZ","S100A12","S100A8","CD14","VCAN","FOS") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("S100A12","VCAN") ,pt.size = 0)

#Non-classical monocytes
VlnPlot(combined_macmod2, assay = "RNA",features = c("CSF1R","FCGR3A","CXCL10","CCL4","CXCL8") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA",features = c("CSF1R","FCGR3A") ,pt.size = 0)

#c1q-hi
VlnPlot(combined_macmod2, assay = "RNA", features = c("C1QA","C1QB","SELENOP","GIPC2","C1QC","CD68","FOLR2","PLTP") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA", features = c("FOLR2","LYVE1","TREM2","SPP1") ,pt.size = 0)

#C1Q+ resident-like
VlnPlot(combined_macmod2, assay = "RNA", features = c("FTL","CD163","HMOX1","MAF","MRC1","SIGLEC1") ,pt.size = 0) #CD169=SIGLEC1,CD206=MRC1
VlnPlot(combined_macmod2, assay = "RNA", features = c("APOC1","APOE","ALOX5AP","INHBA","GPNMB") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA", features = c("C3AR1","CFD","FTH1","ACP5") ,pt.size = 0)
VlnPlot(combined_macmod2, assay = "RNA", features = c("FTL","PLCG2") ,pt.size = 0) #CD169=SIGLEC1,CD206=MRC1

FeaturePlot(combined_macmod2, features=c('IL1B'), min.cutoff=0, max.cutoff='q90')


## Remove neutrophils
combined_macmod3 <- subset(combined_macmod2, idents = c("0","1","2","3","4","5","6"))

DefaultAssay(combined_macmod3) <- "RNA"
combined_macmod3 <- NormalizeData(combined_macmod3, verbose = FALSE)
combined_macmod3 <- FindVariableFeatures(combined_macmod3, selection.method = "vst", nfeatures = 3000)
combined_macmod3 <- ScaleData(combined_macmod3, verbose = FALSE)
#combined_macmod3 <- SCTransform(combined_macmod3, vars.to.regress = c("percent.mt"))
#combined_macmod3 <- RunPCA(combined_macmod3, assay = "SCT", verbose = FALSE)
combined_macmod3 <- RunPCA(combined_macmod3, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_macmod3 <- RunHarmony(combined_macmod3, group.by.vars = "sample", assay.use = "RNA",
                               reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_macmod3, ndims = 50, reduction = "harmony")
combined_macmod3 <- FindNeighbors(combined_macmod3, reduction = "harmony", dims = 1:20)
combined_macmod3 <- FindClusters(combined_macmod3, resolution = 0.45)
combined_macmod3 <- RunUMAP(combined_macmod3, reduction = "harmony", dims = 1:20)

DimPlot(combined_macmod3, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_macmod3, reduction = "umap", split.by = "LoY")
DimPlot(combined_macmod3, reduction = "umap", split.by = "annotation")

table(combined_macmod3$LoY)
table(combined_macmod3$sample)

FeaturePlot(combined_macmod3, features=c('FN1'), min.cutoff=0, max.cutoff='q90')

VlnPlot(combined_macmod3, assay = "RNA",features = c("IL1B","EREG","AREG","HBEGF","CXCL2","VEGFA","CD44","IL1RN","CXCL3") ,pt.size = 0)
VlnPlot(combined_macmod3, assay = "RNA",features = c("LYZ","S100A12","S100A8","CD14","VCAN","FOS") ,pt.size = 0)
VlnPlot(combined_macmod3, assay = "RNA",features = c("TNF","CXCL2","TGFB1","FN1") ,pt.size = 0)
VlnPlot(combined_macmod3, assay = "RNA",features = c("TGFB1", "PDGFB", "OSM", "SPP1", "LGALS3", "THBS1", "FN1", "MMP9", "TIMP1", "CCL2", "IL1B", "IL6") ,pt.size = 0)

Idents(combined_macmod3) <- "LoY"
test <- subset(combined_macmod3, subset = nFeature_RNA > 800)

VlnPlot(test2, assay = "RNA",features = c("FN1","SPP1","TGFB1") ,pt.size = 0)

###Heatmap
markers <- FindAllMarkers(combined_macmod3, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

combined_macmod3_heat <- ScaleData(combined_macmod3, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_macmod3_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=60)


###Proportion calculation
#Idents(combined_macmod3) <- "cell_annotation_final"

# 1. Create a count table of cluster × state
tab <- table(
  cluster = Idents(combined_macmod3),   # or seu$seurat_clusters
  state   = combined_macmod3$LoY      # "A" / "B"
)

tab <- table(
  cluster = combined_macmod3$sample,   # or seu$seurat_clusters
  state   = combined_macmod3$LoY      # "A" / "B"
)

tab

# Row-wise proportions within each cluster
prop_tab <- prop.table(tab, margin = 1)

# Extract only the proportion of A in each cluster
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



combined_macmod3 <- RenameIdents(
  combined_macmod3,
  `0` = "IL1B+ inflammatory mac",
  `1` = "LYVE1+ resident mac",
  `2` = "TNF+ mac",
  `3` = "TREM2+ foamy mac",
  `4` = "classical mono"
)

combined_macmod3 <- AddMetaData(combined_macmod3, Idents(combined_macmod3), "cell_annotation")

##Transfer annotation from combined_macmod3 to combined_macmod2
combined_macmod3$cell_annotation <- as.character(combined_macmod3$cell_annotation)

# 1) Extract only the cell names in the subset
cells_sub <- colnames(combined_macmod3)

# 2) Add a column to the metadata of the original object
combined_macmod2$cell_annotation <- NA

# 3) Transfer the subset annotations to obj
combined_macmod2$cell_annotation[cells_sub] <- combined_macmod3$cell_annotation

# Save the original clusters (coarse clusters)
combined_macmod2$cluster_global <- Idents(combined_macmod2)

# Label column to be used for final visualization and analysis
combined_macmod2$cell_annotation_final <- as.character(combined_macmod2$cluster_global)

# Overwrite only the macrophage subset with subcluster names
idx_mac <- !is.na(combined_macmod2$cell_annotation)
combined_macmod2$cell_annotation_final[idx_mac] <- combined_macmod2$cell_annotation[idx_mac]

# First convert to a factor and check levels
combined_macmod2$cell_annotation_final <- factor(combined_macmod2$cell_annotation_final)
levels(combined_macmod2$cell_annotation_final)

library(dplyr)

combined_macmod2$cell_annotation_final <- recode(
  combined_macmod2$cell_annotation_final,
  `7` = "Neutrophils",
  # Subcluster names remain unchanged unless explicitly specified
)

DimPlot(combined_macmod2, reduction = "umap", split.by = "cell_annotation_final")

##Transfer annotation from combined_macmod2 to combined_macmod
DimPlot(combined_macmod, reduction = "umap", label = TRUE, repel = TRUE)

combined_macmod2$cell_annotation_final <- as.character(combined_macmod2$cell_annotation_final)

# 1) Extract only the cell names in the subset
cells_sub <- colnames(combined_macmod2)

# 2) Add a column to the metadata of the original object
combined_macmod$cell_annotation <- NA

# 3) Transfer the subset annotations to obj
combined_macmod$cell_annotation[cells_sub] <- combined_macmod2$cell_annotation_final

# Save the original clusters (coarse clusters)
combined_macmod$cluster_global <- Idents(combined_macmod)

# Label column to be used for final visualization and analysis
combined_macmod$cell_annotation_final <- as.character(combined_macmod$cluster_global)

# Overwrite only the macrophage subset with subcluster names
idx_mac <- !is.na(combined_macmod$cell_annotation)
combined_macmod$cell_annotation_final[idx_mac] <- combined_macmod$cell_annotation[idx_mac]

# First convert to a factor and check levels
combined_macmod$cell_annotation_final <- factor(combined_macmod$cell_annotation_final)
levels(combined_macmod$cell_annotation_final)

library(dplyr)

combined_macmod$cell_annotation_final <- recode(
  combined_macmod$cell_annotation_final,
  `4` = "cDC2",
  `5` = "Cardiomyocytes"
  # Subcluster names remain unchanged unless explicitly specified
)

DimPlot(combined_macmod, reduction = "umap", split.by = "cell_annotation_final")
DimPlot(combined_macmod, reduction = "umap")

##Transfer annotation from combined_macmod to combined_Mye
DimPlot(combined_Mye, reduction = "umap", label = TRUE, repel = TRUE)

combined_macmod$cell_annotation_final <- as.character(combined_macmod$cell_annotation_final)

# 1) Extract only the cell names in the subset
cells_sub <- colnames(combined_macmod)

# 2) Add a column to the metadata of the original object
combined_Mye$cell_annotation <- NA

# 3) Transfer the subset annotations to obj
combined_Mye$cell_annotation[cells_sub] <- combined_macmod$cell_annotation_final

# Save the original clusters (coarse clusters)
combined_Mye$cluster_global <- Idents(combined_Mye)

# Label column to be used for final visualization and analysis
combined_Mye$cell_annotation_final <- as.character(combined_Mye$cluster_global)

# Overwrite only the macrophage subset with subcluster names
idx_mac <- !is.na(combined_Mye$cell_annotation)
combined_Mye$cell_annotation_final[idx_mac] <- combined_Mye$cell_annotation[idx_mac]

# First convert to a factor and check levels
combined_Mye$cell_annotation_final <- factor(combined_Mye$cell_annotation_final)
levels(combined_Mye$cell_annotation_final)

library(dplyr)
#3: Neutrophils / granulocytes
#4: CD16⁺ non-classical monocytes
#5: Cytotoxic T/NK contamination group
#6: Perivascular fibroblast / SMC / pericyte-like
#7: CLEC9A⁺ cDC1 dendritic cells



combined_Mye$cell_annotation_final <- recode(
  combined_Mye$cell_annotation_final,
  `3` = "Neutrophils",
  `4` = "non classical mono",
  `5` = "Tcell contamination to Mye",
  `6` = "Doublets",
  `7` = "cDC1"
  # Subcluster names remain unchanged unless explicitly specified
)

DimPlot(combined_Mye, reduction = "umap", split.by = "cell_annotation_final")
DimPlot(combined_Mye, reduction = "umap", split.by = "LoY")
DimPlot(combined_Mye, reduction = "umap", split.by = "sample")

combined_Mye <- AddMetaData(combined_Mye, Idents(combined_Mye),col.name = "cluster_Mye")

Idents(combined_Mye) <- "cell_annotation_final"
DimPlot(combined_Mye, reduction = "umap")

combined_Mye <- subset(combined_Mye, idents = c("Tcell contamination to Mye","Doublets", "Cardiomyocytes"), invert=TRUE)

saveRDS(combined_Mye,"~/Desktop/Loss_of_Y_analysis/RDS/With_public_data/All_male_sample_Myeloid_251208_non_doublets.rds")

table(combined_Mye$sample)

DefaultAssay(combined_Mye) <- "RNA"
combined_Mye <- NormalizeData(combined_Mye, verbose = FALSE)
combined_Mye <- FindVariableFeatures(combined_Mye, selection.method = "vst", nfeatures = 12000)
combined_Mye <- ScaleData(combined_Mye, verbose = FALSE)
#combined_Mye <- SCTransform(combined_Mye, vars.to.regress = c("percent.mt"))
#combined_Mye <- RunPCA(combined_Mye, assay = "SCT", verbose = FALSE)
combined_Mye <- RunPCA(combined_Mye, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_Mye <- RunHarmony(combined_Mye, group.by.vars = "sample", assay.use = "RNA",
                           reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_Mye, ndims = 50, reduction = "harmony")
combined_Mye <- FindNeighbors(combined_Mye, reduction = "harmony", dims = 1:7)
combined_Mye <- FindClusters(combined_Mye, resolution = 0.45)
combined_Mye <- RunUMAP(combined_Mye, reduction = "harmony", dims = 1:7)

DimPlot(combined_Mye, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(combined_Mye, reduction = "umap", split.by = "LoY")
DimPlot(combined_Mye, reduction = "umap", split.by = "sample")

Idents(combined_Mye) <- "cell_annotation_final"
DimPlot(combined_Mye, reduction = "umap", label = TRUE, repel = TRUE)

combined_Mye <- AddMetaData(combined_Mye, paste(Idents(combined_Mye),combined_Mye$LoY,sep = "_"),col.name = "cluster_LoY")

table(combined_Mye$sample)

##heatmap
markers <- FindAllMarkers(combined_Mye, assay = "RNA", only.pos = TRUE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

combined_Mye_heat <- ScaleData(combined_Mye, assay = "RNA", features = top10$gene)
p <- DoHeatmap(combined_Mye_heat, assay = "RNA", features = top10$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p

print(top10, n=120)



#Dotplot
combined_Mye_trans <- SCTransform (combined_Mye)

cd_genes <- c("CD14","CD68","CSF3R","FCGR3B","CXCR2","CSF1R","FCGR3A","HLA-DQA1","HLA-DQB1","HLA-DRB5","CLEC9A","CLEC10A","S100A8","S100A12","VCAN","LYZ","IL1B","AREG","EREG","HBEGF","LYVE1","FOLR2","MRC1","TNF","CCL3","CCL4","TREM2","SPP1","APOE")
DotPlot(combined_Mye_trans,features = cd_genes)+RotatedAxis()+coord_flip()

VlnPlot(combined_Mye, assay = "RNA",features = c("CSF3R","FCGR3B","S100A8","CXCR2","NAMPT") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA",features = c("CSF1R","FCGR3A") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA", features = c("HLA-DQA1","HLA-DQB1","HLA-DRB5","HLA-DPB1","CLEC10A","CLEC9A") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA",features = c("S100A12","VCAN") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA",features = c("IL1B","AREG","EREG","HBEGF","TGFB1","CCR2") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA",features = c("LYVE1","FOLR2","IGF1","MRC1") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA",features = c("TNF","CCL3","CCL4","NAMPT") ,pt.size = 0,split.by = "LoY")
VlnPlot(combined_Mye, assay = "RNA",features = c("TREM2","SPP1","LPL","APOE") ,pt.size = 0)
VlnPlot(combined_Mye, assay = "RNA",features = c("TLR4","RAGE") ,pt.size = 0)


DimPlot(combined_Mye,split.by = "LoY")
FeaturePlot(combined_Mye, features=c('RETN'), min.cutoff=0, max.cutoff='q90')
FeaturePlot(combined_Mye, features=c('MS4A7'), min.cutoff=0, max.cutoff='q90')

table(combined_Mye$cell_annotation_final)
table(combined_Mye$sample)



#####Final myeloid volcano plot
DimPlot(combined_Mye, reduction = "umap")
table(combined_Mye$sample)

test <- subset(combined_Mye, idents = "IL1B+ inflammatory mono / mac")
test <- subset(combined_Mye, idents = "TREM2+ foamy mac")
test <- subset(combined_Mye, idents = "Neutrophils")
test <- subset(combined_Mye, idents = "LYVE1+ mac")
test <- combined_Mye
test <- subset(combined_Mye, idents = c("Neutrophils","cDC1","cDC2","non classical mono","classical mono"),invert=TRUE)
test <- subset(combined_Mye, idents = c("IL1B+ inflammatory mono / mac","classical mono","TNF+ inflammatory mac"))
test <- subset(combined_Mye, idents = c("non classical mono"))
test <- subset(combined_Mye, idents = c("LYVE1+ mac","TREM2+ mac"))

DimPlot(test, reduction = "umap")

DefaultAssay(test) <- "RNA"
table(test$LoY)

test_LoY <- subset(test, subset = LoY == "LoY")
test_no_LoY <- subset(test, subset = LoY == "no_LoY")

VlnPlot(test_LoY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(test_no_LoY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


genes_to_keep <- setdiff(rownames(test), c(gene_list))

# Subset in the same way
test[["RNA"]]@counts     <- test[["RNA"]]@counts[genes_to_keep, ]
test[["RNA"]]@data       <- test[["RNA"]]@data[genes_to_keep, ]
#test[["RNA"]]@scale.data <- test[["RNA"]]@scale.data[genes_to_keep, ]

Idents(test) <- "LoY"

## Volcano My0 and My1 Ident1 = right blue My0

Aur.table_chipvsnon <- FindMarkers(test, ident.1 = c("LoY"), ident.2 =c("no_LoY"), verbose = FALSE, logfc.threshold = 0)

Aur.table_chipvsnon$logp <- -log10(Aur.table_chipvsnon$p_val)

Aur.table_chipvsnon_filtered_left = subset(Aur.table_chipvsnon, (logp>=1 & avg_log2FC <= -0.34) | (logp>=5 & avg_log2FC <= -0.25))
Aur.table_chipvsnon_filtered_right = subset(Aur.table_chipvsnon, (logp>=4 & avg_log2FC >= 0.15) | (logp>=2 & avg_log2FC >= 0.35))

Aur.table_chipvsnon_filtered_left = subset(Aur.table_chipvsnon, (logp>=1 & avg_log2FC <= -0.7) | (logp>=2.6 & avg_log2FC <= -0.25))
Aur.table_chipvsnon_filtered_right = subset(Aur.table_chipvsnon, (logp>=1.0 & avg_log2FC >= 0.25) | (logp>=2.5 & avg_log2FC >= 0.15))

genes.to.label.left <- rownames(Aur.table_chipvsnon_filtered_left)
genes.to.label.right <- rownames(Aur.table_chipvsnon_filtered_right)

p1 <- ggplot(Aur.table_chipvsnon, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

## Volcano plot with white background for gene names
# Only needs to be run once
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
    size = 6.5,
    fill  = "white",   # White background
    box.padding   = 0.3,
    point.padding = 0.2,
    label.size    = 0
  ) +
  geom_label_repel(
    data = subset(Aur.table_chipvsnon, gene %in% genes.to.label.left),
    aes(label = gene),
    color = "blue",
    size = 6.5,
    fill  = "white",
    box.padding   = 0.3,
    point.padding = 0.2,
    label.size    = 0
  )
p1

####GSEA
# Aur.table_chipvsnon: FindMarkers result (gene names are row names)
res <- Aur.table_chipvsnon
res$gene <- rownames(res)

# Handle differences in logFC column names
lfc_col <- if ("avg_log2FC" %in% colnames(res)) "avg_log2FC" else "avg_logFC"

# p: adjusted p-values are recommended (clip to avoid zeros)
p <- res$p_val
p[is.na(p)] <- 1
p[p == 0] <- .Machine$double.xmin

# Signed ranking score (for GSEA)
res$rank_score <- res[[lfc_col]] * (-log10(p))

ranks <- res$rank_score
names(ranks) <- res$gene

# If there are duplicated genes, use the one with the largest absolute value
ranks <- tapply(ranks, names(ranks), function(x) x[which.max(abs(x))])

# Sort (fgsea usually expects descending order)
ranks <- sort(ranks, decreasing = TRUE)

#install.packages("msigdbr")
library(msigdbr)
library(fgsea)

# Hallmark (H) example
h <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(h$gene_symbol, h$gs_name)

fg <- fgsea(pathways = pathways, stats = ranks, nperm = 10000)

# Display significant results at the top
fg <- fg[order(fg$padj, -abs(fg$NES)), ]
head(fg, 20)

library(ggplot2)

top <- fg[fg$padj < 0.05, ][1:20, ]
ggplot(top, aes(x=reorder(pathway, NES), y=NES)) +
  geom_col() +
  coord_flip() +
  labs(x=NULL, y="NES", title="GSEA (fgsea) LoY vs noLoY") +
  theme_bw()


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
# Compare clusters
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
    oob = squish   # Values greater than 0.05 are capped at 0.05
  ) +
  theme(axis.text.y = element_text(size = 9, lineheight = 0.7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))


any("organelle fission" == cgBP@compareClusterResult[["Description"]])

# GO:MF

cgMF <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="MF",OrgDb='org.Hs.eg.db')
dotplot(cgMF,showCategory = 10) + 
  theme(axis.text.y = element_text(size = 10, lineheight = 0.9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# GO:CC
cgCC <- compareCluster(geneCluster = genelist, fun = enrichGO, ont="CC",OrgDb='org.Hs.eg.db')
dotplot(cgCC,showCategory = 10) + 
  theme(axis.text.y = element_text(size = 10, lineheight = 0.9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))


###barplot

library(enrichplot)
library(dplyr)

# compareClusterResult -> data.frame
df <- as.data.frame(cgBP)

# Example: type0 only
df0 <- df %>%
  filter(Cluster == "type0") %>%
  arrange(p.adjust) %>%          # Order by significance
  head(15) %>%                    # Equivalent to showCategory
  mutate(Description = factor(Description, levels = rev(Description)),
         mlog10 = -log10(p.adjust))

ggplot(df0, aes(x = Description, y = mlog10)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = NULL, y = "-log10(adj.P)", title = "GO:BP (type0)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, lineheight = 0.7),
        axis.text.x = element_text(size = 8))

-
  
  ### FeaturePlot differences between LoY and no_LoY
  combined_Mye_LoY <- subset(combined_Mye, subset = LoY == c("LoY")) 
combined_Mye_no_LoY <- subset(combined_Mye, subset = LoY == c("no_LoY")) 

VlnPlot(combined_Mye, assay = "RNA",features = c("S100A6","S100A8","S100A9","S100A12","SERPINB2") ,pt.size = 0, split.by = "LoY")

FeaturePlot(combined_Mye_LoY, features=c('S100A12'), min.cutoff=0, max.cutoff='q90')

### Check genes registered to a GO term
library(org.Hs.eg.db)
library(AnnotationDbi)

go_id <- "GO:0030593"  # Example: neutrophil chemotaxis

res <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = go_id,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "GOALL"    # or "GO"
)

# Genes annotated to this GO term (or its descendant terms)
unique(res$SYMBOL)



### Remove neutrophils
combined_Mye_final <- subset(combined_Mye, idents = c("Neutrophils"),invert=TRUE)

DefaultAssay(combined_Mye_final) <- "RNA"
combined_Mye_final <- NormalizeData(combined_Mye_final, verbose = FALSE)
combined_Mye_final <- FindVariableFeatures(combined_Mye_final, selection.method = "vst", nfeatures = 2000)
combined_Mye_final <- ScaleData(combined_Mye_final, verbose = FALSE)
#combined_Mye_final <- SCTransform(combined_Mye_final, vars.to.regress = c("percent.mt"))
#combined_Mye_final <- RunPCA(combined_Mye_final, assay = "SCT", verbose = FALSE)
combined_Mye_final <- RunPCA(combined_Mye_final, assay = "RNA", verbose = FALSE)
library(Rcpp)
library(harmony)
combined_Mye_final <- RunHarmony(combined_Mye_final, group.by.vars = "sample", assay.use = "RNA",
                                 reduction.use = "pca",plot_convergence = TRUE)
ElbowPlot(combined_Mye_final, ndims = 50, reduction = "harmony")
combined_Mye_final <- FindNeighbors(combined_Mye_final, reduction = "harmony", dims = 1:10)
combined_Mye_final <- FindClusters(combined_Mye_final, resolution = 0.1)
combined_Mye_final <- RunUMAP(combined_Mye_final, reduction = "harmony", dims = 1:10)

DimPlot(combined_Mye_final, reduction = "umap", label = FALSE, repel = FALSE,split.by = "LoY")
DimPlot(combined_Mye_final, reduction = "umap", split.by = "annotation")
DimPlot(combined_Mye_final, reduction = "umap", split.by = "sample")

Idents(combined_Mye_final) <- "cell_annotation_final"

combined_Mye_final_trans <- SCTransform (combined_Mye_final)

cd_genes <- c("CD14","CD68","CSF1R","FCGR3A","HLA-DQA1","HLA-DQB1","HLA-DRB5","CLEC9A","CLEC10A","S100A8","S100A12","VCAN","LYZ","IL1B","AREG","EREG","HBEGF","LYVE1","FOLR2","MRC1","TNF","CCL3","CCL4","TREM2","SPP1","APOE")
DotPlot(combined_Mye_final_trans,features = cd_genes)+RotatedAxis()+coord_flip()

saveRDS(combined_Mye_final,"~/Desktop/Loss_of_Y_analysis/RDS/With_public_data/All_male_sample_Myeloid_260305_non_doublets_final4.rds")

## Create LoY proportion plot
table(Idents(combined_Mye))
combined_Mye$
  combined_Mye <- AddMetaData(combined_Mye, paste(Idents(combined_Mye)),col.name = "cluster_Mye")


id <- 1:ncol(combined_Mye)
name <- Idents(combined_Mye)
LoY <- combined_Mye$LoY

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

df_sum$name2 <- factor(df_sum$name, levels = c("TREM2+ mac","TNF+ inflammatory mac","LYVE1+ mac","IL1B+ inflammatory mono / mac","cDC2","cDC1","non classical mono"))  # ← reverse the order


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

## Create LoY proportions by sample

table(Idents(combined_Mye))
table(combined_Mye$sample_LoY)

combined_Mye <- AddMetaData(combined_Mye, paste(combined_Mye$sample, combined_Mye$LoY),col.name = "sample_LoY")


id <- 1:ncol(combined_Mye)
name <- combined_Mye$sample
LoY <- combined_Mye$LoY

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

df_sum$name2 <- factor(df_sum$name, levels = c("TREM2+ mac","TNF+ inflammatory mac","LYVE1+ mac","IL1B+ inflammatory mono / mac","cDC2","cDC1","non classical mono"))  # ← reverse the order


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



###PHATE analysis
# On r-reti2, run pip install phate
library(reticulate)

use_condaenv("r-reti2")

if (!suppressWarnings(require(devtools))) install.packages("devtools")
reticulate::py_install("phate", pip=TRUE)
devtools::install_github("KrishnaswamyLab/phateR")

library(phateR)

Idents(seurat_obj) <- "cell_annotation_final"
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)


seurat_obj <- subset(seurat_obj,idents = c("classical mono","IL1B+ inflammatory mac"), invert=TRUE)

DefaultAssay(seurat_obj) <- "RNA"
set.seed(123)
#seurat_obj <- NormalizeData(seurat_obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_obj <- FindVariableFeatures(seurat_obj, assay = "RNA", selection.method = "vst", nfeatures = 2000)
#seurat_obj <- ScaleData(seurat_obj, assay = "RNA", features = rownames(seurat_obj))

#eml_qc <- RunPCA(seurat_obj, assay = "RNA", features = VariableFeatures(seurat_obj), npcs = 30)
#pca_coords <- seurat_obj@reductions$pca@cell.embeddings
pca_coords <- seurat_obj@reductions$harmony@cell.embeddings
#pca_coords <- seurat_obj@assays$RNA$counts
pca_coords <- as.matrix(pca_coords)

#pca_coords <- t(pca_coords)

dim(pca_coords)
phate_emb <- phate(pca_coords)
#phate_emb <- phate(pca_coords, knn=3, decay=600, t=400, init=phate_emb)
phate_emb <- phate(pca_coords, knn=10,t=10, init=phate_emb)
#phate_emb <- phate(pca_coords, t=19, init=phate_emb)

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
FeaturePlot(seurat_obj,reduction = "phate", features = 'IL1B')

Idents(seurat_obj) <- "LoY"
Idents(seurat_obj) <- "cell_annotation_final"

levels(seurat_obj)


DimPlot(seurat_obj, reduction = "phate",cells.highlight = WhichCells(seurat_obj, idents = c("LoY")))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot")

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(phateR)
library(ggplot2)

# 1) Convert to SCE (cells are columns, genes are rows)
sce <- as.SingleCellExperiment(seurat_obj)

# 2) Store the Harmony embedding in reducedDim
X_harmony <- Embeddings(seurat_obj, "harmony")[, 1:50]

# Match colnames(sce) and rownames(X_harmony) (including order)
X_harmony <- X_harmony[colnames(sce), , drop = FALSE]

# reducedDim must have rows = number of cells (ncol(sce))
reducedDim(sce, "HARMONY") <- X_harmony

# 3) slingshot (use the SCE-side cluster as clusterLabels)
# Use Seurat Idents as they are
colData(sce)$cluster <- Idents(seurat_obj)[colnames(sce)]

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim = "HARMONY",
  start.clus = "IL1B+ inflammatory mono / mac"   # ← change to the less differentiated cluster
)

pt <- slingPseudotime(sce)[, 1]
seurat_obj$SLING <- pt

# 3) PHATE (display)
ph <- phate(X_harmony, knn=10,t=10)
ph_df <- as.data.frame(ph$embedding)
colnames(ph_df) <- c("PHATE1","PHATE2")
ph_df$SLING <- pt

ggplot(ph_df, aes(PHATE1, PHATE2, color = SLING)) +
  geom_point(size = 0.4, alpha = 0.8) +
  theme_classic() +
  labs(title="PHATE colored by Slingshot pseudotime", color="pseudotime")

###Rename labels
Idents(combined_Mye) <- "cell_annotation_final"

combined_Mye$cell_annotation_final <- recode(
  combined_Mye$cell_annotation_final,
  `LYVE1+ resident mac` = "LYVE1+ mac"
  # Subcluster names remain unchanged unless explicitly specified
)

table(combined_Mye$cell_annotation_final)

combined_Mye <- AddMetaData(combined_Mye, paste(Idents(combined_Mye), combined_Mye$LoY),col.name = "cluster_LoY")
table(combined_Mye$cluster_LoY)
