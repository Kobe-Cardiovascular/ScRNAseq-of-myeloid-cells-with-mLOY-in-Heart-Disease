library(Seurat)
library(Matrix)
library(edgeR)
library(dplyr)


library(Seurat)
library(Matrix)
library(edgeR)
library(dplyr)
library(tibble)
library(ggplot2)

#AF
seu <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/With_public_data/All_male_sample_Myeloid_260305_non_doublets_final4.rds")
#DCM
seu <- readRDS("~/Desktop/Loss_of_Y_analysis/RDS/DCM/All_male_DCM_sample_Myeloid_cell_gene3000_260414.rds")

DimPlot(seu)

seu <- subset(seu, subset = nCount_RNA > 3000 & nFeature_RNA < 5000 & percent.mt < 10)

genes_to_keep <- setdiff(rownames(seu), c(gene_list))

seu[["RNA"]]@counts <- seu[["RNA"]]@counts[genes_to_keep, ]
seu[["RNA"]]@data   <- seu[["RNA"]]@data[genes_to_keep, ]

DefaultAssay(seu) <- "RNA"

seu$group <- factor(seu$LoY, levels = c("no_LoY", "LoY"))

seu <- AddMetaData(seu, Idents(seu), "celltype")

seu <- subset(
  seu,
  idents = c(
    "non classical mono",
    "IL1B+ inflammatory mono / mac",
    "TNF+ inflammatory mac",
    "LYVE1+ mac",
    "TREM2+ mac"
  )
)

seu <- subset(
  seu,
  idents = c(
    "IL1B+ inflammatory monos / macs",
    "TNF+ inflammatory macs",
    "LYVE1+ macs",
    "TREM2+ macs"
  )
)

seu <- subset(
  seu,
  idents = c(
    "IL1B+ inflammatory mono / mac",
    "TNF+ inflammatory mac",
    "LYVE1+ mac",
    "TREM2+ mac"
  )
)

seu <- subset(
  seu,
  idents = c(
    "IL1B+ inflammatory monos / macs"
  )
)

seu <- subset(
  seu,
  idents = c(
    "IL1B+ inflammatory mono / mac"
  )
)


seu <- subset(
  seu,
  idents = c(
    "LYVE1+ macs"
  )
)

seu <- subset(
  seu,
  idents = c(
    "LYVE1+ mac"
  )
)

seu <- subset(
  seu,
  idents = c(
    "TREM2+ macs"
  )
)

seu <- subset(
  seu,
  idents = c(
    "TNF+ inflammatory macs"
  )
)



if (!"sample" %in% colnames(seu@meta.data)) {
  if ("Sample" %in% colnames(seu@meta.data)) {
    seu$sample <- seu$Sample
  } else {
    stop("meta.data に sample も Sample もありません。")
  }
}

meta <- seu@meta.data
stopifnot(all(c("sample", "group") %in% colnames(meta)))

library(Seurat)
library(Matrix)
library(edgeR)
library(dplyr)

run_pseudobulk_edger_all <- function(seu,
                                     assay = "RNA",
                                     sample_col = "sample",
                                     group_col = "group",
                                     min_cells_per_sample_group = 20,
                                     min_samples_per_group = 2) {
  
  DefaultAssay(seu) <- assay
  counts <- GetAssayData(seu, assay = assay, slot = "counts")
  meta <- seu@meta.data
  
  stopifnot(all(c(sample_col, group_col) %in% colnames(meta)))
  
  meta$cell_id <- rownames(meta)
  
  sample_group_cell_n <- meta %>%
    dplyr::count(.data[[sample_col]], .data[[group_col]], name = "n_cells")
  
  print(sample_group_cell_n)
  
  valid_pairs <- sample_group_cell_n %>%
    filter(n_cells >= min_cells_per_sample_group)
  
  if (nrow(valid_pairs) == 0) {
    stop("min_cells_per_sample_group を満たす sample×group がありません。")
  }
  
  valid_samples <- valid_pairs %>%
    group_by(.data[[sample_col]]) %>%
    summarise(n_group = n_distinct(.data[[group_col]]), .groups = "drop") %>%
    filter(n_group >= 2) %>%
    pull(.data[[sample_col]])
  
  valid_pairs <- valid_pairs %>%
    filter(.data[[sample_col]] %in% valid_samples)

  group_table <- table(valid_pairs[[group_col]])
  print(group_table)
  
  if (any(group_table < min_samples_per_group)) {
    stop("各群の有効 sample 数が不足しています。")
  }
  
  meta_sub <- meta %>%
    semi_join(valid_pairs, by = c(sample_col, group_col))
  
  if (nrow(meta_sub) == 0) {
    stop("有効な細胞がありません。")
  }
  
  meta_sub$pb_id <- paste(meta_sub[[sample_col]], meta_sub[[group_col]], sep = "__")
  
  pb_ids <- unique(meta_sub$pb_id)
  
  pb_list <- lapply(pb_ids, function(pid) {
    cells_use <- meta_sub$cell_id[meta_sub$pb_id == pid]
    Matrix::rowSums(counts[, cells_use, drop = FALSE])
  })
  
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- pb_ids
  pb_counts <- as.matrix(pb_counts)
  
  # sample metadata
  sample_meta <- meta_sub %>%
    dplyr::select(all_of(c(sample_col, group_col, "pb_id"))) %>%
    dplyr::distinct()
  
  sample_meta <- sample_meta[match(colnames(pb_counts), sample_meta$pb_id), , drop = FALSE]
  rownames(sample_meta) <- sample_meta$pb_id
  
  sample_meta[[sample_col]] <- factor(sample_meta[[sample_col]])
  sample_meta[[group_col]]  <- factor(sample_meta[[group_col]], levels = c("no_LoY", "LoY"))
  
  print(sample_meta)
  
  if (nlevels(sample_meta[[group_col]]) < 2) {
    stop("group が2水準ありません。")
  }
  
  # edgeR
  dge <- DGEList(counts = pb_counts, samples = sample_meta)
  keep <- filterByExpr(dge, group = dge$samples[[group_col]])
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(dge) == 0) {
    stop("filterByExpr を通過した遺伝子がありません。")
  }
  
  dge <- calcNormFactors(dge)
  
  # paired design: sample 
  design <- model.matrix(
    ~ sample_meta[[sample_col]] + sample_meta[[group_col]]
  )
  
  print(design)
  
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design)) 
  
  res <- topTags(qlf, n = Inf)$table
  res$gene <- rownames(res)
  rownames(res) <- NULL
  
  list(
    result = res,
    dge = dge,
    sample_meta = sample_meta,
    sample_group_cell_n = sample_group_cell_n,
    pb_counts = pb_counts,
    design = design
  )
}

pb_all <- run_pseudobulk_edger_all(
  seu = seu,
  assay = "RNA",
  sample_col = "sample",
  group_col = "group",
  min_cells_per_sample_group = 20,
  min_samples_per_group = 2
)

de_table <- pb_all$result %>%
  arrange(FDR)

head(de_table)

pb_FindMarkers_like <- de_table %>%
  transmute(
    gene = gene,
    avg_log2FC = logFC,
    p_val = PValue,
    p_val_adj = FDR
  )

# gene  rownames 
Aur.table_chipvsnon <- pb_FindMarkers_like
rownames(Aur.table_chipvsnon) <- Aur.table_chipvsnon$gene
Aur.table_chipvsnon$gene <- NULL

# =========================
# 7. Volcano plot
# =========================
Aur.table_chipvsnon$logp <- -log10(pmax(Aur.table_chipvsnon$p_val, 1e-300))

Aur.table_chipvsnon_filtered_left <- subset(
  Aur.table_chipvsnon,
  logp >= 1.5 & avg_log2FC <= -0.4
)

Aur.table_chipvsnon_filtered_right <- subset(
  Aur.table_chipvsnon,
  (logp >= 1.5 & avg_log2FC >= 0.4) |
    (logp >= 2.0 & avg_log2FC >= 0.30)
)

genes.to.label.left <- rownames(Aur.table_chipvsnon_filtered_left)
genes.to.label.right <- rownames(Aur.table_chipvsnon_filtered_right)

genes.to.label.right <- c(genes.to.label.right, "LYVE1","RETN","CD47")
genes.to.label.left <- c(genes.to.label.left, "IL1B","NFKB1","NFKBIA","CCL1","CCL3","TNF")
#genes.to.label.right <- c(genes.to.label.right, "LYVE1","S100A8","S100A9","S100A12")
#genes.to.label.right <- c(genes.to.label.right,"S100A8","S100A9","S100A12")
genes.to.label.left <- c(genes.to.label.left, "NFKBIA","CCL8","CCL3","C1QB")


genes.to.label.right <- c(genes.to.label.right, "MS4A7","RETN","TGFB1")

genes.to.label.right <- c("LYVE1","MRC1","FOLR2","RETN")
genes.to.label.left <- c("IL1B","NFKBIA","NFKB1","TNF","CXCL1","CD99","CCL3","C3")

genes.to.label.right <- c("S100A8","S100A9","S100A12")
genes.to.label.left <- c("IL1B","NFKB1","CXCL3","NFKBIA")

genes.to.label.right <- c("LYVE1","FOLR2","RETN")
genes.to.label.left <- c("IL1B","NFKBIA","NFKB1","CCL3","C3")

genes.to.label.right <- c("CXCL1","CXCL3")
genes.to.label.left <- c("CCL8","CCL3","NFKBIA")


Aur.table_chipvsnon$label <- rownames(Aur.table_chipvsnon)

p1 <- ggplot(Aur.table_chipvsnon, aes(avg_log2FC, logp, label = label)) +
  geom_point()

p1 <- LabelPoints(plot = p1, points = genes.to.label.right, color = "red",  repel = TRUE, xnudge = 0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,  color = "blue", repel = TRUE, xnudge = 0)

p1

genes.to.label.left <- c("FABP4","CXCL10","TIMP3","INHBA","CD99","APOE","CD1C","LGALS2","FABP5","TNF","APOC1","CSF2RA","CXCL1","MMP9")

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
    size = 4.5,
    fill  = "white",   
    box.padding   = 0.3,
    point.padding = 0.2,
    label.size    = 0
  ) +
  geom_label_repel(
    data = subset(Aur.table_chipvsnon, gene %in% genes.to.label.left),
    aes(label = gene),
    color = "blue",
    size = 4.5,
    fill  = "white",
    box.padding   = 0.3,
    point.padding = 0.2,
    label.size    = 0
  )
p1

##GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

run_gsea_go <- function(res,
                        gene_col = "gene",
                        logfc_col = "logFC",
                        p_col = "PValue",
                        ont = "BP",
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = 0.5) {
  
  stopifnot(all(c(gene_col, logfc_col, p_col) %in% colnames(res)))
  
  res2 <- res %>%
    filter(!is.na(.data[[gene_col]]),
           !is.na(.data[[logfc_col]]),
           !is.na(.data[[p_col]])) %>%
    mutate(rank_score = .data[[logfc_col]] * (-log10(.data[[p_col]] + 1e-300)))
  
  res2 <- res2 %>%
    group_by(.data[[gene_col]]) %>%
    slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 3. gene symbol -> ENTREZID
  gene_df <- bitr(
    res2[[gene_col]],
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  res3 <- res2 %>%
    inner_join(gene_df, by = setNames("SYMBOL", gene_col))
  
  # 4. ranked vector
  geneList <- res3$rank_score
  names(geneList) <- res3$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  # 5. GSEA
  gsea_res <- gseGO(
    geneList      = geneList,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize,
    pvalueCutoff  = pvalueCutoff,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  
  return(list(
    gsea = gsea_res,
    ranked_table = res3,
    geneList = geneList
  ))
}


gsea_out <- run_gsea_go(pb_out$result, ont = "BP")

head(as.data.frame(gsea_out$gsea))

dotplot(gsea_out$gsea, showCategory = c("neutrophil migration" ,"neutrophil chemotaxis" ,"granulocyte migration","regulation of cell activation","leukocyte chemotaxis","amyloid fibril formation"))

dotplot()

