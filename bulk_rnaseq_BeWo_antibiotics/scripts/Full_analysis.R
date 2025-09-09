----------------------------------------------------------------------------------------------
## Author: Hassan Saei
## Email: hassan.saeiahan@gmail.com 
## Differential expression analysis: edgeR, DESeq2, limma-voom

# Requirements: Install packages once in your environment
#   Bioc: DESeq2, edgeR, limma, EnsDb.Hsapiens.v86, AnnotationDbi, clusterProfiler, msigdbr
#   CRAN: ggplot2, ggrepel, dplyr, RColorBrewer, ggvenn, pheatmap, readr
---------------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(EnsDb.Hsapiens.v86)
  library(AnnotationDbi)
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
  library(RColorBrewer)
  library(ggvenn)
  library(pheatmap)
  library(readr)
})

rm(list = ls())

# inputs
# Provide correct absolute or relative paths for the data
counts_file      <- "HTseq_merged.counts.txt"  # columns: GENEID + sample count columns
group_file       <- "group.txt"                # columns: sample<TAB>condition
out_dir          <- "HTseq"                    # output directory
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "enrichment"), showWarnings = FALSE, recursive = TRUE)

# ID mapping (ENSG -> SYMBOL, ENTREZ optional)
ens2sym <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86,
  keys = keys(EnsDb.Hsapiens.v86),
  columns = c("SYMBOL", "ENTREZID")
) %>% distinct(GENEID = GENEID, .keep_all = TRUE)

# load data
# group.txt expected columns: sample, condition
group_df <- read.table(group_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("sample", "condition") %in% colnames(group_df)))

counts <- read.table(counts_file, header = TRUE, sep = "\t", check.names = FALSE)
stopifnot("GENEID" %in% colnames(counts))

# Map to SYMBOL and deduplicate by SYMBOL
counts <- counts %>%
  inner_join(ens2sym, by = c("GENEID" = "GENEID")) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct(SYMBOL, .keep_all = TRUE)

rownames(counts) <- counts$SYMBOL
sample_cols <- setdiff(colnames(counts), c("GENEID", "SYMBOL", "ENTREZID"))
counts_mat <- as.matrix(counts[, sample_cols, drop = FALSE])
storage.mode(counts_mat) <- "integer"

# Align group to sample columns
stopifnot(setequal(colnames(counts_mat), group_df$sample))
group_df <- group_df[match(colnames(counts_mat), group_df$sample), , drop = FALSE]
condition <- factor(group_df$condition)
# Set reference level (adjust as needed)
condition <- relevel(condition, ref = "G1")

# edgeR pipeline
dge <- DGEList(counts = counts_mat, group = condition)
dge <- calcNormFactors(dge, method = "TMM")

# keep genes with CPM > 1 in at least 2 samples (tune as needed)
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]

design_edger <- model.matrix(~ 0 + condition)
colnames(design_edger) <- gsub("^condition", "", colnames(design_edger))

dge <- estimateDisp(dge, design_edger)
fit_edger <- glmQLFit(dge, design_edger)

# Define contrasts vs reference G1
contrast_levels <- colnames(design_edger)
stopifnot(all(c("G1", "G2", "G3", "G4", "G5") %in% contrast_levels))
contrasts_edger <- makeContrasts(
  SMZvsCTRL  = G2 - G1,
  TMPvsCTRL  = G3 - G1,
  CTXvsCTRL  = G4 - G1,
  DMSOvsCTRL = G5 - G1,
  levels = design_edger
)

qt_SMZ  <- glmQLFTest(fit_edger, contrast = contrasts_edger[, "SMZvsCTRL"])
qt_TMP  <- glmQLFTest(fit_edger, contrast = contrasts_edger[, "TMPvsCTRL"])
qt_CTX  <- glmQLFTest(fit_edger, contrast = contrasts_edger[, "CTXvsCTRL"])
qt_DMSO <- glmQLFTest(fit_edger, contrast = contrasts_edger[, "DMSOvsCTRL"])

edgeR_SMZ  <- topTags(qt_SMZ,  n = Inf)$table
edgeR_TMP  <- topTags(qt_TMP,  n = Inf)$table
edgeR_CTX  <- topTags(qt_CTX,  n = Inf)$table
edgeR_DMSO <- topTags(qt_DMSO, n = Inf)$table

# Normalized counts (edgeR)
write.table(cpm(dge, normalized.lib.sizes = TRUE),
            file = file.path(out_dir, "Normalized_counts_edgeR.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# QC plots
pdf(file.path(out_dir, "edgeR_BCV.pdf")); plotBCV(dge); dev.off()
pdf(file.path(out_dir, "edgeR_MDS.pdf")); plotMDS(dge, gene.selection = "common", main = "MDS"); dev.off()

# DESeq2 pipeline
colData <- data.frame(row.names = colnames(counts_mat), condition = condition)
dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = colData, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# PCA
vst_data <- vst(dds)
pca <- prcomp(t(assay(vst_data)))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- colData$condition
var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  xlab(paste0("PC1: ", var_exp[1], "%")) +
  ylab(paste0("PC2: ", var_exp[2], "%")) +
  labs(color = "Condition")
ggsave(file.path(out_dir, "DESeq2_PCA.pdf"), p_pca, width = 6, height = 5)

# results vs G1
DESeq2_SMZ  <- as.data.frame(results(dds, contrast = c("condition", "G2", "G1"), alpha = 0.05))
DESeq2_TMP  <- as.data.frame(results(dds, contrast = c("condition", "G3", "G1"), alpha = 0.05))
DESeq2_CTX  <- as.data.frame(results(dds, contrast = c("condition", "G4", "G1"), alpha = 0.05))
DESeq2_DMSO <- as.data.frame(results(dds, contrast = c("condition", "G5", "G1"), alpha = 0.05))

# Normalized counts (DESeq2)
normCount <- as.data.frame(counts(dds, normalized = TRUE))
normCount$SYMBOL <- rownames(normCount)
write.table(normCount, file.path(out_dir, "Normalized_DESeq2.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Dispersion plot
pdf(file.path(out_dir, "DESeq2_Dispersion.pdf")); plotDispEsts(dds); dev.off()

# limma-voom pipeline
v <- voom(dge, design_edger, plot = TRUE)
fit_v <- lmFit(v, design_edger)
fit_v2 <- contrasts.fit(fit_v, contrasts_edger)
fit_v2 <- eBayes(fit_v2)

voom_SMZ  <- topTable(fit_v2, coef = "SMZvsCTRL", number = Inf)
voom_TMP  <- topTable(fit_v2, coef = "TMPvsCTRL", number = Inf)
voom_CTX  <- topTable(fit_v2, coef = "CTXvsCTRL", number = Inf)
voom_DMSO <- topTable(fit_v2, coef = "DMSOvsCTRL", number = Inf)

# Normalized expression (limma-voom)
norm_counts_limma <- as.data.frame(v$E)
norm_counts_limma$SYMBOL <- rownames(norm_counts_limma)
write.csv(norm_counts_limma, file = file.path(out_dir, "Normalized_counts_limma.csv"), row.names = FALSE)

# Save results (not filtered)
write.table(edgeR_SMZ,  file.path(out_dir, "edgeR_SMZvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(edgeR_TMP,  file.path(out_dir, "edgeR_TMPvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(edgeR_CTX,  file.path(out_dir, "edgeR_CTXvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(edgeR_DMSO, file.path(out_dir, "edgeR_DMSOvsCTRL_NotFiltered.txt"), sep = "\t", quote = FALSE, col.names = NA)

write.table(DESeq2_SMZ,  file.path(out_dir, "DESeq2_SMZvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(DESeq2_TMP,  file.path(out_dir, "DESeq2_TMPvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(DESeq2_CTX,  file.path(out_dir, "DESeq2_CTXvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(DESeq2_DMSO, file.path(out_dir, "DESeq2_DMSOvsCTRL_NotFiltered.txt"), sep = "\t", quote = FALSE, col.names = NA)

write.table(voom_SMZ,  file.path(out_dir, "limma_voom_SMZvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(voom_TMP,  file.path(out_dir, "limma_voom_TMPvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(voom_CTX,  file.path(out_dir, "limma_voom_CTXvsCTRL_NotFiltered.txt"),  sep = "\t", quote = FALSE, col.names = NA)
write.table(voom_DMSO, file.path(out_dir, "limma_voom_DMSOvsCTRL_NotFiltered.txt"), sep = "\t", quote = FALSE, col.names = NA)

# DE gene lists (FDR < 0.05)
p_cut <- 0.05
lfc_cut <- 0.5

edgeR_DEGs <- list(
  SMZ  = rownames(edgeR_SMZ[edgeR_SMZ$FDR < p_cut, , drop = FALSE]),
  TMP  = rownames(edgeR_TMP[edgeR_TMP$FDR < p_cut, , drop = FALSE]),
  CTX  = rownames(edgeR_CTX[edgeR_CTX$FDR < p_cut, , drop = FALSE]),
  DMSO = rownames(edgeR_DMSO[edgeR_DMSO$FDR < p_cut, , drop = FALSE])
)

DESeq2_DEGs <- list(
  SMZ  = rownames(DESeq2_SMZ[!is.na(DESeq2_SMZ$padj) & DESeq2_SMZ$padj < p_cut, , drop = FALSE]),
  TMP  = rownames(DESeq2_TMP[!is.na(DESeq2_TMP$padj) & DESeq2_TMP$padj < p_cut, , drop = FALSE]),
  CTX  = rownames(DESeq2_CTX[!is.na(DESeq2_CTX$padj) & DESeq2_CTX$padj < p_cut, , drop = FALSE]),
  DMSO = rownames(DESeq2_DMSO[!is.na(DESeq2_DMSO$padj) & DESeq2_DMSO$padj < p_cut, , drop = FALSE])
)

voom_DEGs <- list(
  SMZ  = rownames(voom_SMZ[voom_SMZ$adj.P.Val < p_cut, , drop = FALSE]),
  TMP  = rownames(voom_TMP[voom_TMP$adj.P.Val < p_cut, , drop = FALSE]),
  CTX  = rownames(voom_CTX[voom_CTX$adj.P.Val < p_cut, , drop = FALSE]),
  DMSO = rownames(voom_DMSO[voom_DMSO$adj.P.Val < p_cut, , drop = FALSE])
)

# Venn
deg_lists_DMSO <- list(
  DESeq2 = DESeq2_DEGs$DMSO,
  edgeR  = edgeR_DEGs$DMSO,
  limma  = voom_DEGs$DMSO
)
pdf(file.path(out_dir, "Venn_DEGs_DMSO.pdf"), width = 6, height = 6)
print(ggvenn(deg_lists_DMSO, fill_color = c("#009E73", "#E69F00", "#56B4E9"),
             stroke_size = 0.6, set_name_size = 5))
dev.off()

# Enrichment (Hallmark + GO BP/CC/MF + KEGG + Reactome)
gs_h   <- msigdbr(species = "Homo sapiens", category = "H")                            %>% select(gs_name, gene_symbol)
gs_bp  <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")       %>% select(gs_name, gene_symbol)
gs_cc  <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")       %>% select(gs_name, gene_symbol)
gs_mf  <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")       %>% select(gs_name, gene_symbol)
gs_kegg<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")     %>% select(gs_name, gene_symbol)
gs_rea <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME") %>% select(gs_name, gene_symbol)

run_enrich <- function(sym, TERM2GENE, title_prefix, out_stub) {
  if (length(sym) == 0) return(invisible(NULL))
  eg <- tryCatch(enricher(gene = sym, TERM2GENE = TERM2GENE), error = function(e) NULL)
  if (is.null(eg)) return(invisible(NULL))
  df <- as.data.frame(eg)
  if (!nrow(df)) return(invisible(NULL))
  write.table(df, file = file.path(out_dir, "enrichment", paste0(out_stub, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(file.path(out_dir, "enrichment", paste0(out_stub, "_dotplot.pdf")), width = 8, height = 6)
  print(dotplot(eg, showCategory = 20, title = title_prefix))
  dev.off()
  invisible(df)
}

run_enrich(DESeq2_DEGs$DMSO, gs_h,    "Hallmarks: DMSO vs CTRL (DESeq2)", "Hallmark_DESeq2_DMSO")
run_enrich(DESeq2_DEGs$DMSO, gs_bp,   "GO BP: DMSO vs CTRL (DESeq2)",     "GOBP_DESeq2_DMSO")
run_enrich(DESeq2_DEGs$DMSO, gs_cc,   "GO CC: DMSO vs CTRL (DESeq2)",     "GOCC_DESeq2_DMSO")
run_enrich(DESeq2_DEGs$DMSO, gs_mf,   "GO MF: DMSO vs CTRL (DESeq2)",     "GOMF_DESeq2_DMSO")
run_enrich(DESeq2_DEGs$DMSO, gs_kegg, "KEGG: DMSO vs CTRL (DESeq2)",      "KEGG_DESeq2_DMSO")
run_enrich(DESeq2_DEGs$DMSO, gs_rea,  "Reactome: DMSO vs CTRL (DESeq2)",  "Reactome_DESeq2_DMSO")

# Volcano plots
# DESeq2 CTX vs CTRL
res_ctx <- DESeq2_CTX
res_ctx$category <- ifelse(!is.na(res_ctx$padj) & res_ctx$padj < p_cut & res_ctx$log2FoldChange >  lfc_cut, "up",
                    ifelse(!is.na(res_ctx$padj) & res_ctx$padj < p_cut & res_ctx$log2FoldChange < -lfc_cut, "down", "ns"))
res_ctx$neglog10padj <- -log10(pmax(res_ctx$padj, .Machine$double.xmin))
p_v_ctx <- ggplot(res_ctx, aes(log2FoldChange, neglog10padj, color = category)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(up = "red", down = "blue", ns = "grey70")) +
  theme_minimal() +
  labs(title = "Volcano: CTX vs CTRL (DESeq2)", x = "log2FC", y = "-log10(adj.P)")
ggsave(file.path(out_dir, "Volcano_DESeq2_CTXvsCTRL.pdf"), p_v_ctx, width = 6, height = 5)

# edgeR DMSO vs CTRL
res_dmso <- edgeR_DMSO
res_dmso$category <- ifelse(res_dmso$FDR < p_cut & res_dmso$logFC >  lfc_cut, "up",
                     ifelse(res_dmso$FDR < p_cut & res_dmso$logFC < -lfc_cut, "down", "ns"))
res_dmso$neglog10p <- -log10(pmax(res_dmso$PValue, .Machine$double.xmin))
p_v_dmso <- ggplot(res_dmso, aes(logFC, neglog10p, color = category)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(up = "red", down = "blue", ns = "grey70")) +
  theme_minimal() +
  labs(title = "Volcano: DMSO vs CTRL (edgeR)", x = "log2FC", y = "-log10(P)")
ggsave(file.path(out_dir, "Volcano_edgeR_DMSOvsCTRL.pdf"), p_v_dmso, width = 6, height = 5)

# heatmap
top_sym <- head(rownames(voom_SMZ[order(abs(voom_SMZ$logFC), decreasing = TRUE), , drop = FALSE]), 50)
mat <- v$E[intersect(top_sym, rownames(v$E)), , drop = FALSE]
if (nrow(mat) >= 2) {
  mat_z <- t(scale(t(mat)))
  colnames(mat_z) <- colnames(v$E)
  pdf(file.path(out_dir, "Heatmap_limma_voom_SMZ_top50.pdf"), width = 8, height = 10)
  pheatmap(mat_z,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           main = "limma-voom: SMZ vs CTRL (top |logFC| genes)", fontsize_row = 6)
  dev.off()
}