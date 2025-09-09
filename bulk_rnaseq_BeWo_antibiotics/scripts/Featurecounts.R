# Author: Hassan Saei
# Email: hassan.saeiahan2gmail.com
## Analysis of the gene count matrix (featureCounts)

# ---------------------------
# 0) setup
# ---------------------------
options(stringsAsFactors = FALSE)
rm(list = ls())

# Package management
pkgs <- c(
  "BiocManager", "limma", "edgeR", "DESeq2", "ggplot2",
  "EnsDb.Hsapiens.v86", "dplyr", "ggrepel", "tidyr",
  "RColorBrewer", "clusterProfiler", "msigdbr", "readr"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(to_install)) install.packages(setdiff(to_install, "BiocManager"))
bioc_pkgs <- c("limma", "edgeR", "DESeq2", "EnsDb.Hsapiens.v86", "clusterProfiler")
bioc_missing <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(bioc_missing)) BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(limma); library(edgeR); library(DESeq2); library(ggplot2)
  library(EnsDb.Hsapiens.v86); library(dplyr); library(ggrepel)
  library(tidyr); library(RColorBrewer); library(clusterProfiler); library(msigdbr); library(readr)
})

# ---------------------------
# 1) Parameters and I/O
# ---------------------------
# Set absolute paths
project_dir <- "/Volumes/HassanLaCie/NGSProject_Mehrdad/In vitro"
counts_file <- file.path(project_dir, "FeatureCounts_geneid.txt")
group_file  <- file.path(project_dir, "group.txt")

# Output directories
out_dir <- file.path(project_dir, "HTseq")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "enrichment"), showWarnings = FALSE, recursive = TRUE)

# Read inputs with checks
stopifnot(file.exists(counts_file), file.exists(group_file))
htseq_count <- read.table(counts_file, header = TRUE, sep = "\t", check.names = FALSE)
group_df <- read.table(group_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Expect a column named "condition" with sample groups; ensure factor levels once
stopifnot("condition" %in% colnames(group_df))
group_df$condition <- factor(group_df$condition, levels = c("G1","G2","G3","G4","G5"))

# Ensure sample order alignment between counts columns and group rows
sample_cols <- setdiff(colnames(htseq_count), c("GENEID", "SYMBOL"))
stopifnot(all(sample_cols %in% rownames(group_df)))
group_df <- group_df[sample_cols, , drop = FALSE]

# ---------------------------
# 2) Annotation mapping
# ---------------------------
stopifnot("GENEID" %in% colnames(htseq_count))
gene_ids <- htseq_count$GENEID

ens2sym <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86,
  keys = gene_ids,
  keytype = "GENEID",
  columns = c("SYMBOL")
) %>% distinct(GENEID, .keep_all = TRUE)

# Join and tidy counts
htseq_count_anno <- htseq_count %>%
  left_join(ens2sym, by = "GENEID") %>%
  mutate(GENE = ifelse(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, GENEID)) %>%
  select(GENE, all_of(sample_cols), GENEID)

# De-duplicate genes by keeping the row with highest total counts
count_mat <- htseq_count_anno %>% select(all_of(sample_cols))
dup <- duplicated(htseq_count_anno$GENE)
if (any(dup)) {
  htseq_count_anno <- htseq_count_anno %>%
    group_by(GENE) %>%
    slice_max(order_by = rowSums(across(all_of(sample_cols))), n = 1, with_ties = FALSE) %>%
    ungroup()
  count_mat <- htseq_count_anno %>% select(all_of(sample_cols))
}
rownames(count_mat) <- htseq_count_anno$GENE

# ---------------------------
# 3) edgeR analysis
# ---------------------------
dge <- DGEList(counts = count_mat, group = group_df$condition)
dge <- calcNormFactors(dge, method = "TMM")

# Filtering: keep genes expressed at CPM > 0.1 in at least 1 sample (adjust as needed)
keep <- rowSums(cpm(dge) > 0.1) >= 1
dge <- dge[keep, , keep.lib.sizes = FALSE]

design <- model.matrix(~ 0 + group_df$condition)
colnames(design) <- gsub("group_df\\$condition", "", colnames(design))

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Example QL test (G2 vs G1) using contrasts
qlf_12 <- glmQLFTest(fit, contrast = makeContrasts(G2 - G1, levels = design))
edgeR_12 <- topTags(qlf_12, n = Inf)$table

# Normalized CPM
norm_counts_edgeR <- cpm(dge, normalized.lib.sizes = TRUE)
write.table(norm_counts_edgeR, file = file.path(out_dir, "Normalized_counts_edgeR.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# Diagnostics
pdf(file.path(out_dir, "edgeR_BCV.pdf")); plotBCV(dge); dev.off()
pdf(file.path(out_dir, "edgeR_MDS.pdf"))
mds <- plotMDS(dge, plot = FALSE)
plot(mds$x, mds$y, xlab = "Dim 1", ylab = "Dim 2", main = "MDS Plot", pch = 16, col = as.numeric(group_df$condition))
text(mds$x, mds$y, labels = rownames(group_df), cex = 0.7, pos = 3)
legend("topright", legend = levels(group_df$condition), col = seq_along(levels(group_df$condition)), pch = 16, bty = "n")
dev.off()

# exactTest for pairs vs reference G1
ref <- "G1"
pairs <- c("G2","G3","G4","G5")
edgeR_exact <- lapply(pairs, function(g) exactTest(dge, pair = c(ref, g)))
names(edgeR_exact) <- paste0(pairs, "vs", ref)

edgeR_tables <- lapply(edgeR_exact, function(et) topTags(et, n = Inf, adjust.method = "BH", sort.by = "PValue")$table)
for (nm in names(edgeR_tables)) {
  tab <- edgeR_tables[[nm]]
  tab$Pval_threshold <- tab$PValue < 0.05
  write.table(tab, file = file.path(out_dir, paste0("edgeR_", nm, "_NotFiltered.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(tab[tab$PValue < 0.05, ], file = file.path(out_dir, paste0("edgeR_", nm, "_Filtered_Pval_0.05.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
}

# ---------------------------
# 4) DESeq2 analysis
# ---------------------------
dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = group_df,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# PCA (vst) with labeled samples
vst_mat <- vst(dds)
pca_data <- prcomp(t(assay(vst_mat)))
pca_df <- as.data.frame(pca_data$x)
pca_df$condition <- colData(dds)$condition
pca_df$sample <- rownames(pca_df)
var_expl <- pca_data$sdev^2 / sum(pca_data$sdev^2)

p <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3) +
  theme_bw() +
  labs(title = "PCA (vst)", x = paste0("PC1 (", round(var_expl[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(var_expl[2]*100, 1), "%)"))
ggsave(file.path(out_dir, "DESeq2_PCA.pdf"), p, width = 6, height = 5)

# Normalized counts and dispersion plot
normCount <- counts(dds, normalized = TRUE)
write.table(normCount, file = file.path(out_dir, "Normalized_DESeq2.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
pdf(file.path(out_dir, "DESeq2_dispersion.pdf")); plotDispEsts(dds); dev.off()

# Results vs G1
comparisons <- c("G2","G3","G4","G5")
deseq_results <- lapply(comparisons, function(g) results(dds, contrast = c("condition", g, "G1"), alpha = 0.05))
names(deseq_results) <- paste0(comparisons, "vsG1")

for (nm in names(deseq_results)) {
  res <- as.data.frame(deseq_results[[nm]])
  res$SYMBOL <- rownames(res)
  res$Pval_threshold <- !is.na(res$padj) & res$padj < 0.05
  write.table(res, file = file.path(out_dir, paste0("DESeq2_", nm, "_NotFiltered.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(res[res$Pval_threshold, ], file = file.path(out_dir, paste0("DESeq2_", nm, "_Filtered_Padj_0.05.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
}

# ---------------------------
# 5) limma-voom analysis
# ---------------------------
design_voom <- model.matrix(~ 0 + condition, data = group_df)
colnames(design_voom) <- gsub("condition", "", colnames(design_voom))

v <- voom(dge, design_voom, plot = TRUE)
fit <- lmFit(v, design_voom)
contr <- makeContrasts(
  G2vsG1 = G2 - G1,
  G3vsG1 = G3 - G1,
  G4vsG1 = G4 - G1,
  G5vsG1 = G5 - G1,
  levels = design_voom
)
fit2 <- eBayes(contrasts.fit(fit, contr))

voom_tabs <- lapply(colnames(contr), function(co) topTable(fit2, coef = co, number = Inf))
names(voom_tabs) <- colnames(contr)
for (nm in names(voom_tabs)) {
  tab <- voom_tabs[[nm]]
  tab$Pval_threshold <- tab$adj.P.Val < 0.05
  write.table(tab, file = file.path(out_dir, paste0("limma_voom_", nm, "_NotFiltered.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(tab[tab$adj.P.Val < 0.05, ], file = file.path(out_dir, paste0("limma_voom_", nm, "_Filtered_Padj_0.05.txt")),
              sep = "\t", quote = FALSE, col.names = NA)
}
write.csv(v$E, file = file.path(out_dir, "Normalized_counts_limma.csv"), row.names = TRUE)

# ---------------------------
# 6) DEG sets, Venn, and enrichment
# ---------------------------
# Extract tables for G5 vs G1 across methods
edgeR_15 <- edgeR_tables[["G5vsG1"]]
DESeq2_15 <- as.data.frame(deseq_results[["G5vsG1"]])
voom_15   <- voom_tabs[["G5vsG1"]]

edgeR_DEGs <- rownames(edgeR_15[edgeR_15$PValue < 0.05, ])
DESeq2_DEGs <- rownames(DESeq2_15[!is.na(DESeq2_15$padj) & DESeq2_15$padj < 0.05, ])
voom_DEGs <- rownames(voom_15[voom_15$adj.P.Val < 0.05, ])
DEG_lists <- list(DESeq2 = DESeq2_DEGs, edgeR = edgeR_DEGs, limma_voom = voom_DEGs)

# Venn (CRAN ggVennDiagram preferred over installing from GitHub each run)
if (!"ggVennDiagram" %in% rownames(installed.packages())) install.packages("ggVennDiagram")
suppressPackageStartupMessages(library(ggVennDiagram))
pdf(file.path(out_dir, "Venn_G5vsG1.pdf"), width = 6, height = 6)
print(ggVennDiagram(DEG_lists))
dev.off()

# Enrichment using msigdbr TERM2GENE mappings, with gene symbols
gs_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name, gene_symbol)
gs_go_bp <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>% select(gs_name, gene_symbol)
gs_go_cc <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") %>% select(gs_name, gene_symbol)
gs_go_mf <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF") %>% select(gs_name, gene_symbol)
gs_react <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME") %>% select(gs_name, gene_symbol)

ora_plot <- function(genes, TERM2GENE, title, file_stub) {
  if (length(genes) == 0) return(invisible(NULL))
  egmt <- tryCatch(enricher(gene = genes, TERM2GENE = TERM2GENE), error = function(e) NULL)
  if (is.null(egmt)) return(invisible(NULL))
  egdf <- as.data.frame(egmt)
  write.table(egdf, file = file.path(out_dir, "enrichment", paste0(file_stub, ".txt")),
              sep = "\t", quote = FALSE, col.names = NA)
  pdf(file.path(out_dir, "enrichment", paste0(file_stub, "_dotplot.pdf")), width = 9, height = 7)
  print(dotplot(egmt, showCategory = 20, title = title))
  dev.off()
}

ora_plot(edgeR_DEGs, gs_h, "Hallmarks: G5 vs G1 (edgeR)", "Hallmarks_edgeR_G5vsG1")
ora_plot(edgeR_DEGs, gs_go_bp, "GO-BP: G5 vs G1 (edgeR)", "GO_BP_edgeR_G5vsG1")
ora_plot(edgeR_DEGs, gs_go_cc, "GO-CC: G5 vs G1 (edgeR)", "GO_CC_edgeR_G5vsG1")
ora_plot(edgeR_DEGs, gs_go_mf, "GO-MF: G5 vs G1 (edgeR)", "GO_MF_edgeR_G5vsG1")
ora_plot(edgeR_DEGs, gs_react, "Reactome: G5 vs G1 (edgeR)", "Reactome_edgeR_G5vsG1")

# ---------------------------
# 7) Volcano plot (DESeq2 G5 vs G1)
# ---------------------------
res15 <- as.data.frame(deseq_results[["G5vsG1"]])
res15$SYMBOL <- rownames(res15)
res15$padj_safe <- res15$padj
min_nonzero <- suppressWarnings(min(res15$padj_safe[res15$padj_safe > 0], na.rm = TRUE))
if (!is.finite(min_nonzero)) min_nonzero <- 1
res15$padj_safe[is.na(res15$padj_safe) | res15$padj_safe == 0] <- min_nonzero

threshold_pvalue <- 0.05
threshold_lfc <- 1
res15$category <- with(res15, ifelse(padj_safe < threshold_pvalue & log2FoldChange > threshold_lfc, "up",
                             ifelse(padj_safe < threshold_pvalue & log2FoldChange < -threshold_lfc, "down", "non-significant")))

volcano <- ggplot(res15, aes(x = log2FoldChange, y = -log10(padj_safe), color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "non-significant" = "gray")) +
  labs(title = "Volcano: G5 vs G1 (DESeq2)", x = "Log2 Fold Change", y = "-Log10 adj p-value") +
  geom_hline(yintercept = -log10(threshold_pvalue), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-threshold_lfc, threshold_lfc), linetype = "dashed", color = "blue") +
  ggrepel::geom_text_repel(data = subset(res15, padj_safe < 1e-8 & abs(log2FoldChange) > 1),
                           aes(label = SYMBOL), size = 3, box.padding = 0.5, point.padding = 0.5,
                           segment.color = "gray", max.overlaps = Inf)
ggsave(file.path(out_dir, "Volcano_DESeq2_G5vsG1.pdf"), volcano, width = 7, height = 5)

# ---------------------------
# 8) Common DEGs across methods and merged table
# ---------------------------
edgeR_sig <- edgeR_15[edgeR_15$PValue < 0.05, ]
voom_sig  <- voom_15[voom_15$adj.P.Val < 0.05, ]
deseq_sig <- DESeq2_15[!is.na(DESeq2_15$padj) & DESeq2_15$padj < 0.05, ]

edgeR_sig$SYMBOL <- rownames(edgeR_sig)
voom_sig$SYMBOL  <- rownames(voom_sig)
deseq_sig$SYMBOL <- rownames(deseq_sig)

common_genes <- Reduce(intersect, list(edgeR_sig$SYMBOL, voom_sig$SYMBOL, deseq_sig$SYMBOL))

edgeR_common <- edgeR_sig[edgeR_sig$SYMBOL %in% common_genes, c("SYMBOL", "logFC", "PValue")]
voom_common  <- voom_sig[voom_sig$SYMBOL %in% common_genes, c("SYMBOL", "logFC", "adj.P.Val")]
deseq_common <- deseq_sig[deseq_sig$SYMBOL %in% common_genes, c("SYMBOL", "log2FoldChange", "padj")]

merged_results <- Reduce(function(x, y) merge(x, y, by = "SYMBOL", all = FALSE),
                         list(edgeR_common, voom_common, deseq_common))
colnames(merged_results) <- c("SYMBOL", "logFC_edgeR", "PValue_edgeR", "logFC_limma", "adj.P.Val_limma", "log2FC_DESeq2", "padj_DESeq2")

write.table(merged_results, file = file.path(project_dir, "Merged_G5vsG1.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# ---------------------------
# 9) Session info
# ---------------------------
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))
