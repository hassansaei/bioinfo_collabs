## Author: Hassan Saei
## Email: hassan.saeiahan@gmail.com
## Analysis of Salmon quant files + DESeq2 statistical analysis

#  Libraries (pre-install these outside the script)
library(tximport)
library(ensembldb)
library(AnnotationHub)
library(DESeq2)
library(dplyr)
library(ggplot2)

# Paths (edit these 2 to your environment)
quants_dir <- "/Volumes/HassanLaCie/NGSProject_Mehrdad/In vitro/salmon_out"
sample_table_path <- "/Volumes/HassanLaCie/NGSProject_Mehrdad/In vitro/group_salmon.txt"
out_dir <- "/Volumes/HassanLaCie/NGSProject_Mehrdad/In vitro/results_salmon"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Get EnsDb (Ensembl v109) for tx->gene and gene symbols
get_ensdb_v109 <- function() {
  hub <- AnnotationHub()
  q <- query(hub, c("EnsDb", "Homo sapiens", "109"))
  if (length(q) == 0) stop("No EnsDb v109 found in AnnotationHub.")
  q[[1]]
}
ensdb <- get_ensdb_v109()

tx2gene <- transcripts(ensdb, columns = c("tx_id", "gene_id")) |>
  as.data.frame() |>
  dplyr::select(tx_id, gene_id)

gen2symbol <- genes(ensdb, columns = c("gene_id", "gene_name", "gene_biotype")) |>
  as.data.frame() |>
  dplyr::select(gene_id, gene_name, gene_biotype)

# Collect Salmon quant files and sample names
quant_files <- list.files(quants_dir, pattern = "quant\\.sf$", recursive = TRUE, full.names = TRUE)
if (length(quant_files) == 0) stop("No quant.sf found under quants_dir.")

samples <- basename(dirname(quant_files))
names(quant_files) <- samples

# Import sample table and align with quant files
sample_table <- read.table(sample_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
# Accept either rownames-as-samples or a 'sample' column
if ("sample" %in% tolower(names(sample_table))) {
  sample_col <- names(sample_table)[tolower(names(sample_table)) == "sample"][1]
  rownames(sample_table) <- sample_table[[sample_col]]
}
if (!all(samples %in% rownames(sample_table))) {
  stop("Some samples in Salmon output are missing from sample table.")
}
sample_table <- sample_table[samples, , drop = FALSE]
if (!"condition" %in% names(sample_table)) stop("Sample table must have a 'condition' column.")
sample_table$condition <- factor(sample_table$condition)

# tximport and DESeq2
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Set reference level (edit "G1" if needed)
dds$condition <- relevel(dds$condition, ref = "G1")

dds <- DESeq(dds)

# PCA (VST + plot)
vsd <- vst(dds, blind = TRUE)
pca_df <- as.data.frame(prcomp(t(assay(vsd)))$x)
pca_df$condition <- colData(dds)$condition
p <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal()
ggsave(file.path(out_dir, "PCA_vst.png"), p, width = 6, height = 5, dpi = 300)

# Normalized counts (with gene symbols)
norm_counts <- counts(dds, normalized = TRUE) |> as.data.frame()
norm_counts$gene_id <- rownames(norm_counts)
norm_counts <- norm_counts |>
  dplyr::left_join(gen2symbol, by = "gene_id") |>
  dplyr::relocate(gene_id, gene_name, gene_biotype)
write.table(norm_counts, file.path(out_dir, "DESeq2_normalized_counts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Differential expression vs reference
ref <- levels(dds$condition)[1]
groups <- setdiff(levels(dds$condition), ref)

for (g in groups) {
  res <- results(dds, contrast = c("condition", g, ref), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df |>
    dplyr::left_join(gen2symbol, by = "gene_id") |>
    dplyr::relocate(gene_id, gene_name, gene_biotype)

  comp_name <- paste0("DESeq2_", g, "vs", ref, "_salmon")

  # Unfiltered
  write.table(res_df, file.path(out_dir, paste0(comp_name, "_all.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

  # padj < 0.05
  res_f <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, , drop = FALSE]
  write.table(res_f, file.path(out_dir, paste0(comp_name, "_padj_lt_0.05.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
}

