# ----------------------------------------
# Simple enrichment for bulk RNA-seq
# - ORA (GO BP) using MSigDB via msigdbr + clusterProfiler::enricher
# - GSEA (GO CC) directly with SYMBOLs
# - GSEA (KEGG) after mapping SYMBOL -> ENTREZID
# Outputs saved under RNASeqProject_Mehrdad/Enrichment
# ----------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(DOSE)
})

# Organism settings
organism <- "org.Hs.eg.db"   # human OrgDb
kegg_org <- "hsa"            # KEGG organism code
suppressPackageStartupMessages(library(organism, character.only = TRUE))

# Inputs and outputs
deg_file <- "RNASeqProject_Mehrdad/3_DEG_analysis_limma_edgeR_DESeq2/DESeq2_CTXvsCTRL_result_NotFiltered.txt"
out_dir <- "RNASeqProject_Mehrdad/Enrichment"
contrast_label <- "CTX_vs_CTRL"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load DESeq2 results (rownames expected to be SYMBOL)
DEseq <- read.table(deg_file, sep = "\t", header = TRUE, check.names = FALSE)
DEseq$Symbol <- rownames(DEseq)
if (!"log2FoldChange" %in% colnames(DEseq)) stop("Missing 'log2FoldChange' in DESeq2 table.")

# Named vector of fold-changes for GSEA
fc <- DEseq$log2FoldChange
names(fc) <- DEseq$Symbol
fc <- sort(na.omit(fc), decreasing = TRUE)

# ORA (GO BP)
gs_go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, gene_symbol)

ora <- enricher(
  gene = DEseq$Symbol,
  TERM2GENE = gs_go_bp,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

if (!is.null(ora) && nrow(ora@result) > 0) {
  write.table(ora@result,
              file = file.path(out_dir, paste0("ORA_GO_BP_", contrast_label, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(file.path(out_dir, paste0("ORA_GO_BP_dotplot_", contrast_label, ".png")),
         dotplot(ora, showCategory = 10),
         width = 7, height = 5, dpi = 300)
}

# GSEA (GO CC)
gse_go <- gseGO(
  geneList = fc,
  ont = "CC",
  keyType = "SYMBOL",
  nPerm = 10000,
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  verbose = FALSE,
  OrgDb = organism,
  pAdjustMethod = "BH"
)

if (!is.null(gse_go) && nrow(gse_go@result) > 0) {
  write.table(gse_go@result,
              file = file.path(out_dir, paste0("GSEA_GO_CC_", contrast_label, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(file.path(out_dir, paste0("GSEA_GO_CC_dotplot_", contrast_label, ".png")),
         dotplot(gse_go, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign),
         width = 7, height = 5, dpi = 300)
}

# GSEA (KEGG)
ids <- bitr(names(fc), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism)
ids <- ids[!duplicated(ids$SYMBOL), ]
df_map <- DEseq %>%
  filter(Symbol %in% ids$SYMBOL) %>%
  left_join(ids, by = c("Symbol" = "SYMBOL"))

kegg_fc <- df_map$log2FoldChange
names(kegg_fc) <- df_map$ENTREZID
kegg_fc <- sort(na.omit(kegg_fc), decreasing = TRUE)

gse_kegg <- gseKEGG(
  geneList = kegg_fc,
  organism = kegg_org,
  nPerm = 10000,
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  keyType = "ncbi-geneid",
  verbose = FALSE
)

if (!is.null(gse_kegg) && nrow(gse_kegg@result) > 0) {
  write.table(gse_kegg@result,
              file = file.path(out_dir, paste0("GSEA_KEGG_", contrast_label, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggsave(file.path(out_dir, paste0("GSEA_KEGG_dotplot_", contrast_label, ".png")),
         dotplot(gse_kegg, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign),
         width = 7, height = 5, dpi = 300)
}