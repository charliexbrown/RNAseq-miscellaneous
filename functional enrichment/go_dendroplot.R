#!/usr/bin/env Rscript

# ==========================================================================
# Script: dendrogram_plot.R
# Purpose: Generate dendrogram from GO enrichment output
# Usage:   Rscript dendrogram_plot.R input.csv
# Note: input.csv must contain following columns: Term, genes, pval, ObsExp
# ==========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggdendro)
  library(scales)
  library(readr)
  library(dplyr)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript dendrogram_plot.R <input.csv>")
}

input_file <- args[1]
output_file <- basename(input_file) + "dendroplot.png"

# 1. Parse data
cat("[INFO] Loading input file")
if (!file.exists(input_file)) {
  stop(paste("ERROR: Input file does not exist at: ", input_file))
}
df <- read.csv(input_file, stringsAsFactors = FALSE, sep="\t")

# Deal with column naming
col_mapping <- list(
  Term    = c("Term", "term", "GO_Term", "GO_Name", "NAME", "Description"),
  genes   = c("genes", "Genes", "gene_list", "associated_genes", "core_enrichment"),
  FDR     = c("FDR", "fdr", "adj_pvalue", "padj", "qvalue", "p.adjust"),
  ObsExp  = c("Obs.Exp", "ObsExp", "Obs/Exp", "enrichment_score", "enrichmentRatio")
)

for (ok_cols in names(col_mapping)) {
  match_col <- intersect(colnames(df), col_mapping[[ok_col]])
  if (length(mach_col)>0) {
    colnames(df)[colnames(df) == match_col[1]] <- ok_cols
  } else {
    stop(paste("ERROR: check if columns Term, genes, FDR, ObsExp are in df"))
  }
}

df$Term <- as.character(df$Term)
df$genes <- as.character(df$genes)
df$FDR <- as.numeric(df$FDR)
df$ObsExp <- as.numeric(df$ObsExp)

# Convert genes column to list of gene sets assuming format "gene1,gene2,gene3"
df$gene_set<- lapply(df$genes, function(x) {
  cleaned_str <- gsub("[\\{\\}\\[\\]'\"\\s]", "", x, perl = TRUE)
  trimws(strsplit(cleaned_str, ",")[[1]])
})

# 2. Compute GO_MWU similarity 
## The distance is the n of genes shared among the two GO categories within the analyzed dataset divided by the size of the smaller of the two categories.
n <- nrow(df)
similarity <- matrix(0, n, n)

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    inter <- length(intersect(df$gene_set[[i]], df$gene_set[[j]]))
    minsize <- min(length(df$gene_set[[i]]), length(df$gene_set[[j]]))
    sim <- ifelse(minsize > 0, inter / minsize, 0)
    similarity[i, j] <- sim
    similarity[j, i] <- sim
  }
}

diag(similarity) <- 1
distance <- 1 - similarity

# 3. Hierarchical clustering
hc <- hclust(as.dist(distance), method = "complete")

# Set parameters
cutTreeHeight <- 0.25 # merge if >75% shared genes

# 4. Cut tree and merge clusters
df$cluster <- cutree(hc, h = cutTreeHeight) #h=0.25 (if more than 75% genes are shared the terms are merged)


merged <- df %>%
  group_by(cluster) %>%
  summarise(
    Term = Term[which.max(observed_counts)],
    FDR = min(FDR),
    ObsExp = Obs.Exp[which.max(observed_counts)],
    gene_set = list(unique(unlist(gene_set))), ## gene sets of collapsed terms are merged into one set
    merged_terms = paste(unique(Term), collapse = "; "),
    .groups = "drop"
  )

# 5. Recompute distance for merged set
n2 <- nrow(merged)
similarity2 <- matrix(0, n2, n2)
for (i in 1:(n2-1)) {
  for (j in (i+1):n2) {
    inter <- length(intersect(merged$gene_set[[i]], merged$gene_set[[j]]))
    minsize <- min(length(merged$gene_set[[i]]), length(merged$gene_set[[j]]))
    sim <- ifelse(minsize > 0, inter / minsize, 0)
    similarity2[i, j] <- sim
    similarity2[j, i] <- sim
  }
}
diag(similarity2) <- 1
distance2 <- 1 - similarity2
rownames(distance2)<- merged$Term
hc2 <- hclust(as.dist(distance2), method = "complete")

# 6. Create dendrogram data
dhc <- as.dendrogram(hc2)
ddata <- dendro_data(dhc, type = "rectangle")

# adding statistical parameters for aesthetic plotting
if (!is.null(ddata$labels)) {
  ddata$labels <- merge(ddata$labels, merged[, c("Term", "FDR", "ObsExp")], by.x = "label", by.y = "Term", all.x = TRUE)
} else {
  ddata$labels <- data.frame(x = 1, y = 0, label = merged$Term, FDR = merged$FDR, ObsExp = merged$ObsExp)
}

# log-transform p-values
ddata$labels <- ddata$labels %>% 
  dplyr::mutate(neglogFDR = -log10(FDR + 1e-10))

# --- 8. Plot ---

# dynamic height scaling based on terms in df
image_height <- max(6, min(24, n2 * 0.35))

p <- ggplot(segment(ddata)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = ddata$labels, 
            aes(x = x, y = y, label = label, color=neglogFDR, size = ObsExp),
            vjust = 0, 
            hjust=0, #left-aligned
            position= position_nudge(y=0.015)) + 
  coord_flip(clip="off") + 
  scale_y_reverse(expand = expansion(mult=c(0.2, 1.5))) + #extends space to avoid text truncation
  scale_size(range = c(4, 8)) + # scale font size of labels
  theme_dendro()

ggsave(output_file, plot = p, dpi = 300, width= 13, height = image_height)
