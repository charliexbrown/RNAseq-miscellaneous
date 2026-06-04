#!/usr/bin/env Rscript

# Usage:
#   Rscript degs.R config.yaml 
#
# counts.csv   -> Gene expression matrix (rows=genes, columns=samples)
# metadata.csv -> Sample metadata table with columns: Sample, Group (optional: Covariate)


suppressMessages(library(edgeR))
suppressMessages(library(yaml))
suppressMessages(library("pheatmap"))
suppressMessages(library(RColorBrewer))
suppressMessages(library(readr))

# ----------------------
# 0. Read config file
# ----------------------

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  stop("Usage: Rscript config.yaml")
}

config<- yaml::read_yaml(args[1])

counts_file <- config$counts_file
meta_file   <- config$metadata_file
pval_method <- config$pval_method
pval_thr<- as.numeric(config$pval_thr)
logFC_thr<- as.numeric(config$logFC_thr)
out_file<- paste0(basename(counts_file), "_DE_results.csv")

# --------------------
# 1. Load data
# --------------------
cat("Loading counts...\n")
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE, header=TRUE,sep="\t")

cat("Loading metadata...\n")
meta   <- read.delim(meta_file, header=TRUE, sep="\t") 
meta <- type.convert(meta, as.is = FALSE)

if (length(colnames(meta)) == 2) {
    colnames(meta)<-c("Sample", "Condition")
} else {
    colnames(meta) <- c("Sample", "Condition", "Batch")
}

cat("Counts loaded...\n")
cat("Metadata loaded...\n")

## variable check
cat("Checking if covariates are correlated....\n")
summary(lm(meta$Batch~meta$Condition))

# Ensure sample order matches
if(!all(colnames(counts) == meta$Sample)) {
  stop("Sample names in counts matrix must match 'Sample' column in metadata and be in the same order.")
}

# For final plotting: calculate logCPM
counts_logcpm <- cpm(counts, log=TRUE)

# --------------------
# 2. Build design
# --------------------

### MAKE SURE TO SET THE LEVEL YOU WANT AS CONTROL AS THE FIRST CATEGORY IN YOUR METADATA FILE (ie sample from that category must be in the first row of metadata)
cat("Generating design matrix...\n")
control_group <- as.character(meta$Condition[1])
condition <- factor(meta$Condition, levels=c(control_group, setdiff(unique(meta$Condition), control_group)))

if ( "Batch" %in% colnames(meta)) {
     if (is.numeric(meta$Batch)){
           batch<- meta$Batch
     } else {
           batch <- factor(meta$Batch)
     }     
     design <- model.matrix(~ 0 + condition + batch)
} else {
     design <- model.matrix(~0+ condition)
     colnames(design) <- levels(condition)
}

cat("\nDesign matrix:\n")
print(colnames(design))

# --------------------
# 3. edgeR pipeline
# --------------------
cat("\nStarting edgeR pipeline...\n")
y <- DGEList(counts = counts)

# TMM normalization
y <- calcNormFactors(y, method="TMM")

# Filter lowly expressed genes
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

# Estimate dispersions
y <- estimateDisp(y, design)

png(file="MDS.png", width=600, height=350)
plotMDS(y)
dev.off()

# Sanity check of the dispersion: check biological coefficient of variation and set threshold
bcv <- sqrt(y$common.dispersion)
cat("BCV of dispersion is ", round(bcv,3), "\n")

go4qlf <- TRUE
if (bcv > 0.4) {
    cat("\nWARNING: BCV is higher than 0.4, QLF is unstable with this noise, will switch to LRT\n")
    go4qlf <- FALSE
}

# Fit model: try QLF, if fails due to high dispersion try LRT
test_type <- "QLF"

if (go4qlf) {
    fit <- tryCatch({
         glmQLFit(y, design)
    }, error= function(e) {
         go4qlf <<- FALSE
         return(NULL)
    })
}

if (!go4qlf || is.null(fit)) {
    cat("Forcing LRT...\n") 
    test_type <- "LRT"
    fit <- glmFit(y, design)
}


#pick the first level as control
cat("control:\n")
ctrl <- levels(condition)[1]
print(ctrl)

## get all other groups except ctrl
groups<- setdiff(levels(condition), ctrl)
cat("groups\n")
print(groups)

## loop over groups and run DE
for (grp in groups){
    target_level <- as.character(grp)
    ctrl_level <- as.character(ctrl)
    #conds <- levels(condition)
    
    all_coefs <- colnames(design)
    contrast <- rep(0, length(all_coefs)) 
    names(contrast) <- all_coefs
    
    target_idx <- which(all_coefs == target_level | all_coefs == paste0("condition", target_level)) 
    ctrl_idx <- which(all_coefs == ctrl_level | all_coefs == paste0("condition", ctrl_level))
    
    if (length(target_idx)>0 & length(ctrl_idx)>0) {
          contrast[target_idx]<-1
          contrast[ctrl_idx]<- -1 
          

          de_test <- if (test_type == "QLF") glmQLFTest(fit, contrast = contrast) else glmLRT(fit, contrast= contrast)
    } else {
          stop(paste("Contrast generation failed! Check if", target_level, "or", ctrl_level, "exist\n"))
    }    
    
    # filtering significant degs and saving
    res <- topTags(de_test, n = Inf)$table
    res$adj.p <- p.adjust(res$PValue, method = pval_method)
    print(head(res))
    sig <- subset(res, abs(logFC) > logFC_thr & adj.p < pval_thr)
    fname <- paste0(grp, "VS", ctrl, "_", test_type,".csv")
    plotfile <- paste0(grp, "VS", ctrl, "_", test_type, "_heatmap.png")
    
    # generating heatmap and saving results
    if (dim(sig)[1]>1){
        write.csv(sig, fname, row.names = TRUE)
        deg_list<- rownames(sig)
        counts_degs<- counts_logcpm[rownames(counts_logcpm) %in% deg_list, !colnames(counts_logcpm) %in% c("Row.names")]
        annotation_col <- meta[, c("Sample", "Condition")]
        rownames(annotation_col) <- annotation_col$Sample
        annotation_col<- subset(annotation_col, select= -Sample)
          
    
        if (length(deg_list) <80){
               pheatmap(counts_degs, scale="row", width=10, height=8, show_rownames= TRUE, annotation_col= annotation_col, filename=plotfile)
        } else {
               pheatmap(counts_degs, scale="row", width=10, height=8, show_rownames= FALSE, annotation_col= annotation_col, filename=plotfile)
        }
    
    } else {
        write.csv(res, paste0("not_significant_", fname), row.names = TRUE)
        paste("No significant DEGs for comparison", target_level, "vs", ctrl_level, "!\n")    
        next
    } 

}
cat("Done!\n")

