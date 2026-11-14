#!/usr/bin/env Rscript
# Differential expression with DESeq2 on Salmon output, imported with tximport.
#
#   deseq2.R --quant_dir quant --tx2gene tx2gene.tsv --samplesheet samplesheet.csv \
#            [--design "~ condition"] [--contrast condition,treated,control] \
#            [--gene_names gene_names.tsv] [--alpha 0.05] [--outdir .]
#
# --quant_dir    folder holding one Salmon output folder per sample (<sample>/quant.sf)
# --samplesheet  the pipeline's samplesheet; every column except fastq_1/fastq_2 is
#                sample metadata that the design formula can use
# --design       R model formula. The LAST variable is the one tested by default;
#                put batch/pairing variables first: "~ cell_line + condition"
# --contrast     factor,numerator,denominator -> log2FC = log2(numerator / denominator).
#                Optional when the factor has 2 levels (then the alphabetically first
#                level is the denominator, DESeq2's default); required for 3+ levels.
#
# Uses only DESeq2, tximport and base R (base graphics, no ggplot2), so it runs in a
# container with just those two Bioconductor packages.

suppressPackageStartupMessages({
    library(DESeq2)
    library(tximport)
})

# ---- arguments -------------------------------------------------------------------
parse_args <- function(argv) {
    opts <- list(design = "~ condition", contrast = NA, gene_names = NA, alpha = "0.05", outdir = ".")
    i <- 1
    while (i <= length(argv)) {
        key <- sub("^--", "", argv[i])
        if (i == length(argv)) stop("missing value for --", key)
        opts[[key]] <- argv[i + 1]
        i <- i + 2
    }
    for (req in c("quant_dir", "tx2gene", "samplesheet")) {
        if (is.null(opts[[req]])) stop("--", req, " is required")
    }
    opts$alpha <- as.numeric(opts$alpha)
    opts
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
out <- function(name) file.path(opt$outdir, name)

write_tsv <- function(df, path) {
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# ---- sample table ------------------------------------------------------------------
ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE,
               strip.white = TRUE, colClasses = "character")
coldata <- ss[, setdiff(colnames(ss), c("fastq_1", "fastq_2")), drop = FALSE]
rownames(coldata) <- coldata$sample

design <- as.formula(opt$design)
design_vars <- all.vars(design)
missing_cols <- setdiff(design_vars, colnames(coldata))
if (length(missing_cols)) {
    stop("design uses column(s) not in the samplesheet: ", paste(missing_cols, collapse = ", "))
}
# design variables become factors; DESeq2 treats character columns the same way but warns
for (v in design_vars) coldata[[v]] <- factor(coldata[[v]])

# ---- which comparison ---------------------------------------------------------------
if (!is.na(opt$contrast)) {
    contrast <- trimws(strsplit(opt$contrast, ",")[[1]])
    if (length(contrast) != 3) stop("--contrast must be factor,numerator,denominator; got: ", opt$contrast)
} else {
    factor_name <- tail(design_vars, 1)
    lv <- levels(coldata[[factor_name]])
    if (length(lv) != 2) {
        stop("'", factor_name, "' has ", length(lv), " levels (", paste(lv, collapse = ", "),
             "): give --contrast ", factor_name, ",<numerator>,<denominator>")
    }
    contrast <- c(factor_name, lv[2], lv[1])
    message("No --contrast given: testing ", lv[2], " vs ", lv[1], " (", lv[1],
            " is the reference because it sorts first)")
}
f <- contrast[1]
if (!f %in% design_vars) stop("contrast factor '", f, "' is not in the design ", opt$design)
for (l in contrast[2:3]) {
    if (!l %in% levels(coldata[[f]])) stop("'", l, "' is not a level of ", f, ": ", paste(levels(coldata[[f]]), collapse = ", "))
}
# Make the denominator the reference level, so the model's coefficient matches the contrast
coldata[[f]] <- relevel(coldata[[f]], ref = contrast[3])

# Enough samples to estimate variance? Need more samples than model coefficients.
mm <- model.matrix(design, coldata)
if (nrow(mm) <= qr(mm)$rank) {
    stop("design ", opt$design, " has no residual degrees of freedom with ", nrow(mm),
         " samples: add replicates or simplify the design")
}

# ---- import: transcript-level Salmon -> gene-level counts + length offsets -----------
files <- file.path(opt$quant_dir, coldata$sample, "quant.sf")
names(files) <- coldata$sample
if (!all(file.exists(files))) stop("missing quant.sf for: ", paste(names(files)[!file.exists(files)], collapse = ", "))

tx2gene <- read.delim(opt$tx2gene, header = FALSE, col.names = c("TXNAME", "GENEID"), stringsAsFactors = FALSE)
# tximport sums transcript counts per gene and records each gene's average transcript
# length per sample (weighted by isoform abundance). DESeq2 uses those lengths as offsets,
# so a gene that switches to a longer isoform isn't mistaken for a gene that is more expressed.
# dropInfReps: skip Salmon bootstrap/Gibbs replicates. The pipeline doesn't create them
# (no --numBootstraps) and DESeq2 doesn't use them; without this, tximport tries to read
# meta_info.json for them, which needs the jsonlite package.
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = design)

# Drop genes with almost no reads: they can't be significant and only add to the
# multiple-testing burden. Keep genes with >= 10 reads in at least as many samples
# as the smallest group of the tested factor.
min_group <- min(table(coldata[[f]]))
keep <- rowSums(counts(dds) >= 10) >= min_group
message("Genes kept after low-count filter: ", sum(keep), " of ", length(keep))
dds <- dds[keep, ]

# ---- model ---------------------------------------------------------------------------
# DESeq() = estimateSizeFactors + estimateDispersions + nbinomWaldTest. Run the steps
# separately so the dispersion step can fall back on tiny or noise-free data (the
# synthetic test set has 8 genes), where the trend fit has too few points.
dds <- estimateSizeFactors(dds)
dds <- tryCatch(
    estimateDispersions(dds, quiet = TRUE),
    error = function(e) {
        message("Dispersion trend fit failed (", conditionMessage(e), "); using gene-wise estimates")
        d <- estimateDispersionsGeneEst(dds)
        dispersions(d) <- mcols(d)$dispGeneEst
        d
    })
dds <- nbinomWaldTest(dds, quiet = TRUE)

res <- results(dds, contrast = contrast, alpha = opt$alpha)

# ---- tables ----------------------------------------------------------------------------
tab <- data.frame(gene_id = rownames(res), as.data.frame(res), row.names = NULL, check.names = FALSE)
if (!is.na(opt$gene_names)) {
    gn <- read.delim(opt$gene_names, stringsAsFactors = FALSE)
    tab <- data.frame(gene_id = tab$gene_id, gene_name = gn$gene_name[match(tab$gene_id, gn$gene_id)],
                      tab[, -1], check.names = FALSE)
}
tab <- tab[order(tab$padj, tab$pvalue, na.last = TRUE), ]
write_tsv(tab, out("deseq2_results.tsv"))

norm <- counts(dds, normalized = TRUE)
write_tsv(data.frame(gene_id = rownames(norm), round(norm, 3), check.names = FALSE), out("normalized_counts.tsv"))

# Variance-stabilised values: log2-like scale where variance doesn't depend on the mean,
# which is what PCA and clustering need. blind = TRUE (ignore the design) is the
# recommended setting for QC. vst() subsamples 1000 genes to fit the trend, so tiny
# data falls back to the exact transform, and then to plain log2(normalized + 1).
vsd <- tryCatch({
    if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else varianceStabilizingTransformation(dds, blind = TRUE)
}, error = function(e) {
    message("VST failed (", conditionMessage(e), "); using log2(normalized counts + 1)")
    normTransform(dds)
})
vmat <- assay(vsd)
write_tsv(data.frame(gene_id = rownames(vmat), round(vmat, 4), check.names = FALSE), out("vst_counts.tsv"))

# ---- plots (base graphics -> PDF) --------------------------------------------------------
group <- coldata[[f]]
cols <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(levels(group))],
                 levels(group))

# PCA on the 500 most variable genes (DESeq2's plotPCA default)
top <- head(order(apply(vmat, 1, var), decreasing = TRUE), 500)
pca <- prcomp(t(vmat[top, , drop = FALSE]))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pdf(out("pca.pdf"), width = 6, height = 5)
plot(pca$x[, 1], pca$x[, 2], pch = 19, col = cols[as.character(group)],
     xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
     main = "PCA, variance-stabilised counts")
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.7)
legend("topright", legend = names(cols), col = cols, pch = 19, bty = "n", title = f)
invisible(dev.off())
write_tsv(data.frame(sample = rownames(pca$x), group = group, PC1 = pca$x[, 1], PC2 = pca$x[, 2]), out("pca.tsv"))

pdf(out("ma_plot.pdf"), width = 6, height = 5)
plotMA(res, alpha = opt$alpha, main = paste(contrast[2], "vs", contrast[3]))
invisible(dev.off())

pdf(out("volcano.pdf"), width = 6, height = 5)
ok <- !is.na(tab$padj)
sig <- ok & tab$padj < opt$alpha
plot(tab$log2FoldChange[ok], -log10(tab$pvalue[ok]), pch = 20, cex = 0.5,
     col = ifelse(sig[ok], "#d95f02", "grey60"),
     xlab = paste0("log2 fold change (", contrast[2], " / ", contrast[3], ")"), ylab = "-log10 p-value",
     main = paste0("Volcano: ", sum(sig), " genes with padj < ", opt$alpha))
lab <- head(which(sig), 15)   # tab is sorted by padj, so these are the top hits
if (length(lab)) {
    labels <- if ("gene_name" %in% colnames(tab)) ifelse(is.na(tab$gene_name[lab]), tab$gene_id[lab], tab$gene_name[lab]) else tab$gene_id[lab]
    text(tab$log2FoldChange[lab], -log10(tab$pvalue[lab]), labels = labels, pos = 3, cex = 0.6)
}
invisible(dev.off())

# ---- summary + objects --------------------------------------------------------------------
up   <- sum(sig & tab$log2FoldChange > 0)
down <- sum(sig & tab$log2FoldChange < 0)
summary_lines <- c(
    paste("design:", opt$design),
    paste0("contrast: ", f, " ", contrast[2], " vs ", contrast[3], " (log2FC > 0 = higher in ", contrast[2], ")"),
    paste("samples:", ncol(dds)),
    paste("genes tested:", nrow(dds)),
    paste0("padj < ", opt$alpha, ": ", up + down, " (", up, " up, ", down, " down)"))
writeLines(summary_lines, out("deseq2_summary.txt"))
message(paste(summary_lines, collapse = "\n"))

saveRDS(dds, out("dds.rds"))   # reload in R with readRDS() for further analysis
writeLines(capture.output(sessionInfo()), out("R_sessionInfo.log"))
