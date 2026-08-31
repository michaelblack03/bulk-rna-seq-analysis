# ============================================================
# BULK RNA-SEQ DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================
#
# This script performs downstream analysis of gene- and
# transcript-level count matrices generated from a bulk RNA-seq
# workflow.
#
# Analysis includes:
#   - Count filtering
#   - DESeq2 differential expression analysis
#   - Dispersion and expression QC
#   - rlog transformation
#   - Principal component analysis (PCA)
#   - MA plots
#   - Log2 fold-change shrinkage using ashr
#   - Identification of significantly differentially expressed
#     genes
#   - Exploratory batch-effect correction using limma
#
# Experimental design:
#   9 samples across 3 experimental groups (A, B, C)
#   and 3 batches (1, 2, 3).
#
# Input files:
#   gene_count_matrix.csv
#   transcript_count_matrix.csv
#   Design-Exp.csv
#
# Author: Michael Black
# Institution: The University of Glasgow
# ============================================================


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(DESeq2)
library(ggplot2)
library(limma)
library(vsn)
library(hexbin)
library(svglite)
library(RColorBrewer)
library(genefilter)
library(lattice)
library(ashr)


# ============================================================
# 2. DEFINE PCA PLOTTING FUNCTION
# ============================================================

# Custom PCA plotting function based on the DESeq2 tutorial.
#
# Performs PCA using the most variable features and labels
# individual samples according to their experimental group.

plotPCAWithSampleNames <- function(x,
                                   targets = targets,
                                   intgroup = colnames(targets)[1],
                                   ntop = 500) {
  
  # Calculate variance for each feature.
  rv <- rowVars(x)
  
  # Select the most variable features.
  select <- order(
    rv,
    decreasing = TRUE
  )[seq_len(min(ntop, length(rv)))]
  
  # Perform PCA across samples using the selected features.
  pca <- prcomp(t(x[select, ]))
  
  # Calculate percentage of variance explained by each PC.
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  variance <- round(variance, 3) * 100
  
  # Obtain sample names.
  names <- colnames(x)
  
  # Create factor representing experimental groups.
  fac <- factor(
    apply(
      as.data.frame(
        targets[, intgroup, drop = FALSE]
      ),
      1,
      paste,
      collapse = " : "
    )
  )
  
  # Assign colours according to the number of groups.
  if (nlevels(fac) >= 10) {
    colors <- rainbow(nlevels(fac))
  } else if (nlevels(fac) >= 3) {
    colors <- brewer.pal(nlevels(fac), "Set1")
  } else {
    colors <- c(
      "dodgerblue3",
      "firebrick3"
    )
  }
  
  # Generate PCA plot.
  xyplot(
    PC2 ~ PC1,
    groups = fac,
    data = as.data.frame(pca$x),
    pch = 16,
    cex = 1.5,
    aspect = "fill",
    col = colors,
    xlab = list(
      paste(
        "PC1 (",
        variance[1],
        "%)",
        sep = ""
      ),
      cex = 0.8
    ),
    ylab = list(
      paste(
        "PC2 (",
        variance[2],
        "%)",
        sep = ""
      ),
      cex = 0.8
    ),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...)
      ltext(
        x = x,
        y = y,
        labels = names,
        pos = 1,
        offset = 0.8,
        cex = 0.7
      )
    },
    main = draw.key(
      key = list(
        rect = list(col = colors),
        text = list(levels(fac)),
        rep = FALSE
      )
    )
  )
}


# ============================================================
# 3. IMPORT COUNT MATRICES AND SAMPLE METADATA
# ============================================================

# Gene-level count matrix.
countTableGenes <- read.csv(
  "gene_count_matrix.csv",
  row.names = 1
)

# Transcript-level count matrix.
countTableTranscripts <- read.csv(
  "transcript_count_matrix.csv",
  row.names = 1
)

# Sample metadata containing experimental group and batch.
colTable <- read.csv(
  "Design-Exp.csv",
  row.names = 1
)


# ============================================================
# 4. TRANSCRIPT-LEVEL ANALYSIS
# ============================================================

# Create DESeq2 object using transcript-level counts.
dds <- DESeqDataSetFromMatrix(
  countData = countTableTranscripts,
  colData = colTable,
  design = ~group
)

# Pre-filter low-count transcripts.
#
# A transcript must have at least 10 counts in at least
# three samples to be retained.
smallestGroupSize <- 3

keep <- rowSums(
  counts(dds) >= 10
) >= smallestGroupSize

dds <- dds[keep, ]

summary(dds)

# Run the DESeq2 workflow.
dds <- DESeq(dds)


# ------------------------------------------------------------
# Transcript-level dispersion QC
# ------------------------------------------------------------

plotDispEsts(dds)

svglite("tx.disp.svg")

print(
  plotDispEsts(dds)
)

dev.off()


# ------------------------------------------------------------
# Transcript-level transformations
# ------------------------------------------------------------

# rlog transformation for exploratory analysis.
rld <- rlog(
  dds,
  blind = FALSE
)

# Log2-transformed raw counts.
lgc.raw <- log2(
  counts(
    dds,
    normalized = FALSE
  ) + 1
)

# Log2-transformed normalised counts.
lgc.norm <- log2(
  counts(
    dds,
    normalized = TRUE
  ) + 1
)


# ------------------------------------------------------------
# Transcript-level PCA
# ------------------------------------------------------------

pca <- plotPCAWithSampleNames(
  assay(rld),
  targets = colTable,
  intgroup = "group"
)

pca

svglite("tx.pca.svg")

print(pca)

dev.off()


# ============================================================
# 5. GENE-LEVEL ANALYSIS
# ============================================================

# Create DESeq2 object using gene-level counts.
dds <- DESeqDataSetFromMatrix(
  countData = countTableGenes,
  colData = colTable,
  design = ~group
)

# Pre-filter low-count genes.
#
# A gene must have at least 10 counts in at least
# three samples to be retained.
smallestGroupSize <- 3

keep <- rowSums(
  counts(dds) >= 10
) >= smallestGroupSize

dds <- dds[keep, ]

summary(dds)

# Run the DESeq2 workflow.
dds <- DESeq(dds)


# ------------------------------------------------------------
# Gene-level dispersion QC
# ------------------------------------------------------------

par(mfrow = c(1, 1))

plotDispEsts(dds)

svglite("gene.disp.svg")

print(
  plotDispEsts(dds)
)

dev.off()


# ------------------------------------------------------------
# Gene-level transformations
# ------------------------------------------------------------

# rlog transformation for exploratory analysis.
rld <- rlog(
  dds,
  blind = FALSE
)

# Log2-transformed raw counts.
lgc.raw <- log2(
  counts(
    dds,
    normalized = FALSE
  ) + 1
)

# Log2-transformed normalised counts.
lgc.norm <- log2(
  counts(
    dds,
    normalized = TRUE
  ) + 1
)


# ------------------------------------------------------------
# Gene-level PCA
# ------------------------------------------------------------

pca <- plotPCAWithSampleNames(
  assay(rld),
  targets = colTable,
  intgroup = "group"
)

pca

svglite("gene.pca.svg")

print(pca)

dev.off()


# ============================================================
# 6. GENE-LEVEL TRANSFORMATION QC
# ============================================================

# Examine the relationship between mean expression and
# standard deviation after normalisation.

ntd <- normTransform(dds)

meanSdPlot(
  assay(ntd)
)

svglite("norm.sd_vs_mean.svg")

print(
  meanSdPlot(
    assay(ntd)
  )
)

dev.off()


# Examine the mean-SD relationship after rlog transformation.

meanSdPlot(
  assay(rld)
)

svglite("rlog.sd_vs_mean.svg")

print(
  meanSdPlot(
    assay(rld)
  )
)

dev.off()


# ============================================================
# 7. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================

# Compare group B against reference group A.
#
# Standard null hypothesis: log2 fold-change = 0.

lfc0.B_vs_A <- results(
  dds,
  contrast = c(
    "group",
    "B",
    "A"
  )
)

# Compare group C against reference group A.

lfc0.C_vs_A <- results(
  dds,
  contrast = c(
    "group",
    "C",
    "A"
  )
)

# Summarise differential expression results.

summary(
  lfc0.B_vs_A,
  alpha = 0.05
)

summary(
  lfc0.C_vs_A,
  alpha = 0.05
)

# Inspect the first rows of each result set.

head(lfc0.B_vs_A)

head(lfc0.C_vs_A)


# ------------------------------------------------------------
# Differential expression with LFC threshold
# ------------------------------------------------------------

# Repeat the analysis using a log2 fold-change threshold
# of 1.

lfc1.B_vs_A <- results(
  dds,
  contrast = c(
    "group",
    "B",
    "A"
  ),
  lfcThreshold = 1
)

lfc1.C_vs_A <- results(
  dds,
  contrast = c(
    "group",
    "C",
    "A"
  ),
  lfcThreshold = 1
)

summary(
  lfc1.B_vs_A,
  alpha = 0.05
)

summary(
  lfc1.C_vs_A,
  alpha = 0.05
)


# ============================================================
# 8. MA PLOTS
# ============================================================

# MA plots for the standard null hypothesis (LFC = 0).

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc0.B_vs_A
)

abline(
  h = 0,
  col = "black"
)

DESeq2::plotMA(
  lfc0.C_vs_A
)

abline(
  h = 0,
  col = "black"
)


svglite("ma.lfc0.svg")

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc0.B_vs_A
)

abline(
  h = 0,
  col = "black"
)

DESeq2::plotMA(
  lfc0.C_vs_A
)

abline(
  h = 0,
  col = "black"
)

dev.off()


# ------------------------------------------------------------
# MA plots using LFC threshold of 1
# ------------------------------------------------------------

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc1.B_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)

DESeq2::plotMA(
  lfc1.C_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)


svglite("ma.lfc1.svg")

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc1.B_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)

DESeq2::plotMA(
  lfc1.C_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)

dev.off()


# ============================================================
# 9. LOG2 FOLD-CHANGE SHRINKAGE
# ============================================================

# Estimate shrunken log2 fold-changes using the adaptive
# shrinkage method implemented by the ashr package.

# Standard null hypothesis: LFC = 0.

lfc0.B_vs_A <- lfcShrink(
  dds,
  contrast = c(
    "group",
    "B",
    "A"
  ),
  type = "ashr"
)

lfc0.C_vs_A <- lfcShrink(
  dds,
  contrast = c(
    "group",
    "C",
    "A"
  ),
  type = "ashr"
)


# LFC threshold of 1.

lfc1.B_vs_A <- lfcShrink(
  dds,
  contrast = c(
    "group",
    "B",
    "A"
  ),
  type = "ashr",
  lfcThreshold = 1
)

lfc1.C_vs_A <- lfcShrink(
  dds,
  contrast = c(
    "group",
    "C",
    "A"
  ),
  type = "ashr",
  lfcThreshold = 1
)


# ============================================================
# 10. IDENTIFY SIGNIFICANT DIFFERENTIALLY EXPRESSED GENES
# ============================================================

# B versus A
#
# Sort genes according to shrunken log2 fold-change and
# retain genes with adjusted p-value < 0.05.

B_vs_A.sorted <- lfc1.B_vs_A[
  order(
    lfc1.B_vs_A$log2FoldChange,
    decreasing = TRUE
  ),
]

B_vs_A.sig <- subset(
  B_vs_A.sorted,
  B_vs_A.sorted$padj < 0.05
)

head(B_vs_A.sig)

write.csv(
  B_vs_A.sig,
  file = "Significant.Genes.BvsA.csv",
  quote = FALSE
)


# C versus A

C_vs_A.sorted <- lfc1.C_vs_A[
  order(
    lfc1.C_vs_A$log2FoldChange,
    decreasing = TRUE
  ),
]

C_vs_A.sig <- subset(
  C_vs_A.sorted,
  C_vs_A.sorted$padj < 0.05
)

head(C_vs_A.sig)

write.csv(
  C_vs_A.sig,
  file = "Significant.Genes.CvsA.csv",
  quote = FALSE
)


# ============================================================
# 11. MA PLOTS OF SHRUNKEN RESULTS
# ============================================================

# MA plots for shrunken LFC estimates using the standard
# null hypothesis (LFC = 0).

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc0.B_vs_A
)

abline(
  h = 0,
  col = "black"
)

DESeq2::plotMA(
  lfc0.C_vs_A
)

abline(
  h = 0,
  col = "black"
)


svglite("ma.lfc0.s.svg")

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc0.B_vs_A
)

abline(
  h = 0,
  col = "black"
)

DESeq2::plotMA(
  lfc0.C_vs_A
)

abline(
  h = 0,
  col = "black"
)

dev.off()


# ------------------------------------------------------------
# MA plots of shrunken results using LFC threshold of 1
# ------------------------------------------------------------

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc1.B_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)

DESeq2::plotMA(
  lfc1.C_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)


svglite("ma.lfc1.s.svg")

par(mfrow = c(1, 2))

DESeq2::plotMA(
  lfc1.B_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)

DESeq2::plotMA(
  lfc1.C_vs_A
)

abline(
  h = c(1, -1),
  col = "black"
)

dev.off()


# ============================================================
# 12. EXPLORATORY BATCH-EFFECT CORRECTION
# ============================================================

# Extract the rlog-transformed gene expression matrix.

rld.m <- assay(rld)

# Remove the known batch effect from the rlog-transformed
# expression matrix using limma.
#
# This corrected matrix is used for exploratory visualisation
# and PCA rather than as input to the DESeq2 model.

b.corrected <- removeBatchEffect(
  rld.m,
  batch = colTable$batch
)


# ------------------------------------------------------------
# PCA after batch correction
# ------------------------------------------------------------

pca <- plotPCAWithSampleNames(
  b.corrected,
  targets = colTable,
  intgroup = "group"
)

pca


svglite("gene.pca.b.svg")

print(pca)

dev.off()


# ------------------------------------------------------------
# Export batch-corrected expression matrix
# ------------------------------------------------------------

write.csv(
  b.corrected,
  file = "BatchCorrected.Rlog.csv"
)


# ============================================================
# END OF ANALYSIS
# ============================================================
