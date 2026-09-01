# Bulk RNA-seq Analysis of *Mus musculus*

An end-to-end bulk RNA-seq analysis workflow developed to process and prepare *Mus musculus* chromosome 2 RNA-seq data for differential expression analysis.

The project covers read preprocessing, reference alignment, transcript assembly and quantification, followed by exploratory quality control and differential expression analysis in R using DESeq2. The analysis also investigates the effect of experimental batch on gene expression measurements.

The workflow was originally developed and executed as part of MSc Bioinformatics coursework at the University of Glasgow and has been adapted here as a portfolio project.

## Workflow
1. Raw RNA-seq reads

2. Adapter trimming (Scythe)

3. Quality/length filtering (Sickle)

4. Read alignment (HISAT2)

5. SAM → BAM conversion and sorting (SAMtools)

6. Transcript assembly & quantification (StringTie)

7. Gene/transcript count matrices

8. Exploratory analysis & QC (R)

9. Differential expression analysis (DESeq2): 
   - PCA 
   - Dispersion analysis
   - MA plots
   - Mean-SD analysis
   - Batch-effect investigation (limma)
