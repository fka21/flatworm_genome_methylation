#!/usr/bin/env Rscript

library(GenomicRanges)
library(tidyverse)
library(regioneR)
library(RColorBrewer)
library(ChIPseeker)
library(methylKit)  

# Set working directory
setwd("~/Documents/Projects/CpG_island/")


# Load gene annotations
annot <- GenomicFeatures::makeTxDbFromGFF("./gene_annotations/schMedS3_h1.gff3")

##### Methylation Analysis with methylKit #####
# Read methylation data
myobjDB <- methRead("read_methylation/methylKit_formatted.txt",
                     sample.id = c("ccs"),
                     assembly = "ssmed", mincov = 1)

# Get methylation and coverage statistics
getMethylationStats(myobjDB, plot = TRUE, both.strands = FALSE)
getCoverageStats(myobjDB, plot = TRUE, both.strands = FALSE)

# Read and preprocess methylation data
methyl <- read_tsv("read_methylation/ssmed_5mC-ccsmeth.freq.count.all.bed", col_names = FALSE) %>%
  filter(X10 > 10) %>%
  filter(str_detect(X1, "chr"))

# Create GRanges object for methylation data
methyl_gr <- GRanges(seqnames = methyl$X1,
                     ranges = IRanges(start = methyl$X2, end = methyl$X3),
                     coverage = methyl$X10,
                     methyl_perc = methyl$X11)

# Subset methylation data that overlaps with CpG islands
methyl_gr_subset <- subsetByOverlaps(methyl_gr, island)

# Annotate peaks in methylation data
peakAnno <- annotatePeak(methyl_gr_subset, TxDb = annot,
                         tssRegion = c(-3000, 3000), 
                         addFlankGeneInfo = TRUE, 
                         genomicAnnotationPriority = c("Intergenic", "Downstream", "Promoter", "5UTR", "3UTR", "Intron", "Exon"))

# Plot annotation pie chart
plotAnnoPie(peakAnno)

# Filter high-confidence methylation data
methyl_gr_subset_hconf <- methyl_gr %>%
  plyranges::filter(methyl_perc > 70 & coverage > 10)

# Annotate peaks in high-confidence methylation data
peakAnno <- annotatePeak(methyl_gr_subset_hconf, TxDb = annot,
                         tssRegion = c(-3000, 3000), 
                         addFlankGeneInfo = TRUE, 
                         genomicAnnotationPriority = c("Intergenic", "Downstream", "Promoter", "5UTR", "3UTR", "Intron", "Exon"))

# Plot annotation pie chart for high-confidence data
plotAnnoPie(peakAnno)

# Plot coverage plot for high-confidence methylation data
covplot(methyl_gr_subset_hconf, weightCol = "methyl_perc",
        title = "Coverage of High-confidence Methylation Data", 
        ylab = "Methylation Percentage")

