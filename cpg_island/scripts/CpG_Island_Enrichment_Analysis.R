# Load required libraries
library(GenomicRanges)
library(tidyverse)
library(rtracklayer)
library(ChIPseeker)
library(clusterProfiler)
library(GO.db)

# Set working directory
setwd("~/Documents/Projects/CpG_island/")

# Load CpG island data
island <- read_tsv("./cpg_bed/schMedS3_h1_cpg.bed", col_names = FALSE) %>%
  dplyr::select(-X6) %>%
  filter(str_detect(X1, "chr"))

# Create TxDb object from GFF file
annot <- GenomicFeatures::makeTxDbFromGFF("./gene_annotations/schMedS3_h1.gff3")

# Convert CpG island data to GRanges object
island <- GRanges(seqnames = island$X1, 
                  ranges = IRanges(start = island$X2, end = island$X3),
                  OE_ratio = island$X4,
                  GC_perc = island$X5)

# Peak annotation
peakAnno <- annotatePeak(island, TxDb = annot,
                         tssRegion = c(-3000, 3000), 
                         addFlankGeneInfo = TRUE, 
                         genomicAnnotationPriority = c("Intergenic", "Downstream", "Promoter", "5UTR", "3UTR", "Intron", "Exon"))

# Plot and save annotation pie chart
svg("smed_region_annotation.svg", height = 10, width = 10)
plotAnnoPie(peakAnno)
dev.off()

# Plot and save coverage plot
svg("smed_region_coverage.svg", height = 10, width = 10)
covplot(island, weightCol = "OE_ratio", title = "Observed/expected ratio over the genome", ylab = "Observed/expected ratio")
dev.off()

# Plot and save peak profile plot
svg("smed_region_gene-body.svg", height = 10, width = 10)
plotPeakProf2(island, upstream = 5000, downstream = 5000, conf = 0.95, by = "gene", type = "body", nbin = 100, TxDb = annot, weightCol = "OE_ratio", ignore_strand = TRUE)
dev.off()

# Load ATAC-seq data
atac <- read_tsv("./smed_specific_annotaiton/SMED_accessible_chromatin.narrowPeak", col_names = FALSE) %>%
  filter(str_detect(X4, "primary"))

atac_gr <- GRanges(seqnames = atac$X1, ranges = IRanges(start = atac$X2, end = atac$X3))
genome <- data.frame(island@seqnames@values, island@seqnames@lengths)

# Perform overlap permutation test
pt <- overlapPermTest(A = atac_gr, B = island, ntimes = 500, force.parallel = FALSE, genome = genome)

# Inspect overlaps
numOverlaps(island, atac_gr, count.once = TRUE)

# Subset islands by overlaps with ATAC-seq peaks
atac_island_gr <- subsetByOverlaps(island, atac_gr)

# Peak annotation for ATAC-seq subset
peakAnno <- annotatePeak(atac_island_gr, TxDb = annot,
                         tssRegion = c(-3000, 3000), 
                         addFlankGeneInfo = TRUE, 
                         genomicAnnotationPriority = c("Intergenic", "Downstream", "Promoter", "5UTR", "3UTR", "Intron", "Exon"))

# Plot and save annotation pie chart for ATAC-seq subset
svg("atac_region_annotation.svg", height = 10, width = 10)
plotAnnoPie(peakAnno)
dev.off()

# Plot and save coverage plot for ATAC-seq subset
svg("atac_region_coverage.svg", height = 10, width = 10)
covplot(atac_island_gr, weightCol = "OE_ratio", title = "Observed/expected ratio over the genome", ylab = "Observed/expected ratio")
dev.off()

# Plot and save peak profile plot for ATAC-seq subset
svg("atac_region_gene-body.svg", height = 10, width = 10)
plotPeakProf2(atac_island_gr, upstream = 5000, downstream = 5000, conf = 0.95, by = "gene", type = "body", nbin = 100, TxDb = annot, weightCol = "OE_ratio", ignore_strand = TRUE)
dev.off()

# Load repeat data
repeat1 <- read_tsv("./smed_specific_annotaiton/schMedS3_h1_rRNA.bed", col_names = FALSE) %>%
  filter(str_detect(X1, "chr"))
repeat2 <- read_tsv("./smed_specific_annotaiton/schMedS3_h1_satDNA.bed", col_names = FALSE) %>%
  filter(str_detect(X1, "chr"))
repeat3 <- read_tsv("./smed_specific_annotaiton/schMedS3_h1_barnap.bed", col_names = FALSE) %>%
  filter(str_detect(X1, "chr"))
repeat4 <- read_tsv("./read_methylation/schMedS3_h1.fa.mod.EDTA.TEanno.bed", col_names = FALSE) %>%
  filter(str_detect(X1, "chr"))

repeat_gr <- GRanges(seqnames = c(repeat1$X1, repeat2$X1, repeat3$X1, repeat4$X1),
                     ranges = IRanges(start = c(repeat1$X2, repeat2$X2, repeat3$X2, repeat4$X2),
                                      end = c(repeat1$X3, repeat2$X3, repeat3$X3, repeat4$X3)))

# Perform overlap permutation test with repeats
pt2 <- overlapPermTest(A = repeat_gr, B = island, ntimes = 500, force.parallel = FALSE, genome = genome)

# Functional enrichment analysis with emapper + interproscan functional annotation
interpro <- read_tsv("gene_annotations/schMedS3_h1.tsv", col_names = FALSE) %>%
  separate_rows(X14, sep = "\\|")
desc <- unique(AnnotationDbi::select(GO.db, keys = interpro$X14, columns = c("GOID", "TERM"), keytype = "GOID"))

enrich_bp <- enricher(peakAnno@anno$transcriptId, TERM2GENE = interpro[, c("X14", "X1")], TERM2NAME = desc[, c("GOID", "TERM")])
dotplot(enrich_bp, showCategory = 20) +
  theme(axis.text.y = element_text(size = 1))
ggsave("smed_cpg-island_enriched_terms_interpro.png", width = 10, height = 10, units = 'in', dpi = 320)

# Combine and annotate GO terms from FANTASIA functional annotation pipeline
bpo <- read.table("smed_specific_annotaiton/gopredsim_ssmed_prott5_1_bpo.txt", header = FALSE, sep = "\t") %>%
  dplyr::mutate(V4 = Term(V2)) 
mfo <- read.table("smed_specific_annotaiton/gopredsim_ssmed_prott5_1_mfo.txt", header = FALSE, sep = "\t") %>%
  dplyr::mutate(V4 = Term(V2)) 
bpo <- rbind(bpo, mfo)
desc <- AnnotationDbi::select(GO.db, keys = bpo$V2, columns = c("GOID", "TERM"), keytype = "GOID")

enrich_bp <- enricher(peakAnno@anno$transcriptId, TERM2GENE = bpo[, c("V2", "V1")], TERM2NAME = desc[, c("GOID", "TERM")])
dotplot(enrich_bp)

