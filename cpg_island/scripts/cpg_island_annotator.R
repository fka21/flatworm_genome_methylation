# Load required libraries

library(tidyverse)
library(GenomicFeatures)
library(ChIPseeker)

# Setup path
setwd("~/Documents/Projects/CpG_island/")

# Directories
bed_dir <- "./cpg_bed"
gff3_dir <- "./gene_annotations"
plot_dir <- "./plots"

# Ensure plot directory exists
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Get list of bed files
bed_files <- list.files(bed_dir, pattern = "_cpg.bed$", full.names = TRUE)

# Loop over each bed file
for (bed_file in bed_files) {
  # Extract base name without suffix
  base_name <- sub("_cpg.bed$", "", basename(bed_file))
  
  # Construct associated gff3 file path
  gff3_file <- file.path(gff3_dir, paste0(base_name, ".gff3"))
  
  if (file.exists(gff3_file)) {
    # Read and filter the bed file
    island <- read_tsv(bed_file, col_names = FALSE) %>%
      dplyr::select(-X6) %>%
      filter(str_detect(X1, "chr"))
    
    # Create GRanges object from bed file
    island <- GRanges(seqnames = island$X1, 
                      ranges = IRanges(start = island$X2, end = island$X3),
                      OE_ratio = island$X4,
                      GC_perc = island$X5)
    
    # Create TxDb object from the gff3 file
    annot <- GenomicFeatures::makeTxDbFromGFF(gff3_file)
    
    # Perform peak annotation
    peakAnno <- annotatePeak(island, TxDb = annot,
                             tssRegion = c(-3000, 3000), 
                             addFlankGeneInfo = TRUE, 
                             genomicAnnotationPriority = c("Intergenic", "Promoter", "5UTR", "3UTR", "Intron", "Exon", "Downstream"))
    
    # Plot and save the annotation pie chart
    pie_plot <- plotAnnoPie(peakAnno)
    ggsave(filename = file.path(plot_dir, paste0(base_name, "_anno_pie.png")), plot = pie_plot, width = 6, height = 6)
    
    # Plot and save the coverage plot
    cov_plot <- covplot(island, weightCol = "OE_ratio", title = "Observed/expected ratio over the genome", ylab = "Observed/expected ratio")
    ggsave(filename = file.path(plot_dir, paste0(base_name, "_cov_plot.png")), plot = cov_plot, width = 6, height = 4)
    
    # Plot and save the peak profile plot
    p1 <- plotPeakProf2(island, upstream = 5000, downstream = 5000,
                        conf = 0.95, by = "gene", type = "body", nbin = 100,
                        TxDb = annot, weightCol = "OE_ratio", ignore_strand = TRUE)
    ggsave(filename = file.path(plot_dir, paste0(base_name, "_peak_prof.png")), plot = p1, width = 6, height = 4)
  } else {
    message("GFF3 file not found for: ", base_name)
  }
}
