#!/bin/bash

# Use the kent suit of tools and bedtools

# Directories
input_dir="input"
output_dir="output"

# Ensure output directory exists
mkdir -p $output_dir

# Iterate over all fasta files in the input directory
for fasta_file in $input_dir/*.fa; 
do
    # Get the base name of the file (without directory and extension)
    base_name=$(basename "$fasta_file" .fa)
    
    # Define output file names
    output_file_ucsc="$output_dir/${base_name}_cpg_ucsc.out"
    output_file_bed="$output_dir/${base_name}_cpg.bed"
    genome_size_file="$output_dir/${base_name}.sizes"
    windows_file="$output_dir/${base_name}_windows.bed"
    cpg_counts_file="$output_dir/${base_name}_cpg_counts.bed"
    cpg_counts_sorted_file="$output_dir/${base_name}_cpg_counts_sorted.bed"
    matrix_file_dat="$output_dir/${base_name}_matrix.tss.dat"
    matrix_file_txt="$output_dir/${base_name}_matrix.tss.txt"
    plot_file="$output_dir/${base_name}.pdf"
    
    # Execute the CpG finder
    python3 cpg_island_finder.py --input "$fasta_file" --output "$output_file_bed" -p 100

    # Extract TSS from the GFF3 file
    awk 'BEGIN{OFS="\t"} $3 == "gene" {print $1, $4-1, $5, $9, ".", $7}'  "gene_annotations/${base_name}.gff3" > "$output_dir/${base_name}_tss.bed"

    # Generate genome size file
    faSize -detailed "$fasta_file" > "$genome_size_file"

    # Create windows of 100 bp
    bedtools makewindows -g "$genome_size_file" -w 100 > "$windows_file"

    # Count CpG islands in each window
    bedtools coverage -a "$windows_file" -b "$output_file_bed" -counts > "$cpg_counts_file"

    # Sort the counts file
    bedSort "$cpg_counts_file" "$cpg_counts_sorted_file"
    sort -k1,1 -k2,2n "$cpg_counts_file" > "$cpg_counts_sorted_file"

    # Convert sorted BED to BigWig
    bedGraphToBigWig "$cpg_counts_sorted_file" "$genome_size_file" "$output_dir/${base_name}_cpg_counts.bw"

    # Compute matrix for TSS regions
    computeMatrix scale-regions -S "$output_dir/${base_name}_cpg_counts.bw" -R "$output_dir/${base_name}_tss.bed" --regionBodyLength 5000 -b 5000 -a 5000 --outFileName "$matrix_file_dat" --outFileNameMatrix "$matrix_file_txt" --numberOfProcessors=max

    # Plot heatmap
    plotHeatmap --matrixFile "$matrix_file_dat" --outFileName "$plot_file" --sortRegions descend --sortUsing mean --legendLocation 'none' --heatmapHeight 10 --heatmapWidth 10 --yAxisLabel 'CpG island coverage' --plotType 'se'

done