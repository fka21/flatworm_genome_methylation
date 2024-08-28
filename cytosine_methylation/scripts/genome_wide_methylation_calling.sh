#!/bin/bash

# Script to run CCSMETH commands for methylation calls

# Set variables for file paths and parameters
INPUT_BAM="ssmed_kinteics.bam"
MODEL_FILE="/projects/ferenc.kagan/Tools/ccsmeth/models/model_ccsmeth_5mCpG_call_mods_attbigru2s_b21.v3.ckpt"
OUTPUT_DIR="./ssmed_5mC-ccsmeth"
REFERENCE_GENOME="/projects/ferenc.kagan/Methylation/00_raw_data/genomes/schMedS3_h1.fa"
HIFI_READS="ssmed_5mC-ccsmeth.bam"
MODBAM_OUTPUT="ssmed_5mC-ccsmeth.modbam.bam"
FREQB_OUTPUT="ssmed_5mC-ccsmeth.freq.count.all.bed"
THREADS=60
THREADS_CALL=20

# 1. Call Modifications
ccsmeth call_mods \
    --input "$INPUT_BAM" \
    --model_file "$MODEL_FILE" \
    --output "$OUTPUT_DIR" \
    --threads "$THREADS" \
    --threads_call "$THREADS_CALL" \
    --mode denovo

# 2. Align HiFi Reads
echo "Aligning HiFi reads to reference genome..."
ccsmeth align_hifi \
    --hifireads "$HIFI_READS" \
    --ref "$REFERENCE_GENOME" \
    --output "$MODBAM_OUTPUT" \
    --threads "$THREADS"

# 3. Call Frequency and Counts
echo "Calling frequency and counts of methylation modifications..."
ccsmeth call_freqb \
    --input_bam "$MODBAM_OUTPUT" \
    --ref "$REFERENCE_GENOME" \
    --output "$FREQB_OUTPUT" \
    --threads "$THREADS" \
    --sort \
    --bed \
    --call_mode count
