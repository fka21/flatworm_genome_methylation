#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 09:54:25 2024

@author: Ferenc Kagan
Description: This script identifies CpG islands in multiFASTA sequences using a sliding window approach.
             It calculates CpG island significance using various methods and saves results in BED format.
"""

import argparse
import multiprocessing
import logging
import math
import os

import numpy as np
from Bio import SeqIO
from scipy.stats import binomtest

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_cpg_probability(sequence, window_size, step_size, significance_method):
    """
    Calculate CpG island probability for a given sequence using a sliding window approach.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        window_size (int): The size of the sliding window.
        step_size (int): The step size for the sliding window.
        significance_method (str): Method to calculate significance ('oe_ratio', 'z_score', 'p_value').

    Returns:
        list: A list of tuples containing start, end, score, GC content, and observed/expected ratio.
    """
    results = []
    n = len(sequence)

    for i in range(0, n - window_size + 1, step_size):
        window = sequence[i:i + window_size]
        c_count = window.count('C')
        g_count = window.count('G')
        cg_count = window.count('CG')
        
        observed = cg_count
        expected = (c_count * g_count) / window_size if window_size != 0 else 0
        gc_content = (c_count + g_count) / window_size * 100 if window_size != 0 else 0
        oe_ratio = observed / expected if expected > 0 else 0
        
        if gc_content >= 50 and oe_ratio > 0.6:  # Apply CpG island criteria
            if significance_method == "oe_ratio":
                score = oe_ratio
            elif significance_method == "z_score":
                std_dev = math.sqrt((c_count * g_count * (window_size - c_count) * (window_size - g_count)) / (window_size**2 * (window_size - 1))) if window_size > 1 else 0
                score = (observed - expected) / std_dev if std_dev > 0 else 0
            elif significance_method == "p_value":
                score = calculate_p_value(observed, expected, window_size)
            else:
                raise ValueError("Unknown significance method: " + significance_method)
            
            if score >= 1:  # Adjust this threshold if needed
                results.append((i, i + window_size, score, gc_content, oe_ratio))
    
    return results

def calculate_p_value(observed, expected, window_size):
    """
    Calculate p-value for the observed number of CpG dinucleotides.

    Parameters:
        observed (int): Number of observed CpG dinucleotides.
        expected (float): Expected number of CpG dinucleotides.
        window_size (int): Size of the window.

    Returns:
        float: -log10 of the p-value.
    """
    if expected == 0:
        return 1.0
    prob_cg = expected / window_size if window_size != 0 else 0
    result = binomtest(observed, n=window_size, p=prob_cg, alternative='greater')
    return -math.log10(result.pvalue)

def process_sequence(sequence, window_size, step_size, significance_method):
    """
    Process a single sequence to identify CpG islands and recalculate results.

    Parameters:
        sequence (str): The DNA sequence to analyze.
        window_size (int): The size of the sliding window.
        step_size (int): The step size for the sliding window.
        significance_method (str): Method to calculate significance.

    Returns:
        list: Recalculated results after merging overlapping regions.
    """
    initial_results = calculate_cpg_probability(sequence, window_size, step_size, significance_method)
    return merge_and_recalculate(sequence, initial_results, significance_method)

def merge_and_recalculate(sequence, results, significance_method):
    """
    Merge overlapping CpG island results and recalculate significance.

    Parameters:
        sequence (str): The DNA sequence.
        results (list): Initial list of CpG island results.
        significance_method (str): Method to calculate significance.

    Returns:
        list: Merged and recalculated results.
    """
    if not results:
        return []

    merged_results = []
    current_start, current_end, _, current_gc_content, current_oe_ratio = results[0]
    current_observed = (current_end - current_start) * current_oe_ratio
    current_expected = (current_end - current_start) * (current_gc_content / 100)

    for start, end, score, gc_content, oe_ratio in results[1:]:
        if start <= current_end:
            current_end = end
            current_observed += (end - start) * oe_ratio
            current_expected += (end - start) * (gc_content / 100)
        else:
            merged_results.append((current_start, current_end, current_observed, current_expected))
            current_start, current_end, current_gc_content, current_oe_ratio = start, end, gc_content, oe_ratio
            current_observed = (current_end - current_start) * current_oe_ratio
            current_expected = (current_end - current_start) * (current_gc_content / 100)

    merged_results.append((current_start, current_end, current_observed, current_expected))

    reevaluated_results = []
    for start, end, observed, expected in merged_results:
        window_size = end - start
        gc_content = (sequence[start:end].count('C') + sequence[start:end].count('G')) / window_size * 100 if window_size != 0 else 0
        oe_ratio = observed / expected if expected > 0 else 0
        
        if significance_method == "oe_ratio":
            score = oe_ratio
        elif significance_method == "z_score":
            std_dev = math.sqrt((observed * (window_size - observed)) / (window_size * (window_size - 1))) if window_size > 1 else 0
            score = (observed - expected) / std_dev if std_dev > 0 else 0
        elif significance_method == "p_value":
            score = calculate_p_value(observed, expected, window_size)
        
        if score >= 1 and gc_content >= 50 and oe_ratio > 0.6:
            reevaluated_results.append((start, end, score, gc_content, oe_ratio))
    
    return reevaluated_results

def worker(task):
    """
    Worker function for multiprocessing to process a sequence.

    Parameters:
        task (tuple): Contains sequence ID, sequence, window size, step size, and significance method.

    Returns:
        tuple: Sequence ID and results.
    """
    seq_id, sequence, window_size, step_size, significance_method = task
    results = process_sequence(sequence, window_size, step_size, significance_method)
    return seq_id, results

def main(args):
    """
    Main function to handle argument parsing, sequence processing, and result saving.

    Parameters:
        args (argparse.Namespace): Command-line arguments.
    """
    window_size = args.window_size
    step_size = args.step_size
    significance_method = args.significance_method
    num_processes = args.num_processes

    tasks = []

    try:
        for record in SeqIO.parse(args.input, "fasta"):
            tasks.append((record.id, str(record.seq), window_size, step_size, significance_method))

        logging.info(f"Parsed {len(tasks)} sequences from the input file.")

        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.map(worker, tasks)

        with open(args.output, "w") as bed_file:
            for seq_id, result in results:
                for start, end, score, gc_content, oe_ratio in result:
                    bed_file.write(f"{seq_id}\t{start}\t{end}\t{score:.4f}\t{gc_content:.2f}\t{oe_ratio:.2f}\n")
        logging.info("Results successfully written to the output file.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify CpG islands in multiFASTA sequences using a sliding window approach")
    parser.add_argument("-i", "--input", required=True, help="Input multiFASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument("-w", "--window_size", type=int, default=200, help="Window size for sliding window")
    parser.add_argument("-s", "--step_size", type=int, default=1, help="Step size for sliding window")
    parser.add_argument("-m", "--significance_method", choices=["oe_ratio", "z_score", "p_value"], default="oe_ratio", help="Method to calculate significance")
    parser.add_argument("-p", "--num_processes", type=int, default=1, help="Number of processes to use")

    args = parser.parse_args()
    main(args)
