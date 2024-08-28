#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:21:49 2024

@author: ferenc.kagan
"""

import pandas as pd

# Read the ccsmeth call_freqb file
ccsmeth_file = '/Users/ferenc.kagan/Documents/Projects/CpG_island/read_methylation/ssmed_5mC-ccsmeth.freq.count.all.bed'
output_file = '/Users/ferenc.kagan/Documents/Projects/CpG_island/read_methylation/methylKit_formatted.txt'

# Read the input file into a pandas DataFrame
ccsmeth_df = pd.read_csv(ccsmeth_file, sep='\t', header=None, 
                         names=['chr', 'start', 'end', 'dot', 'score', 'strand', 'start2', 'end2', 'color', 'coverage', 'freq'])

# Create the new DataFrame with the required format
methylation_df = pd.DataFrame({
    'chrBase': ccsmeth_df['chr'] + '.' + ccsmeth_df['start'].astype(str),
    'chr': ccsmeth_df['chr'],
    'base': ccsmeth_df['start'],
    'strand': ccsmeth_df['strand'].apply(lambda x: 'F' if x == '+' else 'R'),
    'coverage': ccsmeth_df['coverage'],
    'freqC': ccsmeth_df['freq'],
    'freqT': 100 - ccsmeth_df['freq']
})

# Save the output to a file
methylation_df.to_csv(output_file, sep='\t', index=False)

print(f"Transformed data saved to {output_file}")

