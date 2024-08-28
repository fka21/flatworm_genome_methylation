# CpG island identification pipeline

Using a custom python script `cpg_island_finder.py` I identify CpG islands in several genomes. Some of these genomes have published data on CpG islands and CpG methylation and were used as a standard during testing the `cpg_island_finder.py` script.

A wrapper script is also present `iterative_cpg_finder_and_plottersh.sh` which iterates through genome files and calls `cpg_island_finder.py` and some other tools for plotting.

Output is analyzed using `CpG_Island_Enrichment_Analysis.R`.

Some of the plots can be found in the `output/`.
