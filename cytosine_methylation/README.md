# CpG methylation identification pipeline

Genome-wide methylation was identified by using [ccsmeth](https://github.com/PengNi/ccsmeth) pipeline. The different steps used are found the `genome_wide_methylation_calling.sh` script.

The input for the analysis first needed to be converted to an appropriate format, this was achieved with `ccsmeth2methylkit.py`.

Analyses were done using the `CpG_Methylation_Analysis.R`.
