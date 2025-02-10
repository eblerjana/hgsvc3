# Pipelines used to generate PanGenie genotyping / phasing results

## Main pipelines:
* ``experiments/genotyping/``: full pipeline to produce genotypes for 1kGP cohort from MC graphs
* ``experiments/phasing/``: full pipeline to phase genotypes and produce consensus haplotypes for 1kGP cohort
* ``experiments/evaluate-consensus-haplotypes/``: full pipeline used to compute variant-based QV estimates for consensus haplotypes

## Preliminary experiments / internal analyses 
* ``experiments/other-analyses/``: scripts to produce supplementary figures comparing MC and PAV
* ``experiments/PAV-merged-set-analysis/``: pipeline used to evaluate PAV calls before and after SV merging
