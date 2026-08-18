# Reproducibility Notes

## Raw data not included

This package does not include raw GWAS summary statistics, large LD/reference files, generated `.mat` files, or large local pleioFDR/conjFDR output directories. Users must obtain the source data and reference resources from the original providers and comply with all applicable data-use terms and restrictions.

## Main raw files expected

The workflow expects the following local raw input filenames:

- `adenomyosis_EUR.txt`
- `endometriosis_EUR_wo_23andMe.txt`
- `endometriosis_wo.adenomyosis_EUR.txt`
- `ECAC_allEC.valid_2.tsv`
- `ECAC_EEC.valid_2.tsv`

These filenames reflect the local analysis environment and are not intended to represent canonical filenames used by the original data providers.

## Key analysis constants

- All-histology endometrial cancer total sample size used for LDSC input preparation and colocalization: `N = 121885`
- Endometrioid endometrial cancer total sample size used for LDSC input preparation and colocalization: `N = 54884`
- Adenomyosis median eligible variant-level sample size used for P-value-based colocalization: `N = 33460.65`
- Adenomyosis case fraction used for colocalization: `8753 / (8753 + 415718) = 0.02062`
- All-histology endometrial cancer case fraction used for colocalization: `12906 / 121885 = 0.10589`
- Endometrioid endometrial cancer case fraction used for colocalization: `8758 / 54884 = 0.15957`
- Six primary cross-trait genetic-correlation tests were evaluated using Bonferroni correction with `alpha = 0.05 / 6`
- Benjamini–Hochberg FDR correction was additionally applied across the same six primary genetic-correlation tests as supportive information
- Matched conjunctional FDR SNP universe: `1,057,700` HapMap3 SNPs
- MHC exclusion window for conjunctional FDR processing: `chr6:26–34 Mb`
- MHC exclusion was not applied to the LDSC analyses
- pleioFDR/conjFDR random-pruning setting: `randprune_n = 20`
- Recurrent lead-SNP records were grouped when the gap from the current region end was `≤500 kb`
- chr12 colocalization region: `rs9668810 ±500 kb`
- Colocalization priors: `p1 = 1e-4`, `p2 = 1e-4`, `p12 = 1e-5`
- Colocalization input preparation is expected upstream of `scripts/06_run_coloc_chr12_rs9668810.R`

## Public-package reproducibility changes

- Removed raw data extraction material, internal conversation-derived text, troubleshooting logs, and other non-public analysis material from the GitHub-facing package.
- Generalized project, software, pleioFDR, MATLAB, and Python paths through `scripts/00_config.sh` and environment-variable overrides.
- Removed user-specific absolute paths from the public scripts.
- Fixed the lambda GC `awk` block in `scripts/03_harmonize_conjfdr_inputs.sh` so empty files do not trigger `next` in an `END` block.
- Made `scripts/04_prepare_pleiofdr_mat_inputs.sh` more robust to different `python_convert` layouts by auto-detecting `sumstats.py` and separating conversion-environment configuration from the script itself.
- Added explicit required-file checks to `scripts/05_prepare_decomposition_common_snps.sh`.
- Made `scripts/06_run_coloc_chr12_rs9668810.R` path-configurable through `$COLOC`, `$PROJECT`, or command-line arguments.
- Extended `scripts/07_make_ldsc_tables_from_logs.py` to reproduce Benjamini–Hochberg FDR correction and approximate between-phenotype genetic-correlation difference tests in addition to LDSC QC and genetic-correlation summary tables.
- Added `scripts/08_collect_decomp_raw_loci.py` to create `table_s9_decomp_conjfdr_loci_raw.tsv` before recurrent-region grouping.
- Retained `scripts/09_group_decomp_loci_table.py` for recurrent-region grouping of pair-specific conjFDR lead-SNP records.
- Added reproducible GTEx V10 scripts for rs9668810 single-tissue eQTL querying and targeted uterus eQTL auditing.

## Interpretive caution

The scripts reproduce or document statistical analyses and should not be interpreted as establishing causality.

- Variation in genetic-correlation point estimates across endometriosis-spectrum phenotypes does not establish statistically distinct or phenotype-exclusive genetic overlap.
- Approximate between-phenotype genetic-correlation difference tests do not model covariance arising from the shared endometrial cancer GWAS or non-independent endometriosis-spectrum estimates.
- Non-detection in a phenotype pair indicates lack of detection under the current GWAS power, phenotype definition, matched SNP universe, and conjunctional FDR threshold; it is not evidence of biological absence.
- Recurrent-locus labels describe observed detection patterns and should not be interpreted as proving phenotype-specific biological effects.
- chr12 colocalization is exploratory because adenomyosis effect estimates and standard errors were unavailable and P-value-based approximate Bayes factors were used.
- Strong colocalization posterior support indicates compatibility with a shared regional signal but does not establish a causal variant or causal gene.
- rs9668810 should be described as an index SNP or regional signal rather than as a causal variant.
- GTEx V10 eQTL associations are exploratory functional annotations and do not establish a causal gene, causal tissue, or biological mechanism.
- Lack of a significant bulk-uterus eQTL association should not be interpreted as evidence that the locus has no regulatory role in disease-relevant tissues or cell types.