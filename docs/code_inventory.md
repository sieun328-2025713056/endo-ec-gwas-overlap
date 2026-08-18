# Code Inventory

## Configuration

- `scripts/00_config.sh`: defines central project paths, tool paths, sample sizes, and output directories using user-configurable environment variables.

## LDSC preprocessing and genetic correlation

- `scripts/01_prepare_project_and_ldsc_inputs.sh`: checks raw input files and converts ECAC beta/SE files into LDSC-compatible `SNP/A1/A2/N/Z/P` files.
- `scripts/02_run_ldsc_qc_and_rg.sh`: runs `munge_sumstats.py`, LDSC SNP-heritability QC, and the six primary genetic-correlation analyses across three endometriosis-spectrum phenotypes and two endometrial cancer outcomes.
- `scripts/07_make_ldsc_tables_from_logs.py`: parses LDSC h² and genetic-correlation logs, generates QC and rg summary tables, applies Benjamini–Hochberg FDR correction across the six primary cross-trait tests, and performs approximate between-phenotype rg difference tests within each endometrial cancer outcome. The difference tests use `sqrt(SE1^2 + SE2^2)` and do not model covariance arising from the shared cancer GWAS or non-independent endometriosis-spectrum estimates.

## conjFDR / pleioFDR input preparation

- `scripts/03_harmonize_conjfdr_inputs.sh`: creates HapMap3-filtered inputs, performs pairwise allele harmonization, flips endometrial cancer Z-scores when necessary to align effect alleles, excludes the MHC region for conjunctional FDR processing, and reports lambda GC sanity checks.
- `scripts/04_prepare_pleiofdr_mat_inputs.sh`: verifies or optionally converts standardized `.std.pval.tsv` files into pleioFDR-compatible `.mat` trait files using `python_convert` and `9545380.ref`.
- `scripts/05_prepare_decomposition_common_snps.sh`: post-processes standardized p-value files to create matched common-SNP inputs for the six phenotype-pair conjunctional FDR analyses. Requires standardized p-value outputs generated upstream.

## Decomposition table generation

- `scripts/08_collect_decomp_raw_loci.py`: collects six-pair conjFDR locus CSVs into normalized `table_s9_decomp_conjfdr_loci_raw.tsv`.
- `scripts/09_group_decomp_loci_table.py`: groups pair-specific conjFDR lead-SNP records into recurrent genomic regions according to chromosome and physical proximity and summarizes their detection patterns across phenotype-pair comparisons.

## Colocalization

- `scripts/06_run_coloc_chr12_rs9668810.R`: runs exploratory P-value-based colocalization for the chr12 rs9668810-indexed region using all-histology and endometrioid endometrial cancer regional inputs.
