# Code Inventory

## Configuration

- `scripts/00_config.sh`: central project paths, tool paths, sample sizes, and output directories.

## LDSC preprocessing and genetic correlation

- `scripts/01_prepare_project_and_ldsc_inputs.sh`: checks raw files and converts ECAC beta/SE files into LDSC-compatible `SNP/A1/A2/N/Z/P` files.
- `scripts/02_run_ldsc_qc_and_rg.sh`: runs `munge_sumstats.py`, LDSC h² QC, and 3 exposure × 2 outcome genetic correlations.
- `scripts/07_make_ldsc_tables_from_logs.py`: parses LDSC h² and rg logs into TSV summary tables.

## conjFDR / pleioFDR input preparation

- `scripts/03_harmonize_conjfdr_inputs.sh`: creates HapMap3-filtered inputs, performs pairwise allele harmonization, flips outcome Z-scores when necessary, excludes MHC, and reports lambda GC sanity checks.
- `scripts/04_prepare_pleiofdr_mat_inputs.sh`: converts harmonized noMHC TSV files into pleioFDR `.mat` trait files using `python_convert` and `9545380.ref`.
- `scripts/05_prepare_decomposition_common_snps.sh`: post-processes standardized p-value files to create matched common-SNP inputs for 6-pair decomposition. Requires standardized p-value outputs generated upstream.

## Decomposition table generation

- `scripts/08_collect_decomp_raw_loci.py`: collects six-pair conjFDR locus CSVs into normalized `table_s9_decomp_conjfdr_loci_raw.tsv`.
- `scripts/09_group_decomp_loci_table.py`: groups raw conjFDR loci into clustered recurrent regions and classifies them as pan-spectrum, ADENO/overall-aligned, adenomyosis-biased candidate, endometriosis-core candidate, or inconclusive.

## Colocalization

- `scripts/06_run_coloc_chr12_rs9668810.R`: runs exploratory p-value-based coloc for chr12 rs9668810-indexed all EC and EEC regional inputs.
