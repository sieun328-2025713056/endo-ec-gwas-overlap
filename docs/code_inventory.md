
# Code Inventory

## Overview

This document describes the public scripts and their roles in the endometriosis-spectrum and endometrial cancer GWAS workflow.

The repository contains original workflow scripts, historical exploratory analyses, and newly added final-analysis reference and verification scripts.

Scripts do not necessarily reproduce the complete study independently. Several require source GWAS data, external reference resources, or previously generated intermediate files that are not distributed in this repository.

## Configuration

- `scripts/00_config.sh`: defines configurable project directories, external tool locations, reference paths, sample sizes, and output directories.

## LDSC preprocessing and genetic correlations

- `scripts/01_prepare_project_and_ldsc_inputs.sh`: checks source GWAS inputs and prepares endometrial cancer summary statistics for LDSC.

- `scripts/02_run_ldsc_qc_and_rg.sh`: performs LDSC summary-statistics munging, SNP-heritability QC, and six bivariate genetic-correlation analyses across three endometriosis-spectrum phenotypes and two endometrial cancer outcomes.

- `scripts/07_make_ldsc_tables_from_logs.py`: extracts LDSC QC and genetic-correlation results from existing logs. It also calculates historical Benjamini-Hochberg corrections and approximate genetic-correlation difference tests.

The six bivariate genetic correlations provide descriptive results. They are not the six final primary hypothesis tests.

The historical approximate difference calculations ignore sampling covariance and are not the final inferential basis.

## Final covariance-aware genetic-correlation inference

- `scripts/12_rgmodel_final_tests.R`: reads saved GenomicSEM `rgmodel()` results for allEC and EEC and calculates covariance-aware omnibus Wald tests and six pairwise contrasts.

The script uses genetic correlations from `R` and their sampling covariance from `V_R`.

Two primary omnibus tests are evaluated using a Bonferroni threshold of P < 0.025.

Six post-hoc pairwise contrasts are corrected jointly using the Holm method.

Required intermediate files:

- `results/rgmodel_allEC.rds`
- `results/rgmodel_EEC.rds`

These intermediate files are not distributed in this repository.

## LDSC block-sensitivity analysis

- `scripts/13_rgmodel_sensitivity.R`: reads existing four-trait LDSC RDS files for the 100-block and 400-block settings, runs GenomicSEM `rgmodel()`, and calculates omnibus Wald statistics.

The modified public script does not save or overwrite GenomicSEM RDS outputs.

It requires previously generated LDSC RDS files and the appropriate R environment.

The public copy passed syntax parsing but has not been independently rerun end-to-end following its repository adjustments.

## ENDO-only MAF sensitivity verification

- `scripts/15_maf05_rgmodel_verification.R`: independently verifies downstream GenomicSEM and Wald calculations using saved LDSC RDS files from the ENDO-only MAF sensitivity analysis.

The upstream ENDO munging used `--maf-min 0.05`.

The verification script:

- Reads saved four-trait LDSC objects.
- Calls GenomicSEM `rgmodel()`.
- Computes two covariance-aware omnibus Wald tests.
- Computes six pairwise contrasts with joint Holm correction.
- Compares selected results against previously reported values.
- Does not rerun GWAS munging or LDSC.
- Does not write or overwrite analysis results.

Execution requires `RUN_MAF05_VERIFICATION=1`.

Optional settings include `GSEM_LIB` and `MAF_LDSC_DIR`.

The downstream verification was executed successfully and matched the reported results within an absolute tolerance of 1e-5.

See `docs/maf05_rgmodel_validation.md`.

## conjFDR / pleioFDR input preparation

- `scripts/03_harmonize_conjfdr_inputs.sh`: prepares HapMap3-based inputs, performs allele harmonization, aligns relevant effect directions, and excludes the MHC region for conjunctional FDR processing.

- `scripts/04_prepare_pleiofdr_mat_inputs.sh`: verifies or optionally generates pleioFDR-compatible MATLAB input files from standardized GWAS inputs.

- `scripts/05_prepare_decomposition_common_snps.sh`: verifies matched common-SNP resources and input files used for the six phenotype-pair conjunctional FDR analyses.

The final matched universe contained 1,057,700 HapMap3 SNPs.

The chr6:26-34 Mb MHC exclusion applies to conjFDR processing, not LDSC.

The public workflow does not independently regenerate every upstream pleioFDR decomposition resource.

## conjFDR locus collection and grouping

- `scripts/08_collect_decomp_raw_loci.py`: collects pair-specific conjFDR locus records into a standardized table.

- `scripts/09_group_decomp_loci_table.py`: groups lead-SNP records by chromosome and genomic proximity and summarizes detection patterns across phenotype-pair comparisons.

The final regional classification comprised:

- 92 lead-SNP records
- 45 grouped genomic regions
- 22 regions detected in more than one phenotype-pair comparison

Multiple detections across dependent comparisons are descriptive, not independent replication.

## Exploratory chr12 colocalization

### Historical implementation

- `scripts/06_run_coloc_chr12_rs9668810.R`: historical exploratory colocalization implementation using P-value-based inputs.

This script is retained for provenance and should not be described as the final colocalization implementation.

### Updated final reference implementation

- `scripts/14_coloc_chr12_final_reference.R`: implements the updated chr12 colocalization input specification.

Final input structure:

Adenomyosis:

- P-values and minor allele frequencies
- Median eligible variant-level `AD_N = 33460.65`
- Assumed case fraction `s = 0.5`

Endometrial cancer:

- Aligned effect estimates
- Squared standard errors as variance estimates
- Outcome-specific total sample sizes and case fractions

Priors:

- `p1 = 1e-4`
- `p2 = 1e-4`
- `p12 = 1e-5`

The derivation of adenomyosis effective sample size and appropriateness of the `s = 0.5` assumption remain unresolved methodological limitations.

The script reads regional inputs from a configurable directory.

Optional SNP-level posterior export is disabled by default and requires both `EXPORT_COLOC_POSTERIOR=1` and an existing `COLOC_POSTERIOR_DIR`.

The script checks whether the destination file already exists before attempting to write.

The public copy passed syntax parsing, but an independent end-to-end rerun after its repository adjustments has not been established.

## GTEx V10 functional annotation

- `scripts/clean/09_query_gtex_v10_rs9668810.py`: queries GTEx V10 for single-tissue eQTL associations involving rs9668810 and targeted uterus eQTL information.

- `scripts/clean/10_make_gtex_uterus_final_audit.py`: summarizes the reproduced uterus-query outputs into an audit table.

These analyses are exploratory functional annotations rather than experimental validation.

## Validation status

The final-analysis reference scripts `12–15` passed R syntax parsing.

The MAF verification script `15` was executed successfully against saved LDSC inputs.

The modified public copies of scripts `13` and `14` have not been independently validated end-to-end following their repository adjustments.

The repository does not provide a fully verified one-command reproduction of the entire analytical workflow.

## Related documentation

- `docs/reproducibility_notes.md`: analysis specifications, computational dependencies, and reproducibility limitations.
- `docs/maf05_rgmodel_validation.md`: verified MAF-sensitivity calculations.
- `docs/final_output_tables.md`: output provenance and historical table names.
- `docs/interpretation_limits.md`: statistical interpretation and methodological limitations.
- `docs/data_sources.md`: input data provenance and restrictions.
