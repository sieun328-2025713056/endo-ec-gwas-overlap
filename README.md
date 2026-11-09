
# Endometriosis-Spectrum and Endometrial Cancer GWAS

## Overview

This repository documents a GWAS summary-statistics workflow investigating shared genetic architecture between endometriosis-spectrum phenotypes and endometrial cancer.

The study combines genome-wide genetic correlation analyses, covariance-aware heterogeneity testing, conjunctional false discovery rate (conjFDR) analyses, and exploratory regional follow-up.

The repository supports analytical transparency and reproducibility. It does not contain all data or intermediate files required for independent end-to-end reproduction.

## Study design

Three endometriosis-spectrum phenotypes were examined:

1. Overall endometriosis, excluding 23andMe
2. Adenomyosis
3. Endometriosis without adenomyosis

Two endometrial cancer outcomes were considered:

1. All-histology endometrial cancer (allEC)
2. Endometrioid endometrial cancer (EEC)

The analysis had three levels:

- **Primary:** LDSC genetic correlations and GenomicSEM covariance-aware heterogeneity tests
- **Secondary:** Six matched conjunctional FDR analyses and recurrent-detection pattern classification
- **Exploratory:** Selected chr12 colocalization and GTEx V10 functional annotation

The exploratory chr12 region was selected following the secondary conjFDR analysis, rather than being an independent primary analysis.

## Data sources

### Endometriosis-spectrum GWAS

European-ancestry summary statistics were obtained from the endometriosis and adenomyosis GWAS resources associated with:

- https://doi.org/10.5281/zenodo.18983492

The overall endometriosis input excluded the 23andMe component.

### Endometrial cancer GWAS

The endometrial cancer GWAS summary statistics correspond to:

- All-histology EC: [GCST006464](https://www.ebi.ac.uk/gwas/studies/GCST006464)
- Endometrioid EC: [GCST006465](https://www.ebi.ac.uk/gwas/studies/GCST006465)

The sample sizes used for the EC analyses were:

| Outcome | Cases | Controls | Total N |
|---|---:|---:|---:|
| allEC | 12,906 | 108,979 | 121,885 |
| EEC | 8,758 | 46,126 | 54,884 |

### Data restrictions

Raw GWAS summary statistics and large reference files are not redistributed here.

Users must obtain source datasets from their original providers and follow applicable access, licensing, and data-use requirements.

Expected local raw filenames:

```text
adenomyosis_EUR.txt
endometriosis_EUR_wo_23andMe.txt
endometriosis_wo.adenomyosis_EUR.txt
ECAC_allEC.valid_2.tsv
ECAC_EEC.valid_2.tsv
```

These are local workflow filenames, not necessarily the providers' original filenames.

## Primary analysis: LDSC and GenomicSEM

Six bivariate LDSC genetic correlations were estimated to characterize genetic overlap.

The final primary inference evaluated whether genetic correlations differed across the three endometriosis-spectrum phenotypes within each EC outcome.

GenomicSEM `rgmodel()` provided:

- `R`: estimated genetic correlations
- `V_R`: sampling covariance matrix

Covariance-aware omnibus Wald tests were performed separately for allEC and EEC.

The primary multiple-testing threshold was:

- Two omnibus tests
- Bonferroni-corrected significance: P < 0.025

Six post-hoc pairwise contrasts were evaluated using Holm correction across all six comparisons.

### Primary omnibus results

| Outcome | Wald ChiSq(2) | P |
|---|---:|---:|
| allEC | 8.155546 | 0.016945159 |
| EEC | 9.632358 | 0.008097669 |

Both primary omnibus tests met the corrected threshold.

### Selected post-hoc comparison

For the overall endometriosis versus adenomyosis contrast:

| Outcome | Raw P | Holm P |
|---|---:|---:|
| allEC | 0.012022315 | 0.06011158 |
| EEC | 0.006044072 | 0.03626443 |

These findings do not establish that every phenotype pair has significantly different genetic correlations.

Earlier approximate tests that ignored covariance are retained for provenance but are not the final inferential basis.

## Sensitivity analyses

The LDSC sampling-covariance estimation was evaluated using different block settings.

| Analysis | allEC omnibus P | EEC omnibus P |
|---|---:|---:|
| 100 blocks | 0.01151326 | 0.007853592 |
| 200 blocks (primary) | 0.016945159 | 0.008097669 |
| 400 blocks | 0.02482712 | 0.01197984 |
| ENDO-only MAF filter | 0.02458479 | 0.01214624 |

For the ENDO-only MAF sensitivity analysis:

- The overall-endometriosis LDSC input was reprocessed with `--maf-min 0.05`.
- The munged ENDO input retained 1,064,328 SNPs.
- Other GWAS inputs in each four-trait LDSC analysis were unchanged.
- Saved LDSC RDS files were used to verify downstream GenomicSEM and Wald calculations.

Both MAF-sensitivity omnibus results met P < 0.025.

However, no individual pairwise contrast survived Holm correction in this sensitivity analysis.

The original overall-endometriosis univariate LDSC intercept was 0.7672 (SE 0.0098). The cause of this deflation remains unresolved.

Detailed validation:

- `docs/maf05_rgmodel_validation.md`
- `docs/reproducibility_notes.md`

## Secondary analysis: conjunctional FDR

Six phenotype-pair comparisons were evaluated using a common matched SNP universe.

Key settings:

- 1,057,700 matched HapMap3 SNPs
- MHC exclusion: chr6:26-34 Mb
- Random pruning: `randprune_n = 20`
- Regional grouping based on a maximum 500-kb gap between lead-SNP records and the current region end

The MHC exclusion applies to conjFDR processing, not to LDSC analyses.

Regional grouping yielded:

- 92 lead-SNP records
- 45 grouped genomic regions
- 22 regions detected in more than one phenotype-pair comparison

The six comparisons are statistically dependent. Repeated detection therefore does not constitute independent replication.

Highlighted recurrent-detection regions include WNT4, SKAP1, chr8 TMEM75/PVT1, and chr12 SSPN-BHLHE41.

These labels represent statistical regions or candidate gene annotations, not confirmed causal genes.

## Exploratory chr12 colocalization

A region indexed by rs9668810 near SSPN-BHLHE41 was selected for exploratory follow-up.

The regional window was rs9668810 +/-500 kb.

The final reference analysis used:

**Adenomyosis**

- P-values and minor allele frequencies
- Median eligible variant-level `AD_N`: 33,460.65
- Assumed case fraction: `s = 0.5`

**Endometrial cancer**

- Aligned beta estimates
- Squared standard errors as `varbeta`
- allEC N = 121,885
- EEC N = 54,884
- Observed outcome-specific case fractions

Colocalization priors were:

- `p1 = 1e-4`
- `p2 = 1e-4`
- `p12 = 1e-5`

### Conditional colocalization results

| Outcome | SNPs | PP.H4 |
|---|---:|---:|
| allEC | 3,159 | 0.9931782 |
| EEC | 3,158 | 0.9880625 |

rs7959150 had the highest conditional SNP.PP.H4 in both analyses.

These posterior probabilities are conditional on the model inputs, priors, and single-causal-variant-per-trait assumption.

**Important unresolved limitation:** the upstream definition of adenomyosis `AD_N` and the appropriateness of `s = 0.5` require further methodological clarification.

High PP.H4 must not be interpreted as proof of an individual causal variant or causal gene.

The historical colocalization script and updated final reference script are intentionally distinguished:

- Historical: `scripts/06_run_coloc_chr12_rs9668810.R`
- Final reference: `scripts/14_coloc_chr12_final_reference.R`

## Functional annotation

GTEx V10 eQTL follow-up was used for exploratory regional functional annotation.

The repository retains scripts for:

- rs9668810 eQTL querying
- Targeted uterus eQTL auditing

eQTL associations do not establish causal mediation, a causal gene, or a causal tissue.

## Repository structure

```text
raw/                       Source GWAS files; not distributed
ref/                       LD/reference resources; not distributed
munged/                    LDSC inputs; not distributed
results/                   Intermediate and final results
scripts/                   Analysis and validation scripts
scripts/clean/             GTEx V10 follow-up scripts
docs/                      Data, provenance and validation documents
README.md                  Project overview
```

Some small GTEx-related outputs are tracked within `results/`. Large inputs and intermediate outputs are not included.

## Scripts

### Original workflow and historical scripts

```text
scripts/00_config.sh
scripts/01_prepare_project_and_ldsc_inputs.sh
scripts/02_run_ldsc_qc_and_rg.sh
scripts/03_harmonize_conjfdr_inputs.sh
scripts/04_prepare_pleiofdr_mat_inputs.sh
scripts/05_prepare_decomposition_common_snps.sh
scripts/06_run_coloc_chr12_rs9668810.R
scripts/07_make_ldsc_tables_from_logs.py
scripts/08_collect_decomp_raw_loci.py
scripts/09_group_decomp_loci_table.py
```

The output of script 07 includes historical approximate genetic-correlation difference tests that are not the final covariance-aware inference.

### Added final-analysis and validation scripts

```text
scripts/12_rgmodel_final_tests.R
scripts/13_rgmodel_sensitivity.R
scripts/14_coloc_chr12_final_reference.R
scripts/15_maf05_rgmodel_verification.R
```

Their purposes are:

- **12:** Final covariance-aware omnibus and pairwise Wald calculations from saved GenomicSEM results
- **13:** LDSC block-sensitivity follow-up from saved LDSC results
- **14:** Updated exploratory chr12 colocalization reference implementation
- **15:** ENDO-only MAF-sensitivity downstream verification

These scripts were added for transparency and should not be represented as unchanged original execution scripts.

All four scripts passed R syntax parsing.

Script 15 was executed against saved LDSC RDS files and reproduced the previously reported validation metrics within the specified tolerance.

Modified public copies of scripts 13 and 14 have not been independently verified end-to-end after their repository adjustments.

### GTEx scripts

```text
scripts/clean/09_query_gtex_v10_rs9668810.py
scripts/clean/10_make_gtex_uterus_final_audit.py
```

## Configuration and execution

The existing workflow uses `scripts/00_config.sh` for local configuration.

Example variables:

```bash
PROJECT=/path/to/endo_ecancer
RAW=/path/to/endo_ecancer/raw
REF=/path/to/endo_ecancer/ref
MUNGED=/path/to/endo_ecancer/munged
RESULTS=/path/to/endo_ecancer/results
CJ=/path/to/endo_ecancer/results/conjfdr
PLEIO=/path/to/pleiofdr
```

Before running any script, review its input requirements, existing outputs, and configuration settings.

Some scripts operate on precomputed results; other scripts prepare or regenerate files and may require substantial computation.

The repository does not currently provide a fully verified, one-command end-to-end reproduction pipeline.

### MAF verification example

The MAF verification script reads existing LDSC RDS objects. It does not rerun GWAS munging or LDSC.

```bash
RUN_MAF05_VERIFICATION=1 \
MAF_LDSC_DIR=/path/to/maf05 \
Rscript scripts/15_maf05_rgmodel_verification.R
```

The local installation may additionally require `GSEM_LIB` to identify the GenomicSEM R library.

The input directory must contain:

```text
ldsc_allEC_4trait.rds
ldsc_EEC_4trait.rds
```

Running this script requires the original intermediate RDS files, which are not distributed in the repository.

## Expected outputs and provenance

The workflow uses or generates the following conjFDR decomposition tables:

```text
results/tables/table_s9_decomp_conjfdr_loci_raw.tsv
results/tables/table_s10_decomp_locus_grouped_from_file.tsv
```

Other retained downstream tables include LDSC QC summaries, phenotype-pair statistics, directionality classifications, and regional summaries.

Internal filenames do not necessarily match final manuscript supplementary-table numbering.

Detailed output descriptions are available in `docs/final_output_tables.md`.

GTEx V10 output examples:

```text
results/reproduced_gtex/gtex_v10_rs9668810_significant_eqtl.tsv
results/reproduced_gtex/gtex_v10_rs9668810_uterus_dynamic_eqtl.tsv
results/reproduced_gtex/gtex_v10_rs9668810_uterus_final_audit.tsv
results/reproduced_gtex/gtex_v10_rs9668810_query_metadata.json
```

These paths describe expected workflow outputs and do not imply that every file is included in the repository.

## Software environment

Software used in the analysis included:

- LDSC v1.0.1 with Python 2.7.18
- pleioFDR/conjFDR with MATLAB R2026a
- R and Python 3
- GenomicSEM 0.0.5
- lavaan 0.6-19
- coloc 5.2.3
- data.table and other R dependencies
- requests and other Python dependencies

An exact cross-language environment lockfile is not currently available.

Package startup warnings observed during validation have not all been independently classified.

Software, dependency, and reference-resource availability must be checked before external reproduction.

## Documentation

Additional files:

- `docs/code_inventory.md`
- `docs/data_sources.md`
- `docs/final_output_tables.md`
- `docs/interpretation_limits.md`
- `docs/maf05_rgmodel_validation.md`
- `docs/reproducibility_notes.md`
- `docs/sample_size_derivation.md`

The documents describe data provenance, output traceability, statistical methods, validation results, and remaining limitations.

## Interpretation limits

The study uses summary-level GWAS data and does not provide individual-level clinical prediction or causal inference.

Important limitations include:

- Phenotype heterogeneity and differences in source GWAS power
- Dependence between genetic-correlation estimates
- Unresolved overall-endometriosis LDSC intercept deflation
- Incomplete provenance of adenomyosis effective sample size for colocalization
- Sensitivity of colocalization posterior support to modeling assumptions
- Non-independent recurrent-locus detections
- Exploratory rather than experimentally validated functional annotations

See `docs/interpretation_limits.md` for a detailed interpretation framework.

## Citation

Users should cite the original GWAS publications and datasets according to the source providers' instructions.

Please cite the associated manuscript once available.

## Contact

For questions about the workflow, contact the corresponding author of the associated manuscript.
