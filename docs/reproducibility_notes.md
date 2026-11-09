
# Reproducibility Notes

## Scope

This repository documents the GWAS summary-statistics workflow used to investigate genetic overlap between endometriosis-spectrum phenotypes and endometrial cancer.

The public repository is not a self-contained computational environment. Reproduction requires the original GWAS summary statistics, reference resources, software dependencies, and intermediate files that are not distributed here.

Historical scripts are retained where useful for provenance. The final inferential methods and exploratory follow-up are distinguished from earlier analyses.

## Raw data not included

The repository does not distribute:

- Raw GWAS summary statistics
- Large LD/reference files
- Generated MATLAB `.mat` files
- Large pleioFDR/conjFDR output directories
- Private working files or unpublished source datasets

Users must obtain source data from the original providers and comply with applicable access and reuse conditions.

## Main raw files expected

Local input filenames:

- `adenomyosis_EUR.txt`
- `endometriosis_EUR_wo_23andMe.txt`
- `endometriosis_wo.adenomyosis_EUR.txt`
- `ECAC_allEC.valid_2.tsv`
- `ECAC_EEC.valid_2.tsv`

These are local analysis filenames, not necessarily the original provider filenames.

## Analysis structure

The investigation includes three endometriosis-spectrum phenotypes:

1. Overall endometriosis excluding 23andMe
2. Adenomyosis
3. Endometriosis without adenomyosis

Each was evaluated against:

1. All-histology endometrial cancer (allEC)
2. Endometrioid endometrial cancer (EEC)

The analysis hierarchy was:

1. PRIMARY: genome-wide LDSC genetic correlations and covariance-aware heterogeneity testing
2. SECONDARY: six matched conjunctional FDR analyses and regional detection-pattern classification
3. EXPLORATORY: selected chr12 colocalization and GTEx V10 annotation

The exploratory chr12 region was selected following the secondary locus-level analysis. It was not an independent prespecified genome-wide primary analysis.

## Primary inference: LDSC and GenomicSEM

Six bivariate LDSC genetic correlations were estimated.

The final primary hypothesis tests were two covariance-aware omnibus Wald tests comparing the three endometriosis-spectrum genetic correlations within each cancer outcome.

GenomicSEM `rgmodel()` supplied:

- `R`: estimated genetic correlations
- `V_R`: sampling covariance matrix of those estimates

The Wald tests accounted for the covariance among genetic-correlation estimates rather than treating them as independent.

The two primary omnibus tests were evaluated using Bonferroni correction:

- Number of primary tests: 2
- Significance threshold: P < 0.025

Six post-hoc pairwise contrasts, three for each cancer outcome, were evaluated with a single Holm correction across all six tests.

### Final primary results

| Outcome | Omnibus ChiSq(2) | P |
|---|---:|---:|
| allEC | 8.155546 | 0.016945159 |
| EEC | 9.632358 | 0.008097669 |

Both primary omnibus tests met the corrected threshold.

For the ENDO-ADENO pairwise comparison:

| Outcome | Raw P | Holm P |
|---|---:|---:|
| allEC | 0.012022315 | 0.06011158 |
| EEC | 0.006044072 | 0.03626443 |

These results should not be described as showing that all individual phenotype pairs differ significantly.

Earlier approximate difference tests based on independent-standard-error assumptions are historical exploratory analyses, not the final inferential results.

## Primary analysis sensitivity checks

Sampling-covariance estimation was examined using different LDSC block settings.

| Setting | allEC omnibus P | EEC omnibus P |
|---|---:|---:|
| 100 blocks | 0.01151326 | 0.007853592 |
| 200 blocks | 0.016945159 | 0.008097669 |
| 400 blocks | 0.02482712 | 0.01197984 |
| ENDO-only MAF filter | 0.02458479 | 0.01214624 |

The 200-block setting represents the principal reported analysis.

### ENDO-only MAF sensitivity

The overall-endometriosis LDSC input was regenerated using `--maf-min 0.05`.

The resulting munged ENDO input retained 1,064,328 SNPs.

The other three GWAS inputs in each four-trait LDSC configuration were unchanged.

The resulting saved LDSC RDS objects were used to recompute GenomicSEM `rgmodel()` outputs, omnibus Wald tests, and Holm-adjusted pairwise contrasts.

The two MAF-sensitivity omnibus results met P < 0.025.

However, no individual pairwise comparison survived Holm correction in this sensitivity analysis.

The ENDO univariate LDSC intercept was 0.7672 (SE 0.0098). The cause of this substantial deflation remains unresolved.

The downstream MAF-sensitivity calculations were independently checked against previously reported values, with absolute differences below 1e-5.

See:

- `scripts/15_maf05_rgmodel_verification.R`
- `docs/maf05_rgmodel_validation.md`

This verification used existing LDSC RDS objects and did not rerun the original GWAS processing or LDSC estimation.

## Secondary analysis: conjunctional FDR

Key settings:

- Six phenotype-pair comparisons
- Common matched SNP universe: 1,057,700 HapMap3 SNPs
- MHC exclusion: chr6:26-34 Mb
- pleioFDR/conjFDR random pruning: `randprune_n = 20`
- Regional grouping: lead-SNP records grouped using a gap of <=500 kb from the current region end

MHC exclusion was applied to conjunctional FDR processing, not to the LDSC analyses.

Regional grouping produced:

- 92 lead-SNP records
- 45 grouped genomic regions
- 22 regions detected in more than one phenotype-pair comparison

These comparisons are statistically dependent. Repeated detection must not be interpreted as independent replication.

The 45 grouped regions should not all be described as recurrent across phenotype pairs.

## Exploratory chr12 colocalization

The selected follow-up region was indexed by rs9668810 in the SSPN-BHLHE41 region.

Regional window:

- rs9668810 +/-500 kb

The current final analysis differs from the historical `scripts/06_run_coloc_chr12_rs9668810.R` workflow.

### Final input specification

Adenomyosis:

- P-values
- Minor allele frequencies
- N = median eligible variant-level `AD_N` = 33460.65
- Assumed case fraction `s = 0.5`

The underlying derivation and interpretation of `AD_N` as an effective sample size have not been fully established.

The value `s = 0.5` was a modeling assumption, not the observed adenomyosis case fraction.

The observed case fraction derived from the reported overall case/control counts is approximately 0.02062. This historical value should not be mistaken for the `s` used in the final colocalization analysis.

Endometrial cancer:

- Aligned effect estimates
- Squared standard errors as `varbeta`
- allEC N = 121885
- EEC N = 54884
- allEC case fraction = 12906 / 121885
- EEC case fraction = 8758 / 54884

Colocalization priors:

- `p1 = 1e-4`
- `p2 = 1e-4`
- `p12 = 1e-5`

### Final conditional results

| Outcome | Eligible SNPs | PP.H4 |
|---|---:|---:|
| allEC | 3159 | 0.9931782 |
| EEC | 3158 | 0.9880625 |

The highest conditional `SNP.PP.H4` was observed for rs7959150:

- allEC: 0.31654469
- EEC: 0.27609080

rs9668810 is the region-indexing SNP, not the highest-posterior SNP.

These results are conditional on the input assumptions, prior probabilities, SNP sets, and the coloc.abf assumption of at most one causal variant per trait in the evaluated region.

They do not establish a causal variant or causal gene.

The updated reference implementation is:

- `scripts/14_coloc_chr12_final_reference.R`

This script can optionally export SNP-level posterior tables when explicitly enabled. The required input files are not distributed in the repository.

## Public scripts and provenance

### Existing workflow scripts

- `scripts/00_config.sh`: configuration and path overrides
- `scripts/03_harmonize_conjfdr_inputs.sh`: harmonization utilities
- `scripts/04_prepare_pleiofdr_mat_inputs.sh`: MATLAB input preparation
- `scripts/05_prepare_decomposition_common_snps.sh`: common-SNP preparation
- `scripts/06_run_coloc_chr12_rs9668810.R`: historical colocalization implementation
- `scripts/07_make_ldsc_tables_from_logs.py`: LDSC summaries and historical approximate comparisons
- `scripts/08_collect_decomp_raw_loci.py`: lead-SNP record collection
- `scripts/09_group_decomp_loci_table.py`: regional grouping

The historical approximate contrast calculations in script 07 are not the final covariance-aware inference.

### Added final-analysis reference scripts

- `scripts/12_rgmodel_final_tests.R`: covariance-aware omnibus and Holm-adjusted pairwise tests from saved GenomicSEM results
- `scripts/13_rgmodel_sensitivity.R`: block-sensitivity calculations from saved LDSC results
- `scripts/14_coloc_chr12_final_reference.R`: updated chr12 colocalization inputs and posterior calculations
- `scripts/15_maf05_rgmodel_verification.R`: independent downstream verification of MAF-sensitivity calculations

Scripts 12-15 were added to improve transparency and reproducibility. They should not be represented as the unchanged original execution scripts.

Scripts 12-15 passed R syntax parsing.

Script 15 was executed and its reported-value comparisons passed. The modified copies of scripts 13 and 14 have not been established as independently rerun end-to-end after their public-repository adjustments.

## Execution environment and file paths

The workflow originally used Ubuntu/WSL and a mixture of R, Python, MATLAB, LDSC, and pleioFDR tools.

The validated GenomicSEM environment included:

- GenomicSEM 0.0.5
- lavaan 0.6-19

An exact complete software lockfile has not been established.

GenomicSEM package startup produced 23 warnings during the MAF verification run; their individual contents have not yet been reviewed.

Scripts may depend on external tools, upstream results, and local configuration. Availability of source data and software must be verified before independent execution.

For the MAF verification script:

- `RUN_MAF05_VERIFICATION=1` enables execution
- `MAF_LDSC_DIR` can override the directory containing saved MAF LDSC RDS files
- `GSEM_LIB` can specify the R package library location

For the final colocalization reference script, relevant environment settings include:

- `PROJECT`
- `COLOC_INPUT_DIR`
- `EXPORT_COLOC_POSTERIOR`
- `COLOC_POSTERIOR_DIR`

Posterior output is optional and must be explicitly enabled.

## Functional annotation

GTEx V10 analyses were exploratory.

The public workflow includes scripts for rs9668810 regional eQTL follow-up and targeted uterus eQTL auditing.

eQTL association does not establish causal mediation, a causal gene, or a disease-relevant tissue of action.

## Interpretation and reproducibility limitations

- The work relies on summary statistics rather than individual-level data.
- Phenotype definitions and GWAS sample sizes differ across source studies.
- Some analyses use statistically dependent phenotype-pair comparisons.
- The overall-endometriosis LDSC intercept deflation remains unresolved.
- Adenomyosis `AD_N` provenance and the appropriateness of `s = 0.5` remain open methodological limitations.
- Colocalization posterior probabilities depend on modeling assumptions and priors.
- Non-detection of a locus does not establish biological absence.
- Genomic regions and genes identified through statistical prioritization are not validated causal mechanisms.
- The public scripts and documentation do not, by themselves, guarantee independent end-to-end reproduction without the missing input data and dependencies.

For detailed inferential cautions, see `docs/interpretation_limits.md`.
