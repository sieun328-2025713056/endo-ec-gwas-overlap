# Endometriosis–Adenomyosis–Endometrial Cancer GWAS Summary-Statistics Workflow

This repository provides a curated reproducibility workflow for a GWAS summary-statistics analysis of genetic overlap among endometriosis, adenomyosis, and endometrial cancer.

The code package is designed to support manuscript-level reproducibility while avoiding redistribution of restricted raw GWAS summary statistics, large LD/reference files, generated `.mat` files, and local pleioFDR/conjFDR output folders.

## Purpose

This workflow supports the following analysis components:

- Preparation of endometrial cancer summary statistics for LD Score Regression (LDSC).
- LDSC munging, SNP-heritability quality control, and genetic correlation analyses.
- Harmonization of GWAS summary statistics for conjFDR/pleioFDR analyses.
- Verification or optional generation of pleioFDR-compatible `.mat` trait files.
- Verification of matched SNP-universe inputs for the six-pair decomposition analysis.
- Exploratory regional colocalization around the chr12 rs9668810-indexed region.
- Reproducible GTEx V10 eQTL follow-up of the chr12 rs9668810-indexed signal.
- Generation or verification of manuscript-related tables from LDSC, pleioFDR, decomposition, colocalization, and functional-annotation outputs.

## Repository status

This is a cleaned code package. It excludes:

- raw GWAS summary statistics;
- large LDSC and pleioFDR reference files;
- generated `.mat` files;
- pleioFDR result folders;
- intermediate large result files;
- local troubleshooting logs;
- WSL installation notes;
- raw conversation extracts.

Several downstream scripts are intentionally written as verify-first scripts. This is because pleioFDR/conjFDR conversion and decomposition output generation can be computationally and I/O intensive, and local installations may differ in output naming conventions.

## Data availability and restrictions

Raw GWAS summary statistics and large reference files are not redistributed in this repository. Users must obtain the original source data from the respective data providers and comply with all applicable data-use terms and restrictions.

The endometriosis-spectrum GWAS summary statistics are available through Zenodo:

- Koller et al. endometriosis-spectrum GWAS summary statistics: DOI `10.5281/zenodo.18983492`

The endometrial cancer GWAS summary-statistics datasets correspond to the following GWAS Catalog accessions:

- All-histology endometrial cancer: `GCST006464`
- Endometrioid endometrial cancer: `GCST006465`

The overall endometriosis analysis in this workflow uses the publicly released European combined GWAS summary-statistics file excluding the 23andMe component. The adenomyosis and endometriosis-without-adenomyosis analyses use the corresponding European-ancestry phenotype-specific summary statistics.

Expected local raw input files are:

```text
raw/
├── adenomyosis_EUR.txt
├── endometriosis_EUR_wo_23andMe.txt
├── endometriosis_wo.adenomyosis_EUR.txt
├── ECAC_allEC.valid_2.tsv
└── ECAC_EEC.valid_2.tsv
```

These filenames reflect the local analysis environment and are provided only to document the expected workflow inputs; the source GWAS files themselves are not included in this repository.

Users should cite the corresponding source GWAS publications and follow the data-use policies of the original data providers when obtaining and using these summary statistics.

## Directory structure

The expected project structure is:

```text
endo_ecancer/
├── raw/                  # restricted raw GWAS summary statistics; not included
├── ref/                  # LD/reference files; not included
├── munged/               # LDSC munged summary statistics
├── results/
│   ├── tables/           # workflow summary tables
│   ├── coloc/            # regional colocalization outputs
│   ├── conjfdr/          # large conjFDR-related intermediate files; not included
│   └── reproduced_gtex/  # reproduced GTEx V10 follow-up outputs
├── scripts/
│   └── clean/            # GTEx V10 reproducibility scripts
├── docs/                 # workflow documentation
└── README.md
```

## Configuration

Before running the workflow, review `scripts/00_config.sh` and either modify the default paths or override them using environment variables.

The configuration file defines local paths such as:

```bash
PROJECT=/path/to/endo_ecancer
RAW=/path/to/endo_ecancer/raw
REF=/path/to/endo_ecancer/ref
MUNGED=/path/to/endo_ecancer/munged
RESULTS=/path/to/endo_ecancer/results
CJ=/path/to/endo_ecancer/results/conjfdr
PLEIO=/path/to/pleiofdr
```

These paths reflect the local analysis environment and can be modified by each user before running the workflow.

## Software environment

The primary analysis workflow used the following software and computational environments:

- LD Score Regression (LDSC) v1.0.1 in a Python 2.7.18 environment for summary-statistics munging, SNP-heritability quality control, and genetic correlation analyses.
- pleioFDR/conjFDR with MATLAB R2026a for conjunctional FDR-related processing and `.mat` file handling.
- R with `coloc` v5.2.3 and `data.table` for exploratory P-value-based regional colocalization.
- Python 3, including the `requests` package, for table processing, recurrent-locus grouping, and GTEx V10 API-based eQTL follow-up.

The workflow also uses externally obtained LD/reference resources and pleioFDR/conjFDR software that are not redistributed in this repository. Local software and resource paths can be configured through `scripts/00_config.sh` and relevant environment variables.

## Main workflow

The analysis workflow is organized as numbered scripts:

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
scripts/clean/09_query_gtex_v10_rs9668810.py
scripts/clean/10_make_gtex_uterus_final_audit.py
```

A typical workflow is:

```bash
source scripts/00_config.sh

bash scripts/01_prepare_project_and_ldsc_inputs.sh
bash scripts/02_run_ldsc_qc_and_rg.sh
bash scripts/03_harmonize_conjfdr_inputs.sh
bash scripts/04_prepare_pleiofdr_mat_inputs.sh
bash scripts/05_prepare_decomposition_common_snps.sh
Rscript scripts/06_run_coloc_chr12_rs9668810.R
python3 scripts/07_make_ldsc_tables_from_logs.py
python3 scripts/08_collect_decomp_raw_loci.py
python3 scripts/09_group_decomp_loci_table.py
python3 scripts/clean/09_query_gtex_v10_rs9668810.py results/reproduced_gtex
python3 scripts/clean/10_make_gtex_uterus_final_audit.py
```

## Verify-first downstream workflow

Several downstream scripts verify existing intermediate files and output tables by default instead of regenerating large files:

```bash
bash scripts/04_prepare_pleiofdr_mat_inputs.sh
bash scripts/05_prepare_decomposition_common_snps.sh
python3 scripts/08_collect_decomp_raw_loci.py
python3 scripts/09_group_decomp_loci_table.py
```

Optional regeneration commands are:

```bash
RUN_CONVERSION=1 bash scripts/04_prepare_pleiofdr_mat_inputs.sh
RUN_REBUILD=1 python3 scripts/08_collect_decomp_raw_loci.py
RUN_REBUILD=1 python3 scripts/09_group_decomp_loci_table.py
```

Scripts 08 and 09 do not overwrite existing output tables if zero rows are detected during rebuild.

## Expected verified outputs

Expected verified pleioFDR and decomposition outputs include:

```text
$PLEIO/traitfiles/adenomyosis_noMHC.mat
$PLEIO/traitfiles/allEC_noMHC.mat
$PLEIO/traitfiles/adenomyosis_EECpair_noMHC.mat
$PLEIO/traitfiles/EEC_noMHC.mat
$PLEIO/traitfiles_decomp_common/

results/tables/table_s9_decomp_conjfdr_loci_raw.tsv
results/tables/table_s10_decomp_locus_grouped_from_file.tsv
```

## Workflow output tables

The workflow generates or uses several intermediate and downstream summary tables. Internal filenames such as `table_s1`, `table_s7`, and `table_s14` reflect workflow-development labels and do not necessarily correspond to the final numbering of manuscript supplementary tables.

LDSC-derived tables generated by `scripts/07_make_ldsc_tables_from_logs.py` include:

```text
results/tables/table_s1_ldsc_h2_qc_from_logs.tsv
results/tables/table_s2_ldsc_rg_from_logs.tsv
results/tables/table_s3_rg_approx_difference_test_from_logs.tsv
results/tables/table_s7_ldsc_rg_multiple_testing_correction.tsv
```

The multiple-testing output applies Benjamini–Hochberg FDR correction across the six primary cross-trait genetic-correlation tests. The approximate between-phenotype difference tests use `sqrt(SE1^2 + SE2^2)` and do not model covariance arising from the shared endometrial cancer GWAS or non-independent endometriosis-spectrum estimates.

Pair-specific and recurrent-locus outputs directly verified or rebuilt by the public decomposition scripts include:

```text
results/tables/table_s9_decomp_conjfdr_loci_raw.tsv
results/tables/table_s10_decomp_locus_grouped_from_file.tsv
```

Additional downstream summary tables retained for workflow traceability include:

```text
results/tables/table_s8_decomp_conjfdr_pair_counts.tsv
results/tables/table_s11_key_recurrent_loci_decomposition.tsv
results/tables/table_s12_pair_level_directionality_summary.tsv
results/tables/table_s13_grouped_locus_classification_direction_summary.tsv
results/tables/table_s14_key_recurrent_loci_compact_summary.tsv
results/tables/table_s15_phenotype_pairwise_rg.tsv
```

Not all downstream summary tables are regenerated by the cleaned public code package. Detailed output provenance is documented in `docs/final_output_tables.md`.

## Regional colocalization outputs

`scripts/06_run_coloc_chr12_rs9668810.R` generates separate outputs for the adenomyosis–all-histology endometrial cancer and adenomyosis–endometrioid endometrial cancer analyses.

Expected outputs include:

```text
results/coloc/rs9668810_region/adeno_allEC_rs9668810_500kb_coloc_summary.tsv
results/coloc/rs9668810_region/adeno_allEC_rs9668810_500kb_coloc_full_results.tsv
results/coloc/rs9668810_region/adeno_allEC_rs9668810_500kb_top_H4_snps.tsv
results/coloc/rs9668810_region/adeno_EEC_rs9668810_500kb_coloc_summary.tsv
results/coloc/rs9668810_region/adeno_EEC_rs9668810_500kb_coloc_full_results.tsv
results/coloc/rs9668810_region/adeno_EEC_rs9668810_500kb_top_H4_snps.tsv
```

## GTEx V10 functional-annotation outputs

The GTEx V10 follow-up scripts generate or summarize:

```text
results/reproduced_gtex/gtex_v10_rs9668810_significant_eqtl.tsv
results/reproduced_gtex/gtex_v10_rs9668810_uterus_dynamic_eqtl.tsv
results/reproduced_gtex/gtex_v10_rs9668810_uterus_final_audit.tsv
results/reproduced_gtex/gtex_v10_rs9668810_query_metadata.json
```

The GTEx V10 follow-up is exploratory. The outputs are used to document functional-annotation results and do not establish a causal gene or tissue-specific mechanism.

## Documentation

Additional documentation is provided in the `docs/` directory:

```text
docs/code_inventory.md
docs/data_sources.md
docs/final_output_tables.md
docs/interpretation_limits.md
docs/reproducibility_notes.md
docs/sample_size_derivation.md
```

These files document source GWAS datasets, redistribution limits, derivation of the 23andMe-excluded overall endometriosis sample size, code provenance, expected workflow outputs, reproducibility notes, and interpretation limits for LDSC, conjFDR, colocalization, and functional annotation.

## Interpretation limits

This workflow supports summary-statistics-based genetic epidemiology analyses. The results should not be interpreted as establishing individual-level risk, clinical prediction, causal variants, causal genes, or biological mechanisms.

Non-detection of a locus in a phenotype pair indicates lack of detection under the current GWAS power, phenotype definition, matched SNP universe, and conjFDR threshold. It should not be interpreted as evidence that the underlying association is absent.

Candidate gene names used in output tables should be interpreted as regional annotations or prioritization labels, not as confirmed causal genes.

Approximate between-phenotype genetic-correlation difference tests do not account for covariance arising from the shared endometrial cancer GWAS or non-independent endometriosis-spectrum estimates and should therefore be interpreted as exploratory.

Regional colocalization and GTEx V10 eQTL follow-up are exploratory and do not establish a causal variant, causal gene, or tissue-specific biological mechanism.

## Citation

Users should cite the original GWAS publications for all source summary statistics and should follow the data-use policies of the original data providers.

If using this workflow, please cite the associated manuscript once available.

## Contact

For questions about the workflow, please contact the corresponding author of the associated manuscript.