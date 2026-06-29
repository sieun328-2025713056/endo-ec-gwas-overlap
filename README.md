# Endometriosis–Adenomyosis–Endometrial Cancer GWAS Summary-Statistics Workflow

This repository provides a curated reproducibility workflow for a GWAS summary-statistics analysis of genetic overlap among endometriosis, adenomyosis, and endometrial cancer.

The code package is designed to support manuscript-level reproducibility while avoiding redistribution of restricted raw GWAS summary statistics, large LD/reference files, generated `.mat` files, and local pleioFDR/conjFDR output folders.

## Purpose

This workflow supports the following analysis components:

* Preparation of endometrial cancer summary statistics for LD Score regression (LDSC).
* LDSC munging, SNP-heritability quality control, and genetic correlation analyses.
* Harmonization of GWAS summary statistics for conjFDR / pleioFDR analyses.
* Verification or optional generation of pleioFDR-compatible `.mat` trait files.
* Verification of matched SNP-universe inputs for the six-pair decomposition analysis.
* Optional regional colocalization around the chr12 rs9668810-indexed region.
* Generation or verification of manuscript tables from LDSC, pleioFDR, and decomposition outputs.

## Repository status

This is a cleaned code package. It excludes:

* raw GWAS summary statistics;
* large LDSC and pleioFDR reference files;
* generated `.mat` files;
* pleioFDR result folders;
* intermediate large result files;
* local troubleshooting logs;
* WSL installation notes;
* raw conversation extracts.

Several downstream scripts are intentionally written as verify-first scripts. This is because pleioFDR/conjFDR conversion and decomposition output generation can be computationally and I/O intensive, and local installations may differ in output naming conventions.

## Data availability and restrictions

Raw GWAS summary statistics and reference files are not redistributed in this repository. Users must obtain the original data from the respective data providers and comply with all applicable data-use restrictions.

Expected raw input files are:

```text
raw/
├── adenomyosis_EUR.txt
├── endometriosis_EUR_wo_23andMe.txt
├── endometriosis_wo.adenomyosis_EUR.txt
├── ECAC_allEC.valid_2.tsv
└── ECAC_EEC.valid_2.tsv
```

The endometriosis-spectrum GWAS summary statistics were obtained from Koller et al. The accompanying data-use terms specify that the data are provided for scientific research purposes and that users should cite the corresponding publication. Endometrial cancer GWAS summary statistics should likewise be obtained from the original source and used according to the relevant data-use policies.

## Directory structure

The expected project structure is:

```text
endo_ecancer/
├── raw/                 # restricted raw GWAS summary statistics; not included
├── ref/                 # LD/reference files; not included
├── munged/              # LDSC munged summary statistics
├── results/
│   ├── tables/          # manuscript and supplementary tables
│   └── conjfdr/         # large conjFDR-related intermediate files; not included
├── scripts/             # reproducibility scripts
├── docs/                # documentation files
└── README.md
```

## Configuration

Before running the workflow, edit:

```text
scripts/00_config.sh
```

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

These paths reflect the local analysis environment and should be modified by each user before running the workflow.

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

Expected verified outputs include:

```text
$PLEIO/traitfiles/adenomyosis_noMHC.mat
$PLEIO/traitfiles/allEC_noMHC.mat
$PLEIO/traitfiles/adenomyosis_EECpair_noMHC.mat
$PLEIO/traitfiles/EEC_noMHC.mat
$PLEIO/traitfiles_decomp_common/

results/tables/table_s9_decomp_conjfdr_loci_raw.tsv
results/tables/table_s10_decomp_locus_grouped_from_file.tsv
```

Expected final manuscript or supplementary tables include:

```text
results/tables/table_s1_ldsc_h2_qc_from_logs.tsv
results/tables/table_s2_ldsc_rg_from_logs.tsv
results/tables/table_s3_rg_approx_difference_test_from_logs.tsv
results/tables/table_s7_ldsc_rg_multiple_testing_correction.tsv
results/tables/table_s8_decomp_conjfdr_pair_counts.tsv
results/tables/table_s9_decomp_conjfdr_loci_raw.tsv
results/tables/table_s10_decomp_locus_grouped_from_file.tsv
results/tables/table_s11_key_recurrent_loci_decomposition.tsv
results/tables/table_s12_pair_level_directionality_summary.tsv
results/tables/table_s13_grouped_locus_classification_direction_summary.tsv
results/tables/table_s14_key_recurrent_loci_compact_summary.tsv
results/tables/table_s15_phenotype_pairwise_rg.tsv
```

## Documentation

Additional documentation is provided in the `docs/` directory:

```text
docs/data_sources.md
docs/sample_size_derivation.md
docs/final_output_tables.md
docs/interpretation_limits.md
docs/code_inventory.md
```

These files document source GWAS datasets, redistribution limits, derivation of the 23andMe-excluded overall endometriosis sample size, expected final tables, and interpretation limits for LDSC, conjFDR, and colocalization analyses.

## Interpretation limits

This workflow supports summary-statistics-based genetic epidemiology analyses. The results should not be interpreted as establishing individual-level risk, clinical prediction, causal variants, causal genes, or biological mechanisms.

Non-detection of a locus in a phenotype pair indicates lack of detection under the current GWAS power, phenotype definition, matched SNP universe, and conjFDR threshold. It should not be interpreted as evidence that the underlying association is absent.

Candidate gene names used in output tables should be interpreted as regional annotations or prioritization labels, not as confirmed causal genes.

## Citation

Users should cite the original GWAS publications for all source summary statistics and should follow the data-use policies of the original data providers.

If using this workflow, please cite the associated manuscript once available.

## Contact

For questions about the workflow, please contact the corresponding author of the associated manuscript.
