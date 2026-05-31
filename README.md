nano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.mdnano README_code_package.md# Endometriosis–Adenomyosis–Endometrial Cancer GWAS Summary Statistics Workflow

This repository provides a curated reproducibility workflow for the GWAS summary-statistics analysis of genetic overlap among endometriosis, adenomyosis, and endometrial cancer.

The code package is organized to support manuscript-level reproducibility while avoiding redistribution of restricted raw GWAS summary statistics, large LD/reference files, and generated pleioFDR outputs.

## Purpose

The workflow supports the following analysis components:

1. Preparation of endometrial cancer summary statistics for LDSC.
2. LDSC munging, SNP-heritability QC, and genetic correlation analyses.
3. Harmonization of GWAS summary statistics for conjFDR / pleioFDR.
4. Verification or optional generation of pleioFDR-compatible `.mat` trait files.
5. Verification of matched SNP-universe inputs for the six-pair decomposition analysis.
6. Optional regional colocalization around chr12 rs9668810.
7. Generation or verification of manuscript tables from LDSC, pleioFDR, and decomposition outputs.

## Repository status

This is a cleaned code package. It excludes:

- raw GWAS summary statistics,
- large LDSC and pleioFDR reference files,
- generated `.mat` files,
- pleioFDR result folders,
- intermediate result tables,
- local troubleshooting logs,
- WSL installation notes,
- raw conversation extracts.

Several downstream scripts are intentionally **verify-first**. This is because pleioFDR/conjFDR conversion and decomposition output generation can be computationally and I/O intensive, and local installations may differ in output naming conventions.

## Data availability

Raw GWAS summary statistics and reference files are not redistributed in this repository. Users must obtain the original data from the respective providers and comply with all data-use restrictions.

Expected raw files:

```text
raw/
├── adenomyosis_EUR.txt
├── endometriosis_EUR_wo_23andMe.txt
├── endometriosis_wo.adenomyosis_EUR.txt
├── ECAC_allEC.valid_2.tsv
└── ECAC_EEC.valid_2.tsv
```

## Important reproducibility notes

This code package uses a verify-first workflow for several downstream pleioFDR/conjFDR steps.

By default, these scripts verify existing intermediate files and output tables instead of regenerating large files:

- `bash scripts/04_prepare_pleiofdr_mat_inputs.sh`
- `bash scripts/05_prepare_decomposition_common_snps.sh`
- `python3 scripts/08_collect_decomp_raw_loci.py`
- `python3 scripts/09_group_decomp_loci_table.py`

Optional regeneration commands:

- Main pleioFDR MAT files: `RUN_CONVERSION=1 bash scripts/04_prepare_pleiofdr_mat_inputs.sh`
- Decomposition raw locus table: `RUN_REBUILD=1 python3 scripts/08_collect_decomp_raw_loci.py`
- Decomposition grouped locus table: `RUN_REBUILD=1 python3 scripts/09_group_decomp_loci_table.py`

Scripts 08 and 09 do not overwrite existing output tables if zero rows are detected during rebuild.

Expected verified outputs include:

- `$PLEIO/traitfiles/adenomyosis_noMHC.mat`
- `$PLEIO/traitfiles/allEC_noMHC.mat`
- `$PLEIO/traitfiles/adenomyosis_EECpair_noMHC.mat`
- `$PLEIO/traitfiles/EEC_noMHC.mat`
- `$PLEIO/traitfiles_decomp_common/`
- `results/tables/table_s9_decomp_conjfdr_loci_raw.tsv`
- `results/tables/table_s10_decomp_locus_grouped_from_file.tsv`

The paths in `scripts/00_config.sh` reflect the local analysis environment and should be edited by each user before running the workflow.

