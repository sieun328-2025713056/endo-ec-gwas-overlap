# Reproducibility Notes

## Raw data not included

This package does not include raw GWAS summary statistics or large reference files. Users must obtain the source data from the original providers and comply with each dataset's data-use terms.

## Main raw files expected

- `adenomyosis_EUR.txt`
- `endometriosis_EUR_wo_23andMe.txt`
- `endometriosis_wo.adenomyosis_EUR.txt`
- `ECAC_allEC.valid_2.tsv`
- `ECAC_EEC.valid_2.tsv`

## Key analysis constants

- All EC sample size used in LDSC/coloc: `N = 121885`
- EEC sample size used in LDSC/coloc: `N = 54884`
- Adenomyosis case-control proportion used in coloc: `8753 / (8753 + 415718)`
- All EC case-control proportion used in coloc: `12906 / (12906 + 108979)`
- EEC case-control proportion used in coloc: `8758 / (8758 + 46126)`
- MHC exclusion window: chr6:26–34 Mb
- chr12 coloc region: rs9668810 ±500 kb; coloc input preparation is expected upstream of `06_run_coloc_chr12_rs9668810.R`

## v2 cleaning changes

- Removed `raw_extraction/` and internal conversation-derived text from the GitHub-facing package.
- Fixed the lambda GC `awk` block in `03_harmonize_conjfdr_inputs.sh` so empty files do not trigger `next` in an `END` block.
- Made `04_prepare_pleiofdr_mat_inputs.sh` more robust to different `python_convert` layouts by auto-detecting `sumstats.py`.
- Added explicit required-file checks to `05_prepare_decomposition_common_snps.sh`.
- Made `06_run_coloc_chr12_rs9668810.R` path-configurable through `$COLOC`, `$PROJECT`, or command-line arguments.
- Added `08_collect_decomp_raw_loci.py` to create `table_s9_decomp_conjfdr_loci_raw.tsv` before grouping loci.

## Interpretive caution

The scripts reproduce statistical analyses. They do not by themselves justify causal interpretation. In particular:

- `adenomyosis-specific` should not be used unless decomposition results support it.
- Non-detection in a phenotype pair should be interpreted as non-detection under current GWAS power and phenotype definition, not as proof of no association.
- chr12 colocalization is exploratory because adenomyosis beta/SE were unavailable and p-value-based coloc was used.
- rs9668810 should be described as an index SNP/regional signal unless fine-mapping identifies a causal variant.
