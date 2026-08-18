# Data sources

This repository contains code and documentation for reproducing the GWAS summary-statistics workflow. Raw GWAS summary statistics, LD/reference files, and large intermediate outputs are not redistributed.

## Endometriosis-spectrum GWAS summary statistics

Endometriosis-spectrum GWAS summary statistics were obtained from Koller et al. (2026) and are available through Zenodo (DOI `10.5281/zenodo.18983492`).

The analyses use the following European-ancestry phenotype definitions:

- Overall endometriosis, using the publicly released European combined GWAS summary-statistics file excluding the 23andMe component.
- Adenomyosis.
- Endometriosis without adenomyosis.

For overall endometriosis, case and control counts were not directly reported for the 23andMe-excluded public release. Approximate counts used for manuscript description were derived from the reported full European combined meta-analysis after subtracting the reported 23andMe contribution. The derivation is documented in `docs/sample_size_derivation.md`.

## Endometrial cancer GWAS summary statistics

Endometrial cancer GWAS summary statistics were obtained from O'Mara et al. (2018).

The analyses use the following datasets:

- All-histology endometrial cancer: GWAS Catalog accession `GCST006464`.
- Endometrioid endometrial cancer: GWAS Catalog accession `GCST006465`.

The total sample sizes used for preparation of the LDSC inputs were:

- All-histology endometrial cancer: N = 121,885, including 12,906 cases and 108,979 controls.
- Endometrioid endometrial cancer: N = 54,884, including 8,758 cases and 46,126 controls.

These sample sizes are also defined in `scripts/00_config.sh`.

## Local workflow filenames

The reproducibility workflow expects the following local raw input filenames:

```text
raw/
├── adenomyosis_EUR.txt
├── endometriosis_EUR_wo_23andMe.txt
├── endometriosis_wo.adenomyosis_EUR.txt
├── ECAC_allEC.valid_2.tsv
└── ECAC_EEC.valid_2.tsv
```

These filenames reflect the local analysis environment and are provided only to document the expected workflow inputs. They should not be interpreted as canonical filenames used by the original data providers.

## Data redistribution and use restrictions

This repository does not redistribute raw GWAS summary statistics, restricted-access data, large LD/reference resources, generated `.mat` files, or large local pleioFDR/conjFDR output directories.

Users must obtain source GWAS summary statistics and reference resources from the respective original providers and comply with all applicable data-use terms and restrictions.

Users should cite the corresponding original GWAS publications when using these data.