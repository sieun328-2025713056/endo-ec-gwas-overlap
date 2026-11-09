# ENDO-only MAF sensitivity: GenomicSEM validation

Validation date: 2026-11-09

## Analysis definition

Only the overall-endometriosis LDSC input was replaced with
the MAF-filtered summary statistics.

The adenomyosis, endometriosis-without-adenomyosis,
and endometrial cancer inputs were unchanged.

LDSC munging used `--maf-min 0.05`.
The resulting ENDO file contained 1,064,328 SNPs with
nonmissing beta values.

## Existing LDSC inputs

- `results/maf05/ldsc_allEC_4trait.rds`
- `results/maf05/ldsc_EEC_4trait.rds`

These paths are relative to the original GenomicSEM analysis directory.

## Verification script

`scripts/15_maf05_rgmodel_verification.R`

This script was newly prepared for reproducibility validation.
It is not claimed to be the original execution script.

It reads existing LDSC RDS files, computes rgmodel results,
performs covariance-aware omnibus and pairwise Wald tests,
and applies Holm correction jointly across six pairwise tests.

No GWAS munging or LDSC rerun is performed.

## Verified results

| Outcome | Omnibus ChiSq(2) | Omnibus P | ENDO-ADENO Holm P |
|---|---:|---:|---:|
| allEC | 7.411255 | 0.02458479 | 0.08779310 |
| EEC | 8.821470 | 0.01214624 | 0.05680075 |

All four values matched the previously reported values
within an absolute tolerance of 1e-5.

Verification status: PASSED.

## Interpretation

Both omnibus tests met the primary two-outcome
Bonferroni threshold of P < 0.025.

No individual pairwise contrast survived Holm correction
across all six comparisons in this sensitivity analysis.

The allEC omnibus result was close to the threshold.
This sensitivity analysis does not resolve the
overall-endometriosis LDSC intercept deflation.

## Software and limitations

Observed software:
- GenomicSEM 0.0.5
- lavaan 0.6-19

Package startup reported 23 warnings whose details
have not yet been independently reviewed.

The verification confirms the downstream calculations
from saved LDSC RDS files, not the correctness of
the original GWAS or LD reference resources.
