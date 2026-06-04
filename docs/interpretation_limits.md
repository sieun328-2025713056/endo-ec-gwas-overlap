# Interpretation limits

This workflow is based on publicly available GWAS summary statistics and should be interpreted accordingly.

## General limitations

- The analysis does not use individual-level genotype or phenotype data.
- Phenotype definition heterogeneity in the source GWAS datasets may influence downstream genetic overlap estimates.
- Differences in sample size and statistical power between endometriosis, adenomyosis, and endometriosis without adenomyosis GWAS datasets may affect observed differences.
- Non-detection of a signal in a given phenotype pair should not be interpreted as evidence of absence.

## LDSC interpretation

- LDSC genetic correlation estimates reflect genome-wide shared genetic architecture.
- Approximate rg difference tests do not account for covariance between partially overlapping GWAS datasets and should be interpreted as exploratory.

## conjFDR interpretation

- conjFDR findings are locus-level statistical prioritizations, not proof of causal variants or causal genes.
- Matched SNP-universe analyses were used to improve comparability across phenotype pairs.

## Colocalization interpretation

- The chr12 SSPN-BHLHE41 follow-up was exploratory.
- P-value-based colocalization can support a shared regional signal but cannot confirm a causal variant or causal gene.
- Gene annotation and eQTL evidence should be treated as prioritization evidence, not functional validation.
