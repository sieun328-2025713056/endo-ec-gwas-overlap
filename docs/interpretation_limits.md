# Interpretation limits

This workflow is based on GWAS summary statistics and should be interpreted within the limitations of summary-statistics-based genetic epidemiology.

## General limitations

- The analysis does not use individual-level genotype or phenotype data.
- Differences in phenotype definitions, sample sizes, and statistical power across the source GWAS datasets may influence the magnitude and detectability of genetic-overlap signals.
- Observed variation across overall endometriosis, adenomyosis, and endometriosis without adenomyosis should not be interpreted as establishing statistically distinct or phenotype-exclusive genetic overlap unless formally supported.
- Non-detection of a signal in a phenotype pair indicates lack of detection under the current GWAS power, phenotype definition, matched SNP universe, and statistical threshold; it should not be interpreted as evidence of biological absence.

## LDSC interpretation

- LDSC genetic-correlation estimates reflect genome-wide shared genetic architecture and do not establish causality.
- The six primary cross-trait genetic-correlation tests were evaluated using Bonferroni correction, with Benjamini–Hochberg FDR q-values reported as supportive information.
- Approximate between-phenotype genetic-correlation difference tests use `sqrt(SE1^2 + SE2^2)` and do not account for covariance arising from the shared endometrial cancer GWAS or non-independent endometriosis-spectrum estimates. These comparisons should therefore be interpreted as exploratory.
- Cross-trait LDSC intercepts may reflect residual covariance, including potential sample overlap or correlated error, and should not be interpreted as direct evidence of participant overlap or bias of a specific direction.

## conjFDR interpretation

- Conjunctional FDR findings are locus-level statistical prioritizations and do not establish causal variants, causal genes, or biological mechanisms.
- The six phenotype-pair comparisons were evaluated within a common matched SNP universe to improve comparability across analyses.
- MHC exclusion at chr6:26–34 Mb was applied to the matched SNP universe used for conjunctional FDR analyses and should not be interpreted as having been applied to the LDSC analyses.
- Recurrent-locus detection patterns are descriptive and power-sensitive. Labels such as pan-spectrum or adenomyosis/overall-endometriosis-aligned describe observed detection patterns rather than phenotype-specific biological effects.

## Colocalization interpretation

- The chr12 rs9668810-indexed SSPN–BHLHE41 regional follow-up was exploratory.
- Because effect estimates and standard errors were unavailable for the adenomyosis GWAS, P-value-based approximate Bayes factors were used for the regional colocalization analyses.
- Strong posterior support for a shared regional signal does not establish that rs9668810, the highest-posterior SNP, or any nearby gene is causal.
- The coloc model assumes at most one causal variant per trait within the evaluated region, and results should be interpreted within that assumption.

## Functional-annotation interpretation

- GTEx V10 eQTL follow-up was exploratory and was used for regional functional annotation rather than functional validation.
- Significant eQTL associations in individual tissues do not establish that the corresponding gene mediates the genetic overlap between the studied phenotypes.
- Lack of a significant bulk-uterus eQTL association should not be interpreted as evidence that the locus has no regulatory role in uterine, endometrial, or other disease-relevant cell types or tissues.
- Gene names reported for recurrent loci should be interpreted as regional annotations or prioritization labels rather than confirmed causal genes.