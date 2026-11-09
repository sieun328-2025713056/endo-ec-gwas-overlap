# Interpretation limits

This workflow uses GWAS summary statistics and should be interpreted within the limitations of summary-statistics-based genetic epidemiology.

## General limitations

- The analysis does not use individual-level genotype or phenotype data.
- Differences in phenotype definitions, sample sizes, and statistical power across source GWAS datasets may influence genetic-overlap estimates and signal detection.
- Differences in genetic-correlation point estimates do not automatically establish statistically significant heterogeneity.
- Non-detection of a genetic signal does not establish biological absence or phenotype specificity.

## LDSC and GenomicSEM interpretation

- Six bivariate LDSC genetic correlations provide descriptive context rather than the primary hypothesis tests.
- Primary inference used two covariance-aware omnibus Wald tests, one for all-histology endometrial cancer and one for endometrioid endometrial cancer.
- The Bonferroni threshold across the two primary omnibus tests was P < 0.025.
- GenomicSEM `rgmodel()` estimates (`R`) and their sampling covariance matrix (`V_R`) were used for omnibus and pairwise contrasts.
- Holm correction was applied jointly across six post-hoc pairwise comparisons.
- Earlier covariance-ignoring approximate difference tests were exploratory and are not the final inferential basis.
- Both primary omnibus tests met the corrected significance threshold, but individual pairwise findings showed limited robustness.
- The ENDO-only `--maf-min 0.05` sensitivity analysis retained omnibus P values below 0.025. No individual pairwise contrast survived Holm correction in this sensitivity analysis.
- The overall-endometriosis univariate LDSC intercept was 0.7672 (SE 0.0098), indicating substantial deflation. Its cause remains unresolved.
- Genetic correlations and differences between genetic correlations do not establish causality or phenotype-exclusive genetic architecture.
- Cross-trait LDSC intercepts do not independently establish participant overlap or a specific source of bias.

## conjFDR interpretation

- Conjunctional FDR findings are locus-level statistical prioritizations and do not establish causal variants, genes, or biological mechanisms.
- Six phenotype-pair comparisons used a common matched universe of 1,057,700 HapMap3 SNPs.
- The chr6:26-34 Mb MHC region was excluded for conjFDR processing, but this exclusion was not applied to LDSC analyses.
- Locus grouping yielded 45 genomic regions from 92 lead-SNP records.
- Of these 45 regions, 22 were detected in more than one phenotype-pair comparison.
- Multiple detections arise from dependent comparisons and should not be interpreted as independent replication.
- Recurrent-locus detection patterns are descriptive and sensitive to statistical power, phenotype definitions, and detection thresholds.
- Labels such as pan-spectrum or adenomyosis/overall-endometriosis-aligned describe observed detection patterns, not established biological specificity.
- Non-detection does not establish the absence of a genetic effect.

## Colocalization interpretation

- The chr12 SSPN-BHLHE41 regional follow-up was exploratory and selected following the secondary conjFDR detection-pattern analysis.
- rs9668810 was used to index the region; it was not the highest-posterior SNP.
- Final colocalization used adenomyosis P-values and minor allele frequencies together with aligned endometrial cancer effect estimates and squared standard errors.
- Adenomyosis N was set to the median eligible variant-level `AD_N` value of 33,460.65.
- The original derivation of adenomyosis effective sample size remains unresolved.
- Adenomyosis case fraction `s = 0.5` was an explicit modeling assumption, not the observed case fraction.
- Colocalization priors were `p1 = 1e-4`, `p2 = 1e-4`, and `p12 = 1e-5`.
- Under these inputs and assumptions, PP.H4 was 0.9931782 for allEC and 0.9880625 for EEC.
- The analyses included 3,159 and 3,158 eligible SNPs, respectively, and therefore used slightly different SNP sets.
- rs7959150 had the highest conditional `SNP.PP.H4` in both analyses (0.3165 for allEC and 0.2761 for EEC).
- High PP.H4 supports the shared-variant hypothesis conditional on the model assumptions and priors; it does not establish an individual causal variant.
- The coloc.abf model assumes at most one causal variant per trait within the evaluated region.
- The results do not establish a causal gene, tissue, or biological mechanism.

## Functional-annotation interpretation

- GTEx V10 eQTL follow-up was exploratory and used for regional functional annotation rather than functional validation.
- Significant eQTL associations in individual tissues do not establish that the corresponding gene mediates genetic overlap between the studied phenotypes.
- Lack of a significant bulk-uterus eQTL association should not be interpreted as evidence that the locus has no regulatory role in uterine, endometrial, or other disease-relevant cell types or tissues.
- Gene names reported for recurrent loci are regional annotations or prioritization labels rather than confirmed causal genes.
