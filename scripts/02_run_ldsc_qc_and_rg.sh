#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

cd "$LDSC_DIR"

# Munge Koller endometriosis-spectrum GWAS files.
./munge_sumstats.py --sumstats "$RAW/endometriosis_EUR_wo_23andMe.txt" \
  --snp SNP --a1 Allele1 --a2 Allele2 --p P.value --N-col N \
  --signed-sumstats Z,0 --merge-alleles "$REF/w_hm3.snplist" \
  --out "$MUNGED/endometriosis_EUR"

./munge_sumstats.py --sumstats "$RAW/adenomyosis_EUR.txt" \
  --snp SNP --a1 Allele1 --a2 Allele2 --p P.value --N-col N \
  --signed-sumstats Z,0 --merge-alleles "$REF/w_hm3.snplist" \
  --out "$MUNGED/adenomyosis_EUR"

./munge_sumstats.py --sumstats "$RAW/endometriosis_wo.adenomyosis_EUR.txt" \
  --snp SNP --a1 Allele1 --a2 Allele2 --p P.value --N-col N \
  --signed-sumstats Z,0 --merge-alleles "$REF/w_hm3.snplist" \
  --out "$MUNGED/endometriosis_wo_adeno_EUR"

# Munge EC outcomes.
./munge_sumstats.py --sumstats "$RAW/ECAC_allEC_for_ldsc.tsv" \
  --snp SNP --a1 A1 --a2 A2 --p P --N-col N \
  --signed-sumstats Z,0 --merge-alleles "$REF/w_hm3.snplist" \
  --out "$MUNGED/endometrial_cancer_allEC"

./munge_sumstats.py --sumstats "$RAW/ECAC_EEC_for_ldsc.tsv" \
  --snp SNP --a1 A1 --a2 A2 --p P --N-col N \
  --signed-sumstats Z,0 --merge-alleles "$REF/w_hm3.snplist" \
  --out "$MUNGED/endometrial_cancer_EEC"

# SNP-heritability QC.
for trait in \
  endometriosis_EUR \
  adenomyosis_EUR \
  endometriosis_wo_adeno_EUR \
  endometrial_cancer_allEC \
  endometrial_cancer_EEC
 do
  ./ldsc.py --h2 "$MUNGED/${trait}.sumstats.gz" \
    --ref-ld-chr "$REF/eur_w_ld_chr/" \
    --w-ld-chr "$REF/eur_w_ld_chr/" \
    --out "$H2/$trait"
done

# Genetic correlations: 3 exposure phenotypes x 2 EC outcomes.
declare -A exposures=(
  [endometriosis_EUR]="endometriosis_EUR"
  [adenomyosis_EUR]="adenomyosis_EUR"
  [endometriosis_wo_adeno_EUR]="endometriosis_wo_adeno_EUR"
)

declare -A outcomes=(
  [allEC]="endometrial_cancer_allEC"
  [EEC]="endometrial_cancer_EEC"
)

for exp_label in "${!exposures[@]}"; do
  for out_label in "${!outcomes[@]}"; do
    ./ldsc.py --rg "$MUNGED/${exposures[$exp_label]}.sumstats.gz,$MUNGED/${outcomes[$out_label]}.sumstats.gz" \
      --ref-ld-chr "$REF/eur_w_ld_chr/" \
      --w-ld-chr "$REF/eur_w_ld_chr/" \
      --out "$RG/rg_${exp_label}_${out_label}"
  done
done

# Log summaries.
grep -E "Total Observed scale h2|Total Observed scale h:|Lambda GC|Mean Chi|Intercept|Ratio" "$H2"/*.log || true
grep -A 20 "Summary of Genetic Correlation Results" "$RG"/*.log || true
