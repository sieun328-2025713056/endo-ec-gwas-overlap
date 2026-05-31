#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

mkdir -p "$CJ/input" "$CJ/qc" "$CJ/harmonized"

# 1. HapMap3-filtered input files.
make_koller_hm3() {
  local infile="$1"
  local outfile="$2"
  awk 'NR==FNR{keep[$1]=1; next}
       FNR==1{print "SNP\tA1\tA2\tZ\tP\tN"; next}
       ($1 in keep){print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$5}' \
       "$REF/w_hm3.snplist" "$infile" > "$outfile"
}

make_ec_hm3() {
  local infile="$1"
  local outfile="$2"
  awk 'NR==FNR{keep[$1]=1; next}
       FNR==1{print "SNP\tA1\tA2\tZ\tP\tN"; next}
       ($1 in keep){print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4}' \
       "$REF/w_hm3.snplist" "$infile" > "$outfile"
}

make_koller_hm3 "$RAW/adenomyosis_EUR.txt" "$CJ/input/adenomyosis_EUR.hm3.tsv"
make_koller_hm3 "$RAW/endometriosis_EUR_wo_23andMe.txt" "$CJ/input/overall_endometriosis.hm3.tsv"
make_koller_hm3 "$RAW/endometriosis_wo.adenomyosis_EUR.txt" "$CJ/input/endo_wo_adeno.hm3.tsv"
make_ec_hm3 "$RAW/ECAC_allEC_for_ldsc.tsv" "$CJ/input/allEC.hm3.tsv"
make_ec_hm3 "$RAW/ECAC_EEC_for_ldsc.tsv" "$CJ/input/EEC.hm3.tsv"

# 2. MHC SNP list from ECAC coordinates. Used for noMHC conjFDR inputs.
awk 'NR>1 && $3==6 && $2>=26000000 && $2<=34000000 {print $1}' "$RAW/ECAC_allEC.valid_2.tsv" > "$CJ/qc/MHC_allEC.tmp"
awk 'NR>1 && $3==6 && $2>=26000000 && $2<=34000000 {print $1}' "$RAW/ECAC_EEC.valid_2.tsv" > "$CJ/qc/MHC_EEC.tmp"
cat "$CJ/qc/MHC_allEC.tmp" "$CJ/qc/MHC_EEC.tmp" | sort -u > "$CJ/qc/MHC_26_34Mb.snplist"

# 3. Pairwise allele harmonization. Outcome Z is flipped when allele orientation is reversed.
harmonize_pair() {
  local exposure_tsv="$1"
  local outcome_tsv="$2"
  local pair_label="$3"
  local exposure_out="$4"
  local outcome_out="$5"

  local exp_sorted="$CJ/qc/${pair_label}.exposure.sorted.tsv"
  local out_sorted="$CJ/qc/${pair_label}.outcome.sorted.tsv"
  local joined="$CJ/qc/${pair_label}.joined.tsv"
  local harmonized="$CJ/harmonized/${pair_label}.harmonized.tsv"

  awk 'NR==1{next} {print $1"\t"toupper($2)"\t"toupper($3)"\t"$4"\t"$5"\t"$6}' "$exposure_tsv" | sort -k1,1 > "$exp_sorted"
  awk 'NR==1{next} {print $1"\t"toupper($2)"\t"toupper($3)"\t"$4"\t"$5"\t"$6}' "$outcome_tsv" | sort -k1,1 > "$out_sorted"
  join -t $'\t' -1 1 -2 1 "$exp_sorted" "$out_sorted" > "$joined"

  awk -F'\t' 'BEGIN{OFS="\t"; print "SNP","A1","A2","Z_EXPOSURE","P_EXPOSURE","N_EXPOSURE","Z_OUTCOME","P_OUTCOME","N_OUTCOME","ALLELE_STATUS"}
    {if($2==$7 && $3==$8) print $1,$2,$3,$4,$5,$6,$9,$10,$11,"same";
     else if($2==$8 && $3==$7) print $1,$2,$3,$4,$5,$6,-$9,$10,$11,"flip"}' "$joined" > "$harmonized"

  awk 'BEGIN{OFS="\t"} NR==1{print "SNP","A1","A2","Z","P","N"; next} {print $1,$2,$3,$4,$5,$6}' "$harmonized" > "$exposure_out"
  awk 'BEGIN{OFS="\t"} NR==1{print "SNP","A1","A2","Z","P","N"; next} {print $1,$2,$3,$7,$8,$9}' "$harmonized" > "$outcome_out"

  # MHC exclusion. FILENAME==ARGV[] is robust even if the MHC file is empty.
  awk 'FILENAME==ARGV[1]{mhc[$1]=1; next} FILENAME==ARGV[2]{if(FNR==1 || !($1 in mhc)) print}' \
    "$CJ/qc/MHC_26_34Mb.snplist" "$exposure_out" > "${exposure_out%.tsv}.noMHC.tsv"
  awk 'FILENAME==ARGV[1]{mhc[$1]=1; next} FILENAME==ARGV[2]{if(FNR==1 || !($1 in mhc)) print}' \
    "$CJ/qc/MHC_26_34Mb.snplist" "$outcome_out" > "${outcome_out%.tsv}.noMHC.tsv"

  echo "### $pair_label"
  wc -l "$joined" "$harmonized" "$exposure_out" "$outcome_out" "${exposure_out%.tsv}.noMHC.tsv" "${outcome_out%.tsv}.noMHC.tsv"
  awk -F'\t' 'BEGIN{same=0; flip=0; other=0} {if($2==$7 && $3==$8) same++; else if($2==$8 && $3==$7) flip++; else other++} END{print "same_allele",same; print "flipped_allele",flip; print "other_mismatch",other; print "total",same+flip+other}' "$joined"

  # Basic p-value sanity check for the two final pair inputs.
  awk 'NR>1 && ($5<=0 || $5>1 || $5=="NA"){bad++} END{print "bad_P_exposure", bad+0}' "$exposure_out"
  awk 'NR>1 && ($5<=0 || $5>1 || $5=="NA"){bad++} END{print "bad_P_outcome", bad+0}' "$outcome_out"
}

# Main ADENO pairs.
harmonize_pair "$CJ/input/adenomyosis_EUR.hm3.tsv" "$CJ/input/allEC.hm3.tsv" \
  "adenomyosis_allEC" "$CJ/harmonized/adenomyosis_for_conjfdr.tsv" "$CJ/harmonized/allEC_for_conjfdr.tsv"

harmonize_pair "$CJ/input/adenomyosis_EUR.hm3.tsv" "$CJ/input/EEC.hm3.tsv" \
  "adenomyosis_EEC" "$CJ/harmonized/adenomyosis_for_conjfdr_EEC.tsv" "$CJ/harmonized/EEC_for_conjfdr.tsv"

# Additional 4 pairs for later 6-pair decomposition.
harmonize_pair "$CJ/input/overall_endometriosis.hm3.tsv" "$CJ/input/allEC.hm3.tsv" \
  "overall_endometriosis_allEC" "$CJ/harmonized/overall_endometriosis_allEC.tsv" "$CJ/harmonized/allEC_overall_endometriosis.tsv"

harmonize_pair "$CJ/input/overall_endometriosis.hm3.tsv" "$CJ/input/EEC.hm3.tsv" \
  "overall_endometriosis_EEC" "$CJ/harmonized/overall_endometriosis_EEC.tsv" "$CJ/harmonized/EEC_overall_endometriosis.tsv"

harmonize_pair "$CJ/input/endo_wo_adeno.hm3.tsv" "$CJ/input/allEC.hm3.tsv" \
  "endo_wo_adeno_allEC" "$CJ/harmonized/endo_wo_adeno_allEC.tsv" "$CJ/harmonized/allEC_endo_wo_adeno.tsv"

harmonize_pair "$CJ/input/endo_wo_adeno.hm3.tsv" "$CJ/input/EEC.hm3.tsv" \
  "endo_wo_adeno_EEC" "$CJ/harmonized/endo_wo_adeno_EEC.tsv" "$CJ/harmonized/EEC_endo_wo_adeno.tsv"

# Lambda GC sanity check. This is not an allele-direction validation because lambda uses Z^2.
for f in "$CJ"/harmonized/*.noMHC.tsv; do
  [[ -s "$f" ]] || continue
  label=$(basename "$f")
  awk 'NR>1{print $4*$4}' "$f" | sort -n | awk -v label="$label" '
    {a[NR]=$1}
    END{
      if(NR==0){printf "%s\tlambda_GC=NA\tN=0\n", label; exit}
      if(NR%2){m=a[(NR+1)/2]} else {m=(a[NR/2]+a[NR/2+1])/2}
      printf "%s\tlambda_GC=%.4f\tN=%d\n", label, m/0.4549364, NR
    }'
done
