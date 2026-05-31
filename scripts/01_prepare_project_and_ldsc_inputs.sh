#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/00_config.sh"

# Required raw files. These are not redistributed with this code package.
required=(
  "adenomyosis_EUR.txt"
  "endometriosis_EUR_wo_23andMe.txt"
  "endometriosis_wo.adenomyosis_EUR.txt"
  "ECAC_allEC.valid_2.tsv"
  "ECAC_EEC.valid_2.tsv"
)

for f in "${required[@]}"; do
  if [[ ! -s "$RAW/$f" ]]; then
    echo "Missing required raw file: $RAW/$f" >&2
    exit 1
  fi
done

rm -f "$RAW/desktop.ini"

# Convert ECAC beta/SE files to LDSC-compatible SNP/A1/A2/N/Z/P format.
awk -v N="$N_ALL_EC" 'BEGIN{OFS="\t"}
  NR==1{print "SNP","A1","A2","N","Z","P"; next}
  NR>1 && $8!="NA" && $9!="NA" && $9!=0 && $10!="NA" {print $1,$6,$7,N,$8/$9,$10}
' "$RAW/ECAC_allEC.valid_2.tsv" > "$RAW/ECAC_allEC_for_ldsc.tsv"

awk -v N="$N_EEC" 'BEGIN{OFS="\t"}
  NR==1{print "SNP","A1","A2","N","Z","P"; next}
  NR>1 && $8!="NA" && $9!="NA" && $9!=0 && $10!="NA" {print $1,$6,$7,N,$8/$9,$10}
' "$RAW/ECAC_EEC.valid_2.tsv" > "$RAW/ECAC_EEC_for_ldsc.tsv"

# Basic checks.
head -n 5 "$RAW/ECAC_allEC_for_ldsc.tsv"
head -n 5 "$RAW/ECAC_EEC_for_ldsc.tsv"
wc -l "$RAW/ECAC_allEC.valid_2.tsv" "$RAW/ECAC_allEC_for_ldsc.tsv"
wc -l "$RAW/ECAC_EEC.valid_2.tsv" "$RAW/ECAC_EEC_for_ldsc.tsv"
