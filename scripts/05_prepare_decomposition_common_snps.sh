#!/usr/bin/env bash
set -euo pipefail

# 05_prepare_decomposition_common_snps.sh
# Verify decomposition common-SNP inputs and pleioFDR MAT files.
#
# Default behavior:
#   - Verify existing common SNP lists.
#   - Verify existing decomposition common standardized TSV files.
#   - Verify existing decomposition common MAT files.
#   - Print MATLAB whos() output for representative MAT files.
#
# This script is intentionally verify-first because the original decomposition
# common MAT files are large and have already been generated successfully.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_config.sh"

QC_DIR="$CJ/qc_decomp_common"
STD_DIR="$CJ/pleio_mat_std_decomp_common"
MAT_DIR="$PLEIO/traitfiles_decomp_common"

mkdir -p "$CJ/logs"

echo "### Checking decomposition common SNP lists"

QC_FILES=(
  "adeno_allEC.snps"
  "overall_allEC.snps"
  "woadeno_allEC.snps"
  "common_allEC_3phenotype.snps"
  "adeno_EEC.snps"
  "overall_EEC.snps"
  "woadeno_EEC.snps"
  "common_EEC_3phenotype.snps"
)

for f in "${QC_FILES[@]}"; do
  path="$QC_DIR/$f"
  if [[ ! -s "$path" ]]; then
    echo "Missing or empty SNP list: $path" >&2
    exit 1
  fi
  echo "OK: $path  lines=$(wc -l < "$path")"
done

echo
echo "### Checking decomposition common standardized TSV files"

STD_FILES=(
  "adeno_allEC_common.std.pval.tsv"
  "allEC_adeno_common.std.pval.tsv"
  "overall_allEC_common.std.pval.tsv"
  "allEC_overall_common.std.pval.tsv"
  "woadeno_allEC_common.std.pval.tsv"
  "allEC_woadeno_common.std.pval.tsv"
  "adeno_EEC_common.std.pval.tsv"
  "EEC_adeno_common.std.pval.tsv"
  "overall_EEC_common.std.pval.tsv"
  "EEC_overall_common.std.pval.tsv"
  "woadeno_EEC_common.std.pval.tsv"
  "EEC_woadeno_common.std.pval.tsv"
)

for f in "${STD_FILES[@]}"; do
  path="$STD_DIR/$f"
  if [[ ! -s "$path" ]]; then
    echo "Missing or empty standardized TSV: $path" >&2
    exit 1
  fi
  header="$(head -n 1 "$path")"
  echo "OK: $path  lines=$(wc -l < "$path")"
  echo "    header: $header"
done

echo
echo "### Checking decomposition common MAT files"

MAT_FILES=(
  "adeno_allEC_common.mat"
  "allEC_adeno_common.mat"
  "overall_allEC_common.mat"
  "allEC_overall_common.mat"
  "woadeno_allEC_common.mat"
  "allEC_woadeno_common.mat"
  "adeno_EEC_common.mat"
  "EEC_adeno_common.mat"
  "overall_EEC_common.mat"
  "EEC_overall_common.mat"
  "woadeno_EEC_common.mat"
  "EEC_woadeno_common.mat"
)

for f in "${MAT_FILES[@]}"; do
  path="$MAT_DIR/$f"
  if [[ ! -s "$path" ]]; then
    echo "Missing or empty MAT file: $path" >&2
    exit 1
  fi
  echo "OK: $path  size=$(du -h "$path" | cut -f1)"
done

echo
echo "### Verifying representative MAT file structure with MATLAB"

# Convert the pleioFDR path for Windows MATLAB when running from WSL.
# Otherwise retain the native path.
if command -v wslpath >/dev/null 2>&1 && [[ "$MATLAB_BIN" == *.exe ]]; then
  PLEIO_MATLAB="$(wslpath -w "$PLEIO")"
else
  PLEIO_MATLAB="$PLEIO"
fi

echo "Using MATLAB executable: $MATLAB_BIN"
echo "Using pleioFDR path in MATLAB: $PLEIO_MATLAB"

"$MATLAB_BIN" -batch "cd('$PLEIO_MATLAB'); whos('-file','traitfiles_decomp_common/adeno_allEC_common.mat'); whos('-file','traitfiles_decomp_common/allEC_adeno_common.mat'); whos('-file','traitfiles_decomp_common/adeno_EEC_common.mat'); whos('-file','traitfiles_decomp_common/EEC_adeno_common.mat')"

echo
echo "### 05_prepare_decomposition_common_snps.sh completed"
