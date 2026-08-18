#!/usr/bin/env bash
set -euo pipefail

# 04_prepare_pleiofdr_mat_inputs.sh
# Verify or optionally generate pleioFDR-compatible .mat trait files.
#
# Default behavior:
#   - Verify that previously generated .mat files exist.
#   - Verify that source standardized .std.pval.tsv files exist.
#   - Print MATLAB whos() output for the .mat files.
#
# Optional conversion:
#   RUN_CONVERSION=1 bash scripts/04_prepare_pleiofdr_mat_inputs.sh
#
# Notes:
#   The successful historical conversion used:
#   $CJ/pleio_mat_std/*.std.pval.tsv
#   not the raw harmonized/*.noMHC.tsv files.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_config.sh"

mkdir -p "$CJ/logs"
mkdir -p "$PLEIO/traitfiles"

SUMSTATS_PY=""
for candidate in \
  "$PYCONVERT_DIR/sumstats.py" \
  "$PYCONVERT_DIR/src/converter/sumstats.py" \
  "$PYCONVERT_DIR/python_convert/sumstats.py"
do
  if [[ -s "$candidate" ]]; then
    SUMSTATS_PY="$candidate"
    break
  fi
done

MAT_NAMES=(
  "adenomyosis_noMHC.mat"
  "allEC_noMHC.mat"
  "adenomyosis_EECpair_noMHC.mat"
  "EEC_noMHC.mat"
)

STD_NAMES=(
  "adenomyosis_noMHC.std.pval.tsv"
  "allEC_noMHC.std.pval.tsv"
  "adenomyosis_EECpair_noMHC.std.pval.tsv"
  "EEC_noMHC.std.pval.tsv"
)

echo "### Checking standardized pleioFDR input TSV files"
for f in "${STD_NAMES[@]}"; do
  path="$CJ/pleio_mat_std/$f"
  if [[ ! -s "$path" ]]; then
    echo "Missing standardized input: $path" >&2
    exit 1
  fi
  echo "OK: $path"
done

echo
echo "### Checking existing pleioFDR MAT trait files"
missing_mat=0

for f in "${MAT_NAMES[@]}"; do
  path="$PLEIO/traitfiles/$f"
  if [[ ! -s "$path" ]]; then
    echo "Missing MAT file: $path" >&2
    missing_mat=1
  else
    echo "OK: $path"
  fi
done

if [[ "${RUN_CONVERSION:-0}" != "1" ]]; then
  if [[ "$missing_mat" -ne 0 ]]; then
    echo
    echo "One or more MAT files are missing."
    echo "To generate them, run:"
    echo "RUN_CONVERSION=1 bash scripts/04_prepare_pleiofdr_mat_inputs.sh"
    exit 1
  fi

  echo
  echo "RUN_CONVERSION is not set to 1."
  echo "Existing MAT files found. Skipping heavy conversion step."
else
  echo
  echo "### RUN_CONVERSION=1: generating MAT files"

  if [[ ! -x "$PYCONVERT_PYTHON" ]]; then
    echo "Missing or non-executable PYCONVERT_PYTHON: $PYCONVERT_PYTHON" >&2
    echo "Set PYCONVERT_PYTHON to the Python executable used for python_convert." >&2
    exit 1
  fi

  if [[ -z "$SUMSTATS_PY" ]]; then
    echo "Could not find sumstats.py under $PYCONVERT_DIR" >&2
    exit 1
  fi

  if [[ ! -s "$PLEIO/9545380.ref" ]]; then
    echo "Missing pleioFDR reference SNP-order file: $PLEIO/9545380.ref" >&2
    exit 1
  fi

  echo "Using PYCONVERT_PYTHON: $PYCONVERT_PYTHON"
  "$PYCONVERT_PYTHON" --version
  echo "Using SUMSTATS_PY: $SUMSTATS_PY"

  for i in "${!MAT_NAMES[@]}"; do
    input="$CJ/pleio_mat_std/${STD_NAMES[$i]}"
    output="$PLEIO/traitfiles/${MAT_NAMES[$i]}"
    screen_log="$PLEIO/traitfiles/${MAT_NAMES[$i]%.mat}_convert_screen.log"

    echo
    echo "Converting:"
    echo "  input:  $input"
    echo "  output: $output"

    "$PYCONVERT_PYTHON" "$SUMSTATS_PY" mat \
      --force \
      --sumstats "$input" \
      --out "$output" \
      --ref "$PLEIO/9545380.ref" \
      --chunksize 1000000 \
      2>&1 | tee "$screen_log"
  done
fi

echo
echo "### Verifying MAT file structure with MATLAB"

# Convert the pleioFDR path for Windows MATLAB when running from WSL.
# Otherwise retain the native path.
if command -v wslpath >/dev/null 2>&1 && [[ "$MATLAB_BIN" == *.exe ]]; then
  PLEIO_MATLAB="$(wslpath -w "$PLEIO")"
else
  PLEIO_MATLAB="$PLEIO"
fi

echo "Using MATLAB executable: $MATLAB_BIN"
echo "Using pleioFDR path in MATLAB: $PLEIO_MATLAB"

"$MATLAB_BIN" -batch "cd('$PLEIO_MATLAB'); whos('-file','traitfiles/adenomyosis_noMHC.mat'); whos('-file','traitfiles/allEC_noMHC.mat'); whos('-file','traitfiles/adenomyosis_EECpair_noMHC.mat'); whos('-file','traitfiles/EEC_noMHC.mat')"

echo
echo "### 04_prepare_pleiofdr_mat_inputs.sh completed"