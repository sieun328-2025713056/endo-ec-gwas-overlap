#!/usr/bin/env bash
set -euo pipefail

# Edit these paths for your environment before running the workflow.
export PROJECT="${PROJECT:-$HOME/projects/endo_ecancer}"
export RAW="${RAW:-$PROJECT/raw}"
export REF="${REF:-$PROJECT/ref}"
export MUNGED="${MUNGED:-$PROJECT/munged}"
export RESULTS="${RESULTS:-$PROJECT/results}"
export H2="${H2:-$RESULTS/h2}"
export RG="${RG:-$RESULTS/rg}"
export CJ="${CJ:-$RESULTS/conjfdr}"
export COLOC="${COLOC:-$RESULTS/coloc/rs9668810_region}"
export TABLES="${TABLES:-$RESULTS/tables}"

# Tool paths.
export LDSC_DIR="${LDSC_DIR:-$HOME/tools/ldsc}"
export PYCONVERT_DIR="${PYCONVERT_DIR:-$HOME/tools/python_convert}"
export PLEIO="${PLEIO:-/mnt/c/GWAS_project/pleiofdr}"
export MATLAB_BIN="${MATLAB_BIN:-/mnt/c/Program Files/MATLAB/R2026a/bin/matlab.exe}"

# Sample sizes used to construct ECAC LDSC inputs.
export N_ALL_EC=121885      # 12,906 cases + 108,979 controls
export N_EEC=54884          # 8,758 cases + 46,126 controls

mkdir -p "$RAW" "$REF" "$MUNGED" "$H2" "$RG" "$CJ" "$COLOC" "$TABLES"
