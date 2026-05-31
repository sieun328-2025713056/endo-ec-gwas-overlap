#!/usr/bin/env python3
"""
08_collect_decomp_raw_loci.py

Verify or optionally rebuild table_s9_decomp_conjfdr_loci_raw.tsv.

Default behavior:
    - Verify that table_s9 already exists and has data rows.
    - Do not touch existing table_s9.

Optional rebuild:
    RUN_REBUILD=1 python3 scripts/08_collect_decomp_raw_loci.py

Safety rule:
    - If zero rows are collected during rebuild, the existing output is not overwritten.
"""

import csv
import glob
import os
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(SCRIPT_DIR)

RESULTS = os.path.join(PROJECT, "results")
TABLES = os.path.join(RESULTS, "tables")
os.makedirs(TABLES, exist_ok=True)

PLEIO = os.environ.get("PLEIO", "/mnt/c/GWAS_project/pleiofdr")

OUTFILE = os.path.join(TABLES, "table_s9_decomp_conjfdr_loci_raw.tsv")
REF_PATH = os.path.join(PLEIO, "9545380.ref")

REQUIRED_COLS = [
    "Analysis",
    "Exposure",
    "Outcome",
    "Locusnum",
    "CHR",
    "BP",
    "Lead_SNP",
    "Ref_A1",
    "Ref_A2",
    "Z_exposure",
    "Z_outcome",
    "P_exposure",
    "P_outcome",
    "conjFDR",
    "Direction",
    "Source_file",
]


RESULT_DIRS = [
    os.path.join(PLEIO, "results_decomp_adeno_allEC_n20"),
    os.path.join(PLEIO, "results_decomp_adeno_EEC_n20"),
    os.path.join(PLEIO, "results_decomp_overall_allEC_n20"),
    os.path.join(PLEIO, "results_decomp_overall_EEC_n20"),
    os.path.join(PLEIO, "results_decomp_woadeno_allEC_n20"),
    os.path.join(PLEIO, "results_decomp_woadeno_EEC_n20"),
]


def fail(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def verify_existing_table(path: str) -> int:
    if not os.path.exists(path):
        fail(f"Missing existing table: {path}")

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"Could not read header from: {path}")

        missing = [c for c in REQUIRED_COLS if c not in reader.fieldnames]
        if missing:
            fail(f"Missing required columns in {path}: {missing}")

        n = sum(1 for _ in reader)

    if n <= 0:
        fail(f"Existing table has zero data rows: {path}")

    print(f"Verified existing table_s9: {path}")
    print(f"Data rows: {n}")
    return n


def infer_exposure_outcome(label: str):
    if label.endswith("_ALL_EC"):
        return label[:-7], "ALL_EC"
    if label.endswith("_EEC"):
        return label[:-4], "EEC"
    return "NA", "NA"


def load_ref_by_pos(ref_path: str):
    if not os.path.exists(ref_path):
        fail(f"Missing reference file: {ref_path}")

    ref_by_pos = {}

    with open(ref_path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split()
            if len(p) >= 6:
                chrom, snp, gp, bp, a1, a2 = p[:6]
                try:
                    ref_by_pos[(str(chrom), int(bp))] = (
                        snp,
                        a1.upper(),
                        a2.upper(),
                    )
                except ValueError:
                    continue

    return ref_by_pos


def collect_locus_files():
    files = []
    for d in RESULT_DIRS:
        pattern = os.path.join(d, "*_zscore_conjfdr_0.05_loci.csv")
        files.extend(glob.glob(pattern))

    files = sorted(files)
    return files


def rebuild_table():
    files = collect_locus_files()

    print(f"Found zscore conjFDR locus CSV files: {len(files)}")
    for f in files:
        print(f"  {f}")

    if not files:
        fail("No locus CSV files found. Existing table_s9 will not be overwritten.")

    print("Loading reference SNP map from 9545380.ref...")
    ref_by_pos = load_ref_by_pos(REF_PATH)

    rows = []

    for f in files:
        analysis = os.path.basename(f).replace("_zscore_conjfdr_0.05_loci.csv", "")
        exposure, outcome = infer_exposure_outcome(analysis)

        with open(f, newline="") as fh:
            reader = csv.DictReader(fh)

            for r in reader:
                chrom = str(r["chrnum"])
                pos = int(float(r["chrpos"]))

                ref_snp, ref_a1, ref_a2 = ref_by_pos.get(
                    (chrom, pos),
                    ("NA", "NA", "NA"),
                )

                zcols = [c for c in r.keys() if c.startswith("zscore_")]
                pcols = [c for c in r.keys() if c.startswith("pval_")]
                conj_cols = [c for c in r.keys() if c.startswith("conjfdr_")]

                if len(zcols) < 2 or len(pcols) < 2 or len(conj_cols) < 1:
                    fail(f"Unexpected columns in {f}: {reader.fieldnames}")

                z_exp = float(r[zcols[0]])
                z_out = float(r[zcols[1]])

                if z_exp * z_out > 0:
                    direction = "concordant"
                elif z_exp * z_out < 0:
                    direction = "discordant"
                else:
                    direction = "zero_or_missing"

                rows.append({
                    "Analysis": analysis,
                    "Exposure": exposure,
                    "Outcome": outcome,
                    "Locusnum": r["locusnum"],
                    "CHR": chrom,
                    "BP": str(pos),
                    "Lead_SNP": ref_snp,
                    "Ref_A1": ref_a1,
                    "Ref_A2": ref_a2,
                    "Z_exposure": "%.6g" % z_exp,
                    "Z_outcome": "%.6g" % z_out,
                    "P_exposure": r[pcols[0]],
                    "P_outcome": r[pcols[1]],
                    "conjFDR": r[conj_cols[0]],
                    "Direction": direction,
                    "Source_file": f,
                })

    if len(rows) == 0:
        fail("Collected zero rows. Existing table_s9 will not be overwritten.")

    fd, tmp_path = tempfile.mkstemp(
        prefix="table_s9_decomp_conjfdr_loci_raw.",
        suffix=".tmp",
        dir=TABLES,
        text=True,
    )
    os.close(fd)

    with open(tmp_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=REQUIRED_COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    os.replace(tmp_path, OUTFILE)

    print(f"Wrote: {OUTFILE}")
    print(f"Data rows: {len(rows)}")


def main():
    if os.environ.get("RUN_REBUILD", "0") != "1":
        print("RUN_REBUILD is not set to 1.")
        print("Verifying existing table_s9 only.")
        verify_existing_table(OUTFILE)
        print("08_collect_decomp_raw_loci.py completed")
        return

    print("RUN_REBUILD=1: rebuilding table_s9.")
    rebuild_table()
    verify_existing_table(OUTFILE)
    print("08_collect_decomp_raw_loci.py completed")


if __name__ == "__main__":
    main()
