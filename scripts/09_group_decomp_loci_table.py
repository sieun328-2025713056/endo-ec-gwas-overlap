#!/usr/bin/env python3
"""
09_group_decomp_loci_table.py

Verify or optionally rebuild table_s10_decomp_locus_grouped_from_file.tsv.

Default behavior:
    - Verify that table_s10 already exists and has data rows.
    - Do not touch existing table_s10.

Optional rebuild:
    RUN_REBUILD=1 python3 scripts/09_group_decomp_loci_table.py

Safety rule:
    - If table_s9 is empty or grouping produces zero rows, existing table_s10 is not overwritten.
"""

import csv
import os
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(SCRIPT_DIR)

TABLES = os.path.join(PROJECT, "results", "tables")
INFILE = os.path.join(TABLES, "table_s9_decomp_conjfdr_loci_raw.tsv")
OUTFILE = os.path.join(TABLES, "table_s10_decomp_locus_grouped_from_file.tsv")

CLUSTER_GAP_BP = 500000

PAIR_ORDER = [
    "ADENO_ALL_EC",
    "ADENO_EEC",
    "OVERALL_ENDO_ALL_EC",
    "OVERALL_ENDO_EEC",
    "ENDO_WO_ADENO_ALL_EC",
    "ENDO_WO_ADENO_EEC",
]

HEADER = [
    "Group_ID",
    "CHR",
    "Observed_region_start",
    "Observed_region_end",
    "N_rows_in_group",
    "Best_lead_SNP",
    "Best_BP",
    "Best_conjFDR",
    "Detected_exposures",
    "Detected_outcomes",
    "Directions_in_group",
    "Classification",
    "ADENO_ALL_EC",
    "ADENO_EEC",
    "OVERALL_ENDO_ALL_EC",
    "OVERALL_ENDO_EEC",
    "ENDO_WO_ADENO_ALL_EC",
    "ENDO_WO_ADENO_EEC",
    "Caution",
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

        missing = [c for c in HEADER if c not in reader.fieldnames]
        if missing:
            fail(f"Missing required columns in {path}: {missing}")

        n = sum(1 for _ in reader)

    if n <= 0:
        fail(f"Existing table has zero data rows: {path}")

    print(f"Verified existing table_s10: {path}")
    print(f"Data rows: {n}")
    return n


def read_s9_rows():
    if not os.path.exists(INFILE):
        fail(f"Missing input table_s9: {INFILE}")

    rows = []

    with open(INFILE, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = [
            "Analysis",
            "Exposure",
            "Outcome",
            "CHR",
            "BP",
            "Lead_SNP",
            "conjFDR",
            "Direction",
            "Z_exposure",
            "Z_outcome",
        ]
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            fail(f"Missing required columns in table_s9: {missing}")

        for r in reader:
            r["BP_int"] = int(r["BP"])
            rows.append(r)

    if len(rows) == 0:
        fail("table_s9 has zero rows. Existing table_s10 will not be overwritten.")

    rows = sorted(rows, key=lambda r: (int(r["CHR"]), r["BP_int"]))
    return rows


def fmt_row(r):
    return (
        f"{r['Lead_SNP']}@{r['BP']}; "
        f"q={r['conjFDR']}; "
        f"dir={r['Direction']}; "
        f"Z={r['Z_exposure']}/{r['Z_outcome']}"
    )


def classify(group_rows):
    exposures = set(r["Exposure"] for r in group_rows)
    outcomes = set(r["Outcome"] for r in group_rows)

    has_adeno = "ADENO" in exposures
    has_overall = "OVERALL_ENDO" in exposures
    has_woadeno = "ENDO_WO_ADENO" in exposures

    both_outcomes = "ALL_EC" in outcomes and "EEC" in outcomes

    if has_adeno and has_overall and has_woadeno:
        cls = "pan-spectrum"
    elif has_adeno and has_overall and not has_woadeno:
        cls = "ADENO/overall-aligned"
    elif has_adeno and not has_overall and not has_woadeno:
        cls = "adenomyosis-biased candidate"
    elif has_woadeno and not has_adeno:
        cls = "endometriosis-core candidate"
    else:
        cls = "inconclusive"

    if not both_outcomes:
        cls = cls + "; single-outcome only"

    return cls


def group_rows(rows):
    groups = []

    for r in rows:
        chrom = r["CHR"]
        bp = r["BP_int"]

        if not groups:
            groups.append({"CHR": chrom, "start": bp, "end": bp, "rows": [r]})
            continue

        g = groups[-1]
        if chrom == g["CHR"] and bp <= g["end"] + CLUSTER_GAP_BP:
            g["rows"].append(r)
            g["end"] = max(g["end"], bp)
        else:
            groups.append({"CHR": chrom, "start": bp, "end": bp, "rows": [r]})

    return groups


def rebuild_table():
    rows = read_s9_rows()
    groups = group_rows(rows)

    if len(groups) == 0:
        fail("Grouped loci count is zero. Existing table_s10 will not be overwritten.")

    out_rows = []

    for i, g in enumerate(groups, start=1):
        group_rows = g["rows"]
        best = sorted(group_rows, key=lambda r: float(r["conjFDR"]))[0]

        pair_cells = {}
        for p in PAIR_ORDER:
            rs = [r for r in group_rows if r["Analysis"] == p]
            if rs:
                rs = sorted(rs, key=lambda r: float(r["conjFDR"]))
                pair_cells[p] = " | ".join(fmt_row(r) for r in rs)
            else:
                pair_cells[p] = "not detected"

        directions = sorted(set(r["Direction"] for r in group_rows))
        outcomes = sorted(set(r["Outcome"] for r in group_rows))
        exposures = sorted(set(r["Exposure"] for r in group_rows))

        out_rows.append({
            "Group_ID": f"G{i}",
            "CHR": g["CHR"],
            "Observed_region_start": str(g["start"]),
            "Observed_region_end": str(g["end"]),
            "N_rows_in_group": str(len(group_rows)),
            "Best_lead_SNP": best["Lead_SNP"],
            "Best_BP": best["BP"],
            "Best_conjFDR": best["conjFDR"],
            "Detected_exposures": ";".join(exposures),
            "Detected_outcomes": ";".join(outcomes),
            "Directions_in_group": ";".join(directions),
            "Classification": classify(group_rows),
            "ADENO_ALL_EC": pair_cells["ADENO_ALL_EC"],
            "ADENO_EEC": pair_cells["ADENO_EEC"],
            "OVERALL_ENDO_ALL_EC": pair_cells["OVERALL_ENDO_ALL_EC"],
            "OVERALL_ENDO_EEC": pair_cells["OVERALL_ENDO_EEC"],
            "ENDO_WO_ADENO_ALL_EC": pair_cells["ENDO_WO_ADENO_ALL_EC"],
            "ENDO_WO_ADENO_EEC": pair_cells["ENDO_WO_ADENO_EEC"],
            "Caution": (
                "Non-detection in a pair indicates lack of detection under current "
                "GWAS power/phenotype definition, not evidence of absence."
            ),
        })

    fd, tmp_path = tempfile.mkstemp(
        prefix="table_s10_decomp_locus_grouped_from_file.",
        suffix=".tmp",
        dir=TABLES,
        text=True,
    )
    os.close(fd)

    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    os.replace(tmp_path, OUTFILE)

    print(f"Wrote: {OUTFILE}")
    print(f"Input rows: {len(rows)}")
    print(f"Grouped loci: {len(out_rows)}")


def main():
    if os.environ.get("RUN_REBUILD", "0") != "1":
        print("RUN_REBUILD is not set to 1.")
        print("Verifying existing table_s10 only.")
        verify_existing_table(OUTFILE)
        print("09_group_decomp_loci_table.py completed")
        return

    print("RUN_REBUILD=1: rebuilding table_s10.")
    rebuild_table()
    verify_existing_table(OUTFILE)
    print("09_group_decomp_loci_table.py completed")


if __name__ == "__main__":
    main()
