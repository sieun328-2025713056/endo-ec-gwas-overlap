#!/usr/bin/env python3
import csv
from pathlib import Path

base = Path("results/reproduced_gtex")
static = base / "gtex_v10_rs9668810_significant_eqtl.tsv"
dynamic = base / "gtex_v10_rs9668810_uterus_dynamic_eqtl.tsv"
output = base / "gtex_v10_rs9668810_uterus_final_audit.tsv"

with static.open(encoding="utf-8") as f:
    significant = {
        (r["gencode_id"].split(".")[0], r["variant_id"], r["tissue"].lower())
        for r in csv.DictReader(f, delimiter="\t")
        if r["tissue"].lower() == "uterus"
    }

with dynamic.open(encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

for r in rows:
    key = (
        r["gencode_id"].split(".")[0],
        r["variant_id"],
        r["tissue"].lower(),
    )
    if r["api_status"] != "CALCULATED":
        status = "UNDETERMINED_API_ERROR"
    elif key in significant:
        status = "SIGNIFICANT_GTEX_V10"
    else:
        status = "NOT_SIGNIFICANT_GTEX_V10"
    r["final_status"] = status
    r["classification_basis"] = (
        "Dynamic uterus eQTL result plus GTEx V10 significant endpoint"
    )

with output.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=list(rows[0].keys()),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote: {output}")
