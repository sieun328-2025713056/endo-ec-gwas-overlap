#!/usr/bin/env python3
import csv
import glob
import math
import os
import re

PROJECT = os.environ.get("PROJECT", os.path.expanduser("~/projects/endo_ecancer"))
H2 = os.path.join(PROJECT, "results", "h2")
RG = os.path.join(PROJECT, "results", "rg")
TABLES = os.path.join(PROJECT, "results", "tables")
os.makedirs(TABLES, exist_ok=True)

trait_labels = {
    "endometriosis_EUR": "Overall endometriosis",
    "adenomyosis_EUR": "Adenomyosis",
    "endometriosis_wo_adeno_EUR": "Endometriosis without adenomyosis",
    "endometrial_cancer_allEC": "Endometrial cancer, all EC",
    "endometrial_cancer_EEC": "Endometrioid endometrial cancer",
}

# h2 QC table.
h2_rows = []
for path in sorted(glob.glob(os.path.join(H2, "*.log"))):
    name = os.path.basename(path).replace(".log", "")
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"Total Observed scale h2?:\s*([\-0-9.eE]+)\s*\(([\-0-9.eE]+)\)", txt)
    if not m:
        continue
    h2 = float(m.group(1)); se = float(m.group(2))
    lam = re.search(r"Lambda GC:\s*([\-0-9.eE]+)", txt)
    mean_chi = re.search(r"Mean Chi\^2:\s*([\-0-9.eE]+)", txt)
    intercept = re.search(r"Intercept:\s*([\-0-9.eE]+)\s*\(([\-0-9.eE]+)\)", txt)
    ratio = re.search(r"Ratio:\s*([\-0-9.eE]+)\s*\(([\-0-9.eE]+)\)", txt)
    ratio_lt0 = "Ratio < 0" in txt
    h2_rows.append({
        "Phenotype": trait_labels.get(name, name),
        "Log_file": os.path.basename(path),
        "h2_observed": f"{h2:.4g}",
        "SE": f"{se:.4g}",
        "h2_Z": f"{h2/se:.4g}" if se else "NA",
        "Lambda_GC": lam.group(1) if lam else "NA",
        "Mean_Chi2": mean_chi.group(1) if mean_chi else "NA",
        "Intercept": intercept.group(1) if intercept else "NA",
        "Intercept_SE": intercept.group(2) if intercept else "NA",
        "Ratio": "<0" if ratio_lt0 else (ratio.group(1) if ratio else "NA"),
        "Ratio_SE": "NA" if ratio_lt0 or not ratio else ratio.group(2),
    })

with open(os.path.join(TABLES, "table_s1_ldsc_h2_qc_from_logs.tsv"), "w", newline="") as f:
    fieldnames = ["Phenotype","Log_file","h2_observed","SE","h2_Z","Lambda_GC","Mean_Chi2","Intercept","Intercept_SE","Ratio","Ratio_SE"]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader(); w.writerows(h2_rows)

# rg table.
rg_rows = []
for path in sorted(glob.glob(os.path.join(RG, "*.log"))):
    txt = open(path, encoding="utf-8", errors="replace").read()
    lines = txt.splitlines()
    header_idx = None
    for i, line in enumerate(lines):
        if "Summary of Genetic Correlation Results" in line:
            header_idx = i
    if header_idx is None:
        continue
    # Find the table row containing .sumstats.gz and numeric rg/se/z/p.
    for line in lines[header_idx:]:
        if ".sumstats.gz" not in line:
            continue
        parts = line.split()
        nums = []
        for token in parts:
            try:
                nums.append(float(token))
            except ValueError:
                pass
        if len(nums) >= 4:
            rg, se, z, p = nums[0], nums[1], nums[2], nums[3]
            break
    else:
        continue
    basename = os.path.basename(path).replace(".log", "")
    rg_rows.append({"Log_file": os.path.basename(path), "Comparison": basename, "rg": rg, "SE": se, "Z": z, "P": p})

with open(os.path.join(TABLES, "table_s2_ldsc_rg_from_logs.tsv"), "w", newline="") as f:
    fieldnames = ["Log_file","Comparison","rg","SE","Z","P"]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader(); w.writerows(rg_rows)

print("Wrote LDSC tables to", TABLES)
