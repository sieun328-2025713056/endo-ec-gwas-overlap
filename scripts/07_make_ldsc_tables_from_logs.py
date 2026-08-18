#!/usr/bin/env python3
import csv
import glob
import math
import os
import re

PROJECT = os.environ.get("PROJECT", os.path.expanduser("~/projects/endo_ecancer"))
H2 = os.environ.get("H2", os.path.join(PROJECT, "results", "h2"))
RG = os.environ.get("RG", os.path.join(PROJECT, "results", "rg"))
TABLES = os.environ.get("TABLES", os.path.join(PROJECT, "results", "tables"))
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

# Multiple-testing correction across the six primary cross-trait rg tests.
# Benjamini-Hochberg FDR is supportive; Bonferroni correction across six
# primary tests (alpha = 0.05 / 6) defines the primary significance threshold.

primary_comparisons = [
    ("rg_adenomyosis_EUR_allEC", "Adenomyosis", "All-histology EC"),
    ("rg_adenomyosis_EUR_EEC", "Adenomyosis", "Endometrioid EC"),
    ("rg_endometriosis_EUR_allEC", "Overall endometriosis", "All-histology EC"),
    ("rg_endometriosis_EUR_EEC", "Overall endometriosis", "Endometrioid EC"),
    ("rg_endometriosis_wo_adeno_EUR_allEC", "Endometriosis without adenomyosis", "All-histology EC"),
    ("rg_endometriosis_wo_adeno_EUR_EEC", "Endometriosis without adenomyosis", "Endometrioid EC"),
]

rg_by_comparison = {row["Comparison"]: row for row in rg_rows}

primary_rows = []
for key, endo_label, ec_label in primary_comparisons:
    if key not in rg_by_comparison:
        raise RuntimeError(f"Missing primary LDSC comparison: {key}")
    row = rg_by_comparison[key]
    primary_rows.append({
        "Comparison": key,
        "Endometriosis_phenotype": endo_label,
        "Endometrial_cancer_outcome": ec_label,
        "rg": float(row["rg"]),
        "SE": float(row["SE"]),
        "Z": float(row["Z"]),
        "P": float(row["P"]),
    })

# Benjamini-Hochberg adjusted P-values.
m = len(primary_rows)
order = sorted(range(m), key=lambda i: primary_rows[i]["P"])
bh = [None] * m
running_min = 1.0

for rank_from_end, idx in enumerate(reversed(order), start=1):
    rank = m - rank_from_end + 1
    adjusted = primary_rows[idx]["P"] * m / rank
    running_min = min(running_min, adjusted)
    bh[idx] = min(running_min, 1.0)

for i, row in enumerate(primary_rows):
    row["BH_FDR_q"] = bh[i]
    row["Bonferroni_significant"] = "Yes" if row["P"] < (0.05 / 6) else "No"

with open(
    os.path.join(TABLES, "table_s7_ldsc_rg_multiple_testing_correction.tsv"),
    "w",
    newline=""
) as f:
    fieldnames = [
        "Comparison",
        "Endometriosis_phenotype",
        "Endometrial_cancer_outcome",
        "rg",
        "SE",
        "Z",
        "P",
        "BH_FDR_q",
        "Bonferroni_significant",
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(primary_rows)


# Approximate between-phenotype rg difference tests within each EC outcome.
# These tests use sqrt(SE1^2 + SE2^2) and therefore do not model covariance
# arising from the shared EC GWAS or non-independent endometriosis estimates.

phenotype_pairs = [
    ("Adenomyosis", "Overall endometriosis"),
    ("Adenomyosis", "Endometriosis without adenomyosis"),
    ("Endometriosis without adenomyosis", "Overall endometriosis"),
]

outcomes = ["All-histology EC", "Endometrioid EC"]

lookup = {
    (row["Endometriosis_phenotype"], row["Endometrial_cancer_outcome"]): row
    for row in primary_rows
}

difference_rows = []

for outcome in outcomes:
    for phenotype1, phenotype2 in phenotype_pairs:
        row1 = lookup[(phenotype1, outcome)]
        row2 = lookup[(phenotype2, outcome)]

        difference = row1["rg"] - row2["rg"]
        se_difference = math.sqrt(row1["SE"] ** 2 + row2["SE"] ** 2)
        z_difference = difference / se_difference

        p_difference = math.erfc(abs(z_difference) / math.sqrt(2.0))

        difference_rows.append({
            "Endometrial_cancer_outcome": outcome,
            "Phenotype_1": phenotype1,
            "Phenotype_2": phenotype2,
            "rg_1": row1["rg"],
            "rg_2": row2["rg"],
            "rg_difference": difference,
            "SE_difference_approx": se_difference,
            "Z_difference_approx": z_difference,
            "P_difference_approx": p_difference,
        })

with open(
    os.path.join(TABLES, "table_s3_rg_approx_difference_test_from_logs.tsv"),
    "w",
    newline=""
) as f:
    fieldnames = [
        "Endometrial_cancer_outcome",
        "Phenotype_1",
        "Phenotype_2",
        "rg_1",
        "rg_2",
        "rg_difference",
        "SE_difference_approx",
        "Z_difference_approx",
        "P_difference_approx",
    ]
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(difference_rows)

print("Wrote LDSC tables to", TABLES)
