#!/usr/bin/env python3
from __future__ import annotations

import csv
import datetime as dt
import json
import sys
from pathlib import Path
from typing import Any

import requests

BASE = "https://gtexportal.org/api/v2"
DATASET = "gtex_v10"
GENCODE_VERSION = "v39"
GENOME_BUILD = "GRCh38/hg38"
VARIANT = "chr12_26273487_T_C_b38"
TISSUE = "Uterus"
GENES = ["BHLHE41", "SSPN", "ENSG00000256894", "ENSG00000255750"]
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "results/reproduced_gtex")
RAW = OUT / "raw_api"
RAW.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.headers.update({"Accept": "application/json", "User-Agent": "endo-ec-gwas-overlap/GTEx-V10-audit"})


def request(path: str, params: dict[str, Any], allow_400: bool = False):
    r = session.get(f"{BASE}/{path}", params=params, timeout=90)
    if allow_400 and r.status_code == 400:
        try:
            payload = r.json()
        except ValueError:
            payload = {"message": r.text}
        return r.url, r.status_code, payload
    r.raise_for_status()
    return r.url, r.status_code, r.json()


def rows(payload: dict[str, Any]):
    data = payload.get("data", [])
    return data if isinstance(data, list) else []


def save_raw(name: str, url: str, status: int, payload: dict[str, Any]):
    (RAW / name).write_text(json.dumps({"url": url, "status_code": status, "payload": payload}, indent=2), encoding="utf-8")


def resolve_gene(query: str):
    url, status, payload = request(
        "reference/gene",
        {
            "geneId": query,
            "gencodeVersion": GENCODE_VERSION,
            "genomeBuild": GENOME_BUILD,
            "itemsPerPage": 250,
        },
    )
    save_raw(f"gene_{query}.json", url, status, payload)
    matches = rows(payload)
    if not matches:
        raise RuntimeError(f"No GTEx V10/GENCODE v39 gene match for {query}")
    q = query.upper()
    for row in matches:
        symbol = str(row.get("geneSymbol", "")).upper()
        gid = str(row.get("gencodeId", ""))
        if symbol == q or gid.split(".")[0].upper() == q:
            return str(row.get("geneSymbol", query)), gid
    row = matches[0]
    return str(row.get("geneSymbol", query)), str(row["gencodeId"])


def pick(obj: dict[str, Any], *keys: str):
    for key in keys:
        value = obj.get(key)
        if value not in (None, ""):
            return value
    return ""


def classify(obj: dict[str, Any]):
    try:
        return "SIGNIFICANT" if float(obj["pValue"]) <= float(obj["pValueThreshold"]) else "NOT_SIGNIFICANT"
    except (KeyError, TypeError, ValueError):
        return "UNDETERMINED"

resolved = [resolve_gene(gene) for gene in GENES]

significant_rows = []
uterus_rows = []

for symbol, gencode_id in resolved:
    url, status, payload = request(
        "association/singleTissueEqtl",
        {
            "gencodeId": gencode_id,
            "variantId": VARIANT,
            "datasetId": DATASET,
            "itemsPerPage": 1000,
        },
    )
    save_raw(f"static_significant_{symbol}.json", url, status, payload)
    for row in rows(payload):
        significant_rows.append({
            "gene_symbol": symbol,
            "gencode_id": gencode_id,
            "variant_id": VARIANT,
            "tissue": row.get("tissueSiteDetailId", ""),
            "p_value": row.get("pValue", ""),
            "nes": row.get("nes", ""),
            "status": "SIGNIFICANT_BY_STATIC_ENDPOINT",
        })

    url, status, payload = request(
        "association/dyneqtl",
        {
            "tissueSiteDetailId": TISSUE,
            "gencodeId": gencode_id,
            "variantId": VARIANT,
            "datasetId": DATASET,
        },
        allow_400=True,
    )
    save_raw(f"dynamic_uterus_{symbol}.json", url, status, payload)

    uterus_rows.append({
        "gene_symbol": symbol,
        "gencode_id": gencode_id,
        "variant_id": VARIANT,
        "tissue": TISSUE,
        "api_status": "CALCULATED" if status == 200 else f"HTTP_{status}",
        "p_value": pick(payload, "pValue", "pvalue") if status == 200 else "",
        "p_value_threshold": pick(payload, "pValueThreshold", "pvalueThreshold") if status == 200 else "",
        "nes": pick(payload, "nes", "normalizedEffectSize") if status == 200 else "",
        "maf": pick(payload, "maf") if status == 200 else "",
        "significance": classify(payload) if status == 200 else "UNDETERMINED",
        "message": pick(payload, "message", "detail") if status != 200 else "",
    })


def write_tsv(path: Path, data):
    if not data:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(data[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(data)

write_tsv(OUT / "gtex_v10_rs9668810_significant_eqtl.tsv", significant_rows)
write_tsv(OUT / "gtex_v10_rs9668810_uterus_dynamic_eqtl.tsv", uterus_rows)

(OUT / "gtex_v10_rs9668810_query_metadata.json").write_text(
    json.dumps({
        "query_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "api_base": BASE,
        "dataset": DATASET,
        "gencode_version": GENCODE_VERSION,
        "genome_build": GENOME_BUILD,
        "variant": VARIANT,
        "target_tissue": TISSUE,
        "genes": [{"gene_symbol": s, "gencode_id": g} for s, g in resolved],
        "interpretation_rule": "Call a uterus association significant only when pValue <= pValueThreshold. HTTP errors or missing values are not null results.",
    }, indent=2),
    encoding="utf-8",
)

print("gene_symbol\tgencode_id\ttissue\tapi_status\tp_value\tp_value_threshold\tnes\tsignificance")
for row in uterus_rows:
    print(row["gene_symbol"], row["gencode_id"], row["tissue"], row["api_status"], row["p_value"], row["p_value_threshold"], row["nes"], row["significance"], sep="\t")
