"""Shared SV-table schema and VCF metadata helpers for pipeline stages."""

from __future__ import annotations

import hashlib


METADATA_COLUMNS = ["sv_id", "chrom", "start", "end", "sv_type", "length", "filter", "imprecise"]


def canonical_sample_id(sample_id: str) -> str:
    if sample_id.startswith("GM") and sample_id[2:].isdigit():
        return f"NA{sample_id[2:]}"
    return sample_id


def normalize_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def parse_info(info_text: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def parse_length(info: dict[str, str | bool], start: int, end: int) -> str:
    value = info.get("SVLEN")
    if isinstance(value, str) and value not in {"", "."}:
        try:
            return str(abs(int(value.split(",")[0])))
        except ValueError:
            pass
    interval_length = end - start
    return str(interval_length) if interval_length > 0 else ""


def simplify_sv_id(sv_id: str, chrom: str, start: int, end: int, sv_type: str, max_length: int = 80) -> str:
    if len(sv_id) <= max_length:
        return sv_id
    digest = hashlib.sha1(sv_id.encode()).hexdigest()[:10]
    return f"SV_{chrom}_{start}_{end}_{sv_type}_{digest}"
