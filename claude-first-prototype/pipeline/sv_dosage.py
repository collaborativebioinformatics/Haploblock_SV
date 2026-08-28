"""Shared conversion of SV genotype columns to numeric alternate-allele dosage."""

from __future__ import annotations

import numpy as np
import pandas as pd

from stage4_classify_af import genotype_counts
from sv_contract import METADATA_COLUMNS


COMMON_GT_DOSAGES = {
    "0|0": 0.0, "0/0": 0.0, "0": 0.0,
    "0|1": 1.0, "1|0": 1.0, "0/1": 1.0, "1/0": 1.0, "1": 1.0,
    "1|1": 2.0, "1/1": 2.0,
    ".": np.nan, ".|.": np.nan, "./.": np.nan,
}


def dosage_matrix_from_genotypes(sv: pd.DataFrame, samples: list[str]) -> np.ndarray:
    """Return variant-by-sample dosages, with a fallback for uncommon GT encodings."""
    matrix = np.empty((len(sv), len(samples)), dtype=float)
    for sample_index, sample in enumerate(samples):
        genotype = sv[sample]
        if pd.api.types.is_numeric_dtype(genotype):
            matrix[:, sample_index] = pd.to_numeric(genotype, errors="coerce")
            continue
        mapped = genotype.map(COMMON_GT_DOSAGES)
        unknown = mapped.isna() & genotype.notna() & ~genotype.isin((".", ".|.", "./."))
        if unknown.any():
            alternate, called = genotype_counts(genotype)
            matrix[:, sample_index] = np.where(called > 0, alternate, np.nan)
        else:
            matrix[:, sample_index] = mapped.to_numpy(dtype=float, na_value=np.nan)
    return matrix


def dosage_table(sv: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    """Retain record metadata beside the numeric sample dosage matrix."""
    return pd.concat(
        [
            sv[METADATA_COLUMNS].reset_index(drop=True),
            pd.DataFrame(dosage_matrix_from_genotypes(sv, samples), columns=samples),
        ],
        axis=1,
    )
