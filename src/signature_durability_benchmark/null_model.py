"""Matched random-signature null generation."""
from __future__ import annotations
from typing import Any
import numpy as np
import pandas as pd
from .scoring import score_signature_in_cohort

def generate_null_signatures(
    gene_count: int,
    universe: list[str],
    n_draws: int,
    seed: int,
) -> list[pd.DataFrame]:
    """Generate random gene signatures matched on size from the cohort universe."""
    rng = np.random.default_rng(seed)
    sigs = []
    for _ in range(n_draws):
        genes = rng.choice(universe, size=min(gene_count, len(universe)), replace=False)
        sig = pd.DataFrame({
            "gene_symbol": list(genes),
            "direction": ["up"] * len(genes),
            "weight": [1.0] * len(genes),
        })
        sigs.append(sig)
    return sigs

def null_separation_p(observed_effect: float, null_effects: list[float]) -> float:
    """Empirical p-value: fraction of null effects >= observed."""
    if not null_effects:
        return 1.0
    null_arr = np.array(null_effects)
    return float(np.mean(np.abs(null_arr) >= abs(observed_effect)))
