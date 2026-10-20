#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Visualise prescribed residual-index densities for base 10.

Run from the repository root:

    python examples/residual_index_lab.py

The density formula is conditional on GRH. The finite prime counts are
unconditional but converge slowly; the plot shows that explicitly.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from modcode.number_theory.residual_index import (
    artin_constant,
    count_residual_indices,
    known_decimal_density_multiple,
    residual_index_convergence,
    residual_index_density,
)
from modcode.viz.number_theory import (
    plot_residual_index_convergence,
    plot_residual_index_density,
    plot_residual_index_histogram,
)

plt.switch_backend("Agg")

OUTPUT_DIR = Path(__file__).resolve().parent / "figures"
RULE = "=" * 78


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    constant = artin_constant(10_000)
    indexes = [1, 2, 3, 4, 6, 8]
    densities = {
        m: residual_index_density(m, terms=100_000) for m in indexes
    }
    predicted = {
        1: constant,
        2: 0.75 * constant,
        4: 29 / 152 * constant,
    }

    limits = [1_000, 5_000, 25_000, 100_000, 500_000]
    convergence = residual_index_convergence(limits)
    counts = count_residual_indices(1_000_000)

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 4.6))
    figure.suptitle(
        "Primes with residual index c_p = (p-1)/ord_p(10)",
        fontsize=12,
    )
    plot_residual_index_density(
        axes[0],
        indexes,
        [densities[m] for m in indexes],
        artin=constant,
        title="predicted density (GRH)",
    )
    plot_residual_index_convergence(
        axes[1],
        limits,
        convergence,
        predicted,
        title="finite counts approach the prediction",
    )
    plot_residual_index_histogram(
        axes[2],
        counts.counts,
        total=counts.total,
        max_index=12,
        title=f"observed at X={counts.limit:,}",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    path = OUTPUT_DIR / "residual_index_density.png"
    figure.savefig(path, dpi=140)
    plt.close(figure)

    print(RULE)
    print("RESIDUAL-INDEX DENSITIES FOR BASE 10")
    print(RULE)
    print(f"Artin constant A            : {constant:.10f}")
    print(f"prime count up to 1,000,000 : {counts.total:,}")
    print()
    print(f"{'m':>4}{'density':>12}{'multiple of A':>16}{'observed':>12}")
    print("-" * 78)
    for m in indexes:
        multiple = known_decimal_density_multiple(m)
        multiple_text = (
            str(multiple) if multiple is not None else "not tabulated"
        )
        print(
            f"{m:>4}{densities[m]:>12.6f}{multiple_text:>16}"
            f"{counts.ratio(m):>12.6f}"
        )
    print()
    print("  m=1 is Artin's full-period case: density A.")
    print("  m=2 is the half-period case: density 3A/4, not a coding")
    print("  'semi-cyclic' class; the word is overloaded across fields.")
    print("  m=4 is 29A/152; the frequently quoted 29A/52 does not match")
    print("  the field-degree formula or the direct count.")
    print(f"  figure -> {path}")
    print(RULE)


if __name__ == "__main__":
    main()
