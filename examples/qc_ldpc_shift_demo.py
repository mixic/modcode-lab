# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Print a cyclic identity block and an illustrative lifted QC matrix."""

import numpy as np

from modcode.fec.qc_ldpc import circulant_identity, expand_qc_base_matrix


def main() -> None:
    lifting_factor = 5
    shift = 2
    block = circulant_identity(lifting_factor, shift)
    print(f"{lifting_factor}x{lifting_factor} identity shifted right by {shift}:")
    print(block)

    base_matrix = np.array(
        [
            [2, -1, 0],
            [-1, 1, 3],
        ],
        dtype=np.int64,
    )
    expanded = expand_qc_base_matrix(base_matrix, lifting_factor)
    print("\nIllustrative base matrix (-1 denotes a zero block):")
    print(base_matrix)
    print("Expanded binary parity-check matrix:")
    print(expanded)


if __name__ == "__main__":
    main()
