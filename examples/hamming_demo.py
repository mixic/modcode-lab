# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Milan Nikolic
#
# This file is part of modcode-lab.
# modcode-lab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Encode and correct a one-bit error with cyclic Hamming(7,4)."""

import numpy as np

from modcode.fec.hamming import Hamming74


def main() -> None:
    code = Hamming74()
    data = np.array([1, 0, 1, 0], dtype=np.uint8)
    encoded = code.encode(data)
    corrupted = encoded.copy()
    corrupted[2] ^= 1
    decoded = code.decode(corrupted)
    print(f"data:      {''.join(map(str, data))}")
    print(f"encoded:   {''.join(map(str, encoded))}")
    print(f"corrupted: {''.join(map(str, corrupted))}")
    print(f"decoded:   {''.join(map(str, decoded))}")


if __name__ == "__main__":
    main()
