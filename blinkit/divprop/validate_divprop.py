"""
Soundness check for blink_divprop.

For cubes small enough to sum exhaustively, compare what the SAT model claims is
balanced against what the real cipher actually does. A model-claimed balanced bit
that comes out 1 in the experiment means the model is unsound and any
distinguisher built on it is worthless.

Run:  python3 validate_divprop.py
"""

import random

from blinkit.cipher.blink_ref import set_variant
from blinkit.divprop.blink_divprop import Distinguisher, verify_experimentally, NBITS

KEY = int("d6a102d888a467e4d1d7dec33a246943e07c1dc6f302c57e762c2df9de6f0d21"
          "6dd387874a0b52ce3022e0ad78c78a0697779021b38e7fa1", 16)
TWEAK = 0x0123456789abcdef


def cell_cube(cells):
    """All four bits of each listed cell."""
    return [4 * c + b for c in cells for b in range(4)]


def run_case(n_sbox, active, key, tweak, label):
    d = Distinguisher(n_sbox)
    bal = d.balanced_bits(active)
    d.close()
    if not bal:
        print(f"  {label}: model claims nothing balanced — skipped")
        return None
    acc, bad = verify_experimentally(n_sbox, active, bal, key, tweak)
    status = "OK" if not bad else f"UNSOUND, bits {bad} are 1"
    print(f"  {label}: model says {len(bal):2d} balanced, "
          f"sum = {acc:#018x} -> {status}")
    return not bad


def main():
    set_variant(64, 64, 2, 3)
    random.seed(20260805)
    all_ok = True

    print("exhaustive cube checks (real cipher, test-vector key)\n")
    for n_sbox in range(2, 7):
        print(f"S layers = {n_sbox}")
        for cells in ([0, 1, 2], [0, 4, 8], [0, 1, 2, 3], [0, 5, 10, 15]):
            active = cell_cube(cells)
            r = run_case(n_sbox, active, KEY, TWEAK, f"cells {cells}")
            if r is False:
                all_ok = False
        print()

    print("random keys and tweaks, 16-dim cube, S layers = 5\n")
    for trial in range(4):
        k = random.getrandbits(448)
        t = random.getrandbits(64) | 1          # spec requires t != 0
        active = cell_cube([0, 1, 2, 3])
        r = run_case(5, active, k, t, f"trial {trial}")
        if r is False:
            all_ok = False

    print("\nVERDICT:", "model sound on all tested cubes" if all_ok
          else "MODEL IS UNSOUND — do not use")


if __name__ == "__main__":
    main()
