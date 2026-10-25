"""
End-to-end key recovery, at a scale small enough to run in seconds.

Structure of the attack, using k=1 so the cube is only dimension 4:

    S#4   distinguisher:  XOR X4 = 0                      (cube dim 4)
    S#5   free extension: XOR S(X5) = 0                   (no key needed)
    S#6   KEY RECOVERY:   guess 16 bits of K4, invert one round, check

Between S#5 and S#6 the cipher applies  M, then AK with K4 = rk4 ^ rc4, then P,
then S.  So

    X6 = S( P( M(X5) ^ K4 ) )    =>    X5 = M( Pbar( S(X6) ) ^ K4 )

using that both S and M are involutions.

One nibble of X5 depends on one MixColumn column, i.e. 4 nibbles of K4 = 16 bits.
And that same column yields 4 nibbles of X5, so a single 16-bit guess is checked
by a 16-bit condition -- enough to pin it down, with a second cube to break ties.

Scaled up to the real attack this is 2^16 guesses x 2^37 data ~= 2^53 per column.
The point here is to confirm the technique actually recovers the right key bits
before trusting that arithmetic.
"""

import random

from blinkit.cipher.blink_ref import (SBOX, P_64_INV, set_variant, S_layer, P_layer,
                       to_state, to_int, split_master, RC)
from blinkit.cipher.blink_ops import build_ops, nonlinear_prefix, evaluate, CFG

NCOL = 4


def collect(n_sbox, tweak_bits, m, key, base):
    """Observed states after n_sbox S layers, over the tweak cube."""
    ops = nonlinear_prefix(build_ops(), n_sbox)
    out = []
    for combo in range(1 << len(tweak_bits)):
        t = base
        for idx, bit in enumerate(tweak_bits):
            if (combo >> idx) & 1:
                t ^= 1 << bit
        out.append(evaluate(ops, m, t, key, CFG))
    return out


def column_samples(x6_list, col):
    """Pbar(S(X6)) restricted to one MixColumn column, per sample."""
    rows = []
    for x6 in x6_list:
        y = P_layer(S_layer(to_state(x6, True)), P_64_INV)
        rows.append((y[col], y[col + NCOL], y[col + 2 * NCOL],
                     y[col + 3 * NCOL]))
    return rows


def check_guess(rows, g):
    """XOR of S(X5) over the cube, for the four nibbles of this column."""
    g0, g1, g2, g3 = g
    a0 = a1 = a2 = a3 = 0
    for y0, y1, y2, y3 in rows:
        i0 = y0 ^ g0
        i1 = y1 ^ g1
        i2 = y2 ^ g2
        i3 = y3 ^ g3
        a0 ^= SBOX[i1 ^ i2 ^ i3]          # M row 0: XOR of the other three
        a1 ^= SBOX[i0 ^ i2 ^ i3]
        a2 ^= SBOX[i0 ^ i1 ^ i3]
        a3 ^= SBOX[i0 ^ i1 ^ i2]
    return a0 | a1 | a2 | a3


def recover_column(rows_a, rows_b, col, rows_c=None):
    """2^16 exhaustive search over one column of K4, filtered by two cubes."""
    cands = []
    for g in range(1 << 16):
        gg = (g & 0xf, (g >> 4) & 0xf, (g >> 8) & 0xf, (g >> 12) & 0xf)
        if (check_guess(rows_a, gg) == 0 and check_guess(rows_b, gg) == 0
                and (rows_c is None or check_guess(rows_c, gg) == 0)):
            cands.append(gg)
    return cands


def main():
    set_variant(64, 64, 2, 3)
    random.seed(11)

    key = random.getrandbits(448)
    m = random.getrandbits(64)
    _, _, rk = split_master(key)

    # the round key added between S#5 and S#6 is  rk4 ^ rc4
    K4 = to_state(rk[3] ^ RC[3], True)
    print(f"target round key K4 = rk4 ^ rc4 = "
          f"{to_int(K4, True):#018x}\n")

    # two independent cubes, both dimension 4 (deg = 3 at this depth)
    cubes = []
    for _ in range(3):
        base = random.getrandbits(64) | (1 << 63)
        bits = random.sample(range(63), 4)
        cubes.append((base, bits))

    # sanity: the distinguisher and its free extension must hold first
    for r, label in ((4, "distinguisher  XOR X4"),
                     (5, "free extension XOR S(X5)")):
        base, bits = cubes[0]
        xs = collect(r, bits, m, key, base)
        if r == 4:
            s = 0
            for x in xs:
                s ^= x
        else:
            s = 0
            for x in xs:
                s ^= to_int(S_layer(to_state(x, True)), True)
        print(f"  {label:26s} = {s:#018x}  {'OK' if s == 0 else 'FAIL'}")
    print()

    # the attack proper, at 6 S layers
    obs = [collect(6, b, m, key, base) for base, b in cubes]

    total_ok = 0
    for col in range(NCOL):
        rows_a = column_samples(obs[0], col)
        rows_b = column_samples(obs[1], col)
        cands = recover_column(rows_a, rows_b, col, column_samples(obs[2], col))
        truth = (K4[col], K4[col + NCOL], K4[col + 2 * NCOL],
                 K4[col + 3 * NCOL])
        hit = truth in cands
        total_ok += hit
        print(f"  column {col}: {len(cands):3d} candidate(s) from 2^16   "
              f"true key {truth} {'RECOVERED' if hit else 'MISSED'}")

    print(f"\n{total_ok}/4 columns recovered "
          f"({'16' if total_ok else '0'} key bits each)")


if __name__ == "__main__":
    main()
