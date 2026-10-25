"""
Experimental confirmation of the multi-tweak degree bounds.

For small k the predicted cube dimension is tiny, so the tweak cube can be summed
exhaustively against the real cipher. Two things are checked:

  soundness  a cube of dimension deg_bound + 1 really does sum to zero
  tightness  a cube of dimension deg_bound generally does NOT

Tightness matters: a bound that is merely correct but wildly loose would make the
2^56 claim at k = 4 worthless.
"""

import random

from blinkit.cipher.blink_ref import set_variant
from blinkit.cipher.blink_ops import build_ops, nonlinear_prefix, evaluate, CFG
from blinkit.divprop.blink_mt import ops_from_seam, truncate_sboxes, degree_bound

PI1_SBOXES = 3           # pi1 is constant under a fixed plaintext


def tweak_cube_sum(n_sbox_total, tweak_bits, m_fixed, key, base_tweak):
    """XOR of ciphertexts over an affine tweak subspace."""
    ops = nonlinear_prefix(build_ops(), n_sbox_total)
    acc = 0
    for combo in range(1 << len(tweak_bits)):
        t = base_tweak
        for idx, bit in enumerate(tweak_bits):
            if (combo >> idx) & 1:
                t ^= 1 << bit
        acc ^= evaluate(ops, m_fixed, t, key, CFG)
    return acc


def main():
    set_variant(64, 64, 2, 3)
    random.seed(20260805)
    full = ops_from_seam('h1')

    print("multi-tweak cube sums against the real cipher\n")
    for k in (1, 2):
        D = degree_bound(truncate_sboxes(full, k))
        total = PI1_SBOXES + k
        print(f"k = {k}  (cipher truncated to {total} S layers)  deg <= {D}")

        for trial in range(3):
            key = random.getrandbits(448)
            m_fixed = random.getrandbits(64)
            # keep the whole cube away from t = 0, which the spec forbids
            base = random.getrandbits(64) | (1 << 63)

            bits_ok = random.sample(range(63), D + 1)
            s_ok = tweak_cube_sum(total, bits_ok, m_fixed, key, base)

            bits_tight = bits_ok[:D]
            s_tight = tweak_cube_sum(total, bits_tight, m_fixed, key, base)

            print(f"  trial {trial}: dim {D+1:2d} -> {s_ok:#018x} "
                  f"{'OK' if s_ok == 0 else 'FAILED'}   "
                  f"dim {D:2d} -> {s_tight:#018x} "
                  f"{'(tight)' if s_tight else '(also zero)'}")
        print()


if __name__ == "__main__":
    main()
