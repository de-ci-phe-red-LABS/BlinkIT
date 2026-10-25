"""
Does the distinguisher extend one S layer for free?

Claim. Let X_r be the state after r S layers, with the cube sum XOR X_r = 0.
One layer further the state is X_{r+1} = S( L(X_r) + const ), where L is linear
and `const` is a round key, round constant, or a secret tweak hash h(t).

Summing over the cube:

    XOR S(X_{r+1}) = XOR ( L(X_r) + const )        (S is an INVOLUTION)
                   = L( XOR X_r ) + XOR const
                   = L(0) + 0                      (2^d copies of a constant)
                   = 0

The secret tweak hash vanishes too: h is LINEAR in t and the cube is an affine
subspace of dimension >= 2, so XOR t = 0 and hence XOR h(t) = h(0) = 0.

So a distinguisher on r S layers should also be a distinguisher on r+1 layers,
with the check "XOR S(ciphertext) = 0" instead of "XOR ciphertext = 0" -- and it
needs no key material at all.

This script tests that at small k, where the cubes are tiny.
"""

import random

from blinkit.cipher.blink_ref import set_variant, S_layer, to_state, to_int
from blinkit.cipher.blink_ops import build_ops, nonlinear_prefix, evaluate, CFG

PI1 = 3          # pi1 is constant under a fixed plaintext


def cube_sums(n_sbox, tweak_bits, m_fixed, key, base_tweak):
    """Return (plain sum, sum after applying S to each output)."""
    ops = nonlinear_prefix(build_ops(), n_sbox)
    acc_plain = 0
    acc_sinv = 0
    for combo in range(1 << len(tweak_bits)):
        t = base_tweak
        for idx, bit in enumerate(tweak_bits):
            if (combo >> idx) & 1:
                t ^= 1 << bit
        out = evaluate(ops, m_fixed, t, key, CFG)
        acc_plain ^= out
        acc_sinv ^= to_int(S_layer(to_state(out, True)), True)
    return acc_plain, acc_sinv


def main():
    set_variant(64, 64, 2, 3)
    random.seed(7)
    print("k  = S layers traversed after the h1 injection")
    print("r  = PI1 + k = total S layers covered\n")

    for k, dim in ((1, 4), (2, 8), (3, 17)):
        r = PI1 + k
        print(f"k={k}  cube dim {dim}")
        for trial in range(3):
            key = random.getrandbits(448)
            m = random.getrandbits(64)
            base = random.getrandbits(64) | (1 << 63)
            bits = random.sample(range(63), dim)

            p_r, _ = cube_sums(r, bits, m, key, base)
            p_r1, s_r1 = cube_sums(r + 1, bits, m, key, base)

            print(f"  trial {trial}: "
                  f"r={r} plain {p_r:#018x} {'OK' if p_r == 0 else 'FAIL'}  |  "
                  f"r={r+1} plain {p_r1:#018x} {'zero' if p_r1 == 0 else 'nonzero'}"
                  f"  S-applied {s_r1:#018x} "
                  f"{'OK' if s_r1 == 0 else 'FAIL'}")
        print()


if __name__ == "__main__":
    main()
