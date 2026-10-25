"""
Soundness cross-check: conventional BDP S-box table vs exact-ANF table.

Everything reported so far used the exact-ANF (monomial-prediction) transition
set while modelling key/constant addition as the identity. Those come from two
different frameworks: under monomial prediction a constant addition admits
u -> v with v >= u, and dropping those transitions can produce a FALSE UNSAT,
i.e. a false claim of balancedness.

The conventional BDP table has no such gap -- the existential "there exists
w >= k" in its definition is exactly what makes constant addition sound to skip.
It is a superset of the ANF table, so it yields fewer UNSAT results and weaker
(but trustworthy) conclusions.

If the 6/7 boundary and the multi-tweak degree bounds survive the swap, the
framework question is moot and every number stands on standard ground.

Run: python3 compare_tables.py
"""

import time

from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153

import blinkit.divprop.blink_divprop as dp
from blinkit.divprop.blink_divprop import NBITS
from blinkit.cipher.blink_ops import build_ops, nonlinear_prefix
import blinkit.divprop.blink_mt as blink_mt


def st_balanced(r, active):
    d = dp.Distinguisher(r)
    bal = d.balanced_bits(active)
    d.close()
    return bal


def mt_trail_exists(k, bound):
    """True (SAT) iff a division trail of total input weight >= bound exists."""
    seg = blink_mt.truncate_sboxes(blink_mt.ops_from_seam('h1'), k)
    m, inputs, out = blink_mt.build_mt_cnf(seg)
    cnf = list(m.cnf)
    cnf.append(list(out))
    for i in range(NBITS):
        for j in range(i + 1, NBITS):
            cnf.append([-out[i], -out[j]])
    card = CardEnc.atleast(lits=inputs, bound=bound, vpool=m.pool,
                           encoding=EncType.seqcounter)
    s = Cadical153(bootstrap_with=cnf + card.clauses)
    sat = s.solve()
    s.delete()
    return sat, len(inputs)


def mt_degree(k, hi):
    """Largest weight admitting a trail = deg(F). Binary search, fast queries."""
    lo = 0
    while lo < hi:
        mid = (lo + hi + 1) // 2
        sat, _ = mt_trail_exists(k, mid)
        if sat:
            lo = mid
        else:
            hi = mid - 1
    return lo


def main():
    for kind in ('anf', 'bdp'):
        n = dp.set_sbox_table(kind)
        label = ('exact ANF / monomial prediction' if kind == 'anf'
                 else 'conventional BDP')
        print(f"\n{'='*66}")
        print(f"{label}: {n}/256 valid transitions")
        print('=' * 66, flush=True)

        print("\nsingle-tweak calibration (63-bit plaintext cube)")
        active = list(range(1, NBITS))
        for r in range(4, 9):
            t0 = time.time()
            bal = st_balanced(r, active)
            print(f"  S layers {r}: {len(bal):2d} balanced bits "
                  f"({time.time()-t0:.1f}s)", flush=True)

        print("\nmulti-tweak degree bounds (3 free S layers from pi1)")
        for k in range(1, 5):
            t0 = time.time()
            _, nin = mt_trail_exists(k, 1)
            D = mt_degree(k, nin)
            print(f"  k={k}  covers {3+k} S layers  deg(F)={D:3d}  "
                  f"cube dim {D+1:3d} -> 2^{D+1} data  ({time.time()-t0:.1f}s)",
                  flush=True)


if __name__ == "__main__":
    main()
