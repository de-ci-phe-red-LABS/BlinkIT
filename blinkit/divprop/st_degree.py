"""
Control experiment: the single-tweak degree bound.

The multi-tweak claim is only interesting if the single-tweak setting genuinely
needs a big cube. Establishing "63 active bits works" is not enough -- a referee
will ask whether a much smaller plaintext cube also reaches 6 S layers.

Here the plaintext is the cube and the tweak is fixed, so every tweak injection
is a constant. The degree bound over the first r S layers gives the cube
dimension that is guaranteed to work: any cube of dimension > deg works.

Usage: python3 st_degree.py <r> <bound>
"""

import sys
import time

from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153

from blinkit.divprop.blink_divprop import build_cnf, NBITS
from blinkit.cipher.blink_ops import build_ops, nonlinear_prefix


def decide(r, bound):
    ops = nonlinear_prefix(build_ops(), r)
    m, inp, out = build_cnf(ops)

    cnf = list(m.cnf)
    cnf.append(list(out))
    for i in range(NBITS):
        for j in range(i + 1, NBITS):
            cnf.append([-out[i], -out[j]])
    card = CardEnc.atleast(lits=inp, bound=bound, vpool=m.pool,
                           encoding=EncType.seqcounter)
    cnf += card.clauses

    s = Cadical153(bootstrap_with=cnf)
    t0 = time.time()
    sat = s.solve()
    dt = time.time() - t0
    s.delete()
    print(f"  ST r={r} weight>={bound}: {'SAT' if sat else 'UNSAT'} ({dt:.1f}s)",
          flush=True)
    return sat


if __name__ == "__main__":
    r = int(sys.argv[1])
    lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 \
        else (int(sys.argv[2]), int(sys.argv[2]))
    for b in range(lo, hi + 1):
        decide(r, b)
