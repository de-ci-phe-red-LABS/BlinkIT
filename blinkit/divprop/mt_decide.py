"""
The decisive multi-tweak query, asked directly.

We do not need deg(F) exactly -- only whether it clears the Blink-64 data limit.
A tweak cube of dimension d sums to zero when d > deg(F), and the spec allows
d <= 56. So the single question that matters at k S layers is:

    is there a division trail of total input weight >= BOUND ?

    UNSAT -> deg(F) < BOUND -> a 2^BOUND tweak cube is a distinguisher
    SAT   -> no distinguisher at this weight

One SAT call per (k, BOUND) instead of a binary search.

Usage: python3 mt_decide.py <k> [bound]
"""

import sys
import time

from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153

from blinkit.divprop.blink_mt import ops_from_seam, truncate_sboxes, build_mt_cnf, DATA_LIMIT
from blinkit.divprop.blink_divprop import NBITS

PI1_SBOXES = 3


def decide(k, bound):
    seg = truncate_sboxes(ops_from_seam('h1'), k)
    m, inputs, out = build_mt_cnf(seg)

    cnf = list(m.cnf)
    cnf.append(list(out))
    for i in range(NBITS):
        for j in range(i + 1, NBITS):
            cnf.append([-out[i], -out[j]])

    card = CardEnc.atleast(lits=inputs, bound=bound, vpool=m.pool,
                           encoding=EncType.seqcounter)
    cnf += card.clauses

    print(f"k = {k}: {len(inputs)} input vars, {len(cnf)} clauses, "
          f"asking for a trail of weight >= {bound}", flush=True)

    s = Cadical153(bootstrap_with=cnf)
    t0 = time.time()
    sat = s.solve()
    dt = time.time() - t0
    s.delete()

    print(f"  -> {'SAT' if sat else 'UNSAT'}  ({dt:.1f}s)", flush=True)
    if not sat:
        print(f"  deg(F) < {bound}; a dimension-{bound} tweak cube zeroes "
              f"{PI1_SBOXES + k} S layers at 2^{bound} data "
              f"(limit 2^{DATA_LIMIT})", flush=True)
    else:
        print(f"  a weight-{bound} trail exists; no distinguisher at this "
              f"cube dimension", flush=True)
    return sat


if __name__ == "__main__":
    k = int(sys.argv[1])
    bound = int(sys.argv[2]) if len(sys.argv) > 2 else DATA_LIMIT
    decide(k, bound)
