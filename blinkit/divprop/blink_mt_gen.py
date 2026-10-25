"""
Multi-tweak degree bounds for every Blink variant.

blink_mt.py is hard-wired to the 64-bit block; this generalises the same model
to (n, tau, a, b) so the 128-bit variants can be measured.

Model, unchanged in substance:
  S      conventional bit-based division property table (superset of exact ANF)
  M      COPY + XOR over each 4-bit slice of each column
  P/Pb   re-indexing
  AK/X   constants -> identity on trails
  seam injections inside the segment -> fresh input variables, XOR rule

deg(F) = largest total input weight admitting a trail to one output bit.
A tweak cube of dimension deg(F)+1 then zeroes the output, and with the
free-round extension (S is an involution) the distinguisher reaches one S layer
further still.

Usage:  python3 blink_mt_gen.py [variant ...]
"""

import itertools
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153

from blinkit.cipher.blink_gen import Blink, VARIANTS, SBOX
from blinkit.divprop.blink_divprop import sbox_anf_transitions, sbox_bdp_transitions

BDP = sbox_bdp_transitions(SBOX)
INVALID = [(u, v) for u in range(16) for v in range(16) if (u, v) not in BDP]


class Model:
    def __init__(self):
        self.pool = IDPool()
        self.cnf = []
        self._n = 0

    def fresh(self, k=1):
        out = []
        for _ in range(k):
            self._n += 1
            out.append(self.pool.id(('v', self._n)))
        return out

    def add(self, *cl):
        self.cnf.extend(cl)

    def sum_eq(self, total, parts):
        for p in parts:
            self.add([-p, total])
        self.add([-total] + list(parts))
        for x, y in itertools.combinations(parts, 2):
            self.add([-x, -y])


def build_cnf(B, ops):
    """CNF for the op segment. Returns (model, all inputs, outputs)."""
    m = Model()
    nb = B.n
    cur = m.fresh(nb)
    inputs = list(cur)

    for op in ops:
        if op == 'S':
            nxt = list(cur)
            for cell in range(B.ncell):
                inv = [cur[4 * cell + b] for b in range(4)]
                outv = m.fresh(4)
                for u, v in INVALID:
                    cl = []
                    for b in range(4):
                        cl.append(-inv[b] if (u >> b) & 1 else inv[b])
                    for b in range(4):
                        cl.append(-outv[b] if (v >> b) & 1 else outv[b])
                    m.add(cl)
                for b in range(4):
                    nxt[4 * cell + b] = outv[b]
            cur = nxt
        elif op == 'M':
            nxt = list(cur)
            c = B.ncol
            for col in range(c):
                cells = [col, col + c, col + 2 * c, col + 3 * c]
                for b in range(4):
                    a = [cur[4 * x + b] for x in cells]
                    out = m.fresh(4)
                    cp = {}
                    for i in range(4):
                        br = m.fresh(3)
                        cp[i] = dict(zip([j for j in range(4) if j != i], br))
                        m.sum_eq(a[i], br)
                    for j in range(4):
                        m.sum_eq(out[j], [cp[i][j] for i in range(4) if i != j])
                    for i, x in enumerate(cells):
                        nxt[4 * x + b] = out[i]
            cur = nxt
        elif op in ('P', 'Pb'):
            p = B.P if op == 'P' else B.Pinv
            nxt = [0] * nb
            for i in range(B.ncell):
                for b in range(4):
                    nxt[4 * i + b] = cur[4 * p[i] + b]
            cur = nxt
        elif isinstance(op, tuple) and op[0] == 'X' and op[1] in ('h', 'h2'):
            inj = m.fresh(nb)
            inputs += inj
            nxt = m.fresh(nb)
            for i in range(nb):
                m.sum_eq(nxt[i], [cur[i], inj[i]])
            cur = nxt
        # w1 / w2 / AK / AKh are constants -> identity
    return m, inputs, cur


def segment(B, k):
    """Ops from just after the h1 injection, truncated to k S layers."""
    ops = B.build_ops()
    i = ops.index(('X', 'h1'))
    rest = ops[i + 1:]
    out, seen = [], 0
    for op in rest:
        if op == 'S':
            if seen == k:
                break
            seen += 1
        out.append(op)
    while out and out[-1] != 'S':
        out.pop()
    return out


def trail_exists(B, k, bound):
    m, inputs, out = build_cnf(B, segment(B, k))
    cnf = list(m.cnf)
    cnf.append(list(out))
    for i in range(B.n):
        for j in range(i + 1, B.n):
            cnf.append([-out[i], -out[j]])
    card = CardEnc.atleast(lits=inputs, bound=bound, vpool=m.pool,
                           encoding=EncType.seqcounter)
    s = Cadical153(bootstrap_with=cnf + card.clauses)
    sat = s.solve()
    s.delete()
    return sat, len(inputs)


def degree(B, k, hi=None):
    _, nin = trail_exists(B, k, 1)
    lo, hi = 0, hi or nin
    while lo < hi:
        mid = (lo + hi + 1) // 2
        sat, _ = trail_exists(B, k, mid)
        if sat:
            lo = mid
        else:
            hi = mid - 1
    return lo


def sweep(name, kmax=5):
    n, tau, a, b, _, _, _ = VARIANTS[name]
    B = Blink(n, tau, a, b)
    total = B.n_sboxes()
    free = a + 1
    limit = 56 if n == 64 else 80
    print(f"\n{name}: {total} S layers, {free} free, tweak {tau} bits, "
          f"data limit 2^{limit}")
    print(f"  {'k':>2} {'S layers':>9} {'deg(F)':>7} {'cube':>6} "
          f"{'+free rnd':>10}  verdict")
    for k in range(1, kmax + 1):
        t0 = time.time()
        d = degree(B, k)
        cube = d + 1
        covered = free + k
        ok = cube <= limit
        print(f"  {k:2d} {covered:9d} {d:7d} {cube:6d} {covered+1:10d}  "
              f"{'2^%d — OK' % cube if ok else 'exceeds limit'}"
              f"   ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    names = sys.argv[1:] or ['Blink-64a', 'Blink-128a', 'Blink-128A']
    for nm in names:
        sweep(nm)
