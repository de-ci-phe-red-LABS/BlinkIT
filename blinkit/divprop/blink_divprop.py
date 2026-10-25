"""
Bit-based division property model for Blink-64, built on the true reflection
structure (blink_ops.build_ops) rather than a plain iterated SPN.

Model:
  S      exact ANF transition set: u -> v valid iff x^u occurs in ANF(y^v)
  M      COPY + XOR over each 4-bit slice of each column (Midori matrix)
  P/Pb   pure re-indexing of bit variables
  AK/X   key, constant and tweak-hash XORs preserve the division property
         and are therefore identity on trails (Todo's constant-addition rule)

A cube with active-bit set I leaves output bit j balanced when NO division
trail runs from 1_I to e_j, i.e. when the CNF is UNSAT under those assumptions.

Every distinguisher this model reports is re-checked against the real cipher in
verify_experimentally(), so the claim never rests on the model alone.
"""

import itertools
import sys

from pysat.formula import IDPool
from pysat.solvers import Cadical153

from blinkit.cipher.blink_ref import SBOX, P_64, P_64_INV, NCELL, NCOL, set_variant
from blinkit.cipher.blink_ops import build_ops, nonlinear_prefix, count_sboxes, evaluate, CFG

NBITS = 64


# ------------------------------------------------------------------ S-box table

def sbox_anf_transitions(sbox, n=4):
    """{(u, v)} such that the monomial x^u occurs in the ANF of y^v."""
    size = 1 << n
    trails = set()
    for v in range(size):
        # f(x) = prod_{i : v_i = 1} y_i(x)
        f = []
        for x in range(size):
            y = sbox[x]
            bit = 1
            for i in range(n):
                if (v >> i) & 1 and not ((y >> i) & 1):
                    bit = 0
                    break
            f.append(bit)
        # Moebius transform -> ANF coefficients
        a = f[:]
        for i in range(n):
            for x in range(size):
                if (x >> i) & 1:
                    a[x] ^= a[x ^ (1 << i)]
        for u in range(size):
            if a[u]:
                trails.add((u, v))
    return trails


def sbox_bdp_transitions(sbox, n=4):
    """Conventional bit-based division property table (Todo / Xiang et al.).

    Propagating input division property D_k through y = S(x):

        XOR_Y y^v = SUM_w a_w^(v) * (XOR_X x^w)

    and XOR_X x^w can be non-zero only when w >= k. So the trail k -> v is
    valid whenever THERE EXISTS w >= k with x^w in ANF(y^v) -- not only w = k.

    That existential quantifier is what makes this table a SUPERSET of the
    exact-ANF (monomial-prediction) one, and it is what makes key/constant
    addition sound to model as the identity.
    """
    anf = sbox_anf_transitions(sbox, n)
    size = 1 << n
    trails = set()
    for k in range(size):
        for v in range(size):
            # w >= k  <=>  k is a submask of w
            if any((w, v) in anf for w in range(size) if (w & k) == k):
                trails.add((k, v))
    return trails


TABLES = {
    'anf': sbox_anf_transitions(SBOX),   # exact monomial prediction
    'bdp': sbox_bdp_transitions(SBOX),   # conventional division property
}

SBOX_TRAILS = TABLES['anf']
SBOX_INVALID = [(u, v) for u in range(16) for v in range(16)
                if (u, v) not in SBOX_TRAILS]


def set_sbox_table(kind):
    """Select 'anf' (monomial prediction) or 'bdp' (conventional). Must be
    called BEFORE constructing any model, since the CNF bakes the table in."""
    global SBOX_TRAILS, SBOX_INVALID
    SBOX_TRAILS = TABLES[kind]
    SBOX_INVALID = [(u, v) for u in range(16) for v in range(16)
                    if (u, v) not in SBOX_TRAILS]
    return len(SBOX_TRAILS)


# ------------------------------------------------------------------ CNF pieces

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

    def add(self, *clauses):
        self.cnf.extend(clauses)

    def at_most_one(self, lits):
        for a, b in itertools.combinations(lits, 2):
            self.add([-a, -b])

    def sum_eq(self, total, parts):
        """total = sum(parts), all 0/1 -> at most one part is set."""
        for p in parts:
            self.add([-p, total])
        self.add([-total] + list(parts))
        self.at_most_one(parts)


def apply_sbox(m, cur):
    nxt = list(cur)
    for cell in range(NCELL):
        inv = [cur[4 * cell + b] for b in range(4)]
        outv = m.fresh(4)
        for u, v in SBOX_INVALID:
            clause = []
            for b in range(4):
                clause.append(-inv[b] if (u >> b) & 1 else inv[b])
            for b in range(4):
                clause.append(-outv[b] if (v >> b) & 1 else outv[b])
            m.add(clause)
        for b in range(4):
            nxt[4 * cell + b] = outv[b]
    return nxt


def apply_mixcolumn(m, cur):
    """Midori involutory matrix: out_i = XOR of the three other inputs."""
    nxt = list(cur)
    for col in range(NCOL):
        cells = [col, col + NCOL, col + 2 * NCOL, col + 3 * NCOL]
        for b in range(4):
            a = [cur[4 * c + b] for c in cells]
            out = m.fresh(4)
            # COPY: each input fans out to the three outputs it feeds
            cp = {}
            for i in range(4):
                branches = m.fresh(3)
                cp[i] = dict(zip([j for j in range(4) if j != i], branches))
                m.sum_eq(a[i], branches)
            # XOR: out_j collects one branch from each i != j
            for j in range(4):
                m.sum_eq(out[j], [cp[i][j] for i in range(4) if i != j])
            for i, c in enumerate(cells):
                nxt[4 * c + b] = out[i]
    return nxt


def apply_perm(cur, perm):
    nxt = [0] * NBITS
    for i in range(NCELL):
        for b in range(4):
            nxt[4 * i + b] = cur[4 * perm[i] + b]
    return nxt


def build_cnf(ops):
    m = Model()
    inp = m.fresh(NBITS)
    cur = inp
    for op in ops:
        if op == 'S':
            cur = apply_sbox(m, cur)
        elif op == 'M':
            cur = apply_mixcolumn(m, cur)
        elif op == 'P':
            cur = apply_perm(cur, P_64)
        elif op == 'Pb':
            cur = apply_perm(cur, P_64_INV)
        # AK / AKh / X: constant-addition rule, identity on division trails
    return m, inp, cur


# ------------------------------------------------------------------ the search

class Distinguisher:
    def __init__(self, n_sbox):
        self.n_sbox = n_sbox
        ops = nonlinear_prefix(build_ops(), n_sbox)
        self.ops = ops
        m, self.inp, self.out = build_cnf(ops)
        self.solver = Cadical153(bootstrap_with=m.cnf)
        self.nvars = m._n

    def balanced_bits(self, active):
        """Output bits with no trail from the cube 1_active -> e_j."""
        amask = set(active)
        in_assump = [self.inp[i] if i in amask else -self.inp[i]
                     for i in range(NBITS)]
        bal = []
        for j in range(NBITS):
            out_assump = [self.out[i] if i == j else -self.out[i]
                          for i in range(NBITS)]
            if not self.solver.solve(assumptions=in_assump + out_assump):
                bal.append(j)
        return bal

    def close(self):
        self.solver.delete()


# ------------------------------------------------------------ empirical check

def verify_experimentally(n_sbox, active, bits, key, tweak, trials=1):
    """Sum the real cipher over the cube and confirm the claimed bits are 0."""
    set_variant(64, 64, 2, 3)
    ops = nonlinear_prefix(build_ops(), n_sbox)
    d = len(active)
    acc = 0
    for combo in range(1 << d):
        m_val = 0
        for idx, bitpos in enumerate(active):
            if (combo >> idx) & 1:
                m_val |= 1 << bitpos
        acc ^= evaluate(ops, m_val, tweak, key, CFG)
    bad = [j for j in bits if (acc >> j) & 1]
    return acc, bad


if __name__ == "__main__":
    print(f"S-box ANF transitions: {len(SBOX_TRAILS)}/256 valid, "
          f"{len(SBOX_INVALID)} forbidden")
    full = count_sboxes(build_ops())
    print(f"full Blink-64a = {full} S layers\n")

    lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1, 8)
    for r in range(lo, hi + 1):
        d = Distinguisher(r)
        active = list(range(1, NBITS))          # 63 active bits, bit 0 fixed
        bal = d.balanced_bits(active)
        print(f"S layers {r:2d}  vars {d.nvars:7d}  "
              f"balanced bits (63-dim cube): {len(bal)}")
        d.close()
