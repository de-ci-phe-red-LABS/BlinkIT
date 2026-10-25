"""
Multi-tweak integral analysis of Blink-64 --- the case the designers assert but
never search ("strong resistance against ST / MT integral attacks", THF 7.4).

Setting
-------
Fix the plaintext and vary only the tweak t over a d-dimensional affine subspace.
Then pi1 is applied to a constant, so its 3 S layers contribute nothing: the cube
is injected straight into the pi2 input as c1 ^ h1(t).

The tweak reaches the state at three seams, as h1(t), h(t) = h1(t)^h2(t) and
h2(t) -- every one of them LINEAR in t (key-dependent Toeplitz matrices). So

    c = F( L1(t), L2(t), L3(t) ),   L_i linear

and deg_t(c) <= deg(F), the total degree of the nonlinear composition F in its
three joint 64-bit arguments. If d > deg(F) the cube sum vanishes.

This bound is KEY-INDEPENDENT and needs no rank assumption on T1, T2: a linear
substitution cannot raise degree, whatever the secret matrices are.

Degree is bounded with the validated division-property model: deg(F_j) <= the
largest total input weight w admitting a division trail to output bit e_j. A
cardinality constraint on the input variables turns that into one SAT query.
"""

import sys

from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153

from blinkit.cipher.blink_ref import P_64, P_64_INV
from blinkit.cipher.blink_ops import build_ops
from blinkit.divprop.blink_divprop import (Model, apply_sbox, apply_mixcolumn, apply_perm,
                           NBITS)

DATA_LIMIT = 56          # THF Sec. 5.5: at most 2^56 blocks per key for Blink-64


def ops_from_seam(seam='h1'):
    """Op list from just after the named tweak injection to the ciphertext."""
    ops = build_ops()
    idx = ops.index(('X', seam))
    return ops[idx + 1:]


def truncate_sboxes(ops, k):
    out, seen = [], 0
    for op in ops:
        if op == 'S':
            if seen == k:
                break
            seen += 1
        out.append(op)
    # drop trailing linear ops so the segment ends on an S layer
    while out and out[-1] != 'S':
        out.pop()
    return out


def build_mt_cnf(ops):
    """CNF for the segment; returns (model, all_input_vars, output_vars).

    Seam injections inside the segment contribute FRESH input variables joined
    by the XOR rule, because there the tweak re-enters as a variable, not a
    constant.
    """
    m = Model()
    cur = m.fresh(NBITS)
    inputs = list(cur)

    for op in ops:
        if op == 'S':
            cur = apply_sbox(m, cur)
        elif op == 'M':
            cur = apply_mixcolumn(m, cur)
        elif op == 'P':
            cur = apply_perm(cur, P_64)
        elif op == 'Pb':
            cur = apply_perm(cur, P_64_INV)
        elif op[0] == 'X' and op[1] in ('h', 'h2'):
            inj = m.fresh(NBITS)
            inputs += inj
            nxt = m.fresh(NBITS)
            for i in range(NBITS):
                m.sum_eq(nxt[i], [cur[i], inj[i]])
            cur = nxt
        # ('X','w1'/'w2') and AK/AKh are constants -> identity on trails
    return m, inputs, cur


def degree_bound(ops, verbose=False):
    """Largest total input weight admitting a trail to a single output bit.

    'Exactly one output bit set' is encoded in the CNF, so the solver picks the
    best output bit itself -- one query per weight instead of one per bit.
    """
    m, inputs, out = build_mt_cnf(ops)

    base = list(m.cnf)
    base.append(list(out))                       # at least one output bit
    for i in range(NBITS):                       # at most one output bit
        for j in range(i + 1, NBITS):
            base.append([-out[i], -out[j]])

    lo, hi = 0, len(inputs)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        card = CardEnc.atleast(lits=inputs, bound=mid, vpool=m.pool,
                               encoding=EncType.seqcounter)
        s = Cadical153(bootstrap_with=base + card.clauses)
        sat = s.solve()
        s.delete()
        if verbose:
            print(f"    weight >= {mid}: {'SAT' if sat else 'UNSAT'}")
        if sat:
            lo = mid
        else:
            hi = mid - 1
    return lo


def main():
    seam = sys.argv[1] if len(sys.argv) > 1 else 'h1'
    full = ops_from_seam(seam)
    n_s = sum(1 for o in full if o == 'S')
    print(f"multi-tweak analysis, cube injected at seam '{seam}'")
    print(f"segment to ciphertext has {n_s} S layers")
    print(f"data limit for Blink-64: 2^{DATA_LIMIT}\n")
    print(f"{'S layers':>9}  {'deg bound':>9}  {'cube dim needed':>15}  verdict")

    for k in range(1, n_s + 1):
        seg = truncate_sboxes(full, k)
        D = degree_bound(seg)
        need = D + 1
        ok = need <= DATA_LIMIT
        # 3 free S layers from pi1, which the fixed plaintext makes constant
        covered = 3 + k
        verdict = (f"distinguisher on {covered} S layers at 2^{need} tweaks"
                   if ok else "exceeds data limit")
        print(f"{k:>9}  {D:>9}  {need:>15}  {verdict}")


if __name__ == "__main__":
    main()
