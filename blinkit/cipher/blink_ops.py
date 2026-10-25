"""
Blink-64 operation sequence — single source of truth.

The same op list drives (a) a bit-exact evaluator checked against the ToSC test
vectors and (b) the division-property CNF in blink_divprop.py, so the SAT model
provably describes the real cipher rather than a simplified SPN.

Op tokens:
    'S'            S-box layer
    'M'            MixColumn (Midori involutory matrix, per column)
    'P' / 'Pb'     cell shuffle / its inverse
    ('AK', i)      XOR round key rk[i] and constant rc[i]
    ('AKh', i)     XOR M(rk[i]) and M(rc[i])   (backward rounds)
    ('X', name)    XOR a tweak/whitening value: w1, h1, h, h2, w2

Structure (Blink-64a, a=2, b=3), from THF Sec. 5.1:
    pi1 = M . S . R^a(. ^ w1)      R    = P . AC . AK . M . S
    pi2 = M . S . R^b . P          Rbar = S . MK . Pbar
    pi3 = Pbar . Rbar^b . S
    pi4 = Rbar^a . S . M(.) ^ w2
Round-key order is the Figure-2 literal reading recovered in blink_ref.py:
pi3 uses rk_1..rk_b and pi4 uses rk_{b+1}..rk_{a+b} (NOT the mirror order).
"""

from blinkit.cipher.blink_ref import (SBOX, P_64, P_64_INV, RC, RCP, NCELL, NCOL,
                       to_state, to_int, S_layer, M_layer, P_layer,
                       split_master, derive_hash_keys, toeplitz, key_orders,
                       set_variant, M_TV, T_TV, K_TV, C_TV)

# the configuration uniquely determined by the ToSC test vectors
CFG = (True, True, False, False, False, False, False, False, True, 1)


def build_ops(a=2, b=3, scheme=1):
    """Full Blink-64a op sequence, plaintext -> ciphertext."""
    f_a, f_b, b_b, b_a = key_orders(scheme)
    ops = [('X', 'w1')]

    # pi1 = M . S . R^a(. ^ w1)
    for i in f_a:
        ops += ['S', 'M', ('AK', i), 'P']
    ops += ['S', 'M', ('X', 'h1')]

    # pi2 = M . S . R^b . P
    ops += ['P']
    for i in f_b:
        ops += ['S', 'M', ('AK', i), 'P']
    ops += ['S', 'M', ('X', 'h')]

    # pi3 = Pbar . Rbar^b . S
    ops += ['S']
    for i in b_b:
        ops += ['Pb', 'M', ('AKh', i), 'S']
    ops += ['Pb', ('X', 'h2')]

    # pi4 = Rbar^a . S . M(.) ^ w2
    ops += ['M', 'S']
    for i in b_a:
        ops += ['Pb', 'M', ('AKh', i), 'S']
    ops += [('X', 'w2')]
    return ops


def nonlinear_prefix(ops, n_sbox):
    """Truncate ops after the n_sbox-th S layer."""
    out, seen = [], 0
    for op in ops:
        out.append(op)
        if op == 'S':
            seen += 1
            if seen == n_sbox:
                return out
    return out


def count_sboxes(ops):
    return sum(1 for op in ops if op == 'S')


# ------------------------------------------------------------------ evaluator

def evaluate(ops, m, t, k_int, cfg=CFG):
    """Bit-exact evaluation of an op list — used to validate against ToSC F.1."""
    lsb, out_row0, key_msb, sub_tail, k1_last, slice_rev, \
        hat_h1, hat_h2, hat_rk, _ = cfg

    w1, w2, rk_words = split_master(k_int)
    k1b, k2b = derive_hash_keys(k_int, key_msb, sub_tail, k1_last, slice_rev)
    h1 = toeplitz(t, k1b, out_row0)
    h2 = toeplitz(t, k2b, out_row0)

    rk = [to_state(x, lsb) for x in rk_words]
    rc = [to_state(x, lsb) for x in RC]
    rcp = [to_state(x, lsb) for x in RCP]
    vals = {'w1': to_state(w1, lsb), 'w2': to_state(w2, lsb),
            'h1': to_state(h1, lsb), 'h2': to_state(h2, lsb),
            'h': to_state(h1 ^ h2, lsb)}
    if hat_h1:
        vals['h1'] = M_layer(vals['h1'])
    if hat_h2:
        vals['h2'] = M_layer(vals['h2'])

    s = to_state(m, lsb)
    for op in ops:
        if op == 'S':
            s = S_layer(s)
        elif op == 'M':
            s = M_layer(s)
        elif op == 'P':
            s = P_layer(s, P_64)
        elif op == 'Pb':
            s = P_layer(s, P_64_INV)
        elif op[0] == 'AK':
            k, c = rk[op[1]], rc[op[1]]
            s = [x ^ kk ^ cc for x, kk, cc in zip(s, k, c)]
        elif op[0] == 'AKh':
            k, c = rk[op[1]], rcp[op[1]]
            if hat_rk:
                k, c = M_layer(k), M_layer(c)
            s = [x ^ kk ^ cc for x, kk, cc in zip(s, k, c)]
        elif op[0] == 'X':
            s = [x ^ y for x, y in zip(s, vals[op[1]])]
        else:
            raise ValueError(op)
    return to_int(s, lsb)


def self_test():
    set_variant(64, 64, 2, 3)
    ops = build_ops()
    got = evaluate(ops, M_TV, T_TV, K_TV)
    ok = got == C_TV
    print(f"op-list evaluator vs ToSC F.1: {'PASS' if ok else 'FAIL'} "
          f"({got:#018x})")
    print(f"total S layers in full Blink-64a: {count_sboxes(ops)}")
    return ok


if __name__ == "__main__":
    self_test()
