"""
Blink-64a reference implementation + convention search.

The THF/Blink paper (ePrint 2025/1314, ToSC 2025) specifies Blink only up to a
number of bit-order conventions that Figure 2 and Sec. 5.3-5.4 leave implicit.
This script parameterises every such choice and brute-forces the combination
that reproduces the ToSC Appendix F.1 test vector.

Spec references:
  Sec. 5.1  pi1 = M . S . R^a(. ^ w1)
            pi2 = M . S . R^b . P
            pi3 = Pbar . Rbar^b . S
            pi4 = Rbar^a . S . M(.) ^ w2
            R    = P . AC . AK . M . S
            Rbar = S . MK_{rk^ ^ rc^} . Pbar,  MK_k(x) = M(x) ^ k,  z^ = M(z)
  Sec. 5.2  S involutive; M = Midori involutory matrix; P as below
  Sec. 5.3  T[i][j] = k[n-1-i+j];  h(t) = h1(t) ^ h2(t);  t_0 = LSB of t
  Sec. 5.4  k = rk_{a+b} || ... || rk_1 || w2 || w1
            k'_i = k_{11i mod (a+b+2)n}
            k2||k1 = least significant 2n + 2tau - 2 bits of k'
"""

import itertools

SBOX = [1, 0, 9, 3, 8, 5, 0xe, 7, 4, 2, 0xc, 0xb, 0xa, 0xf, 6, 0xd]
P_64 = [0, 5, 11, 10, 1, 6, 4, 13, 2, 12, 9, 15, 3, 7, 14, 8]
P_64_INV = [P_64.index(i) for i in range(16)]

RC = [0x13198a2e03707344, 0x082efa98ec4e6c89, 0xbe5466cf34e90c6c,
      0x3f84d5b5b5470917, 0xd1310ba698dfb5ac]
RCP = [0x0d95748f728eb658, 0x7b54a41dc25a59b5, 0xc5d1b023286085f0,
       0x8e79dcb0603a180e, 0xd71577c1bd314b27]

N, TAU, A, B = 64, 64, 2, 3
NCELL = N // 4          # 16 nibbles
NCOL = N // 16          # 4 columns
KLEN = (A + B + 2) * N  # 448


def set_variant(n, tau, a, b):
    """Blink-64a = (64,64,2,3); Blink-64b = (64,128,2,3)."""
    global N, TAU, A, B, NCELL, NCOL, KLEN
    N, TAU, A, B = n, tau, a, b
    NCELL, NCOL, KLEN = n // 4, n // 16, (a + b + 2) * n


# ---------------------------------------------------------------- state <-> int

def to_state(v, lsb_first):
    if lsb_first:
        return [(v >> (4 * i)) & 0xf for i in range(NCELL)]
    return [(v >> (4 * (NCELL - 1 - i))) & 0xf for i in range(NCELL)]


def to_int(s, lsb_first):
    v = 0
    for i in range(NCELL):
        sh = 4 * i if lsb_first else 4 * (NCELL - 1 - i)
        v |= (s[i] & 0xf) << sh
    return v


# ---------------------------------------------------------------- round layers

def S_layer(s):
    return [SBOX[x] for x in s]


def M_layer(s):
    out = [0] * NCELL
    for j in range(NCOL):
        a0, a1, a2, a3 = s[j], s[j + NCOL], s[j + 2 * NCOL], s[j + 3 * NCOL]
        out[j] = a1 ^ a2 ^ a3
        out[j + NCOL] = a0 ^ a2 ^ a3
        out[j + 2 * NCOL] = a0 ^ a1 ^ a3
        out[j + 3 * NCOL] = a0 ^ a1 ^ a2
    return out


def P_layer(s, perm):
    return [s[perm[i]] for i in range(NCELL)]


def R_fwd(s, rk, rc):
    """R = P . AC . AK . M . S"""
    s = M_layer(S_layer(s))
    s = [x ^ k ^ c for x, k, c in zip(s, rk, rc)]
    return P_layer(s, P_64)


def R_bwd(s, rk_hat, rc_hat):
    """Rbar = S . MK_{rk^ ^ rc^} . Pbar,  MK_k(x) = M(x) ^ k"""
    s = M_layer(P_layer(s, P_64_INV))
    s = [x ^ k ^ c for x, k, c in zip(s, rk_hat, rc_hat)]
    return S_layer(s)


# ---------------------------------------------------------------- toeplitz hash

def toeplitz(t, kbits, out_lsb_row0):
    """h_T(t) = T.t with T[i][j] = kbits[n-1-i+j]; t_0 is the LSB of t."""
    tb = [(t >> j) & 1 for j in range(TAU)]
    out = 0
    for i in range(N):
        v = 0
        base = N - 1 - i
        for j in range(TAU):
            if tb[j] and kbits[base + j]:
                v ^= 1
        out |= v << (i if out_lsb_row0 else N - 1 - i)
    return out


# ---------------------------------------------------------------- key schedule

def split_master(k_int):
    """k = rk5||rk4||rk3||rk2||rk1||w2||w1, written MSB-first."""
    words = [(k_int >> (64 * i)) & ((1 << 64) - 1) for i in range(7)]
    # words[0] is the rightmost 64 bits = w1
    w1, w2 = words[0], words[1]
    rk = [words[2 + i] for i in range(A + B)]  # rk[0]=rk_1 ... rk[4]=rk_5
    return w1, w2, rk


def derive_hash_keys(k_int, key_idx_msb, subset_tail, k1_last, slice_rev):
    """k'_i = k_{11i mod 448}; k2||k1 = least significant 2n+2tau-2 bits of k'."""
    if key_idx_msb:
        kb = [(k_int >> (KLEN - 1 - i)) & 1 for i in range(KLEN)]
    else:
        kb = [(k_int >> i) & 1 for i in range(KLEN)]

    kp = [kb[(11 * i) % KLEN] for i in range(KLEN)]

    m = 2 * N + 2 * TAU - 2  # 254
    subset = kp[-m:] if subset_tail else kp[:m]

    half = N + TAU - 1  # 127
    if k1_last:
        k2b, k1b = subset[:half], subset[half:]
    else:
        k1b, k2b = subset[:half], subset[half:]

    if slice_rev:
        k1b, k2b = k1b[::-1], k2b[::-1]
    return k1b, k2b


# ---------------------------------------------------------------- key orderings

def key_orders(scheme):
    """(fwd_a, fwd_b, bwd_b, bwd_a) as index lists into rk[0..4] / rc[0..4].

    Figure 2 renders the bottom-row labels ambiguously, so the mirror reading
    and the literal left-to-right reading are both tried.
    """
    f_a = [0, 1]           # pi1: rk_1, rk_2
    f_b = [2, 3, 4]        # pi2: rk_3, rk_4, rk_5
    if scheme == 0:        # mirror: pi3 = rk_5..rk_3, pi4 = rk_2..rk_1
        return f_a, f_b, [4, 3, 2], [1, 0]
    if scheme == 1:        # figure-literal: pi3 = rk_1..rk_3, pi4 = rk_4..rk_5
        return f_a, f_b, [0, 1, 2], [3, 4]
    if scheme == 2:        # pi3 = rk_3..rk_1, pi4 = rk_5..rk_4
        return f_a, f_b, [2, 1, 0], [4, 3]
    return f_a, f_b, [3, 4, 0], [1, 2]  # scheme 3: shifted


# ---------------------------------------------------------------- encryption

def encrypt(m, t, k_int, cfg):
    lsb, out_row0, key_msb, sub_tail, k1_last, slice_rev, \
        hat_h1, hat_h2, hat_rk, scheme = cfg

    w1, w2, rk_words = split_master(k_int)
    k1b, k2b = derive_hash_keys(k_int, key_msb, sub_tail, k1_last, slice_rev)

    h1 = toeplitz(t, k1b, out_row0)
    h2 = toeplitz(t, k2b, out_row0)
    h = h1 ^ h2

    rk = [to_state(x, lsb) for x in rk_words]
    rc = [to_state(x, lsb) for x in RC]
    rcp = [to_state(x, lsb) for x in RCP]

    f_a, f_b, b_b, b_a = key_orders(scheme)

    def hat(s):
        return M_layer(s) if hat_rk else s

    s = to_state(m ^ w1, lsb)

    # pi1 = M . S . R^a(. ^ w1)
    for i in f_a:
        s = R_fwd(s, rk[i], rc[i])
    s = M_layer(S_layer(s))

    # ^ h1(t)
    h1s = to_state(h1, lsb)
    s = [x ^ y for x, y in zip(s, M_layer(h1s) if hat_h1 else h1s)]

    # pi2 = M . S . R^b . P
    s = P_layer(s, P_64)
    for i in f_b:
        s = R_fwd(s, rk[i], rc[i])
    s = M_layer(S_layer(s))

    # ^ h(t)
    s = [x ^ y for x, y in zip(s, to_state(h, lsb))]

    # pi3 = Pbar . Rbar^b . S
    s = S_layer(s)
    for i in b_b:
        s = R_bwd(s, hat(rk[i]), hat(rcp[i]))
    s = P_layer(s, P_64_INV)

    # ^ h2^(t)
    h2s = to_state(h2, lsb)
    s = [x ^ y for x, y in zip(s, M_layer(h2s) if hat_h2 else h2s)]

    # pi4 = Rbar^a . S . M(.) ^ w2
    s = S_layer(M_layer(s))
    for i in b_a:
        s = R_bwd(s, hat(rk[i]), hat(rcp[i]))

    return to_int(s, lsb) ^ w2


# ---------------------------------------------------------------- search

FIELDS = ["lsb_first", "toeplitz_row0_is_lsb", "key_idx_msb", "subset_tail",
          "k1_is_last", "slice_reversed", "hat_h1", "hat_h2", "hat_rk",
          "key_order_scheme"]

# ToSC Appendix F.1
M_TV = 0x0
T_TV = 0x0123456789abcdef
K_TV = int("d6a102d888a467e4d1d7dec33a246943e07c1dc6f302c57e762c2df9de6f0d21"
           "6dd387874a0b52ce3022e0ad78c78a0697779021b38e7fa1", 16)
C_TV = 0xa4a0d10502be846e


def search():
    bools = [(False, True)] * 9
    space = itertools.product(*bools, range(4))
    hits, tried = [], 0
    for cfg in space:
        tried += 1
        try:
            if encrypt(M_TV, T_TV, K_TV, cfg) == C_TV:
                hits.append(cfg)
                print("MATCH:", dict(zip(FIELDS, cfg)))
        except Exception:
            continue
    print(f"\ntried {tried} configurations, {len(hits)} match(es)")
    return hits


# Blink-64b: same master key, tau = 128, different tweak (ToSC App. F.2)
T_TV_B = 0x0123456789abcdef0123456789abcdef
C_TV_B = 0x743e142f17caaae1


def verify(cfg):
    """Confirm cfg on 64a, then cross-check the SAME cfg on 64b (tau=128)."""
    set_variant(64, 64, 2, 3)
    ok_a = encrypt(M_TV, T_TV, K_TV, cfg) == C_TV
    set_variant(64, 128, 2, 3)
    got_b = encrypt(M_TV, T_TV_B, K_TV, cfg)
    ok_b = got_b == C_TV_B
    set_variant(64, 64, 2, 3)
    print(f"Blink-64a: {'PASS' if ok_a else 'FAIL'}")
    print(f"Blink-64b: {'PASS' if ok_b else 'FAIL'}  got {got_b:#018x} "
          f"want {C_TV_B:#018x}")
    return ok_a and ok_b


if __name__ == "__main__":
    assert [SBOX[SBOX[x]] for x in range(16)] == list(range(16)), "S not involutive"
    hits = search()
    for cfg in hits:
        print()
        verify(cfg)
