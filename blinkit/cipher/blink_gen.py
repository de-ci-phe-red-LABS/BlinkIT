"""
Generalised Blink — all six variants from one implementation.

blink_ref.py is hard-wired to the 64-bit block. This module parameterises
(n, tau, a, b) so Blink-128a/b/A/B can be reached, and validates every variant
against its ToSC Appendix F test vector.

The bit-order conventions were recovered for Blink-64a by exhaustive search
(blink_ref.py). This module assumes they carry over to the 128-bit variants and
CHECKS that assumption against the published vectors -- if a variant fails, the
convention does not generalise and must be searched separately.

Variants (THF Table 2):

    name         n    tau   a  b   key bits   grid
    Blink-64a    64    64   2  3      448     4 x 4
    Blink-64b    64   128   2  3      448     4 x 4
    Blink-128a  128   128   3  3     1024     4 x 8
    Blink-128b  128   256   3  3     1024     4 x 8
    Blink-128A  128   128   3  5     1280     4 x 8
    Blink-128B  128   256   3  5     1280     4 x 8

Key bits = (a + b + 2) * n.
"""

SBOX = [1, 0, 9, 3, 8, 5, 0xe, 7, 4, 2, 0xc, 0xb, 0xa, 0xf, 6, 0xd]

P_64 = [0, 5, 11, 10, 1, 6, 4, 13, 2, 12, 9, 15, 3, 7, 14, 8]

P_128 = [5, 12, 4, 1, 17, 9, 10, 16, 28, 14, 21, 22, 11, 27, 8, 13,
         2, 25, 18, 3, 30, 6, 19, 20, 0, 23, 24, 31, 7, 15, 29, 26]

RC_64 = [0x13198a2e03707344, 0x082efa98ec4e6c89, 0xbe5466cf34e90c6c,
         0x3f84d5b5b5470917, 0xd1310ba698dfb5ac]
RCP_64 = [0x0d95748f728eb658, 0x7b54a41dc25a59b5, 0xc5d1b023286085f0,
          0x8e79dcb0603a180e, 0xd71577c1bd314b27]

RC_128 = [0x243f6a8885a308d313198a2e03707344,
          0xa4093822299f31d0082efa98ec4e6c89,
          0x452821e638d01377be5466cf34e90c6c,
          0xc0ac29b7c97c50dd3f84d5b5b5470917,
          0x9216d5d98979fb1bd1310ba698dfb5ac,
          0x2ffd72dbd01adfb7b8e1afed6a267e96,
          0xba7c9045f12c7f9924a19947b3916cf7,
          0x0801f2e2858efc16636920d871574e69]
RCP_128 = [0xa458fea3f4933d7e0d95748f728eb658,
           0x718bcd5882154aee7b54a41dc25a59b5,
           0x9c30d5392af26013c5d1b023286085f0,
           0xca417918b8db38ef8e79dcb0603a180e,
           0x6c9e0e8bb01e8a3ed71577c1bd314b27,
           0x78af2fda55605c60e65525f3aa55ab94,
           0x5748986263e8144055ca396a2aab10b6,
           0xb4cc5c341141e8cea15486af7c72e993]


class Blink:
    """One Blink variant."""

    def __init__(self, n, tau, a, b):
        self.n, self.tau, self.a, self.b = n, tau, a, b
        self.ncell = n // 4              # 16 or 32 nibbles
        self.ncol = n // 16              # 4 or 8 columns
        self.nwords = a + b + 2
        self.klen = self.nwords * n
        self.P = P_64 if n == 64 else P_128
        self.Pinv = [self.P.index(i) for i in range(self.ncell)]
        self.RC = (RC_64 if n == 64 else RC_128)[:a + b]
        self.RCP = (RCP_64 if n == 64 else RCP_128)[:a + b]
        self.half = n + tau - 1          # bits per Toeplitz key
        self.subset = 2 * n + 2 * tau - 2

    # ------------------------------------------------------------- state

    def to_state(self, v):
        return [(v >> (4 * i)) & 0xf for i in range(self.ncell)]

    def to_int(self, s):
        v = 0
        for i in range(self.ncell):
            v |= (s[i] & 0xf) << (4 * i)
        return v

    # ------------------------------------------------------------- layers

    def S(self, s):
        return [SBOX[x] for x in s]

    def M(self, s):
        c = self.ncol
        out = [0] * self.ncell
        for j in range(c):
            a0, a1, a2, a3 = s[j], s[j + c], s[j + 2 * c], s[j + 3 * c]
            out[j] = a1 ^ a2 ^ a3
            out[j + c] = a0 ^ a2 ^ a3
            out[j + 2 * c] = a0 ^ a1 ^ a3
            out[j + 3 * c] = a0 ^ a1 ^ a2
        return out

    def perm(self, s, p):
        return [s[p[i]] for i in range(self.ncell)]

    # ------------------------------------------------------- key schedule

    def split_key(self, k):
        """k = rk_{a+b} || ... || rk_1 || w2 || w1, written MSB-first."""
        mask = (1 << self.n) - 1
        w = [(k >> (self.n * i)) & mask for i in range(self.nwords)]
        return w[0], w[1], [w[2 + i] for i in range(self.a + self.b)]

    def hash_keys(self, k):
        """k'_i = k_{11i mod klen}; k2||k1 from the low 2n+2tau-2 bits of k'."""
        kb = [(k >> i) & 1 for i in range(self.klen)]
        kp = [kb[(11 * i) % self.klen] for i in range(self.klen)]
        sub = kp[:self.subset]
        return sub[:self.half], sub[self.half:]

    def toeplitz(self, t, kb):
        """h_T(t) = T.t, T[i][j] = kb[n-1-i+j]; t_0 and out bit 0 are LSBs."""
        tb = [(t >> j) & 1 for j in range(self.tau)]
        out = 0
        for i in range(self.n):
            v = 0
            base = self.n - 1 - i
            for j in range(self.tau):
                if tb[j] and kb[base + j]:
                    v ^= 1
            out |= v << i
        return out

    # ------------------------------------------------------------- ops

    def key_order(self):
        """Figure-2 literal reading: pi3 takes rk_1..rk_b, pi4 rk_{b+1}..rk_{a+b}."""
        a, b = self.a, self.b
        return (list(range(a)), list(range(a, a + b)),
                list(range(b)), list(range(b, a + b)))

    def build_ops(self):
        f_a, f_b, b_b, b_a = self.key_order()
        ops = [('X', 'w1')]
        for i in f_a:
            ops += ['S', 'M', ('AK', i), 'P']
        ops += ['S', 'M', ('X', 'h1')]
        ops += ['P']
        for i in f_b:
            ops += ['S', 'M', ('AK', i), 'P']
        ops += ['S', 'M', ('X', 'h')]
        ops += ['S']
        for i in b_b:
            ops += ['Pb', 'M', ('AKh', i), 'S']
        ops += ['Pb', ('X', 'h2')]
        ops += ['M', 'S']
        for i in b_a:
            ops += ['Pb', 'M', ('AKh', i), 'S']
        ops += [('X', 'w2')]
        return ops

    def n_sboxes(self, ops=None):
        return sum(1 for o in (ops or self.build_ops()) if o == 'S')

    def prefix(self, ops, n_sbox):
        out, seen = [], 0
        for op in ops:
            out.append(op)
            if op == 'S':
                seen += 1
                if seen == n_sbox:
                    return out
        return out

    # ------------------------------------------------------- evaluation

    def evaluate(self, ops, m, t, k, h1=None, h2=None):
        w1, w2, rkw = self.split_key(k)
        if h1 is None:
            k1b, k2b = self.hash_keys(k)
            h1 = self.toeplitz(t, k1b)
            h2 = self.toeplitz(t, k2b)
        rk = [self.to_state(x) for x in rkw]
        rc = [self.to_state(x) for x in self.RC]
        rcp = [self.to_state(x) for x in self.RCP]
        val = {'w1': self.to_state(w1), 'w2': self.to_state(w2),
               'h1': self.to_state(h1), 'h2': self.to_state(h2),
               'h': self.to_state(h1 ^ h2)}

        s = self.to_state(m)
        for op in ops:
            if op == 'S':
                s = self.S(s)
            elif op == 'M':
                s = self.M(s)
            elif op == 'P':
                s = self.perm(s, self.P)
            elif op == 'Pb':
                s = self.perm(s, self.Pinv)
            elif op[0] == 'AK':
                kk, cc = rk[op[1]], rc[op[1]]
                s = [x ^ y ^ z for x, y, z in zip(s, kk, cc)]
            elif op[0] == 'AKh':
                kk = self.M(rk[op[1]])
                cc = self.M(rcp[op[1]])
                s = [x ^ y ^ z for x, y, z in zip(s, kk, cc)]
            elif op[0] == 'X':
                s = [x ^ y for x, y in zip(s, val[op[1]])]
            else:
                raise ValueError(op)
        return self.to_int(s)


# --------------------------------------------------------------- vectors

K64 = int("d6a102d888a467e4d1d7dec33a246943e07c1dc6f302c57e762c2df9de6f0d21"
          "6dd387874a0b52ce3022e0ad78c78a0697779021b38e7fa1", 16)
K128 = int("d6a102d888a467e4d1d7dec33a246943e07c1dc6f302c57e762c2df9de6f0d21"
           "6dd387874a0b52ce3022e0ad78c78a0697779021b38e7fa15e2b66350517f80f"
           "2961c648d578bae174d70cb769c30a45cc40300fe8a342ca57a0bd0251ae39b6"
           "21b8f104904374bbd6a102e234a664e421b8f104904374bbd6a102d888a666e4", 16)
K128L = int("d6a102d888a467e4d1d7dec33a246943e07c1dc6f302c57e762c2df9de6f0d21"
            "6dd387874a0b52ce3022e0ad78c78a0697779021b38e7fa15e2b66350517f80f"
            "2961c648d578bae174d70cb769c30a45cc40300fe8a342ca57a0bd0251ae39b6"
            "21b8f104904374bbd6a102e234a664e421b8f104904374bbd6a102d888a666e4"
            "28962a4c96893eda752c17026a6395c2d6963be43b2fc10813d73f5a4a48d28d", 16)

T64 = 0x0123456789abcdef
T128 = 0x0123456789abcdef0123456789abcdef
T256 = 0x0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef

VARIANTS = {
    'Blink-64a':  (64, 64, 2, 3, K64, T64, 0xa4a0d10502be846e),
    'Blink-64b':  (64, 128, 2, 3, K64, T128, 0x743e142f17caaae1),
    'Blink-128a': (128, 128, 3, 3, K128, T128,
                   0xb722eef350bb182074a6ff13c967a593),
    'Blink-128b': (128, 256, 3, 3, K128, T256,
                   0x20705a38e00412165bdabcac1dcbdec2),
    'Blink-128A': (128, 128, 3, 5, K128L, T128,
                   0x82449f141c183601195b5046eac2b026),
    'Blink-128B': (128, 256, 3, 5, K128L, T256,
                   0x8dc41b223bc8cd9923b1297dd27583fc),
}


def make(name):
    n, tau, a, b, _, _, _ = VARIANTS[name]
    return Blink(n, tau, a, b)


def self_test():
    print(f"{'variant':12s} {'n':>4} {'tau':>4} {'a':>2} {'b':>2} "
          f"{'key':>5} {'S-lay':>6} {'free':>5}  result")
    ok = True
    for name, (n, tau, a, b, k, t, c) in VARIANTS.items():
        B = Blink(n, tau, a, b)
        assert k.bit_length() <= B.klen, f"{name}: key too long"
        got = B.evaluate(B.build_ops(), 0, t, k)
        good = got == c
        ok &= good
        free = a + 1          # pi1 = a rounds + S,M -> a+1 S layers
        print(f"{name:12s} {n:4d} {tau:4d} {a:2d} {b:2d} {B.klen:5d} "
              f"{B.n_sboxes():6d} {free:5d}  "
              f"{'PASS' if good else f'FAIL got {got:#034x}'}")
    return ok


if __name__ == "__main__":
    print("Blink — all variants vs ToSC Appendix F\n")
    print("PASS" if self_test() else "\nsome variants FAILED")
