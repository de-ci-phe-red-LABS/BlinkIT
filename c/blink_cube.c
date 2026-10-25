/*
 * blink_cube.c -- multi-tweak cube evaluator for Blink-64a.
 *
 * Purpose: settle experimentally what the division-property model only predicts.
 * The model says a tweak cube of dimension deg(F)+1 zeroes the output; k=1 and
 * k=2 were checkable in Python, but k=3 needs 2^28 encryptions.
 *
 * Two things make that cheap:
 *
 *   1. h1 and h2 are LINEAR in t, so walking the cube in Gray-code order lets
 *      each step update them with a single XOR -- no Toeplitz hash per encryption.
 *   2. With nibbles stored LSB-first, the four MixColumn rows land in the four
 *      16-bit chunks of the state word, so the whole Midori matrix is three XORs.
 *
 * The op sequence is interpreted from a table that mirrors blink_ops.py exactly,
 * rather than hand-inlined, so this provably runs the same cipher.
 *
 * Build: make -f Makefile.cube      (or: gcc -O3 -march=native -fopenmp)
 * Run:   ./blink_cube               (self-test, then k=1,2,3)
 *        ./blink_cube <k> <dim>     (single case)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* ------------------------------------------------------------------ cipher */

static const uint8_t SBOX[16] = {1,0,9,3,8,5,0xe,7,4,2,0xc,0xb,0xa,0xf,6,0xd};
static const uint8_t P64[16]  = {0,5,11,10,1,6,4,13,2,12,9,15,3,7,14,8};
static uint8_t P64I[16];

static const uint64_t RC[5] = {
    0x13198a2e03707344ULL, 0x082efa98ec4e6c89ULL, 0xbe5466cf34e90c6cULL,
    0x3f84d5b5b5470917ULL, 0xd1310ba698dfb5acULL};
static const uint64_t RCP[5] = {
    0x0d95748f728eb658ULL, 0x7b54a41dc25a59b5ULL, 0xc5d1b023286085f0ULL,
    0x8e79dcb0603a180eULL, 0xd71577c1bd314b27ULL};

static uint16_t SBOX16[65536];   /* S applied to four nibbles at once */

static void init_tables(void) {
    for (int i = 0; i < 16; i++) P64I[P64[i]] = (uint8_t)i;
    for (uint32_t x = 0; x < 65536; x++) {
        uint16_t y = 0;
        for (int n = 0; n < 4; n++) y |= (uint16_t)SBOX[(x >> (4*n)) & 0xf] << (4*n);
        SBOX16[x] = y;
    }
}

static inline uint64_t S_layer(uint64_t s) {
    return (uint64_t)SBOX16[s & 0xffff]
         | ((uint64_t)SBOX16[(s >> 16) & 0xffff] << 16)
         | ((uint64_t)SBOX16[(s >> 32) & 0xffff] << 32)
         | ((uint64_t)SBOX16[(s >> 48) & 0xffff] << 48);
}

/* Midori involutory matrix. Column j is the j-th nibble of each 16-bit chunk,
 * so out_r = XOR of the three other rows, applied chunk-wise. */
static inline uint64_t M_layer(uint64_t s) {
    uint64_t r0 = s & 0xffffULL, r1 = (s >> 16) & 0xffffULL;
    uint64_t r2 = (s >> 32) & 0xffffULL, r3 = (s >> 48) & 0xffffULL;
    uint64_t n0 = r1 ^ r2 ^ r3, n1 = r0 ^ r2 ^ r3;
    uint64_t n2 = r0 ^ r1 ^ r3, n3 = r0 ^ r1 ^ r2;
    return n0 | (n1 << 16) | (n2 << 32) | (n3 << 48);
}

static inline uint64_t perm(uint64_t s, const uint8_t *p) {
    uint64_t o = 0;
    for (int i = 0; i < 16; i++) o |= ((s >> (4*p[i])) & 0xfULL) << (4*i);
    return o;
}

/* ------------------------------------------------------------- op sequence */

enum { OP_S, OP_M, OP_P, OP_PB, OP_AK, OP_AKH, OP_XW1, OP_XW2, OP_XH1, OP_XH, OP_XH2 };

typedef struct { uint8_t op, idx; } Op;

static Op OPS[64];
static int NOPS;

/* Mirrors blink_ops.build_ops() with key_order_scheme = 1:
 * pi3 uses rk_1..rk_3 and pi4 uses rk_4..rk_5 -- the literal Figure 2 reading,
 * NOT the mirrored order the reflection structure suggests. */
static void build_ops(void) {
    int n = 0;
    const int f_a[2] = {0,1}, f_b[3] = {2,3,4}, b_b[3] = {0,1,2}, b_a[2] = {3,4};
    OPS[n++] = (Op){OP_XW1,0};
    for (int i = 0; i < 2; i++) {
        OPS[n++] = (Op){OP_S,0}; OPS[n++] = (Op){OP_M,0};
        OPS[n++] = (Op){OP_AK,(uint8_t)f_a[i]}; OPS[n++] = (Op){OP_P,0};
    }
    OPS[n++] = (Op){OP_S,0}; OPS[n++] = (Op){OP_M,0}; OPS[n++] = (Op){OP_XH1,0};
    OPS[n++] = (Op){OP_P,0};
    for (int i = 0; i < 3; i++) {
        OPS[n++] = (Op){OP_S,0}; OPS[n++] = (Op){OP_M,0};
        OPS[n++] = (Op){OP_AK,(uint8_t)f_b[i]}; OPS[n++] = (Op){OP_P,0};
    }
    OPS[n++] = (Op){OP_S,0}; OPS[n++] = (Op){OP_M,0}; OPS[n++] = (Op){OP_XH,0};
    OPS[n++] = (Op){OP_S,0};
    for (int i = 0; i < 3; i++) {
        OPS[n++] = (Op){OP_PB,0}; OPS[n++] = (Op){OP_M,0};
        OPS[n++] = (Op){OP_AKH,(uint8_t)b_b[i]}; OPS[n++] = (Op){OP_S,0};
    }
    OPS[n++] = (Op){OP_PB,0}; OPS[n++] = (Op){OP_XH2,0};
    OPS[n++] = (Op){OP_M,0}; OPS[n++] = (Op){OP_S,0};
    for (int i = 0; i < 2; i++) {
        OPS[n++] = (Op){OP_PB,0}; OPS[n++] = (Op){OP_M,0};
        OPS[n++] = (Op){OP_AKH,(uint8_t)b_a[i]}; OPS[n++] = (Op){OP_S,0};
    }
    OPS[n++] = (Op){OP_XW2,0};
    NOPS = n;
}

/* number of ops to run so that exactly `nsbox` S layers are applied */
static int prefix_len(int nsbox) {
    int seen = 0;
    for (int i = 0; i < NOPS; i++) {
        if (OPS[i].op == OP_S && ++seen == nsbox) return i + 1;
    }
    return NOPS;
}

/* ---------------------------------------------------------------- key data */

typedef struct {
    uint64_t w1, w2;
    uint64_t rk[5], Mrk[5], Mrcp[5];
    uint8_t k1b[127], k2b[127];
} Key;

static void key_expand(const uint64_t kw[7], Key *K) {
    K->w1 = kw[0];
    K->w2 = kw[1];
    for (int i = 0; i < 5; i++) {
        K->rk[i]   = kw[2+i];
        K->Mrk[i]  = M_layer(kw[2+i]);
        K->Mrcp[i] = M_layer(RCP[i]);
    }
    /* k'_i = k_{11i mod 448}, k_0 the LSB; hash keys from the FIRST 254 bits,
     * k1 in the low half. */
    uint8_t kp[448];
    for (int i = 0; i < 448; i++) {
        int b = (11 * i) % 448;
        kp[i] = (uint8_t)((kw[b >> 6] >> (b & 63)) & 1ULL);
    }
    memcpy(K->k1b, kp,       127);
    memcpy(K->k2b, kp + 127, 127);
}

/* h_T(t) = T.t with T[i][j] = kbits[63-i+j]; output row 0 is the LSB */
static uint64_t toeplitz(uint64_t t, const uint8_t *kb) {
    uint64_t out = 0;
    for (int i = 0; i < 64; i++) {
        uint64_t v = 0;
        for (int j = 0; j < 64; j++)
            if ((t >> j) & 1ULL) v ^= kb[63 - i + j];
        out |= (v & 1ULL) << i;
    }
    return out;
}

/* ------------------------------------------------------------- encryption */

static inline uint64_t encrypt_pref(int plen, uint64_t m, uint64_t h1,
                                    uint64_t h2, const Key *K) {
    uint64_t s = m;
    for (int i = 0; i < plen; i++) {
        switch (OPS[i].op) {
        case OP_S:   s = S_layer(s); break;
        case OP_M:   s = M_layer(s); break;
        case OP_P:   s = perm(s, P64); break;
        case OP_PB:  s = perm(s, P64I); break;
        case OP_AK:  s ^= K->rk[OPS[i].idx] ^ RC[OPS[i].idx]; break;
        case OP_AKH: s ^= K->Mrk[OPS[i].idx] ^ K->Mrcp[OPS[i].idx]; break;
        case OP_XW1: s ^= K->w1; break;
        case OP_XW2: s ^= K->w2; break;
        case OP_XH1: s ^= h1; break;
        case OP_XH:  s ^= h1 ^ h2; break;
        case OP_XH2: s ^= h2; break;
        }
    }
    return s;
}

/* ------------------------------------------------------------- self-test */

static const char *KEY_HEX =
    "d6a102d888a467e4d1d7dec33a246943e07c1dc6f302c57e762c2df9de6f0d21"
    "6dd387874a0b52ce3022e0ad78c78a0697779021b38e7fa1";

static void parse_key(const char *hex, uint64_t kw[7]) {
    size_t len = strlen(hex);
    char buf[17];
    for (int w = 0; w < 7; w++) {                 /* w=0 is the LOW 64 bits */
        size_t end = len - (size_t)w * 16;
        memcpy(buf, hex + end - 16, 16);
        buf[16] = 0;
        kw[w] = strtoull(buf, NULL, 16);
    }
}

static int self_test(void) {
    uint64_t kw[7]; Key K;
    parse_key(KEY_HEX, kw);
    key_expand(kw, &K);
    uint64_t t = 0x0123456789abcdefULL;
    uint64_t h1 = toeplitz(t, K.k1b), h2 = toeplitz(t, K.k2b);
    uint64_t c = encrypt_pref(NOPS, 0, h1, h2, &K);
    uint64_t want = 0xa4a0d10502be846eULL;
    printf("ToSC F.1 test vector: %016llx %s\n",
           (unsigned long long)c, c == want ? "PASS" : "FAIL");
    return c == want;
}

/* ------------------------------------------------------------- cube sum */

/* XOR of the truncated cipher over a d-dimensional affine tweak subspace.
 * Uses the linearity of h1,h2: each Gray-code step is one XOR per hash. */
static uint64_t cube_sum(int nsbox, int d, const int *bits, uint64_t m,
                         uint64_t base, const Key *K) {
    int plen = prefix_len(nsbox);
    uint64_t d1[64], d2[64];
    for (int j = 0; j < d; j++) {
        uint64_t e = 1ULL << bits[j];
        d1[j] = toeplitz(e, K->k1b);
        d2[j] = toeplitz(e, K->k2b);
    }
    uint64_t h1_0 = toeplitz(base, K->k1b), h2_0 = toeplitz(base, K->k2b);

    int dlow = d < 16 ? d : 16;
    int dhigh = d - dlow;
    uint64_t acc = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(^:acc) schedule(static)
#endif
    for (long hi = 0; hi < (1L << dhigh); hi++) {
        uint64_t h1 = h1_0, h2 = h2_0;
        for (int j = 0; j < dhigh; j++)
            if ((hi >> j) & 1L) { h1 ^= d1[dlow + j]; h2 ^= d2[dlow + j]; }
        uint64_t local = encrypt_pref(plen, m, h1, h2, K);
        for (uint64_t g = 1; g < (1ULL << dlow); g++) {
            int j = __builtin_ctzll(g);
            h1 ^= d1[j]; h2 ^= d2[j];
            local ^= encrypt_pref(plen, m, h1, h2, K);
        }
        acc ^= local;
    }
    return acc;
}

/* Same cube walk, but S is applied to each output before XORing. Verifies the
 * free-round extension one S layer past the distinguisher. */
static uint64_t cube_sum_sapplied(int nsbox, int d, const int *bits, uint64_t m,
                                  uint64_t base, const Key *K) {
    int plen = prefix_len(nsbox);
    uint64_t d1[64], d2[64];
    for (int j = 0; j < d; j++) {
        uint64_t e = 1ULL << bits[j];
        d1[j] = toeplitz(e, K->k1b);
        d2[j] = toeplitz(e, K->k2b);
    }
    uint64_t h1_0 = toeplitz(base, K->k1b), h2_0 = toeplitz(base, K->k2b);
    int dlow = d < 16 ? d : 16;
    int dhigh = d - dlow;
    uint64_t acc = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(^:acc) schedule(static)
#endif
    for (long hi = 0; hi < (1L << dhigh); hi++) {
        uint64_t h1 = h1_0, h2 = h2_0;
        for (int j = 0; j < dhigh; j++)
            if ((hi >> j) & 1L) { h1 ^= d1[dlow + j]; h2 ^= d2[dlow + j]; }
        uint64_t local = S_layer(encrypt_pref(plen, m, h1, h2, K));
        for (uint64_t g = 1; g < (1ULL << dlow); g++) {
            int j = __builtin_ctzll(g);
            h1 ^= d1[j]; h2 ^= d2[j];
            local ^= S_layer(encrypt_pref(plen, m, h1, h2, K));
        }
        acc ^= local;
    }
    return acc;
}

/* --------------------------------------------------- single-tweak cube sum */

/* The control experiment. The claim "single-tweak needs a cube of dimension
 * >= 62" came from a SAT answer, and SAT proves nothing: a trail may exist while
 * the monomial cancels. Measuring the ST threshold directly, at the same S-layer
 * counts as the multi-tweak runs, puts the comparison on solid ground.
 *
 * Here the tweak is fixed and the PLAINTEXT is the cube, so h1 and h2 are
 * constants and only m varies. */
static uint64_t st_cube_sum(int nsbox, int d, const int *bits, uint64_t base_m,
                            uint64_t h1, uint64_t h2, const Key *K) {
    int plen = prefix_len(nsbox);
    int dlow = d < 16 ? d : 16;
    int dhigh = d - dlow;
    uint64_t acc = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(^:acc) schedule(static)
#endif
    for (long hi = 0; hi < (1L << dhigh); hi++) {
        uint64_t m = base_m;
        for (int j = 0; j < dhigh; j++)
            if ((hi >> j) & 1L) m ^= 1ULL << bits[dlow + j];
        uint64_t local = encrypt_pref(plen, m, h1, h2, K);
        for (uint64_t g = 1; g < (1ULL << dlow); g++) {
            int j = __builtin_ctzll(g);
            m ^= 1ULL << bits[j];
            local ^= encrypt_pref(plen, m, h1, h2, K);
        }
        acc ^= local;
    }
    return acc;
}

/* ------------------------------------------------------------------- main */

static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;
static uint64_t rnd(void) {
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17; return rng_state;
}

static void run_case(int k, int dim, int trials) {
    int nsbox = 3 + k;                 /* pi1's 3 S layers are free */
    printf("\nk = %d  ->  %d S layers,  cube dimension %d (2^%d tweaks)\n",
           k, nsbox, dim, dim);
    for (int tr = 0; tr < trials; tr++) {
        uint64_t kw[7]; Key K;
        for (int i = 0; i < 7; i++) kw[i] = rnd();
        key_expand(kw, &K);
        uint64_t m = rnd();
        uint64_t base = rnd() | (1ULL << 63);      /* spec requires t != 0 */

        int bits[64];
        for (int i = 0; i < 63; i++) bits[i] = i;
        for (int i = 62; i > 0; i--) {             /* shuffle, pick dim of them */
            int j = (int)(rnd() % (uint64_t)(i + 1));
            int tmp = bits[i]; bits[i] = bits[j]; bits[j] = tmp;
        }
        uint64_t s_ok  = cube_sum(nsbox, dim,     bits, m, base, &K);
        uint64_t s_low = cube_sum(nsbox, dim - 1, bits, m, base, &K);
        printf("  trial %d: dim %2d -> %016llx %-4s | dim %2d -> %016llx %s\n",
               tr, dim, (unsigned long long)s_ok, s_ok ? "FAIL" : "OK",
               dim - 1, (unsigned long long)s_low,
               s_low ? "(tight)" : "(also zero)");
    }
}

/* Single-tweak control: cube the plaintext, at the SAME total S-layer count the
 * multi-tweak runs cover, so the two are directly comparable. */
static void run_st_case(int nsbox, int dim, int trials) {
    printf("\nST  %d S layers,  plaintext cube dimension %d (2^%d data)\n",
           nsbox, dim, dim);
    for (int tr = 0; tr < trials; tr++) {
        uint64_t kw[7]; Key K;
        for (int i = 0; i < 7; i++) kw[i] = rnd();
        key_expand(kw, &K);
        uint64_t t = rnd() | (1ULL << 63);
        uint64_t h1 = toeplitz(t, K.k1b), h2 = toeplitz(t, K.k2b);
        uint64_t base_m = rnd();

        int bits[64];
        for (int i = 0; i < 64; i++) bits[i] = i;
        for (int i = 63; i > 0; i--) {
            int j = (int)(rnd() % (uint64_t)(i + 1));
            int tmp = bits[i]; bits[i] = bits[j]; bits[j] = tmp;
        }
        uint64_t s = st_cube_sum(nsbox, dim, bits, base_m, h1, h2, &K);
        printf("  trial %d: dim %2d -> %016llx %s\n",
               tr, dim, (unsigned long long)s, s ? "FAIL" : "OK");
    }
}

/* A single zero bit could be luck (1/2). A bit that is zero across every trial,
 * each with an independent key/tweak/plaintext/bit-selection, is balanced with
 * confidence 2^-trials. Reports that count rather than all-or-nothing. */
static void run_balanced_bits(int mode, int nsbox, int dim, int trials) {
    uint64_t always_zero = ~0ULL;
    for (int tr = 0; tr < trials; tr++) {
        uint64_t kw[7]; Key K;
        for (int i = 0; i < 7; i++) kw[i] = rnd();
        key_expand(kw, &K);
        uint64_t m = rnd(), base = rnd() | (1ULL << 63);
        int bits[64];
        int navail = (mode == 0) ? 63 : 64;      /* MT cubes tweak bits 0..62 */
        for (int i = 0; i < navail; i++) bits[i] = i;
        for (int i = navail - 1; i > 0; i--) {
            int j = (int)(rnd() % (uint64_t)(i + 1));
            int tmp = bits[i]; bits[i] = bits[j]; bits[j] = tmp;
        }
        uint64_t s;
        if (mode == 2) {
            /* free-round extension: S is an involution, so applying it back
             * costs no key. The round key and the secret h(t) both cancel in
             * the sum, so XOR S(C) = 0 one layer past the distinguisher. */
            s = cube_sum_sapplied(nsbox, dim, bits, m, base, &K);
        } else if (mode == 0) {
            s = cube_sum(nsbox, dim, bits, m, base, &K);
        } else {
            uint64_t h1 = toeplitz(base, K.k1b), h2 = toeplitz(base, K.k2b);
            s = st_cube_sum(nsbox, dim, bits, m, h1, h2, &K);
        }
        always_zero &= ~s;
    }
    printf("  %s  %d S layers  dim %2d : %2d / 64 bits balanced in all %d trials\n",
           mode == 0 ? "MT" : (mode == 2 ? "MT+S" : "ST"), nsbox, dim,
           __builtin_popcountll(always_zero), trials);
}

int main(int argc, char **argv) {
    init_tables();
    build_ops();
    setvbuf(stdout, NULL, _IOLBF, 0);        /* line-buffer even when redirected */
    if (!self_test()) { fprintf(stderr, "self-test failed, aborting\n"); return 1; }
    if (argc > 3 && argv[1][0] == 's') {     /* st <nsbox> <dim> [trials] */
        run_st_case(atoi(argv[2]), atoi(argv[3]),
                    argc > 4 ? atoi(argv[4]) : 3);
        return 0;
    }
    if (argc > 4 && argv[1][0] == 'b') {     /* bal <mt|st> <nsbox> <dim> [tr] */
        run_balanced_bits(argv[2][0] == 'm' ? 0 : (argv[2][0] == 'f' ? 2 : 1), atoi(argv[3]),
                          atoi(argv[4]), argc > 5 ? atoi(argv[5]) : 3);
        return 0;
    }
#ifdef _OPENMP
    printf("OpenMP threads: %d\n", omp_get_max_threads());
#endif
    if (argc > 2) {
        run_case(atoi(argv[1]), atoi(argv[2]), argc > 3 ? atoi(argv[3]) : 3);
        return 0;
    }

    run_case(1, 4,  3);      /* deg 3  -> predicted tight */
    run_case(2, 10, 3);      /* deg 9  */
    run_case(3, 28, 2);      /* deg 27 -- the one Python could not reach */
    return 0;
}
