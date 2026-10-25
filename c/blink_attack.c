/*
 * blink_attack.c -- the complete multi-tweak key-recovery attack on Blink-64a.
 *
 * Chain of the attack (see ATTACK_STRATEGY.md):
 *
 *   r = 3 + k     distinguisher   XOR X_r        = 0     (tweak cube, dim > deg F)
 *   r + 1         free round      XOR S(X_{r+1}) = 0     (S is an involution)
 *   r + 2         KEY RECOVERY    guess 16 bits per MixColumn column
 *
 * At layer r+2 the cipher applies one round, so
 *
 *   forward transition :  X_{r+2} = S( P( M(X_{r+1}) ^ K ) )
 *                         X_{r+1} = M( Pbar(S(X_{r+2})) ^ K )
 *   backward transition:  X_{r+2} = S( M( Pbar(X_{r+1}) ) ^ K )
 *                         X_{r+1} = P( M( S(X_{r+2}) ^ K ) )
 *
 * Both have the shape  X_{r+1} = OUTER( M( INNER(S(X_{r+2})) ^ K ) ), and OUTER
 * is a cell permutation, so it commutes with S and with the cube sum and can be
 * dropped. The test is therefore
 *
 *      XOR over the cube of  S( M( INNER(S(X_{r+2})) ^ K ) )  ==  0
 *
 * M acts column-wise, so one column of the result depends on 4 nibbles (16 bits)
 * of K and yields 4 nibbles -- a 16-bit guess tested by a 16-bit condition.
 *
 * PARTIAL SUMS. The 16-bit column value is all that matters, so values occurring
 * an even number of times cancel. One pass builds a 2^16 parity table per column
 * (32 KB total); then each guess costs 2^16 instead of 2^d. Time collapses from
 * 2^16 * 2^d to 2^d + 2^32.
 *
 * Build: gcc -O3 -march=native -fopenmp -o blink_attack blink_attack.c
 * Run:   ./blink_attack            k=2, dim 8   (instant smoke test)
 *        ./blink_attack 4 37       the real attack (~15-20 min on 14 threads)
 *        ./blink_attack <k> <dim> <trials> [cubeid]   (cubeid>0: same key, cube #cubeid)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define PI1_SBOXES 3            /* pi1 contributes a+1 = 3 free S layers */
#define NCOL       4
#define TABW      (1 << 16)     /* 16-bit column value */
#define TABWORDS  (TABW / 64)

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

static uint16_t SBOX16[65536];

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

static inline uint64_t M_layer(uint64_t s) {
    uint64_t r0 = s & 0xffffULL, r1 = (s >> 16) & 0xffffULL;
    uint64_t r2 = (s >> 32) & 0xffffULL, r3 = (s >> 48) & 0xffffULL;
    return (r1^r2^r3) | ((r0^r2^r3) << 16) | ((r0^r1^r3) << 32) | ((r0^r1^r2) << 48);
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

static int prefix_len(int nsbox) {
    int seen = 0;
    for (int i = 0; i < NOPS; i++)
        if (OPS[i].op == OP_S && ++seen == nsbox) return i + 1;
    return NOPS;
}

/* ---------------------------------------------------------------- key data */

typedef struct {
    uint64_t w1, w2, rk[5], Mrk[5], Mrcp[5];
    uint8_t k1b[127], k2b[127];
} Key;

static void key_expand(const uint64_t kw[7], Key *K) {
    K->w1 = kw[0]; K->w2 = kw[1];
    for (int i = 0; i < 5; i++) {
        K->rk[i] = kw[2+i];
        K->Mrk[i] = M_layer(kw[2+i]);
        K->Mrcp[i] = M_layer(RCP[i]);
    }
    uint8_t kp[448];
    for (int i = 0; i < 448; i++) {
        int b = (11 * i) % 448;
        kp[i] = (uint8_t)((kw[b >> 6] >> (b & 63)) & 1ULL);
    }
    memcpy(K->k1b, kp, 127);
    memcpy(K->k2b, kp + 127, 127);
}

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

/* ------------------------------------------------------ the target round */

/* Ops strictly between S layer `a` and S layer `a+1` determine whether the
 * transition is a forward or a backward round, and which round key it uses. */
typedef struct {
    int backward;        /* 1 = backward round (Pbar,M,AKh,S) */
    int keyidx;
    uint64_t K;          /* the value XORed: rk^rc, or M(rk)^M(rc') */
} Target;

static void find_target(int a, const Key *Kk, Target *T) {
    int lo = prefix_len(a), hi = prefix_len(a + 1);
    T->backward = 0; T->keyidx = -1;
    for (int i = lo; i < hi; i++) {
        if (OPS[i].op == OP_PB) T->backward = 1;
        if (OPS[i].op == OP_AK)  { T->keyidx = OPS[i].idx; T->backward = 0; }
        if (OPS[i].op == OP_AKH) { T->keyidx = OPS[i].idx; T->backward = 1; }
    }
    /* Three transitions carry a tweak seam instead of a round key: S#3->S#4,
     * S#7->S#8 and S#11->S#12. There is nothing to guess there, and inverting
     * would need the per-sample h(t), which is secret and varies with t. So k
     * must be chosen so the target round is a keyed one. */
    if (T->keyidx < 0) {
        fprintf(stderr,
            "target round (S#%d -> S#%d) is a TWEAK SEAM, not a keyed round.\n"
            "  Nothing to guess, and h(t) varies per sample so it cannot be\n"
            "  inverted. Choose k = 1, 2 or 4 (k=3 lands on the seam).\n",
            a, a + 1);
        exit(1);
    }
    T->K = T->backward ? (Kk->Mrk[T->keyidx] ^ Kk->Mrcp[T->keyidx])
                       : (Kk->rk[T->keyidx] ^ RC[T->keyidx]);
}

/* INNER(S(X)) -- the value the guess is XORed into */
static inline uint64_t inner(uint64_t x, int backward) {
    uint64_t s = S_layer(x);
    return backward ? s : perm(s, P64I);      /* forward case needs Pbar */
}

/* four nibbles of column c, packed into 16 bits */
static inline uint32_t col16(uint64_t w, int c) {
    return (uint32_t)(((w >> (4*c)) & 0xf)
                   | (((w >> (4*(c+NCOL))) & 0xf) << 4)
                   | (((w >> (4*(c+2*NCOL))) & 0xf) << 8)
                   | (((w >> (4*(c+3*NCOL))) & 0xf) << 12));
}

/* ------------------------------------------------- pass 1: parity tables */

/* One walk over the cube. For each column, toggle parity[col][value]. */
static void collect(int nsbox_obs, int d, const int *bits, uint64_t m,
                    uint64_t base, const Key *K, int backward,
                    uint64_t par[NCOL][TABWORDS]) {
    int plen = prefix_len(nsbox_obs);
    uint64_t d1[64], d2[64];
    for (int j = 0; j < d; j++) {
        uint64_t e = 1ULL << bits[j];
        d1[j] = toeplitz(e, K->k1b);
        d2[j] = toeplitz(e, K->k2b);
    }
    uint64_t h1_0 = toeplitz(base, K->k1b), h2_0 = toeplitz(base, K->k2b);
    int dlow = d < 16 ? d : 16, dhigh = d - dlow;

    for (int c = 0; c < NCOL; c++) memset(par[c], 0, sizeof(uint64_t)*TABWORDS);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        uint64_t loc[NCOL][TABWORDS];
        for (int c = 0; c < NCOL; c++) memset(loc[c], 0, sizeof(uint64_t)*TABWORDS);
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
        for (long hi = 0; hi < (1L << dhigh); hi++) {
            uint64_t h1 = h1_0, h2 = h2_0;
            for (int j = 0; j < dhigh; j++)
                if ((hi >> j) & 1L) { h1 ^= d1[dlow+j]; h2 ^= d2[dlow+j]; }
            for (uint64_t g = 0; g < (1ULL << dlow); g++) {
                if (g) {
                    int j = __builtin_ctzll(g);
                    h1 ^= d1[j]; h2 ^= d2[j];
                }
                uint64_t w = inner(encrypt_pref(plen, m, h1, h2, K), backward);
                for (int c = 0; c < NCOL; c++) {
                    uint32_t v = col16(w, c);
                    loc[c][v >> 6] ^= 1ULL << (v & 63);
                }
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        for (int c = 0; c < NCOL; c++)
            for (int i = 0; i < TABWORDS; i++) par[c][i] ^= loc[c][i];
    }
}

/* ------------------------------------------ pass 2: 2^16 guesses per column */

/* For guess g, XOR S(M(v ^ g)) over every v of odd parity. Returns the four
 * output nibbles packed; 0 means the guess survives. */
static uint32_t test_guess(const uint64_t *par, uint32_t g) {
    uint32_t a0 = 0, a1 = 0, a2 = 0, a3 = 0;
    for (int wi = 0; wi < TABWORDS; wi++) {
        uint64_t word = par[wi];
        while (word) {
            int b = __builtin_ctzll(word);
            word &= word - 1;
            uint32_t v = ((uint32_t)wi << 6 | (uint32_t)b) ^ g;
            uint32_t i0 = v & 0xf, i1 = (v >> 4) & 0xf,
                     i2 = (v >> 8) & 0xf, i3 = (v >> 12) & 0xf;
            a0 ^= SBOX[i1 ^ i2 ^ i3];       /* M row 0 = XOR of other three */
            a1 ^= SBOX[i0 ^ i2 ^ i3];
            a2 ^= SBOX[i0 ^ i1 ^ i3];
            a3 ^= SBOX[i0 ^ i1 ^ i2];
        }
    }
    return a0 | a1 | a2 | a3;
}

static int recover_column(const uint64_t *par, uint32_t truth,
                          uint32_t *cands, int maxc) {
    int n = 0;
    int hit = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long g = 0; g < TABW; g++) {
        if (test_guess(par, (uint32_t)g) == 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                if (n < maxc) cands[n] = (uint32_t)g;
                n++;
                if ((uint32_t)g == truth) hit = 1;
            }
        }
    }
    return hit ? n : -n;      /* negative n signals the true key was missed */
}

/* ------------------------------------------------------------------- main */

static uint64_t rng = 0x9e3779b97f4a7c15ULL;
static uint64_t rnd(void) {
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng;
}

static const char *KEY_HEX =
    "d6a102d888a467e4d1d7dec33a246943e07c1dc6f302c57e762c2df9de6f0d21"
    "6dd387874a0b52ce3022e0ad78c78a0697779021b38e7fa1";

static void parse_key(const char *hex, uint64_t kw[7]) {
    size_t len = strlen(hex); char buf[17];
    for (int w = 0; w < 7; w++) {
        memcpy(buf, hex + (len - (size_t)w*16) - 16, 16); buf[16] = 0;
        kw[w] = strtoull(buf, NULL, 16);
    }
}

static int self_test(void) {
    uint64_t kw[7]; Key K;
    parse_key(KEY_HEX, kw); key_expand(kw, &K);
    uint64_t t = 0x0123456789abcdefULL;
    uint64_t c = encrypt_pref(NOPS, 0, toeplitz(t, K.k1b), toeplitz(t, K.k2b), &K);
    int ok = c == 0xa4a0d10502be846eULL;
    printf("ToSC F.1 test vector: %016llx %s\n\n",
           (unsigned long long)c, ok ? "PASS" : "FAIL");
    return ok;
}

/* cubeid > 0: SAME key for every run (the key of trial 0), cube `cubeid`
 * (plaintext, base tweak, positions) drawn from a seed derived from cubeid;
 * cubeid == 1 reproduces trial 0 exactly. Full candidate lists are printed so
 * that independent runs can be intersected offline. */
static void attack(int k, int dim, int trials, int cubeid) {
    int r        = PI1_SBOXES + k;   /* distinguisher depth */
    int obs      = r + 2;            /* we observe two layers further */
    printf("=== attack: k=%d  distinguisher %d S layers, free round -> %d, "
           "key recovery at %d ===\n", k, r, r+1, obs);
    printf("cube dimension %d (2^%d chosen tweaks)\n\n", dim, dim);

    for (int tr = 0; tr < trials; tr++) {
        uint64_t kw[7]; Key K;
        for (int i = 0; i < 7; i++) kw[i] = rnd();
        key_expand(kw, &K);
        if (cubeid > 1) { rng = 0x9e3779b97f4a7c15ULL ^ ((uint64_t)cubeid * 0xd1b54a32d192ed03ULL); for (int i = 0; i < 7; i++) rnd(); }
        uint64_t m = rnd(), base = rnd() | (1ULL << 63);   /* spec: t != 0 */

        int bits[64];
        for (int i = 0; i < 63; i++) bits[i] = i;
        for (int i = 62; i > 0; i--) {
            int j = (int)(rnd() % (uint64_t)(i+1));
            int tmp = bits[i]; bits[i] = bits[j]; bits[j] = tmp;
        }

        Target T; find_target(r + 1, &K, &T);
        printf("trial %d (cube %d): target round is %s, round key index %d, K = %016llx\n",
               tr, cubeid, T.backward ? "backward" : "forward", T.keyidx,
               (unsigned long long)T.K);

        static uint64_t par[NCOL][TABWORDS];
        collect(obs, dim, bits, m, base, &K, T.backward, par);

        /* Diagnostic: how many of the 2^16 column values have odd parity.
         * If this is 0 the filter is vacuous -- every guess passes. Too large a
         * cube drives it there, so the attack wants the SMALLEST cube that still
         * gives the distinguisher. */
        for (int c = 0; c < NCOL; c++) {
            long pop = 0;
            for (int i = 0; i < TABWORDS; i++) pop += __builtin_popcountll(par[c][i]);
            printf("  column %d parity table: %ld / %d odd entries\n", c, pop, TABW);
        }

        int recovered = 0;
        for (int c = 0; c < NCOL; c++) {
            uint32_t truth = col16(T.K, c);
            static uint32_t cands[65536];
            int n = recover_column(par[c], truth, cands, 65536);
            int found = n >= 0;
            if (!found) n = -n;
            recovered += found;
            printf("  column %d: %3d candidate(s) of 2^16   true %04x  %s\n",
                   c, n, truth, found ? "RECOVERED" : "MISSED");
            if (cubeid > 0) {
                printf("  CANDS %d:", c);
                for (int i = 0; i < n && i < 65536; i++) printf(" %04x", cands[i]);
                printf("\n");
            }
        }
        printf("  -> %d/4 columns recovered\n\n", recovered);
    }
}

int main(int argc, char **argv) {
    init_tables();
    build_ops();
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (!self_test()) return 1;
#ifdef _OPENMP
    printf("OpenMP threads: %d\n\n", omp_get_max_threads());
#endif
    int k      = argc > 1 ? atoi(argv[1]) : 2;
    int dim    = argc > 2 ? atoi(argv[2]) : 8;
    int trials = argc > 3 ? atoi(argv[3]) : 1;
    int cubeid = argc > 4 ? atoi(argv[4]) : 0;   /* same-key multi-cube mode */
    attack(k, dim, trials, cubeid);
    return 0;
}
