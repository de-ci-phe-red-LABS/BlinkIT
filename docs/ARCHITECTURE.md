# Blink-64 architecture — AES-style component view

Every diagram below is for **Blink-64** (64-bit block, 16 nibbles). Blink-128 uses
the same components on a 4 x 8 grid.

---

## 1. The state

64 bits = 16 nibbles, arranged in **4 rows x 4 columns**, filled **row-major**
(THF Sec. 5.2). Cell `s0` is the least significant nibble.

```
                 col0    col1    col2    col3
              ┌───────┬───────┬───────┬───────┐
       row0   │  s0   │  s1   │  s2   │  s3   │
              ├───────┼───────┼───────┼───────┤
       row1   │  s4   │  s5   │  s6   │  s7   │
              ├───────┼───────┼───────┼───────┤
       row2   │  s8   │  s9   │  s10  │  s11  │
              ├───────┼───────┼───────┼───────┤
       row3   │  s12  │  s13  │  s14  │  s15  │
              └───────┴───────┴───────┴───────┘

       each cell = one nibble (4 bits)
```

---

## 2. Component ↔ AES correspondence

| AES | Blink | Difference that matters |
|---|---|---|
| SubBytes | **S** (SubCells) | 4-bit S-box, and it is an **involution** (`S∘S = id`) |
| ShiftRows | **P** (ShuffleCell) | **not** a row rotation — an arbitrary cell permutation |
| MixColumns | **M** | **binary** matrix (XOR only), no GF multiplication — this is the low-latency choice |
| AddRoundKey | **AK** + **AC** | round key and round constant, added separately |

Both `S` and `M` are involutions. That is deliberate: it is what lets Blink run the
second half of the cipher **backwards through the same hardware** (the reflection
structure in §6).

---

## 3. S — SubCells

The 4-bit S-box applied to all 16 cells in parallel.

```
   ┌────┬────┬────┬────┐              ┌────┬────┬────┬────┐
   │ s0 │ s1 │ s2 │ s3 │              │S(s0)│S(s1)│S(s2)│S(s3)│
   ├────┼────┼────┼────┤              ├────┼────┼────┼────┤
   │ s4 │ s5 │ s6 │ s7 │    ──S──▶    │S(s4)│S(s5)│S(s6)│S(s7)│
   ├────┼────┼────┼────┤              ├────┼────┼────┼────┤
   │ s8 │ s9 │s10 │s11 │              │S(s8)│S(s9)│S(s10)│S(s11)│
   ├────┼────┼────┼────┤              ├────┼────┼────┼────┤
   │s12 │s13 │s14 │s15 │              │S(s12)│S(s13)│S(s14)│S(s15)│
   └────┴────┴────┴────┘              └────┴────┴────┴────┘
```

```
    x     0  1  2  3  4  5  6  7  8  9  a  b  c  d  e  f
  S(x)    1  0  9  3  8  5  e  7  4  2  c  b  a  f  6  d
```

**Algebraic degree 3.** Involutive: `S(S(x)) = x`.

---

## 4. M — MixColumn

Applied **column by column**, exactly like AES MixColumns, but the matrix is binary
so it is pure XOR.

```
         ⎡ 0 1 1 1 ⎤   ⎡ s_j    ⎤       ⎡ s_{j+4} ⊕ s_{j+8} ⊕ s_{j+12} ⎤
         ⎢ 1 0 1 1 ⎥ · ⎢ s_{j+4}⎥   =   ⎢ s_j     ⊕ s_{j+8} ⊕ s_{j+12} ⎥
         ⎢ 1 1 0 1 ⎥   ⎢ s_{j+8}⎥       ⎢ s_j     ⊕ s_{j+4} ⊕ s_{j+12} ⎥
         ⎣ 1 1 1 0 ⎦   ⎣s_{j+12}⎦       ⎣ s_j     ⊕ s_{j+4} ⊕ s_{j+8}  ⎦
```

**Each output cell is the XOR of the other three cells in its column.**

```
              col0    col1    col2    col3
           ┌───────┬───────┬───────┬───────┐
     row0  │  s0   │  s1   │  s2   │  s3   │  ─┐
           ├───────┼───────┼───────┼───────┤   │
     row1  │  s4   │  s5   │  s6   │  s7   │   │  M acts DOWN
           ├───────┼───────┼───────┼───────┤   │  each column
     row2  │  s8   │  s9   │  s10  │  s11  │   │  independently
           ├───────┼───────┼───────┼───────┤   │
     row3  │  s12  │  s13  │  s14  │  s15  │  ─┘
           └───────┴───────┴───────┴───────┘
               │       │       │       │
               ▼       ▼       ▼       ▼
              col0    col1    col2    col3
```

This is the **Midori** involutory near-MDS matrix: branch number 4 (AES MixColumns
has 5), but it costs only XOR gates, which is the whole point for latency.

`M ∘ M = identity`.

---

## 5. AK / AC — AddRoundKey and AddRoundConstant

Cell-wise XOR over the whole state.

```
   ┌────┬────┬────┬────┐   ┌────┬────┬────┬────┐   ┌────┬────┬────┬────┐
   │ s0 │ s1 │ s2 │ s3 │   │ k0 │ k1 │ k2 │ k3 │   │ c0 │ c1 │ c2 │ c3 │
   ├────┼────┼────┼────┤   ├────┼────┼────┼────┤   ├────┼────┼────┼────┤
   │ s4 │ s5 │ s6 │ s7 │ ⊕ │ k4 │ k5 │ k6 │ k7 │ ⊕ │ c4 │ c5 │ c6 │ c7 │
   ├────┼────┼────┼────┤   ├────┼────┼────┼────┤   ├────┼────┼────┼────┤
   │ s8 │ s9 │s10 │s11 │   │ k8 │ k9 │k10 │k11 │   │ c8 │ c9 │c10 │c11 │
   ├────┼────┼────┼────┤   ├────┼────┼────┼────┤   ├────┼────┼────┼────┤
   │s12 │s13 │s14 │s15 │   │k12 │k13 │k14 │k15 │   │c12 │c13 │c14 │c15 │
   └────┴────┴────┴────┘   └────┴────┴────┴────┘   └────┴────┴────┴────┘
          state                round key            round constant
```

**No key schedule.** The 448-bit master key is simply sliced:

```
      k  =  rk5 ‖ rk4 ‖ rk3 ‖ rk2 ‖ rk1 ‖ w2 ‖ w1      (7 x 64 bits)
```

so any four consecutive round keys are independent — which is exactly why key
recovery against Blink is expensive.

---

## 6. P — ShuffleCell

AES rotates each row. **Blink does not** — `P` is an arbitrary cell permutation:

```
    P = [0, 5, 11, 10, 1, 6, 4, 13, 2, 12, 9, 15, 3, 7, 14, 8]

    new_s[i] = old_s[P[i]]
```

Written on the grid — each destination cell shows **where it pulls from**:

```
            BEFORE                                AFTER
    ┌────┬────┬────┬────┐              ┌────┬────┬────┬────┐
    │ s0 │ s1 │ s2 │ s3 │              │ s0 │ s5 │s11 │s10 │
    ├────┼────┼────┼────┤              ├────┼────┼────┼────┤
    │ s4 │ s5 │ s6 │ s7 │   ──P──▶     │ s1 │ s6 │ s4 │s13 │
    ├────┼────┼────┼────┤              ├────┼────┼────┼────┤
    │ s8 │ s9 │s10 │s11 │              │ s2 │s12 │ s9 │s15 │
    ├────┼────┼────┼────┤              ├────┼────┼────┼────┤
    │s12 │s13 │s14 │s15 │              │ s3 │ s7 │s14 │ s8 │
    └────┴────┴────┴────┘              └────┴────┴────┴────┘
```

Note column 0 of the result is `s0, s1, s2, s3` — the **first row** of the input.
So `P` mixes rows into columns, which is what lets `M` (a column operation) spread
differences across the whole state. Together `P` and `M` give **full diffusion in
3 rounds**.

Cells `s0` and `s14` are fixed points of `P`.

---

## 7. The round function

```
        R   =   P  ∘  AC  ∘  AK  ∘  M  ∘  S              (forward)
```

Read right to left — the order of application is:

```
    ┌───┐     ┌───┐     ┌────┐     ┌────┐     ┌───┐
 ──▶│ S │────▶│ M │────▶│ AK │────▶│ AC │────▶│ P │──▶
    └───┘     └───┘     └────┘     └────┘     └───┘
   SubCells  MixCol   AddRoundKey  AddConst  ShuffleCell
```

The backward round, used in the reflection halves:

```
        R̄   =   S  ∘  MK  ∘  P̄        where  MK_k(x) = M(x) ⊕ k

    ┌────┐     ┌───┐     ┌──────────────┐     ┌───┐
 ──▶│ P̄  │────▶│ M │────▶│ ⊕ M(rk)⊕M(rc)│────▶│ S │──▶
    └────┘     └───┘     └──────────────┘     └───┘
```

Note the round key and constant are themselves passed through `M` in the backward
direction (`r̂k = M(rk)`), which is what keeps `R̄` a true inverse-shaped round while
reusing the same `M` hardware.

---

## 8. Full cipher — Blink-64a

Reflection structure, `a = 2`, `b = 3`, **14 S-box layers** total.

```
                              m   (plaintext)
                              │
                            ⊕ w1
                              │
        ╔═════════════════════▼═════════════════════════════════╗
        ║  pi1        R(rk1) ─▶ R(rk2) ─▶ S ─▶ M                 ║  3 S
        ╚═════════════════════╤═════════════════════════════════╝
                              │
                          ⊕ h1(t)   ◄═══════ TWEAK (1st injection)
                              │
        ╔═════════════════════▼═════════════════════════════════╗
        ║  pi2   P ─▶ R(rk3) ─▶ R(rk4) ─▶ R(rk5) ─▶ S ─▶ M      ║  4 S
        ╚═════════════════════╤═════════════════════════════════╝
                              │
                        ⊕ h(t) = h1(t) ⊕ h2(t)   ◄═══ TWEAK (2nd)
                              │
        ╔═════════════════════▼═════════════════════════════════╗
        ║  pi3   S ─▶ R̄(rk1) ─▶ R̄(rk2) ─▶ R̄(rk3) ─▶ P̄          ║  4 S
        ╚═════════════════════╤═════════════════════════════════╝
                              │
                          ⊕ h2(t)   ◄═══════ TWEAK (3rd injection)
                              │
        ╔═════════════════════▼═════════════════════════════════╗
        ║  pi4   M ─▶ S ─▶ R̄(rk4) ─▶ R̄(rk5)                     ║  3 S
        ╚═════════════════════╤═════════════════════════════════╝
                              │
                            ⊕ w2
                              │
                              c   (ciphertext)
```

**Round-key order.** `pi3` uses `rk1, rk2, rk3` and `pi4` uses `rk4, rk5` — the
literal reading of the specification's Figure 2, **not** the PRINCE-style reversal
the reflection shape suggests. This is easy to get wrong; assuming the mirrored
order is why an earlier implementation never matched the official test vectors.

---

## 9. The tweak path

This is what distinguishes Blink from SKINNY-style designs.

```
                    t   (tweak, public)
                    │
        ┌───────────┴───────────┐
        │                       │
        ▼                       ▼
   ┌─────────┐            ┌─────────┐
   │ h1 = T1·t│           │ h2 = T2·t│      T1, T2 = Toeplitz matrices
   └────┬────┘            └────┬────┘       whose entries come FROM THE KEY
        │                      │
        │      ┌───────────────┤
        │      │               │
        ▼      ▼               ▼
     ⊕ h1    ⊕ (h1⊕h2)      ⊕ h2
    after pi1  after pi2    after pi3
```

Two properties, both load-bearing:

1. **The tweak schedule is SECRET** — `T1`, `T2` depend on the key. In SKINNY,
   Deoxys or QARMAv2 the tweakey schedule is public and can be modelled directly.
   Here it cannot. That looks like a defensive advantage.

2. **The injections are LINEAR in `t`.** And that is what defeats the first point:
   a linear substitution cannot raise algebraic degree, so
   `deg_t(c) ≤ deg(F)` **whatever `T1` and `T2` are**. The secret never enters the
   bound. See `ATTACK_STRATEGY.md` §4.2.

The tweak reaches the state at only **three seams**, and the first is *after*
`pi1` — so with the plaintext held fixed, `pi1`'s three S-box layers are traversed
by a constant and cost the attacker nothing.

---

## 10. Implementation note

With nibbles packed least-significant-first, **each row of the grid lands in one
16-bit chunk** of the 64-bit state word:

```
    bits  0..15   →  row0  (s0  s1  s2  s3 )
    bits 16..31   →  row1  (s4  s5  s6  s7 )
    bits 32..47   →  row2  (s8  s9  s10 s11)
    bits 48..63   →  row3  (s12 s13 s14 s15)
```

Since `M` acts down the columns — i.e. **across** the rows — the entire MixColumn
layer for all four columns at once is three XORs on 16-bit chunks:

```c
    n0 = r1 ^ r2 ^ r3;   n1 = r0 ^ r2 ^ r3;
    n2 = r0 ^ r1 ^ r3;   n3 = r0 ^ r1 ^ r2;
```

This, plus Gray-code enumeration of the tweak cube (legitimate precisely because
`h1` and `h2` are linear), is what makes 2^38 encryptions feasible in `blink_cube.c`.

---

## 11. Parameters

| | Blink-64a | Blink-64b | Blink-128a/b | Blink-128A/B |
|---|---|---|---|---|
| Block `n` | 64 | 64 | 128 | 128 |
| Tweak `tau` | 64 | 128 | 128 / 256 | 128 / 256 |
| `a` | 2 | 2 | 3 | 3 |
| `b` | 3 | 3 | 3 | 5 |
| Key | 448 | 448 | 1024 | 1280 |
| Grid | 4 x 4 | 4 x 4 | 4 x 8 | 4 x 8 |
| Data limit | 2^56 | 2^56 | 2^80 | 2^80 |

Key length is `(a + b + 2) · n`.
