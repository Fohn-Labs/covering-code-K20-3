# Research notes: structure of the 1536-word code and searches below it

Working notes from a search for K(20,3) < 1536 (2026-11-03). No smaller code has been
found; the best upper bound is still K(20,3) <= 1536. These notes record what was
checked, so later searches can build on it.

Labels: **exact** means an exhaustive computation in the tools below; **solver** means
a CP-SAT status (not an independently checked certificate); **heuristic** means a search
that found nothing, which proves nothing.

## 1. The record lives in one coset of a [20,15] code (exact)

For C = { x : Hx in S }, let V = span(S - S) in F_2^r. If dim V = d < r, change
coordinates in F_2^r so that S lies in F_2^d x {0}. Then H splits as
H_top (d rows) over H_bot (r - d rows), C lies in a single coset of the linear code
M = ker H_bot of dimension 20 - (r - d), and the covering condition splits into
2^(r-d) fibers: for each y in F_2^(r-d),

    T + B(y) = F_2^d,   B(y) = { H_top e : wt(e) <= 3, H_bot e = y },

where T is S in the new coordinates. Since the problem is translation invariant, the LP
bound for this H is |T| >= 2^d / min_y |B(y)|.

| Code | r | dim span(S - S) | Fibers | Fiber sizes | LP bound on \|S\| | \|S\| |
|--:|--:|--:|--:|:--|--:|--:|
| 1792 | 12 | 6 | 64 | 13 to 24 | 4.92 | 7 |
| 1664 | 13 | 10 | 8 | 96 to 228 | 10.67 | 13 |
| 1600 | 14 | 11 | 8 | 96 to 228 | 21.33 | 25 |
| 1536 | 14 | 9 | 32 | 24 to 60 | 21.33 | 24 |

For the 1536-word code the 12 smallest fibers have 24 points each and the 24 translates
of T cover F_2^9 (512 points) with 448 points once and 64 twice: a near tiling. The
fibers fall into classes with identical multiplicity profiles (sizes 24 x 12, 31 x 6,
36 x 6, 48, 55, 57 x 4, 58, 60), which led to item 2.

Reproduce: `python3 research/tools/fibers.py`.

## 2. Automorphisms of the record (exact)

- ker H, the [20,6] code, has 55296 coordinate automorphisms.
- Exactly 384 coordinate permutations sigma map C to a translate of itself
  (sigma(C) = C + c). Element orders: 1 (1), 2 (159), 3 (8), 4 (96), 6 (120).
  Coordinate orbits (bit positions 0..19 of the integer form): {0,14}, {1,5,7,16},
  {3,4,6,8,12,15}, {13,18}, and six fixed coordinates {2,9,10,11,17,19}.
- With the translations by ker H this gives an affine automorphism group of order
  384 * 64 = 24576.

The search enumerates all permutations that preserve single and pairwise weight
invariants of ker H, so it is complete. The 384 permutations are in
`research/data/aut1536_perms.json`.

Reproduce: `python3 research/tools/autgroup.py construction-1536.json out.json` (about
1 minute).

## 3. Codes with prescribed symmetry from the record's group

`research/tools/symcover.py` takes a subgroup K of the affine group above. It computes
the K-orbits on F_2^20 and solves the exact covering problem over orbits with CP-SAT:
choose orbits whose union has covering radius 3, minimising the total size. The record
is always feasible, since it is K-invariant.

| K | \|K\| | Orbits | Result |
|:--|--:|--:|:--|
| full affine group | 24576 | 1200 | **solver: OPTIMAL 1536** in 4 s. No smaller code has the record's full symmetry. |
| 8 random elements + 2-dim kernel translations | 6144 | 3840 | solver: best 1536, bound 576 after 600 s |
| 8 random elements, no added translations | 3072 | 7680 | solver: best 1536, bound 512 after 600 s |
| 3 random elements + 3-dim kernel translations | 768 | 3584 | solver: best 1536, bound 1024 after 600 s |
| 3 random elements + full ker H | 1024 | 3840 | solver: best 1536, bound 1472 after 120 s |

Reproduce the first row:
`python3 research/tools/symcover.py --random-gens 8 --trans-dim 4 --seed 13 --time 120 --out /tmp/full`
(needs OR-Tools; tested with ortools 9.12.4544).

## 4. Other checks

- **Golay kernels (exact).** Shortening the [23,12,7] Golay code on 3 coordinates gives a
  [20,9,7] kernel. Every sum of at most 3 columns is distinct (1351 of 2048 syndromes),
  but no second coset covers the remaining 697 syndromes: the best covers 551. So this
  kernel gives no 1024-word code. `python3 research/tools/golay.py`.
- **Golay [20,8,8] kernel (heuristic).** Shortening the extended Golay code on 4
  coordinates gives all 1351 sums distinct in F_2^12. Tabu search over S still leaves
  70 syndromes uncovered with 6 cosets and 181 with 5. Distinct sums are not enough.
- **Record H, fewer syndromes (heuristic).** With H of the 1536-word code fixed, tabu
  search over S finds a 24-syndrome cover from a random start in about 1 s. With 23
  syndromes, 4 runs of 90 s bottom out at 253 uncovered syndromes, the record minus one
  coset.
- **Balanced fibers with random H_top (heuristic).** Choosing H_bot with 20 distinct
  nonzero columns in F_2^5 makes the fiber subset counts nearly equal (about 42 each).
  With random H_top, about 32 syndromes are needed (2048 words). The record's H_top is
  far from random.
- **Joint (H, S) tabu from random starts (heuristic).** `research/tools/jtabu.c` (exact
  Walsh-Hadamard scoring of every column replacement, optional restriction of S to a
  subspace) does not reach 1536 from random starts in 5 minutes per run (best 112 to 533
  uncovered over 8 runs). The known matrix-method records came from descending
  1792 -> 1664 -> 1600 -> 1536, not from random starts.

## Tools

| File | Purpose |
|:--|:--|
| `tools/jtabu.c` | joint weighted tabu over (H, S); `gcc -O3 -o jtabu jtabu.c`; options in the header comment |
| `tools/fibers.py` | span(S - S), fiber decomposition and LP bound for the construction files |
| `tools/autgroup.py` | coordinate automorphisms of a matrix-method code, up to translation |
| `tools/affgroup.py` | affine maps, group closure and orbits on F_2^20 |
| `tools/symcover.py` | exact orbit covering with CP-SAT for a subgroup of the record's group |
| `tools/golay.py` | Golay-derived kernels and the two-coset check |
| `tools/export.py` | writes `K20_3_<N>.txt/.hex` and `construction-<N>.json` from a search result, with an independent dilation check; then run `verify.py` |
