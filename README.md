# K(20,3) ≤ 1536

**A binary code of length 20 with 1536 codewords and covering radius 3, so K(20,3) ≤ 1536.**

This improves the best published upper bound, K(20,3) ≤ 1905 (K. Özdemir, Zenodo,
2026-10-26); the standard tables of Kéri and Litsyn list 2048. It is an upper
bound only: the best lower bound is 892, so the exact value of K(20,3) remains
open.

K(n, R) is the smallest number of binary words of length n such that every binary
word of length n is within Hamming distance R of at least one of them. This
repository gives an explicit 1536-word code for n = 20, R = 3, a short algebraic
description of it, the earlier codes found on the way (1904, 1792, 1664 and 1600
words), and a standard-library Python script that checks every one of them against
all 2<sup>20</sup> = 1,048,576 words exactly and rebuilds each matrix-method code
from its construction.

| Bound on K(20,3) | Value | Source |
|:--|--:|:--|
| Upper bound, this code | **1536** | this repository, v1.1 (found 2026-10-29) |
| Previous published upper bound | 1905 | K. Özdemir, Zenodo, 2026-10-26, [doi:10.5281/zenodo.23676952](https://doi.org/10.5281/zenodo.23676952) |
| Upper bound in the standard tables | 2048 | G. Kéri, [Tables for bounds on covering codes](https://old.sztaki.hu/~keri/codes/) (binary table [2_tables.pdf](https://old.sztaki.hu/~keri/codes/2_tables.pdf), dated 2009-10-15); S. Litsyn, [Table of covering codes](https://web.archive.org/web/20160306160151/http://www.eng.tau.ac.il:80/~litsyn/tablecr/) (last modified 1999-08-23). This is 4 · K(18,3) ≤ 4 · 512, with K(18,3) ≤ 512 from Brualdi and Pless (1990). |
| Lower bound | 892 | W. Haas (2007–2008), as credited in Kéri's table (Litsyn's 1999 table lists 889) |
| Sphere-covering bound | 777 | ⌈2<sup>20</sup> / 1351⌉ |

## History

The standard tables list 2048 (Litsyn, 1999; Kéri, 2009); the first published
improvement we know of is Özdemir's 1905-word code of 2026-10-26. Every code from
2026-10-29 below was accepted by the arena of [fohn.app](https://fohn.app) task #10
([fohn.app/t/10](https://fohn.app/t/10)), whose server checks each submitted code
exactly against all 2<sup>20</sup> words, and was checked again for this repository
(see [Verify](#verify)). Times are the arena's acceptance times, US Eastern
(EDT, UTC−4); in UTC they fall on 2026-10-30. The 1904- and 1792-word codes were
found earlier that evening and published here as v1.0 at 22:32 EDT, before their
arena submissions.

| Date | Words | Density | Construction | Found by | Code |
|:--|--:|--:|:--|:--|:--|
| 1990 | 2048 | 2.639 | 4 · K(18,3) ≤ 4 · 512 (direct sum with F<sub>2</sub><sup>2</sup>) | Brualdi and Pless (K(18,3) ≤ 512); listed in Litsyn's table (1999) and Kéri's table (2009) | |
| 2026-10-26 | 1905 | 2.454 | local search (ruin-and-recreate) seeded from a 1988-word code: a 1792-word matrix-method base plus 196 patch words | K. Özdemir, with Claude (Anthropic) | [doi:10.5281/zenodo.23676952](https://doi.org/10.5281/zenodo.23676952) (version 2; version 1 had 1917 words) |
| 2026-10-29 22:41 | 1904 | 2.453 | local search from the 1905-word code (removes 5 words, adds 4); not a union of cosets of a nontrivial linear code | GPT-6.1 Sol (OpenAI) in Codex; fohn candidate [#15](https://api.fohn.app/api/tasks/10/candidates/15) | [intermediate/K20_3_1904.txt](intermediate/K20_3_1904.txt) |
| 2026-10-29 22:42 | 1792 | 2.309 | matrix method, H 12 × 20, \|S\| = 7 (7 cosets of a [20,8] code) | Claude Opus 5.5 (Anthropic) in Claude Code; fohn candidate [#17](https://api.fohn.app/api/tasks/10/candidates/17) | [previous/K20_3_1792.txt](previous/K20_3_1792.txt) (v1.0 of this repository) |
| 2026-10-29 23:03 | 1664 | 2.144 | matrix method, H 13 × 20, \|S\| = 13 (13 cosets of a [20,7] code) | Claude Opus 5.5 in Claude Code; fohn candidate [#19](https://api.fohn.app/api/tasks/10/candidates/19) | [previous/K20_3_1664.txt](previous/K20_3_1664.txt) |
| 2026-10-29 23:21 | 1600 | 2.061 | matrix method, H 14 × 20, \|S\| = 25 (25 cosets of a [20,6] code) | Claude Opus 5.5 in Claude Code; fohn candidate [#21](https://api.fohn.app/api/tasks/10/candidates/21) | [previous/K20_3_1600.txt](previous/K20_3_1600.txt) |
| 2026-10-29 23:24 | **1536** | **1.979** | matrix method, H 14 × 20, \|S\| = 24 (24 cosets of a [20,6] code) | Claude Opus 5.5 in Claude Code; fohn candidate [#22](https://api.fohn.app/api/tasks/10/candidates/22) | [K20_3_1536.txt](K20_3_1536.txt) |

Density is N · 1351 / 2<sup>20</sup>, the average number of codewords within
distance 3 of a word (1 would be a perfect code). The candidate links return the
words as accepted by the arena (`object.words`, 5-digit hex); they are the same
sets as the files here.

## Construction

The code is built with the matrix method (Blokhuis and Lam for radius 1,
generalised by Östergård, see the references). Let H be the 14 × 20 binary
matrix below and S the set of 24 syndromes below. Then

C = { x ∈ F<sub>2</sub><sup>20</sup> : H x<sup>T</sup> ∈ S }.

H has rank 14, so its kernel is a linear [20, 6] code with 64 words, and C is the
union of 24 of its cosets: |C| = 24 · 2<sup>6</sup> = 1536.

C has covering radius 3 because every one of the 16384 syndromes in
F<sub>2</sub><sup>14</sup> can be written as s + h<sub>1</sub> + … + h<sub>t</sub>
with s ∈ S and t ≤ 3 columns h<sub>1</sub>, …, h<sub>t</sub> of H. (The sums of at
most 3 columns give 1139 distinct syndromes, and their 24 translates by S cover all
16384.) Given any word x, pick such an s and columns for the syndrome H x<sup>T</sup>;
flipping the matching t ≤ 3 coordinates of x gives a word with syndrome s, which is
a codeword.

Columns are in the order of the characters of the strings in
[K20_3_1536.txt](K20_3_1536.txt), left to right; syndromes are read top to bottom.

```
H =
0 0 1 1 0 0 1 1 0 1 0 0 0 1 0 1 1 1 1 0
0 1 0 0 1 0 1 0 1 0 0 0 0 0 0 1 1 0 0 0
0 0 1 0 1 1 1 0 1 0 1 0 0 1 1 0 1 1 1 0
1 1 1 0 1 0 0 0 1 1 0 0 0 1 1 0 1 0 1 1
1 1 0 0 0 0 0 1 1 0 0 1 0 1 1 0 0 0 1 1
1 1 0 0 0 1 0 1 1 1 1 1 1 0 0 1 0 1 1 0
1 0 1 1 0 0 0 1 1 1 1 0 1 0 1 0 1 0 1 0
0 0 1 0 1 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0
1 1 1 1 1 0 1 1 1 0 0 1 1 0 1 0 1 0 1 1
1 1 0 0 0 0 1 0 1 0 0 0 1 0 0 0 0 1 1 1
1 1 0 0 0 1 1 1 0 1 0 1 0 0 0 1 0 0 0 0
1 1 0 1 1 0 0 0 0 1 0 1 1 0 1 0 0 0 1 0
1 0 0 1 1 0 0 1 0 1 0 0 0 1 0 0 0 1 1 0
0 1 1 0 1 0 0 1 0 1 0 0 1 0 1 1 0 1 0 0

S = { 10001111010010, 11001111010111, 10010111110111, 00010111110111,
      01001111010111, 00101011001111, 01010110001010, 11010110001010,
      10011111011000, 01110011010010, 11110011010010, 10000011110101,
      01100001101101, 00001111010010, 10101011001111, 11100001101101,
      10111000001000, 00111000001000, 10100011011101, 00111111001101,
      00011111100101, 10111111110000, 00100011100000, 00000011001000 }
```

[construction-1536.json](construction-1536.json) gives the same H and S as a 0/1
matrix and as 14-bit integers, with the conventions spelled out. In the integer
form a codeword is the integer w whose 20-digit binary expansion (most significant
bit first) is its string, and bit j of w (value 2<sup>j</sup>) selects
`H_columns[j]`; so character k of the string (1-based, left to right) selects
column k of the matrix, which is `H_columns[20 − k]`, and a column read top to
bottom is a binary number with the top row as its most significant bit:

```
H_columns = [1584, 12214, 10515, 15520, 12553, 3749, 11842, 437, 812, 2432,
             9679, 8176, 9195, 14392, 2312, 7271, 8358, 11489, 5949, 1982]
S = [9170, 13271, 9719, 1527, 5079, 2767, 5514, 13706, 10200, 7378, 15570, 8437,
     6253, 978, 10959, 14445, 11784, 3592, 10461, 4045, 2021, 12272, 2272, 200]
```

The coset structure can also be read off the words alone: the stabiliser
L = { d : C + d = C } of the 1536 words is a linear code of dimension 6, and it is
exactly ker H. Any parity-check matrix H′ of L gives the same code as
{ x : H′x<sup>T</sup> ∈ S′ }, with S′ the 24 syndromes of the codewords under H′.
The same holds for the 1792-, 1664- and 1600-word codes, whose constructions are in
[previous/](previous/) in the same format. The 1904-word code has trivial stabiliser
{0}, so it is not a union of cosets of a nontrivial linear code.

Other properties of the 1536-word code, all checked exhaustively: every codeword is
needed (deleting any one leaves between 234 and 343 words uncovered); the minimum
distance is 2; each word of F<sub>2</sub><sup>20</sup> is within distance 3 of
between 1 and 8 codewords; the covering density is 1536 · 1351 / 2<sup>20</sup> ≈
1.979. With the direct sum and K(n+1, R+1) ≤ K(n, R), the new bound does not
improve any other entry of Kéri's binary table; it gives K(21,3) ≤ 3072, equal to
the table, and the amalgamated direct sum would at best give K(22,4) ≤ 1536, also
equal to the table (and it needs normality, which we have not checked).

## How it was found

All of the 2026-10-29 codes were found on [fohn.app](https://fohn.app), task #10,
by AI agents lent by Eric Encinger.

The matrix-method codes come from one Claude Opus 5.5 agent (Anthropic) running in
Claude Code, which searched over pairs (H, S) with a weighted tabu search. Each
iteration scores every replacement of one column of H or one element of S exactly,
using Walsh–Hadamard transforms to evaluate all 2<sup>r</sup> replacement values at
once, and the weights of uncovered syndromes grow until they are covered.

| Words | r × 20, \|S\| | Start | Seed | Iteration with 0 uncovered syndromes | Time |
|--:|:--|:--|--:|--:|--:|
| 1792 | 12 × 20, 7 | random (H, S) | 107 | 17977 | 108 s |
| 1664 | 13 × 20, 13 | a saved state of an earlier 13 × 20 run with 56 uncovered syndromes | 4302 | 1465 | 26 s |
| 1600 | 14 × 20, 25 | the 1664-word code lifted to 14 rows (one row appended to H, S replaced by S × {0, 1}: the same 1664 words) with one syndrome dropped, leaving 80 uncovered syndromes | 9513 | 108 | 5.4 s |
| 1536 | 14 × 20, 24 | the 1600-word code with its least harmful syndrome dropped, leaving 64 uncovered syndromes | 10401 | 15 | 0.8 s |

For the 1600- and 1536-word codes, moves that would lower the rank of H or create a
weight-1 word in its row space (a coloop) were not allowed. Each construction file
records the seed, iteration and time.

Özdemir's Zenodo deposit describes its starting point as a 1988-word code, a
1792-word matrix-method base plus 196 patch words; the 1792-word code here is a
complete matrix-method code of that size, with nothing added.

Independently, a GPT-6.1 Sol agent (OpenAI) running in Codex found the 1904-word
code by local search from Özdemir's 1905-word code:
[intermediate/K20_3_1904.txt](intermediate/K20_3_1904.txt). It shares 1900 words
with the 1905-word code.

## Verify

Python 3.6 or later, no packages needed:

```
python3 verify.py
```

It checks every code file in the repository (about 5 seconds). For each code it
marks the radius-3 ball (1351 words) around every codeword in a table of all
2<sup>20</sup> words, counting multiplicities, and requires every word to be
covered; it checks that radius-2 balls leave words uncovered, so the covering
radius is exactly 3; it checks that every codeword is needed; and it computes the
stabiliser L. When a `construction-<N>.json` sits next to the code file, it
rebuilds C = { x : Hx ∈ S } by computing the syndrome of all 2<sup>20</sup> words,
compares it with the file as sets, checks that ker H = L and that the 0/1 matrix
form agrees with the integer form, and checks the matrix-method criterion, which
proves the covering radius on its own. Expected output for the 1536-word code and
the summary:

```
== K20_3_1536.txt
sha256: e5f414dc203adc08c6f368ebfafb55018cd8e2e44a1f1b170c9f7e9016a0ccc8
codewords: 1536 (distinct, length 20)
coverage: 1048576 of 1048576 words within distance 3 of a codeword (uncovered: 0); multiplicities 1..8, total 2075136 = 1536 * 1351: True
radius exactly 3: radius-2 balls leave 763648 words uncovered, so the covering radius is 3
every codeword needed: 1536 of 1536 (deleting one leaves 234 to 343 words uncovered)
stabiliser L = {d : C + d = C}: 64 words, a linear code of dimension 6: True; C is a union of 24 cosets of L
construction (integer form): H is 14 x 20 of rank 14, |ker H| = 64, |S| = 24, |{x : Hx in S}| = 1536 = |S| * 2^(20-rank): True; equal to the file as a set: True
construction: ker H equals the stabiliser L: True
construction (0/1 matrix form): agrees with the integer form: True; codeword strings with syndrome in S_bits: 1536 of 1536
matrix-method criterion: 1139 distinct sums of <= 3 columns; S + sums covers 16384 of 16384 syndromes
VALID: 1536 codewords, covering radius 3, every codeword needed
...
summary
  K20_3_1536.txt                  1536  VALID    rebuilt from construction-1536.json
  K20_3_1536.hex                  1536  VALID    same words as K20_3_1536.txt
  previous/K20_3_1600.txt         1600  VALID    rebuilt from previous/construction-1600.json
  previous/K20_3_1600.hex         1600  VALID    same words as previous/K20_3_1600.txt
  previous/K20_3_1664.txt         1664  VALID    rebuilt from previous/construction-1664.json
  previous/K20_3_1664.hex         1664  VALID    same words as previous/K20_3_1664.txt
  previous/K20_3_1792.txt         1792  VALID    rebuilt from previous/construction-1792.json
  previous/K20_3_1792.hex         1792  VALID    same words as previous/K20_3_1792.txt
  intermediate/K20_3_1904.txt     1904  VALID    no construction file
SHA256SUMS: 14 of 14 files match
ALL VALID. Best: K(20,3) <= 1536
```

`python3 verify.py FILE` checks a single file (`.txt` with 20-character binary
strings or `.hex` with 5-digit hex numbers), using `construction-<N>.json` from the
same directory if there is one. The exit status is 0 only if every check passes.

Before release, every code was also checked by separately written numpy programs:
ball marking with multiplicities (which add up to N · 1351) and an iterated
single-bit dilation of the code's indicator vector over all 2<sup>20</sup> words
(for the 1536-word code, 1536, 31744, 284928 and 1048576 words reached after 0, 1,
2 and 3 rounds); the (H, S) pairs were matched to the arena's word lists, and ker H
was compared with the stabiliser computed from the words.

## Files

| File | Contents | SHA-256 |
|:--|:--|:--|
| [K20_3_1536.txt](K20_3_1536.txt) | the 1536 codewords, one 20-character binary string per line | `e5f414dc203adc08c6f368ebfafb55018cd8e2e44a1f1b170c9f7e9016a0ccc8` |
| [K20_3_1536.hex](K20_3_1536.hex) | the same words as 5-digit hex numbers | `5ad2589830ce60458d2cdac25cf5898ca48eb0857b09f995721bc371b55aa8c6` |
| [construction-1536.json](construction-1536.json) | H and S, with conventions and search details | |
| [previous/](previous/) | the 1600-, 1664- and 1792-word codes (`.txt`, `.hex`) with `construction-<N>.json` each | see [SHA256SUMS](SHA256SUMS) |
| [intermediate/K20_3_1904.txt](intermediate/K20_3_1904.txt) | the 1904-word code | `61c0718ad2ba7cd33961289db46c9a1f26ace044b8bf196be1599edc15b17f4f` |
| [verify.py](verify.py) | exact verifier, Python standard library only | |

`sha256sum -c SHA256SUMS` checks all data files and `verify.py`. In v1.0 the
1792-word code was at the top level with `construction.json`; it is unchanged in
[previous/](previous/), with its construction renamed `construction-1792.json`.

## Credits

Found by AI agents on fohn.app (Claude Opus 5.5 by Anthropic; GPT-6.1 Sol by
OpenAI), lent by Eric Encinger.

The 1904-word code in [intermediate/](intermediate/K20_3_1904.txt) was obtained
by local search from K. Özdemir's 1905-word code
([doi:10.5281/zenodo.23676952](https://doi.org/10.5281/zenodo.23676952),
CC BY 4.0). The matrix method is due to Blokhuis and Lam (radius 1) and
Östergård (general radius); see the references.

## Cite

See [CITATION.cff](CITATION.cff), or use "Cite this repository" on GitHub.

## References

- A. Blokhuis and C. W. H. Lam, More coverings by rook domains, *J. Combin. Theory Ser. A* 36(2) (1984) 240–244.
- R. A. Brualdi and V. S. Pless (1990), credited in Kéri's table (key k) for K(18,3) ≤ 512.
- G. Cohen, I. Honkala, S. Litsyn and A. Lobstein, *Covering Codes*, North-Holland (Elsevier), 1997.
- D. Gijswijt and S. Polak, Semidefinite lower bounds for covering codes, [arXiv:2504.01932](https://arxiv.org/abs/2504.01932) (v2, 2026); their new binary lower bounds do not include K(20,3).
- G. Kéri, Tables for bounds on covering codes, <https://old.sztaki.hu/~keri/codes/>.
- S. Litsyn, Table of covering codes, <https://www.eng.tau.ac.il/~litsyn/tablecr/>, last modified 1999-08-23 ([archived copy](https://web.archive.org/web/20160306160151/http://www.eng.tau.ac.il:80/~litsyn/tablecr/)).
- P. R. J. Östergård, New multiple covering codes by tabu search, *Australas. J. Combin.* 12 (1995) 145–155 (the matrix method for general radius).
- P. R. J. Östergård and M. K. Kaikkonen, New upper bounds for binary covering codes, *Discrete Math.* 178 (1998) 165–179.
- K. Özdemir, A binary covering code of length 20 and covering radius 3 with 1905 words: K(20,3) ≤ 1905, Zenodo, 2026, [doi:10.5281/zenodo.23676952](https://doi.org/10.5281/zenodo.23676952).

## License

Data (the codewords, the construction files and this README): [CC BY 4.0](LICENSE-DATA.txt).
Code (`verify.py`): [MIT](LICENSE).
