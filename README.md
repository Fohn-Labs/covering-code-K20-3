# K(20,3) ≤ 1792

**A binary code of length 20 with 1792 codewords and covering radius 3, so K(20,3) ≤ 1792.**

K(n, R) is the smallest number of binary words of length n such that every binary
word of length n is within Hamming distance R of at least one of them. This
repository gives an explicit 1792-word code for n = 20, R = 3, a short algebraic
description of it, and a standard-library Python script that checks all
2<sup>20</sup> = 1,048,576 words exactly.

| Bound on K(20,3) | Value | Source |
|:--|--:|:--|
| Upper bound, this code | **1792** | this repository (found 2026-10-29) |
| Previous upper bound | 1905 | K. Özdemir, Zenodo, 2026-10-26, [doi:10.5281/zenodo.23676952](https://doi.org/10.5281/zenodo.23676952) |
| Upper bound in the standard tables | 2048 | G. Kéri, [Tables for bounds on covering codes](https://old.sztaki.hu/~keri/codes/) (binary table [2_tables.pdf](https://old.sztaki.hu/~keri/codes/2_tables.pdf), dated 2009-10-15); S. Litsyn, [Table of covering codes](https://web.archive.org/web/20160306160151/http://www.eng.tau.ac.il:80/~litsyn/tablecr/) (last modified 1999-08-23) |
| Lower bound | 892 | W. Haas (2007–2008), as credited in Kéri's table (Litsyn's 1999 table lists 889) |
| Sphere-covering bound | 777 | ⌈2<sup>20</sup> / 1351⌉ |

## Construction

The code is built with the matrix method (Blokhuis and Lam for radius 1,
generalised by Östergård, see the references). Let H be the 12 × 20 binary
matrix below and S the set of 7 syndromes below. Then

C = { x ∈ F<sub>2</sub><sup>20</sup> : H x<sup>T</sup> ∈ S }.

H has rank 12, so its kernel is a linear [20, 8] code with 256 words, and C is
the union of 7 of its cosets: |C| = 7 · 2<sup>8</sup> = 1792.

C has covering radius 3 because every one of the 4096 syndromes in
F<sub>2</sub><sup>12</sup> can be written as s + h<sub>1</sub> + … + h<sub>t</sub>
with s ∈ S and t ≤ 3 columns h<sub>1</sub>, …, h<sub>t</sub> of H. (The sums of at most 3
columns give 1164 distinct syndromes, and their 7 translates by S cover all 4096.)
Given any word x, pick such an s and columns for the syndrome H x<sup>T</sup>;
flipping the matching t ≤ 3 coordinates of x gives a word with syndrome s, which
is a codeword.

Columns are in the order of the characters of the strings in
[K20_3_1792.txt](K20_3_1792.txt), left to right; syndromes are read top to bottom.

```
H =
0 0 0 0 1 0 1 0 1 0 1 0 0 1 0 1 0 1 1 0
1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 0
1 0 1 0 0 0 0 1 0 0 1 0 0 0 0 0 0 1 0 0
0 0 0 1 0 0 1 1 0 1 1 0 1 0 1 1 1 0 1 1
1 0 0 1 1 1 1 0 0 0 1 0 0 1 0 1 0 0 0 0
1 1 1 0 0 0 1 1 1 0 1 0 0 0 0 0 0 0 0 0
0 0 0 1 1 0 0 1 0 1 1 0 1 1 0 0 1 1 0 0
1 0 1 0 1 1 0 0 0 0 0 0 1 0 1 0 0 0 0 0
1 0 1 0 1 0 1 0 1 1 1 0 1 1 1 1 0 1 0 0
0 1 1 0 1 0 0 0 1 0 1 1 0 0 1 0 0 1 1 1
0 1 1 1 1 1 1 1 1 1 0 1 0 0 1 0 0 0 1 0
0 1 0 0 1 0 1 1 0 0 0 0 1 1 0 0 1 1 1 1

S = { 010001111001, 011011000110, 111001001010, 100010110111,
      101111010100, 100111110100, 010001111110 }
```

[construction.json](construction.json) gives the same H and S as a 0/1 matrix and
as 12-bit integers (columns `[261, 3335, 3629, 1313, 2440, 1310, 3241, 1337, 6,
4076, 1322, 3150, 1891, 3531, 1170, 3263, 1442, 1630, 1095, 1752]`, where bit j of a
codeword's integer value selects column j, and S = `[1145, 1734, 3658, 2231, 3028,
2548, 1150]`), with the conventions spelled out.

Other properties, all checked exhaustively: every codeword is needed (deleting any
one leaves between 88 and 222 words uncovered); the minimum distance is 2; each
word of F<sub>2</sub><sup>20</sup> is within distance 3 of between 1 and 8
codewords. The covering density is 1792 · 1351 / 2<sup>20</sup> ≈ 2.309; for
comparison, the best known code of length 21 and radius 3 (3072 words,
Östergård and Kaikkonen 1998) has density ≈ 2.288. With the direct sum and
K(n+1, R+1) ≤ K(n, R), the new bound does not improve any other entry of Kéri's
binary table.

## How it was found

The code was found on [fohn.app](https://fohn.app), task #10, part #58 (matrix
method), by a Claude Opus 5.5 agent (Anthropic) running in Claude Code. It ran a
tabu search over pairs (H, S) with 12 rows and 7 syndromes; the run with seed 107
reached full coverage at iteration 17977, after 108 seconds.

Özdemir's record mentions a 1792-word matrix-method base that needed 196 extra
words; the code here is a complete matrix-method code of that size, with nothing
added.

Independently, on part #32 of the same task, a GPT-6.1 Sol agent (OpenAI)
running in Codex found a 1904-word code by local search from Özdemir's
1905-word code: [intermediate/K20_3_1904.txt](intermediate/K20_3_1904.txt), also
verified by `verify.py`.

## Verify

Python 3.6 or later, no packages needed:

```
python3 verify.py
```

Expected output (about a second):

```
file: K20_3_1792.txt
sha256: 19ab0796dcc3244e718da3aa48736179930eb2e76817218eef1e7771a3f838ea
codewords: 1792 (distinct, length 20)
uncovered words at radius 3: 0 of 1048576
negative controls: radius 2 leaves 729600 uncovered; deleting the first codeword leaves 125 uncovered at radius 3
construction (integer form): rank(H) = 12, |ker H| = 256, |{x : Hx in S}| = 1792, equal to the file as a set: True
construction (0/1 matrix form): rank = 12, codewords with syndrome in S: 1792 of 1792, |S| * 2^(20-rank) = 1792
matrix-method criterion: 1164 distinct sums of <= 3 columns; S + sums covers 4096 of 4096 syndromes
VALID: K(20,3) <= 1792 (1792 codewords, covering radius <= 3)
```

The first check marks the radius-3 ball (1351 words) around every codeword in a
table of all 2<sup>20</sup> words and counts what is left. The second rebuilds the
code from H and S, compares it with the file as sets, and checks the syndrome
criterion above, which proves the covering radius on its own. `python3 verify.py
intermediate/K20_3_1904.txt` checks the 1904-word code; `.hex` files are also
accepted.

Before release, the code was also checked by separately written programs: an
iterated single-bit dilation of the code's indicator vector over all 2<sup>20</sup>
words (full coverage after 3 rounds, 729600 words still uncovered after 2), and
ball-marking checkers in pure Python and in numpy, one of which also confirms that
the ball multiplicities add up to 1792 · 1351.

## Files

| File | Contents | SHA-256 |
|:--|:--|:--|
| [K20_3_1792.txt](K20_3_1792.txt) | the 1792 codewords, one 20-character binary string per line | `19ab0796dcc3244e718da3aa48736179930eb2e76817218eef1e7771a3f838ea` |
| [K20_3_1792.hex](K20_3_1792.hex) | the same words as 5-digit hex numbers | `6142242e7b0384404b13836bef801dd08da843a31517d846e755120b65643d7d` |
| [construction.json](construction.json) | H and S, with conventions | |
| [verify.py](verify.py) | exact verifier, Python standard library only | |
| [intermediate/K20_3_1904.txt](intermediate/K20_3_1904.txt) | the independent 1904-word code | `61c0718ad2ba7cd33961289db46c9a1f26ace044b8bf196be1599edc15b17f4f` |

`sha256sum -c SHA256SUMS` checks all files.

## Credits

Found by AI agents on fohn.app (Claude Opus 5.5 by Anthropic; GPT-6.1 Sol by
OpenAI) lent by Eric Encinger.

## Cite

See [CITATION.cff](CITATION.cff), or use "Cite this repository" on GitHub.

## References

- A. Blokhuis and C. W. H. Lam, More coverings by rook domains, *J. Combin. Theory Ser. A* 36(2) (1984) 240–244.
- G. Cohen, I. Honkala, S. Litsyn and A. Lobstein, *Covering Codes*, North-Holland (Elsevier), 1997.
- D. Gijswijt and S. Polak, Semidefinite lower bounds for covering codes, [arXiv:2504.01932](https://arxiv.org/abs/2504.01932) (v2, 2026); their new binary lower bounds do not include K(20,3).
- G. Kéri, Tables for bounds on covering codes, <https://old.sztaki.hu/~keri/codes/>.
- S. Litsyn, Table of covering codes, <https://www.eng.tau.ac.il/~litsyn/tablecr/>, last modified 1999-08-23 ([archived copy](https://web.archive.org/web/20160306160151/http://www.eng.tau.ac.il:80/~litsyn/tablecr/)).
- P. R. J. Östergård, New multiple covering codes by tabu search, *Australas. J. Combin.* 12 (1995) 145–155 (the matrix method for general radius).
- P. R. J. Östergård and M. K. Kaikkonen, New upper bounds for binary covering codes, *Discrete Math.* 178 (1998) 165–179.
- K. Özdemir, A binary covering code of length 20 and covering radius 3 with 1905 words: K(20,3) ≤ 1905, Zenodo, 2026, [doi:10.5281/zenodo.23676952](https://doi.org/10.5281/zenodo.23676952).

## License

Data (the codewords, construction.json and this README): [CC BY 4.0](LICENSE-DATA.txt).
Code (`verify.py`): [MIT](LICENSE).
