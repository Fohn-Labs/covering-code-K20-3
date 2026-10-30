#!/usr/bin/env python3
"""Exact check that a binary code of length 20 has covering radius at most 3.

Usage:
    python3 verify.py                          # checks K20_3_1792.txt and its construction
    python3 verify.py intermediate/K20_3_1904.txt
    python3 verify.py K20_3_1792.hex           # hex lines (5 digits) also accepted

Python 3.6 or later, standard library only.

Check 1 (exhaustive): every one of the 2^20 = 1,048,576 binary words of length 20
is within Hamming distance 3 of some codeword. The radius-3 ball of each codeword
(1 + 20 + 190 + 1140 = 1351 words) is marked in a table of all 2^20 words, and the
unmarked words are counted. Negative controls follow: the same procedure must report
uncovered words for radius 2, and it reports how many words become uncovered when the
first codeword is deleted.

Check 2 (construction, run when construction.json sits next to this script and the
checked file is the 1792-word code): rebuild C = { x : H x in S } from the 12 x 20
parity-check matrix H and the 7 syndromes S, compare it with the file as sets, and
check the matrix-method criterion: every syndrome in F_2^12 is the sum of an element
of S and at most 3 columns of H. That criterion alone also proves covering radius
at most 3, independently of Check 1.
"""
import hashlib
import itertools
import json
import os
import sys

N = 20
R = 3
SPACE = 1 << N
HERE = os.path.dirname(os.path.abspath(__file__))


def read_code(path):
    """Return (sorted list of distinct 20-bit integers, number of lines, sha256)."""
    raw = open(path, 'rb').read()
    lines = [ln.strip() for ln in raw.decode('ascii').splitlines()]
    lines = [ln for ln in lines if ln and not ln.startswith('#')]
    words = []
    for ln in lines:
        if len(ln) == N and set(ln) <= set('01'):
            words.append(int(ln, 2))
        elif len(ln) == 5 and all(c in '0123456789abcdefABCDEF' for c in ln):
            words.append(int(ln, 16))
        else:
            sys.exit('unrecognised line: %r' % ln)
    if len(set(words)) != len(words):
        sys.exit('FAIL: %d duplicate lines' % (len(words) - len(set(words))))
    if any(not 0 <= w < SPACE for w in words):
        sys.exit('FAIL: word out of range')
    return sorted(words), len(lines), hashlib.sha256(raw).hexdigest()


def ball_masks(radius):
    """All error patterns of weight <= radius, as integers."""
    masks = []
    for k in range(radius + 1):
        for pos in itertools.combinations(range(N), k):
            m = 0
            for p in pos:
                m |= 1 << p
            masks.append(m)
    return masks


def uncovered_count(words, radius):
    """Number of words of {0,1}^20 at distance > radius from every codeword."""
    masks = ball_masks(radius)
    covered = bytearray(SPACE)
    for c in words:
        for m in masks:
            covered[c ^ m] = 1
    return covered.count(0)


def gf2_rank(vectors):
    basis = []  # kept reduced: each new vector is reduced by all previous ones
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def check_construction(words):
    cons = json.load(open(os.path.join(HERE, 'construction.json')))
    ok = True

    # Integer form: bit j of a word selects column H[j].
    H = cons['integer_form']['H_columns']
    S = cons['integer_form']['S']
    rank = gf2_rank(H)
    syn = [0] * SPACE
    for x in range(1, SPACE):
        low = x & -x
        syn[x] = syn[x ^ low] ^ H[low.bit_length() - 1]
    S_set = set(S)
    rebuilt = [x for x in range(SPACE) if syn[x] in S_set]
    kernel = syn.count(0)
    same = rebuilt == words
    print('construction (integer form): rank(H) = %d, |ker H| = %d, |{x : Hx in S}| = %d, '
          'equal to the file as a set: %s' % (rank, kernel, len(rebuilt), same))
    ok &= rank == 12 and kernel == 256 and len(rebuilt) == 1792 and same

    # Matrix form: 12 x 20 0/1 matrix acting on the 20-character strings, left to right.
    M = cons['matrix_form']['H_matrix']
    S_bits = set(cons['matrix_form']['S_bits'])
    assert len(M) == 12 and all(len(row) == N for row in M)
    def syndrome_bits(w):
        x = format(w, '020b')
        return ''.join(str(sum(row[k] * (x[k] == '1') for k in range(N)) % 2) for row in M)
    in_S = sum(1 for w in words if syndrome_bits(w) in S_bits)
    cols = [int(''.join(str(M[i][k]) for i in range(12)), 2) for k in range(N)]
    rank_m = gf2_rank(cols)
    # rank 12 => each syndrome has exactly 2^(20-12) = 256 preimages, so |{x : Hx in S}| = 7 * 256
    print('construction (0/1 matrix form): rank = %d, codewords with syndrome in S: %d of %d, '
          '|S| * 2^(20-rank) = %d' % (rank_m, in_S, len(words), len(S_bits) * 2 ** (N - rank_m)))
    ok &= rank_m == 12 and in_S == len(words) == len(S_bits) * 2 ** (N - rank_m)

    # Matrix-method criterion: S + {sums of <= 3 columns} = F_2^12.
    sums = set()
    for k in range(R + 1):
        for pos in itertools.combinations(range(N), k):
            s = 0
            for p in pos:
                s ^= H[p]
            sums.add(s)
    reach = {s ^ t for s in S for t in sums}
    print('matrix-method criterion: %d distinct sums of <= 3 columns; S + sums covers %d of 4096 '
          'syndromes' % (len(sums), len(reach)))
    ok &= len(reach) == 4096
    return ok


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'K20_3_1792.txt')
    words, nlines, digest = read_code(path)
    print('file: %s' % os.path.relpath(path))
    print('sha256: %s' % digest)
    print('codewords: %d (distinct, length %d)' % (len(words), N))

    bad = uncovered_count(words, R)
    print('uncovered words at radius 3: %d of %d' % (bad, SPACE))
    valid = bad == 0

    # Negative controls: the checker must find uncovered words in these cases.
    bad2 = uncovered_count(words, 2)
    bad_del = uncovered_count(words[1:], R)
    print('negative controls: radius 2 leaves %d uncovered; deleting the first codeword leaves %d '
          'uncovered at radius 3' % (bad2, bad_del))
    controls_ok = bad2 > 0  # the deletion count is informative only (a codeword may be redundant)

    cons_ok = True
    if len(words) == 1792 and os.path.exists(os.path.join(HERE, 'construction.json')):
        cons_ok = check_construction(words)

    if valid and controls_ok and cons_ok:
        print('VALID: K(20,3) <= %d (%d codewords, covering radius <= 3)' % (len(words), len(words)))
        return 0
    print('INVALID')
    return 1


if __name__ == '__main__':
    sys.exit(main())
