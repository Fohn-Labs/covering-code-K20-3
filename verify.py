#!/usr/bin/env python3
"""Exact checks for the binary covering codes of length 20 and covering radius 3 in this
repository.

Usage:
    python3 verify.py                     # every code file in the repository (about 15 s)
    python3 verify.py FILE [FILE ...]     # the given files (.txt: 20-character binary
                                          # strings; .hex: 5-digit hex numbers)

Python 3.6 or later, standard library only.

For each code file C (M codewords):

1. Coverage, exhaustively: the radius-3 ball of every codeword (1 + 20 + 190 + 1140 = 1351
   words) is marked in a table of all 2^20 = 1,048,576 words, counting multiplicities, and
   every word must be marked. The multiplicities must add up to M * 1351.
2. Radius exactly 3: the same marking with radius-2 balls must leave words unmarked.
3. Every codeword is needed: each codeword's ball must contain words marked once (covered by
   no other codeword), so deleting any codeword leaves words uncovered.
4. Stabiliser: L = { d : C + d = C } is computed from the words. C is a union of cosets of the
   linear code L (for the 1904-word code, L = {0}).
5. Construction, when construction-<M>.json sits next to the file: C = { x : H x in S } is
   rebuilt by computing the syndrome of all 2^20 words and compared with the file as a set;
   rank(H), |ker H| and |S| * 2^(20 - rank) = M are checked, ker H must equal L, the 0/1
   matrix form must agree with the integer form and send every codeword string into S, and
   the matrix-method criterion is checked: every syndrome in F_2^r is s + (sum of at most 3
   columns of H) for some s in S. That criterion proves covering radius <= 3 on its own,
   independently of check 1 (flip the at most 3 coordinates of x whose columns, added to
   H x, give an element of S).

Bit convention (also in each construction file): a codeword is the integer w whose 20-digit
binary expansion, most significant bit first, is the string in the .txt file; bit j of w
(value 2^j) selects H_columns[j], so character k (1-based, left to right) selects column k
of H_matrix, which is H_columns[20 - k].
"""
import glob
import hashlib
import itertools
import json
import os
import sys
from array import array

N = 20
R = 3
SPACE = 1 << N
HERE = os.path.dirname(os.path.abspath(__file__))


def read_code(path):
    """Return (sorted list of distinct 20-bit integers, list of 20-character strings, sha256)."""
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
            raise SystemExit('unrecognised line in %s: %r' % (path, ln))
    if len(set(words)) != len(words):
        raise SystemExit('FAIL: %s has %d duplicate lines' % (path, len(words) - len(set(words))))
    if any(not 0 <= w < SPACE for w in words):
        raise SystemExit('FAIL: %s has a word out of range' % path)
    words.sort()
    return words, [format(w, '020b') for w in words], hashlib.sha256(raw).hexdigest()


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


MASKS3 = ball_masks(3)
MASKS2 = ball_masks(2)


def gf2_rank(vectors):
    basis = []  # each new vector is reduced by all previous basis vectors
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def check_words(words):
    """Checks 1-4. Returns (ok, report lines, facts)."""
    out = []
    count = array('H', bytes(2 * SPACE))  # multiplicities; at most 1351 < 65536
    for c in words:
        for m in MASKS3:
            count[c ^ m] += 1
    uncovered = count.count(0)
    total = sum(count)
    mult = [x for x in set(count)]
    out.append('coverage: %d of %d words within distance 3 of a codeword (uncovered: %d); '
               'multiplicities %d..%d, total %d = %d * 1351: %s'
               % (SPACE - uncovered, SPACE, uncovered, min(mult), max(mult), total, len(words),
                  total == len(words) * 1351))

    marked2 = bytearray(SPACE)
    for c in words:
        for m in MASKS2:
            marked2[c ^ m] = 1
    unc2 = marked2.count(0)
    radius = ('more than 3' if uncovered else '3' if unc2 else 'at most 2')
    out.append('radius exactly 3: radius-2 balls leave %d words uncovered, so the covering '
               'radius is %s' % (unc2, radius))

    private = [sum(1 for m in MASKS3 if count[c ^ m] == 1) for c in words]
    needed = sum(1 for p in private if p > 0)
    out.append('every codeword needed: %d of %d (deleting one leaves %d to %d words uncovered)'
               % (needed, len(words), min(private), max(private)))

    member = bytearray(SPACE)
    for c in words:
        member[c] = 1
    c0 = words[0]
    stab = sorted(d for d in (c ^ c0 for c in words) if all(member[c ^ d] for c in words))
    k = gf2_rank(stab)
    linear = len(stab) == 1 << k
    out.append('stabiliser L = {d : C + d = C}: %d word%s, a linear code of dimension %d: %s; '
               'C is a union of %d cosets of L' % (len(stab), '' if len(stab) == 1 else 's', k, linear,
                                                   len(words) // len(stab)))

    ok = (uncovered == 0 and total == len(words) * 1351 and unc2 > 0 and needed == len(words)
          and linear)
    return ok, out, {'stabiliser': stab}


def check_construction(path, words, strings, stab):
    """Check 5 against construction-<M>.json. Returns (ok, report lines)."""
    cons = json.load(open(path, encoding='utf-8'))
    out = []
    H = cons['integer_form']['H_columns']
    S = cons['integer_form']['S']
    r = cons['redundancy_r']
    if len(H) != N or len(set(S)) != len(S) or any(not 0 <= v < 1 << r for v in H + S):
        return False, ['construction: malformed H or S']

    # Rebuild from H and S over all 2^20 words.
    syn = [0] * SPACE
    for x in range(1, SPACE):
        low = x & -x
        syn[x] = syn[x ^ low] ^ H[low.bit_length() - 1]
    S_set = set(S)
    rebuilt = [x for x in range(SPACE) if syn[x] in S_set]
    kernel = [x for x in range(SPACE) if syn[x] == 0]
    rank = gf2_rank(H)
    same = rebuilt == words
    out.append('construction (integer form): H is %d x 20 of rank %d, |ker H| = %d, |S| = %d, '
               '|{x : Hx in S}| = %d = |S| * 2^(20-rank): %s; equal to the file as a set: %s'
               % (r, rank, len(kernel), len(S), len(rebuilt),
                  len(rebuilt) == len(S) << (N - rank), same))
    ker_is_stab = kernel == stab
    out.append('construction: ker H equals the stabiliser L: %s' % ker_is_stab)
    ok = (rank == r and len(kernel) == 1 << (N - r) and same and ker_is_stab
          and len(rebuilt) == len(S) << (N - rank) == cons['size'] == len(words))

    # 0/1 matrix form, acting on the 20-character strings left to right.
    Mx = cons['matrix_form']['H_matrix']
    S_bits = cons['matrix_form']['S_bits']
    shape_ok = len(Mx) == r and all(len(row) == N for row in Mx) and len(S_bits) == len(S)
    cols = [int(''.join(str(Mx[i][k]) for i in range(r)), 2) for k in range(N)] if shape_ok else []
    forms_agree = (shape_ok and cols == [H[N - 1 - k] for k in range(N)]
                   and [int(s, 2) for s in S_bits] == S)
    in_S = 0
    S_bits_set = set(S_bits)
    for x in strings:
        ones = [k for k in range(N) if x[k] == '1']
        s = ''.join(str(sum(Mx[i][k] for k in ones) % 2) for i in range(r))
        in_S += s in S_bits_set
    out.append('construction (0/1 matrix form): agrees with the integer form: %s; codeword strings '
               'with syndrome in S_bits: %d of %d' % (forms_agree, in_S, len(strings)))
    ok &= forms_agree and in_S == len(strings)

    # Matrix-method criterion: S + {sums of <= 3 columns} = F_2^r.
    sums = set()
    for t in range(R + 1):
        for pos in itertools.combinations(range(N), t):
            v = 0
            for p in pos:
                v ^= H[p]
            sums.add(v)
    reach = {s ^ t for s in S for t in sums}
    out.append('matrix-method criterion: %d distinct sums of <= 3 columns; S + sums covers %d of %d '
               'syndromes' % (len(sums), len(reach), 1 << r))
    ok &= len(reach) == 1 << r
    return ok, out


def default_files():
    files = []
    for sub in ('', 'previous', 'intermediate'):
        found = [f for f in glob.glob(os.path.join(HERE, sub, 'K20_3_*')) if f.endswith(('.txt', '.hex'))]
        files += sorted(found, key=lambda f: (os.path.basename(f)[:-4], f.endswith('.hex')))
    return files


def check_sums():
    path = os.path.join(HERE, 'SHA256SUMS')
    if not os.path.exists(path):
        return
    good, bad = 0, []
    for ln in open(path):
        if not ln.strip():
            continue
        digest, name = ln.split(None, 1)
        name = name.strip().lstrip('*')
        f = os.path.join(HERE, name)
        if os.path.exists(f) and hashlib.sha256(open(f, 'rb').read()).hexdigest() == digest:
            good += 1
        else:
            bad.append(name)
    print('SHA256SUMS: %d of %d files match%s' % (good, good + len(bad),
          '' if not bad else ' (differ: %s; line-ending conversion on checkout changes hashes, '
          'not codewords)' % ', '.join(bad)))


def main():
    files = sys.argv[1:] or default_files()
    if not files:
        raise SystemExit('no code files found')
    seen = {}
    summary = []
    all_ok = True
    for path in files:
        words, strings, digest = read_code(path)
        rel = os.path.relpath(os.path.abspath(path), HERE)
        print('== %s' % rel)
        print('sha256: %s' % digest)
        print('codewords: %d (distinct, length %d)' % (len(words), N))
        key = tuple(words)
        if key in seen:
            ok, note = seen[key]
            print('same codeword set as %s, which is checked above: %s' % (note, 'VALID' if ok else 'INVALID'))
            summary.append((rel, len(words), ok, 'same words as ' + note))
            all_ok &= ok
            print()
            continue
        ok, lines, facts = check_words(words)
        for ln in lines:
            print(ln)
        cons = os.path.join(os.path.dirname(os.path.abspath(path)), 'construction-%d.json' % len(words))
        if os.path.exists(cons):
            c_ok, c_lines = check_construction(cons, words, strings, facts['stabiliser'])
            for ln in c_lines:
                print(ln)
            ok &= c_ok
            how = 'rebuilt from ' + os.path.relpath(cons, HERE)
        else:
            how = 'no construction file'
            print('construction: none (no construction-%d.json next to this file)' % len(words))
        print('VALID: %d codewords, covering radius 3, every codeword needed' % len(words) if ok
              else 'INVALID: %d codewords' % len(words))
        print()
        seen[key] = (ok, rel)
        summary.append((rel, len(words), ok, how))
        all_ok &= ok
    print('summary')
    for rel, m, ok, how in summary:
        print('  %-30s %5d  %-7s  %s' % (rel, m, 'VALID' if ok else 'INVALID', how))
    if not sys.argv[1:]:
        check_sums()
    if all_ok:
        best = min(m for _, m, _, _ in summary)
        print('ALL VALID. Best: K(20,3) <= %d' % best)
        return 0
    print('INVALID')
    return 1


if __name__ == '__main__':
    sys.exit(main())
