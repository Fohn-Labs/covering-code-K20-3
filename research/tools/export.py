"""Export a search result {r, H_columns, S} as K20_3_<N>.txt/.hex and construction-<N>.json
in the conventions of Fohn-Labs/covering-code-K20-3, then an independent coverage check.

usage: python3 export.py result.json outdir [method note]
"""
import json, os, sys, hashlib

N = 20


def gf2_rank(vs):
    basis = []
    for v in vs:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def main():
    res = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    note = sys.argv[3] if len(sys.argv) > 3 else ''
    r, H, S = res['r'], res['H_columns'], sorted(res['S'])
    assert len(H) == N and gf2_rank(H) == r and len(set(S)) == len(S)
    syn = [0] * (1 << N)
    for x in range(1, 1 << N):
        low = x & -x
        syn[x] = syn[x ^ low] ^ H[low.bit_length() - 1]
    Sset = set(S)
    words = [x for x in range(1 << N) if syn[x] in Sset]
    size = len(words)
    assert size == len(S) << (N - r)

    # independent check: three rounds of single-bit dilation of the indicator reach every word
    reach = bytearray(1 << N)
    for w in words:
        reach[w] = 1
    frontier = list(words)
    counts = [len(words)]
    for _ in range(3):
        nxt = []
        for x in frontier:
            for j in range(N):
                y = x ^ (1 << j)
                if not reach[y]:
                    reach[y] = 1
                    nxt.append(y)
        frontier = nxt
        counts.append(counts[-1] + len(nxt))
    print('dilation counts after 0..3 rounds:', counts, 'all covered:', counts[-1] == 1 << N)
    if counts[-1] != 1 << N:
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    strings = [format(w, '020b') for w in words]
    txt = ''.join(s + '\n' for s in strings)
    hx = ''.join(format(w, '05x') + '\n' for w in words)
    base = 'K20_3_%d' % size
    open(os.path.join(outdir, base + '.txt'), 'w').write(txt)
    open(os.path.join(outdir, base + '.hex'), 'w').write(hx)
    Hm = [[(H[N - 1 - k] >> (r - 1 - i)) & 1 for k in range(N)] for i in range(r)]
    S_bits = [format(s, '0%db' % r) for s in S]
    sums = set()
    import itertools
    for t in range(4):
        for pos in itertools.combinations(range(N), t):
            v = 0
            for p in pos:
                v ^= H[p]
            sums.add(v)
    cons = {
        'description': 'Matrix-method construction of a binary covering code of length 20, covering radius 3, '
                       'with %d words: C = { x in F_2^20 : H x in S }, a union of %d cosets of the [20,%d] linear '
                       'code ker H.' % (size, len(S), N - r),
        'n': N, 'R': 3, 'redundancy_r': r, 'kernel_dimension_k': N - r, 'number_of_syndromes': len(S),
        'size': size,
        'matrix_form': {
            'convention': 'Read a codeword as the 20-character string x_1..x_20 of %s.txt, left to right. '
                          'H_matrix is %d rows by 20 columns over F_2, column k (left to right) belonging to x_k. '
                          'The syndrome H x^T is a %d-bit column read top to bottom, written as a %d-character '
                          'string. C = { x : H x^T in S_bits }.' % (base, r, r, r),
            'H_matrix': Hm, 'S_bits': S_bits},
        'integer_form': {
            'convention': 'A codeword is the integer w whose 20-digit binary expansion (most significant bit '
                          'first) is the string in %s.txt and whose 5-digit hex form is the line in %s.hex. '
                          'Bit j of w (value 2^j, j = 0..19) selects H_columns[j]; the syndrome of w is the XOR '
                          'of the selected columns, a %d-bit integer. C = { w : syndrome(w) in S }. H_columns[j] '
                          'is column 20-j (1-based, left to right) of H_matrix read top to bottom as a binary '
                          'number; S[i] is S_bits[i] read as a binary number.' % (base, base, r),
            'H_columns': H, 'S': S},
        'facts_checked_by_verify_py': {
            'rank_H_over_F2': r, 'kernel_size': 1 << (N - r),
            'distinct_syndromes_from_at_most_3_columns': len(sums),
            'syndromes_covered_by_S_plus_at_most_3_columns': len({s ^ t for s in S for t in sums}),
        },
        'found_by': {'method': note},
        'files': {base + '.txt': hashlib.sha256(txt.encode()).hexdigest(),
                  base + '.hex': hashlib.sha256(hx.encode()).hexdigest()},
    }
    json.dump(cons, open(os.path.join(outdir, 'construction-%d.json' % size), 'w'), indent=1)
    print('wrote', outdir, base, 'size', size)


if __name__ == '__main__':
    main()
