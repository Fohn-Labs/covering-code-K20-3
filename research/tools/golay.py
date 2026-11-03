"""Golay-derived kernels for the matrix method at n=20, R=3 (stdlib only)."""
import itertools, sys

def golay23():
    # cyclic [23,12,7] code, generator g(x) = x^11+x^10+x^6+x^5+x^4+x^2+1
    g = (1<<11)|(1<<10)|(1<<6)|(1<<5)|(1<<4)|(1<<2)|1
    rows = [g << i for i in range(12)]
    return rows  # ints with 23 bits

def span(rows):
    out = [0]
    for r in rows:
        out += [x ^ r for x in out]
    return out

def golay24():
    rows = golay23()
    return [(r << 1) | (bin(r).count('1') & 1) for r in rows]  # parity bit at bit 0

def rref(rows, nbits):
    rows = list(rows); piv = []; out = []
    for bit in reversed(range(nbits)):
        idx = next((i for i, r in enumerate(rows) if (r >> bit) & 1), None)
        if idx is None: continue
        p = rows.pop(idx)
        rows = [r ^ p if (r >> bit) & 1 else r for r in rows]
        out = [o ^ p if (o >> bit) & 1 else o for o in out]
        out.append(p); piv.append(bit)
    return out, piv

def kernel_code(rows, n, drop):
    """shorten on bit positions in drop (codewords zero there), then delete those positions."""
    words = [w for w in span(rows) if all(not (w >> d) & 1 for d in drop)]
    keep = [b for b in range(n) if b not in drop]
    def proj(w): return sum(((w >> b) & 1) << i for i, b in enumerate(keep))
    return sorted(set(proj(w) for w in words))

def parity_check(code_words, n):
    """columns of a parity-check matrix H (r x n) with ker H = code; returns list of n ints (r bits)."""
    basis, _ = rref([w for w in code_words if w], n)
    k = len(basis)
    # dual space: all v with <v, c> = 0 for c in basis
    # solve via brute force on structure: build dual basis by rref on generator
    # use standard method: put generator into systematic form on pivot columns
    gen, piv = rref(basis, n)
    nonpiv = [b for b in range(n) if b not in piv]
    dual = []
    for j in nonpiv:
        v = 1 << j
        for g, p in zip(gen, piv):
            if (g >> j) & 1: v |= 1 << p
        dual.append(v)
    r = len(dual)
    cols = [sum(((d >> b) & 1) << i for i, d in enumerate(dual)) for b in range(n)]
    for c in code_words:
        assert all(bin(c & d).count('1') % 2 == 0 for d in dual)
    return cols, r

def ball3(cols, r):
    s = set([0])
    for t in range(1, 4):
        for comb in itertools.combinations(cols, t):
            x = 0
            for c in comb: x ^= c
            s.add(x)
    return s

if __name__ == '__main__':
    g23 = golay23(); w23 = span(g23)
    md = min(bin(w).count('1') for w in w23 if w)
    print('G23 words', len(w23), 'min dist', md)
    g24 = golay24(); w24 = span(g24)
    print('G24 min dist', min(bin(w).count('1') for w in w24 if w))
    # r=11: shorten G23 on 3 coords
    for name, rows, n, drop in [('G23 shortened x3 [20,9]', g23, 23, [0, 1, 2]),
                                ('G24 shortened x4 [20,8]', g24, 24, [0, 1, 2, 3])]:
        code = kernel_code(rows, n, drop)
        cols, r = parity_check(code, 20)
        B = ball3(cols, r)
        print(name, 'dim', len(code).bit_length()-1, 'r', r, 'min dist', min(bin(w).count('1') for w in code if w), '|B3|', len(B), 'of', 1 << r)
        if r == 11:
            comp = set(range(1 << r)) - B
            best = max(len(comp & set(s ^ b for b in B)) for s in range(1 << r))
            print('  |S|=2: best coverage of complement', best, 'of', len(comp))
