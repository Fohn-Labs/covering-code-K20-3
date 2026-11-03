"""Coordinate permutations sigma with sigma(C) = C + c for the matrix-method code C = {x : Hx in S}."""
import json, sys, itertools
from collections import Counter
def load(path):
    f = json.load(open(path))['integer_form']; return f['H_columns'], f['S']
def kernel(H, n=20):
    syn = [0] * (1 << n)
    for x in range(1, 1 << n):
        low = x & -x; syn[x] = syn[x ^ low] ^ H[low.bit_length() - 1]
    return syn, [x for x in range(1 << n) if syn[x] == 0]
def perm_word(x, p):
    y = 0
    for i in range(20):
        if (x >> i) & 1: y |= 1 << p[i]
    return y
def aut_of_code(words, n=20):
    W = set(words)
    sup = [[i for i in range(n) if (w >> i) & 1] for w in words]
    inv1 = [tuple(sorted(Counter(len(s) for s in sup if i in s).items())) for i in range(n)]
    inv2 = [[tuple(sorted(Counter(len(s) for s in sup if i in s and j in s).items())) for j in range(n)] for i in range(n)]
    res = []
    p = [-1] * n; used = [False] * n
    def bt(i):
        if i == n:
            if all(perm_word(w, p) in W for w in words): res.append(tuple(p))
            return
        for t in range(n):
            if used[t] or inv1[t] != inv1[i]: continue
            if any(inv2[i][j] != inv2[t][p[j]] for j in range(i)): continue
            p[i] = t; used[t] = True; bt(i + 1); used[t] = False; p[i] = -1
    bt(0)
    return res
if __name__ == '__main__':
    H, S = load(sys.argv[1])
    syn, L = kernel(H)
    autL = aut_of_code(L)
    print('|ker H| =', len(L), ' |Aut(ker H)| =', len(autL))
    Sset = set(S)
    C = [x for x in range(1 << 20) if syn[x] in Sset]
    Cset = set(C)
    good = []
    for p in autL:
        img = [perm_word(x, p) for x in C]
        c = img[0] ^ C[0]
        # sigma(C) = C + c' for some c'; test all c' = img[0] ^ y for y in C would be slow; use syndromes:
        # sigma(C) + c' = C  <=>  {syn(sigma x)} = S + syn(c')
        simg = set(syn[y] for y in img)
        ok = any(set(s ^ d for s in simg) == Sset for d in set(s0 ^ S[0] for s0 in simg))
        if ok: good.append(p)
    print('|Aut(C)| restricted to coordinates (with translations) =', len(good))
    from collections import Counter as Ct
    def order(p):
        q = list(range(20)); k = 0
        while True:
            q = [p[i] for i in q]; k += 1
            if q == list(range(20)): return k
    print('element orders', sorted(Ct(order(p) for p in good).items()))
    json.dump(good, open(sys.argv[2], 'w'))
