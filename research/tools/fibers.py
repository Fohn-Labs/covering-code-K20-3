"""Change basis of F_2^r so that span(S - S) = first coordinates; report fiber sizes."""
import json, itertools, sys
from collections import Counter, defaultdict
def reduce_basis(vs):
    basis = []
    for v in vs:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis
def fibers(H, S, r):
    s0 = S[0]; V = reduce_basis([s ^ s0 for s in S]); d = len(V)
    # complete V to a basis of F_2^r, then express coordinates: new coords = (V-part, complement-part)
    full = list(V)
    for e in range(r):
        if len(reduce_basis(full + [1 << e])) > len(full): full.append(1 << e)
    # solve x = sum a_i full_i ; build inverse by gaussian elimination on 2^r? r <= 16 so brute table
    coord = {}
    for m in range(1 << r):
        x = 0
        for i in range(r):
            if (m >> i) & 1: x ^= full[i]
        coord[x] = m
    lowmask = (1 << d) - 1
    Ht = [coord[h] for h in H]
    fib = defaultdict(set)
    for t in range(4):
        for J in itertools.combinations(range(20), t):
            v = 0
            for j in J: v ^= Ht[j]
            fib[v >> d].add(v & lowmask)
    sizes = sorted(len(fib[y]) for y in range(1 << (r - d)))
    Tset = sorted(coord[s ^ s0] for s in S)
    assert all(t >> d == 0 for t in Tset)
    return d, sizes, [h >> d for h in Ht], Tset
if __name__ == '__main__':
    import os
    R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..') + '/'
    for name in ['previous/construction-1792.json', 'previous/construction-1664.json', 'previous/construction-1600.json', 'construction-1536.json']:
        f = json.load(open(R + name))['integer_form']; H, S = f['H_columns'], f['S']
        r = {'1792': 12, '1664': 13, '1600': 14, '1536': 14}[name.split('-')[1][:4]]
        d, sizes, bot, T = fibers(H, S, r)
        print(name, 'r', r, 'dim span(S-S)', d, 'fibers', 1 << (r - d), 'fiber sizes min/max', sizes[0], sizes[-1], 'LP bound |T| >=', round((1 << d) / sizes[0], 2), '|S|', len(S))
        print('   sizes', sizes)
        print('   H_bot columns', Counter(bot))
