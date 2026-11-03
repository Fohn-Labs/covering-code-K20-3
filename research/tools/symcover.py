"""Exact covering over orbits of a prescribed group K (a subgroup of the affine automorphism group of the
record 1536-word code).  Minimise |C| over K-invariant codes C with covering radius 3.

usage: symcover.py --perm-gens I,J,... [--trans-dim t] [--seed s] [--time T] [--target M] --out prefix
  --perm-gens: indices into the record's 384 coordinate automorphisms (each with its translation)
  --random-gens m: instead pick m random elements
  --trans-dim t: add a random t-dimensional subspace of ker H as pure translations
"""
import argparse, itertools, json, os, sys, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from affgroup import record_affine_group, orbits, closure, SPACE, N

HERE = os.path.dirname(os.path.abspath(__file__))
REC = os.path.join(HERE, '..', '..', 'construction-1536.json')
PERMS = os.path.join(HERE, '..', 'data', 'aut1536_perms.json')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--perm-gens', default='')
    ap.add_argument('--random-gens', type=int, default=0)
    ap.add_argument('--trans-dim', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--time', type=float, default=60)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--target', type=int, default=0, help='stop at first solution with size <= target')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    perms = json.load(open(PERMS))
    aff, L, C = record_affine_group(REC, perms)
    gens = []
    if a.perm_gens:
        gens += [aff[int(i)] for i in a.perm_gens.split(',')]
    if a.random_gens:
        gens += rng.sample(aff, a.random_gens)
    if a.trans_dim:
        # random t-dim subspace of L
        basis = []
        def rank(vs):
            bs = []
            for v in vs:
                for b in bs: v = min(v, v ^ b)
                if v: bs.append(v)
            return len(bs)
        Lnz = [x for x in L if x]
        while len(basis) < a.trans_dim:
            v = rng.choice(Lnz)
            if rank(basis + [v]) > len(basis): basis.append(v)
        gens += [(tuple(range(N)), v) for v in basis]
    t0 = time.time()
    if gens:
        grp = closure(gens)
        oid, norb = orbits(gens)
    else:
        grp = {(tuple(range(N)), 0)}
        oid, norb = list(range(SPACE)), SPACE
    size = [0] * norb
    rep = [-1] * norb
    for x in range(SPACE):
        size[oid[x]] += 1
        if rep[oid[x]] < 0: rep[oid[x]] = x
    print('group order', len(grp), 'orbits', norb, 'time %.1fs' % (time.time() - t0), flush=True)
    # record's orbits (must be a union of orbits if K <= Aut)
    Cset = set(C)
    rec_orbs = sorted(set(oid[x] for x in C))
    assert sum(size[o] for o in rec_orbs) == len(C), 'record is not K-invariant'
    # ball patterns
    pats = [0]
    for t in range(1, 4):
        for J in itertools.combinations(range(N), t):
            v = 0
            for j in J: v |= 1 << j
            pats.append(v)
    cover = []
    for o in range(norb):
        p = rep[o]
        cover.append(sorted(set(oid[p ^ e] for e in pats)))
    print('built constraints %.1fs' % (time.time() - t0), flush=True)
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    xv = [m.NewBoolVar('x%d' % i) for i in range(norb)]
    for cv in cover:
        m.AddBoolOr([xv[i] for i in cv])
    obj = sum(size[i] * xv[i] for i in range(norb))
    m.Minimize(obj)
    if a.target:
        m.Add(obj <= a.target)
    recset = set(rec_orbs)
    if not a.target:
        for i in range(norb):
            m.AddHint(xv[i], i in recset)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = a.time
    s.parameters.num_workers = a.workers
    s.parameters.random_seed = a.seed

    class CB(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.best = None
        def on_solution_callback(self):
            val = int(self.ObjectiveValue())
            print('  solution %d (bound %d) at %.1fs' % (val, int(self.BestObjectiveBound()), time.time() - t0), flush=True)
            if self.best is None or val < self.best[0]:
                self.best = (val, [i for i in range(norb) if self.Value(xv[i])])
                words = sorted(x for x in range(SPACE) if oid[x] in set(self.best[1]))
                open(a.out + '_%d.hex' % val, 'w').write(''.join(format(w, '05x') + '\n' for w in words))

    cb = CB()
    st = s.Solve(m, cb)
    print('status', s.StatusName(st), 'best', cb.best[0] if cb.best else None, 'bound', s.BestObjectiveBound(),
          'time %.1fs' % (time.time() - t0), flush=True)
    json.dump({'gens': [[list(g[0]), g[1]] for g in gens], 'group_order': len(grp), 'orbits': norb,
               'status': s.StatusName(st), 'best': cb.best[0] if cb.best else None,
               'bound': s.BestObjectiveBound()}, open(a.out + '.json', 'w'))


if __name__ == '__main__':
    main()
