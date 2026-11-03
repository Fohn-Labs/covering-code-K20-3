/* Joint weighted tabu search over (H, S) for the matrix method, n = 20, R = 3.
 * C = { x in F_2^20 : H x in S }, H is r x 20 of rank r, |C| = |S| 2^(20-r).
 * Covering radius <= 3 iff S + B = F_2^r, B = sums of at most 3 columns of H.
 *
 * Moves: replace one element of S (all pairs scored exactly), or replace one
 * column of H by any value (all 2^r values scored exactly by a Walsh-Hadamard
 * correlation).  Scores use weights on uncovered syndromes that grow each
 * iteration.  Rank of H is kept equal to r.
 *
 * usage: jtabu -r r -k k [-s seed] [-t seconds] [-o outprefix] [-h hcols_per_iter]
 *              [-f startfile] [-d d] [-b]
 *   -f: text file "r c0 .. c19 [s0 ..]" to start from (S filled up at random)
 *   -d: restrict S to F_2^d x {0} (the low d bits); H moves then keep the high bits
 *       of each column unless -b is given
 *   -h 0: S-only search with H fixed
 * Writes outprefix.json with the best (H, S) found.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>

#define NC 20
static int r, N, k, nb;
static int cols[NC];
static uint8_t *inB;
static int *B;
static int *S;
static uint8_t *inS;
static uint16_t *cnt;
static int64_t *w;
static int *U, *posU, nU;
static int *tabuS;           /* per syndrome value */
static long coltabu[NC];

static uint64_t rs = 88172645463325252ULL;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static void addU(int x) { posU[x] = nU; U[nU++] = x; }
static void delU(int x) { int p = posU[x], y = U[--nU]; U[p] = y; posU[y] = p; posU[x] = -1; }
static void add_s(int s) { for (int i = 0; i < nb; i++) { int x = s ^ B[i]; if (cnt[x]++ == 0) delU(x); } inS[s] = 1; }
static void rem_s(int s) { for (int i = 0; i < nb; i++) { int x = s ^ B[i]; if (--cnt[x] == 0) addU(x); } inS[s] = 0; }

static void build_B(void) {
    memset(inB, 0, N); inB[0] = 1;
    for (int a = 0; a < NC; a++) { inB[cols[a]] = 1;
        for (int b = a + 1; b < NC; b++) { int ab = cols[a] ^ cols[b]; inB[ab] = 1;
            for (int c = b + 1; c < NC; c++) inB[ab ^ cols[c]] = 1; } }
    nb = 0; for (int x = 0; x < N; x++) if (inB[x]) B[nb++] = x;
}
static void recompute(void) {
    build_B();
    memset(cnt, 0, N * sizeof(uint16_t));
    for (int i = 0; i < k; i++) for (int t = 0; t < nb; t++) cnt[S[i] ^ B[t]]++;
    nU = 0; for (int x = 0; x < N; x++) { posU[x] = -1; if (!cnt[x]) addU(x); }
}
static int rank_of(const int *v, int m, int *normal) {
    /* rank of m vectors in F_2^r; if rank == r-1 and normal != NULL, return a normal vector h with h.v = 0 */
    int basis[32]; int piv[32]; int nbasis = 0;
    for (int i = 0; i < m; i++) { int x = v[i];
        for (int j = 0; j < nbasis; j++) if ((x >> piv[j]) & 1) x ^= basis[j];
        if (x) { int p = 31 - __builtin_clz(x);
            for (int j = 0; j < nbasis; j++) if ((basis[j] >> p) & 1) basis[j] ^= x;
            basis[nbasis] = x; piv[nbasis++] = p; } }
    if (normal && nbasis == r - 1) {
        /* free coordinate: the bit not used as pivot */
        int used = 0; for (int j = 0; j < nbasis; j++) used |= 1 << piv[j];
        int f = -1; for (int b = 0; b < r; b++) if (!((used >> b) & 1)) { f = b; break; }
        /* h has h_f = 1 and h_{piv_j} = bit f of basis_j (reduced form) */
        int h = 1 << f; for (int j = 0; j < nbasis; j++) if ((basis[j] >> f) & 1) h |= 1 << piv[j];
        *normal = h;
    }
    return nbasis;
}
static void wht(int64_t *a) {
    for (int len = 1; len < N; len <<= 1)
        for (int i = 0; i < N; i += len << 1)
            for (int j = i; j < i + len; j++) { int64_t u = a[j], v = a[j + len]; a[j] = u + v; a[j + len] = u - v; }
}
static void save(const char *path, int unc) {
    FILE *f = fopen(path, "w");
    fprintf(f, "{\"r\": %d, \"H_columns\": [", r);
    for (int j = 0; j < NC; j++) fprintf(f, "%s%d", j ? ", " : "", cols[j]);
    fprintf(f, "], \"S\": [");
    for (int i = 0; i < k; i++) fprintf(f, "%s%d", i ? ", " : "", S[i]);
    fprintf(f, "], \"uncovered_syndromes\": %d, \"words\": %d}\n", unc, k << (20 - r));
    fclose(f);
}

int main(int argc, char **argv) {
    uint64_t seed = 1; double secs = 60; const char *out = "out"; const char *startf = NULL;
    int hper = 4, d = 0, botmoves = 0, opt;
    r = 14; k = 24;
    while ((opt = getopt(argc, argv, "r:k:s:t:o:h:f:d:b")) != -1) switch (opt) {
        case 'r': r = atoi(optarg); break; case 'k': k = atoi(optarg); break;
        case 's': seed = strtoull(optarg, 0, 10); break; case 't': secs = atof(optarg); break;
        case 'o': out = optarg; break; case 'h': hper = atoi(optarg); break;
        case 'f': startf = optarg; break; case 'd': d = atoi(optarg); break; case 'b': botmoves = 1; break;
        default: fprintf(stderr, "usage: jtabu -r r -k k -s seed -t secs -o out [-h hcols] [-f start] [-d d] [-b]\n"); return 2; }
    rs ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 50; i++) rnd();
    if (d <= 0 || d > r) d = r;
    int SN = 1 << d, lowmask = SN - 1;
    N = 1 << r;
    inB = calloc(N, 1); B = malloc(sizeof(int) * N); S = malloc(sizeof(int) * k); inS = calloc(N, 1);
    cnt = calloc(N, sizeof(uint16_t)); w = malloc(sizeof(int64_t) * N); U = malloc(sizeof(int) * N);
    posU = malloc(sizeof(int) * N); tabuS = calloc(N, sizeof(int));
    for (int x = 0; x < N; x++) w[x] = 1;
    if (startf) {
        FILE *f = fopen(startf, "r"); int rr;
        if (fscanf(f, "%d", &rr) != 1 || rr != r) { fprintf(stderr, "start r mismatch\n"); return 2; }
        for (int j = 0; j < NC; j++) if (fscanf(f, "%d", &cols[j]) != 1) return 2;
        int m = 0, v; while (m < k && fscanf(f, "%d", &v) == 1) if (!inS[v]) { S[m++] = v; inS[v] = 1; }
        while (m < k) { v = rnd() % SN; if (!inS[v]) { S[m++] = v; inS[v] = 1; } }
        fclose(f);
    } else {
        do { for (int j = 0; j < NC; j++) cols[j] = 1 + rnd() % (N - 1); } while (rank_of(cols, NC, NULL) < r);
        for (int m = 0; m < k;) { int v = rnd() % SN; if (!inS[v]) { S[m++] = v; inS[v] = 1; } }
    }
    recompute();
    fprintf(stderr, "r=%d k=%d d=%d words=%d |B|=%d start uncovered=%d\n", r, k, d, k << (20 - r), nb, nU);
    int64_t *cw = calloc(N, sizeof(int64_t)), *pw = calloc(N, sizeof(int64_t));
    int64_t *fa = malloc(sizeof(int64_t) * N), *fd = malloc(sizeof(int64_t) * N);
    uint8_t *inA = malloc(N), *inD = malloc(N), *b2 = malloc(N), *b3 = malloc(N);
    int *b2l = malloc(sizeof(int) * N), *b3l = malloc(sizeof(int) * N);
    int *touched = malloc(sizeof(int) * N), *P = malloc(sizeof(int) * N);
    int bestU = nU, *bestS = malloc(sizeof(int) * k), bestcols[NC];
    memcpy(bestS, S, sizeof(int) * k); memcpy(bestcols, cols, sizeof cols);
    double t0 = now(), tlast = t0; long it = 0;
    char path[1024];
    int ten_in = 2 + k / 4, ten_out = 1 + k / 8, ten_col = 3;
    while (nU > 0) {
        it++;
        int64_t Wcur = 0; for (int u = 0; u < nU; u++) Wcur += w[U[u]];
        int64_t bestDelta = INT64_MIN; int mtype = -1, bo = -1, bi = -1; int nties = 0;
        /* ---- S moves ---- */
        memset(cw, 0, sizeof(int64_t) * N);
        for (int u = 0; u < nU; u++) { int x = U[u]; int64_t wx = w[x]; for (int i = 0; i < nb; i++) cw[x ^ B[i]] += wx; }
        for (int o = 0; o < k; o++) {
            int so = S[o]; int np = 0; int64_t cost = 0;
            for (int i = 0; i < nb; i++) { int x = so ^ B[i]; if (cnt[x] == 1) { P[np++] = x; cost += w[x]; } }
            int nt = 0;
            for (int q = 0; q < np; q++) { int x = P[q]; int64_t wx = w[x];
                for (int i = 0; i < nb; i++) { int y = x ^ B[i]; if (pw[y] == 0) touched[nt++] = y; pw[y] += wx; } }
            int tab_o = tabuS[so] > it;
            for (int s = 0; s < SN; s++) {
                if (inS[s]) continue;
                int64_t delta = cw[s] + pw[s] - cost;
                if ((tab_o || tabuS[s] > it) && delta <= 0) continue;
                if (delta > bestDelta) { bestDelta = delta; mtype = 0; bo = o; bi = s; nties = 1; }
                else if (delta == bestDelta && rnd() % (++nties) == 0) { mtype = 0; bo = o; bi = s; }
            }
            for (int q = 0; q < nt; q++) pw[touched[q]] = 0;
        }
        /* ---- H column moves ---- */
        for (int h = 0; h < hper; h++) {
            int j = rnd() % NC;
            if (coltabu[j] > it) continue;
            int oth[NC], m = 0; for (int t = 0; t < NC; t++) if (t != j) oth[m++] = cols[t];
            int normal = 0; int rk = rank_of(oth, m, &normal);
            if (rk < r - 1) continue;
            memset(b2, 0, N); memset(b3, 0, N); b2[0] = 1;
            for (int a = 0; a < m; a++) { b2[oth[a]] = 1;
                for (int b = a + 1; b < m; b++) { int ab = oth[a] ^ oth[b]; b2[ab] = 1;
                    for (int c = b + 1; c < m; c++) b3[ab ^ oth[c]] = 1; } }
            int n2 = 0, n3 = 0;
            for (int x = 0; x < N; x++) { if (b2[x]) { b2l[n2++] = x; b3[x] = 1; } if (b3[x]) b3l[n3++] = x; }
            memset(inA, 0, N); memset(inD, 0, N);
            for (int i = 0; i < k; i++) { int s = S[i];
                for (int t = 0; t < n3; t++) inA[s ^ b3l[t]] = 1;
                for (int t = 0; t < n2; t++) inD[s ^ b2l[t]] = 1; }
            int64_t sumA = 0;
            for (int x = 0; x < N; x++) { fa[x] = inA[x] ? 0 : w[x]; sumA += fa[x]; fd[x] = inD[x]; }
            wht(fa); wht(fd);
            for (int x = 0; x < N; x++) fa[x] *= fd[x];
            wht(fa);
            /* W(v) = sumA - corr(v); corr(v) = fa[v] / N */
            int hi = cols[j] & ~lowmask;
            for (int vv = botmoves ? 1 : 0; vv < (botmoves ? N : SN); vv++) {
                int v = botmoves ? vv : (hi | vv);
                if (v == 0 || v == cols[j]) continue;
                if (rk == r - 1 && !(__builtin_popcount(v & normal) & 1)) continue;
                int64_t Wv = sumA - fa[v] / N;
                int64_t delta = Wcur - Wv;
                if (delta > bestDelta) { bestDelta = delta; mtype = 1; bo = j; bi = v; nties = 1; }
                else if (delta == bestDelta && rnd() % (++nties) == 0) { mtype = 1; bo = j; bi = v; }
            }
        }
        if (mtype == 0) {
            int so = S[bo]; rem_s(so); add_s(bi); S[bo] = bi;
            tabuS[so] = it + ten_in + rnd() % 3; tabuS[bi] = it + ten_out + rnd() % 2;
        } else if (mtype == 1) {
            cols[bo] = bi; coltabu[bo] = it + ten_col + rnd() % 3; recompute();
        }
        for (int u = 0; u < nU; u++) w[U[u]] += 1;
        if (nU < bestU) { bestU = nU; memcpy(bestS, S, sizeof(int) * k); memcpy(bestcols, cols, sizeof cols);
            if (bestU < 64) fprintf(stderr, "it %ld t %.1fs uncovered %d |B|=%d\n", it, now() - t0, nU, nb); }
        double tn = now();
        if (tn - tlast > 60) { tlast = tn; fprintf(stderr, "it %ld t %.0fs cur %d best %d |B|=%d\n", it, tn - t0, nU, bestU, nb); }
        if (tn - t0 > secs) break;
    }
    memcpy(S, bestS, sizeof(int) * k); memcpy(cols, bestcols, sizeof cols);
    snprintf(path, sizeof path, "%s.json", out);
    save(path, bestU);
    printf("%s r=%d k=%d words=%d best_uncovered=%d iterations=%ld time=%.1f\n", out, r, k, k << (20 - r), bestU, it, now() - t0);
    return 0;
}
