/*
 * ECM stage 2 (standard baby-step/giant-step, x-only "product of point
 * differences"). Depends on mont128.cl and ecm.cl (point, mdbl, madd,
 * mont_*, from_mont) being prepended by the host.
 *
 * Input per work item: the stage-1 result point Q = [E]P0 in Montgomery form
 * (xz_in) and its curve constant b. The host chooses the giant step 2D and
 * passes, per prime B1<p<=B2, a pair (j,k) with p = 2D*j +/- k, 1<=k<=D.
 * We build baby points B[k]=[k]Q (k=1..D) and giant points G[j]=[2D*j]Q,
 * then accumulate prod (X(G_j)*Z(B_k) - X(B_k)*Z(G_j)); if p | ord(Q) the
 * corresponding factor vanishes mod the hidden prime, so gcd(prod,n) reveals
 * it. One gcd per work item.
 *
 * MAXD / MAXG bound the per-work-item tables (private memory).
 */
#ifndef MAXD
#define MAXD 64
#endif
#ifndef MAXG
#define MAXG 512
#endif

__kernel void ecm_stage2(__global const u128* n_g,
                         __global const ulong* invm_g,
                         __global const u128* b_g,
                         __global const u128* xz_in,     /* Q=[E]P0, mont */
                         const uint D,                    /* baby 1..D */
                         const uint n_giant,              /* giants 1..n_giant */
                         __global const uint* pk_j,       /* per prime */
                         __global const uint* pk_k,
                         const uint n_primes,
                         __global u128* g_out) {
    uint i = get_global_id(0);
    u128 n = n_g[i]; ulong invm = invm_g[i]; u128 b = b_g[i];
    point Q; Q.x = xz_in[2*i]; Q.z = xz_in[2*i+1];

    point B[MAXD + 1];
    B[1] = Q;
    if (D >= 2) B[2] = mdbl(Q, n, invm, b);
    for (uint k = 3; k <= D; k++)
        B[k] = madd(B[k-1], Q, B[k-2], n, invm);   /* [k]Q, diff [k-2]Q */

    /* giant step is [2D]Q; G[j] = [2D*j]Q */
    point step = mdbl(B[D], n, invm, b);           /* [2D]Q */
    point G[MAXG + 1];
    G[1] = step;
    if (n_giant >= 2) G[2] = mdbl(step, n, invm, b);
    for (uint j = 3; j <= n_giant; j++)
        G[j] = madd(G[j-1], step, G[j-2], n, invm);

    u128 acc; int have = 0;
    for (uint p = 0; p < n_primes; p++) {
        uint j = pk_j[p], k = pk_k[p];
        u128 t1 = mont_mul(G[j].x, B[k].z, n, invm);
        u128 t2 = mont_mul(B[k].x, G[j].z, n, invm);
        u128 diff = mont_sub(t1, t2, n);
        if (!have) { acc = diff; have = 1; }
        else acc = mont_mul(acc, diff, n, invm);
    }
    u128 a; a.x = 0; a.y = 0;
    if (have) a = from_mont(acc, n, invm);
    g_out[i] = gcd128(a, n);
}
