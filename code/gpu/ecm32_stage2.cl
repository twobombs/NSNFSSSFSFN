/*
 * ECM stage 2 with 32-bit-limb Montgomery arithmetic (mont32.cl + ecm32.cl),
 * the GPU-friendly counterpart of ecm_stage2.cl. Same baby-step/giant-step
 * "product of point differences" over primes in (B1,B2]; only the field type
 * changes (u128_32 / point32). mont32.cl and ecm32.cl are prepended by the
 * host.
 */
#ifndef MAXD
#define MAXD 64
#endif
#ifndef MAXG
#define MAXG 512
#endif

__kernel void ecm32_stage2(__global const u128_32* n_g,
                           __global const uint* ninv_g,
                           __global const u128_32* b_g,
                           __global const u128_32* xz_in,   /* Q=[E]P0, mont */
                           const uint D,
                           const uint n_giant,
                           __global const uint* pk_j,
                           __global const uint* pk_k,
                           const uint n_primes,
                           __global u128_32* g_out) {
    uint i = get_global_id(0);
    u128_32 n = n_g[i]; uint ninv = ninv_g[i]; u128_32 b = b_g[i];
    point32 Q; Q.x = xz_in[2*i]; Q.z = xz_in[2*i+1];

    point32 B[MAXD + 1];
    B[1] = Q;
    if (D >= 2) B[2] = mdbl32(Q, n, ninv, b);
    for (uint k = 3; k <= D; k++)
        B[k] = madd32(B[k-1], Q, B[k-2], n, ninv);

    point32 step = mdbl32(B[D], n, ninv, b);          /* [2D]Q */
    point32 G[MAXG + 1];
    G[1] = step;
    if (n_giant >= 2) G[2] = mdbl32(step, n, ninv, b);
    for (uint j = 3; j <= n_giant; j++)
        G[j] = madd32(G[j-1], step, G[j-2], n, ninv);

    u128_32 acc; int have = 0;
    for (uint p = 0; p < n_primes; p++) {
        uint j = pk_j[p], k = pk_k[p];
        u128_32 t1 = m32_mul(G[j].x, B[k].z, n, ninv);
        u128_32 t2 = m32_mul(B[k].x, G[j].z, n, ninv);
        u128_32 diff = m32_sub(t1, t2, n);
        if (!have) { acc = diff; have = 1; }
        else acc = m32_mul(acc, diff, n, ninv);
    }
    u128_32 a = (u128_32)(0, 0, 0, 0);
    if (have) a = from_mont32(acc, n, ninv);
    g_out[i] = m32_gcd(a, n);
}
