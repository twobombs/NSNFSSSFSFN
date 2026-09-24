/*
 * ECM stage 1 (and stage 2) on 128-bit cofactors, OpenCL, one work-item per
 * (cofactor, curve). Montgomery curve, matching CADO's ec_arith_Montgomery
 * differential add/double. The host supplies, per work item and in Montgomery
 * form: the modulus n, invm, and the Brent-Suyama starting point (x0,z0) and
 * curve constant b = (A+2)/4 (computed exactly on the host so curves match
 * CADO bit-for-bit). Stage 1 multiplies P0 by the scalar E = prod p^k <= B1,
 * whose bits are shared by all work items. Stage 2 walks primes B1<p<=B2 via
 * a precomputed baby-step table and accumulates a product of z-differences,
 * then one gcd. Output: gcd(result, n); 1<g<n is a factor.
 */

/* mont128.cl is prepended by the host at build time (shared with the tests). */

typedef struct { u128 x, z; } point;

/* Montgomery differential double: Q = 2P, given curve constant b=(A+2)/4. */
inline point mdbl(point P, u128 n, ulong invm, u128 b) {
    u128 u = mont_add(P.x, P.z, n); u = mont_sqr(u, n, invm);   /* (x+z)^2 */
    u128 v = mont_sub(P.x, P.z, n); v = mont_sqr(v, n, invm);   /* (x-z)^2 */
    point Q;
    Q.x = mont_mul(u, v, n, invm);                              /* x2 */
    u128 w = mont_sub(u, v, n);                                 /* 4xz */
    u = mont_mul(w, b, n, invm);
    u = mont_add(u, v, n);
    Q.z = mont_mul(w, u, n, invm);
    return Q;
}

/* Montgomery differential add: R = P+Q given D=P-Q. */
inline point madd(point P, point Q, point D, u128 n, ulong invm) {
    u128 u = mont_mul(mont_sub(P.x, P.z, n), mont_add(Q.x, Q.z, n), n, invm);
    u128 v = mont_mul(mont_add(P.x, P.z, n), mont_sub(Q.x, Q.z, n), n, invm);
    u128 w = mont_add(u, v, n);
    v = mont_sub(u, v, n);
    w = mont_sqr(w, n, invm);
    v = mont_sqr(v, n, invm);
    point R;
    R.x = mont_mul(w, D.z, n, invm);
    R.z = mont_mul(v, D.x, n, invm);
    return R;
}

/* [E]P0 by the Montgomery ladder; E given MSB-first in ebits[0..ebitlen-1]. */
inline point ladder(point P0, __global const uchar* ebits, uint ebitlen,
                    u128 n, ulong invm, u128 b) {
    point R0 = P0;
    point R1 = mdbl(P0, n, invm, b);
    for (uint i = 1; i < ebitlen; i++) {   /* skip the leading 1 */
        if (ebits[i]) {
            R0 = madd(R0, R1, P0, n, invm);
            R1 = mdbl(R1, n, invm, b);
        } else {
            R1 = madd(R0, R1, P0, n, invm);
            R0 = mdbl(R0, n, invm, b);
        }
    }
    return R0;
}

inline u128 from_mont(u128 a, u128 n, ulong invm) {
    u128 one; one.x = 1; one.y = 0;
    return mont_mul(a, one, n, invm);
}

/*
 * Stage 1 only. Output g = gcd(z_stage1, n) per work item.
 * n_g, invm_g, x0_g, z0_g, b_g are per-work-item arrays (Montgomery form for
 * x0,z0,b). Stage-1 result point is left in Montgomery form; we convert z and
 * gcd with n.
 */
__kernel void ecm_stage1(__global const u128* n_g,
                         __global const ulong* invm_g,
                         __global const u128* x0_g,
                         __global const u128* z0_g,
                         __global const u128* b_g,
                         __global const uchar* ebits, const uint ebitlen,
                         __global u128* xz_out,    /* 2 per item: x,z (mont) */
                         __global u128* g_out) {
    uint i = get_global_id(0);
    u128 n = n_g[i]; ulong invm = invm_g[i];
    point P0; P0.x = x0_g[i]; P0.z = z0_g[i];
    point R = ladder(P0, ebits, ebitlen, n, invm, b_g[i]);
    xz_out[2*i] = R.x; xz_out[2*i+1] = R.z;
    u128 z = from_mont(R.z, n, invm);
    g_out[i] = gcd128(z, n);
}

/* Stage 2 lives in ecm_stage2.cl (added and validated separately). */
