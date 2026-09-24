/*
 * ECM stage 1 with 32-bit-limb Montgomery arithmetic (mont32.cl), the
 * GPU-friendly variant for cards with slow 64-bit integer multiply.
 * Same algorithm and formulas as ecm.cl; only the field type changes
 * (u128_32 = uint4 instead of ulong2). One work item per (cofactor, curve).
 *
 * mont32.cl is prepended by the host.
 */

typedef struct { u128_32 x, z; } point32;

inline point32 mdbl32(point32 P, u128_32 n, uint ninv, u128_32 b) {
    u128_32 u = m32_add(P.x, P.z, n); u = m32_sqr(u, n, ninv);
    u128_32 v = m32_sub(P.x, P.z, n); v = m32_sqr(v, n, ninv);
    point32 Q;
    Q.x = m32_mul(u, v, n, ninv);
    u128_32 w = m32_sub(u, v, n);
    u = m32_mul(w, b, n, ninv);
    u = m32_add(u, v, n);
    Q.z = m32_mul(w, u, n, ninv);
    return Q;
}

inline point32 madd32(point32 P, point32 Q, point32 D, u128_32 n, uint ninv) {
    u128_32 u = m32_mul(m32_sub(P.x, P.z, n), m32_add(Q.x, Q.z, n), n, ninv);
    u128_32 v = m32_mul(m32_add(P.x, P.z, n), m32_sub(Q.x, Q.z, n), n, ninv);
    u128_32 w = m32_add(u, v, n);
    v = m32_sub(u, v, n);
    w = m32_sqr(w, n, ninv);
    v = m32_sqr(v, n, ninv);
    point32 R;
    R.x = m32_mul(w, D.z, n, ninv);
    R.z = m32_mul(v, D.x, n, ninv);
    return R;
}

inline u128_32 from_mont32(u128_32 a, u128_32 n, uint ninv) {
    u128_32 one = (u128_32)(1, 0, 0, 0);
    return m32_mul(a, one, n, ninv);
}

__kernel void ecm32_stage1(__global const u128_32* n_g,
                           __global const uint* ninv_g,
                           __global const u128_32* x0_g,
                           __global const u128_32* z0_g,
                           __global const u128_32* b_g,
                           __global const uchar* ebits, const uint ebitlen,
                           __global u128_32* xz_out,
                           __global u128_32* g_out) {
    uint i = get_global_id(0);
    u128_32 n = n_g[i]; uint ninv = ninv_g[i]; u128_32 b = b_g[i];
    point32 P0; P0.x = x0_g[i]; P0.z = z0_g[i];
    point32 R0 = P0;
    point32 R1 = mdbl32(P0, n, ninv, b);
    for (uint k = 1; k < ebitlen; k++) {
        if (ebits[k]) { R0 = madd32(R0, R1, P0, n, ninv); R1 = mdbl32(R1, n, ninv, b); }
        else          { R1 = madd32(R0, R1, P0, n, ninv); R0 = mdbl32(R0, n, ninv, b); }
    }
    xz_out[2*i] = R0.x; xz_out[2*i+1] = R0.z;
    u128_32 z = from_mont32(R0.z, n, ninv);
    g_out[i] = m32_gcd(z, n);
}
