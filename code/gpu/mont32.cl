/*
 * 128-bit Montgomery arithmetic with 32-bit limbs, for GPUs whose 64-bit
 * integer multiply is synthesized (AMD GCN / Vega, most consumer/pro cards).
 * Every partial product here is a 32x32 -> 64 multiply (native `uint*uint`
 * widened to `ulong`), so it needs no `mul_hi(ulong,...)`. Same math as
 * mont128.cl (CADO modredc_2ul2), different limb width.
 *
 * A 128-bit value is uint4 = (limb0..limb3), little-endian. Modulus n < 2^127
 * odd; R = 2^128; ninv = -n^{-1} mod 2^32. Residues are in Montgomery form.
 *
 * NLIMBS is fixed at 4 (128-bit). The CIOS loops are written generically over
 * NLIMBS so the same code can be widened later.
 */
#define NLIMBS 4

typedef uint4 u128_32;

inline void ld(uint *l, u128_32 a) { l[0]=a.x; l[1]=a.y; l[2]=a.z; l[3]=a.w; }
inline u128_32 st(const uint *l) { return (u128_32)(l[0], l[1], l[2], l[3]); }

inline int geq(const uint *a, const uint *b) {
    for (int i = NLIMBS - 1; i >= 0; i--) {
        if (a[i] != b[i]) return a[i] > b[i];
    }
    return 1;
}

inline void sub_in(uint *a, const uint *b) {   /* a -= b (assumes a>=b) */
    ulong borrow = 0;
    for (int i = 0; i < NLIMBS; i++) {
        ulong d = (ulong)a[i] - b[i] - borrow;
        a[i] = (uint)d;
        borrow = (d >> 63) & 1;   /* set if underflow */
    }
}

inline u128_32 m32_add(u128_32 A, u128_32 B, u128_32 N) {
    uint a[NLIMBS], b[NLIMBS], n[NLIMBS];
    ld(a, A); ld(b, B); ld(n, N);
    ulong c = 0;
    for (int i = 0; i < NLIMBS; i++) { ulong s = (ulong)a[i] + b[i] + c; a[i] = (uint)s; c = s >> 32; }
    if (c || geq(a, n)) sub_in(a, n);
    return st(a);
}

inline u128_32 m32_sub(u128_32 A, u128_32 B, u128_32 N) {
    uint a[NLIMBS], b[NLIMBS], n[NLIMBS];
    ld(a, A); ld(b, B); ld(n, N);
    if (geq(a, b)) { sub_in(a, b); return st(a); }
    /* a + n - b */
    ulong c = 0;
    for (int i = 0; i < NLIMBS; i++) { ulong s = (ulong)a[i] + n[i] + c; a[i] = (uint)s; c = s >> 32; }
    sub_in(a, b);
    return st(a);
}

/* CIOS Montgomery multiply, 32-bit limbs, 64-bit accumulator. */
inline u128_32 m32_mul(u128_32 A, u128_32 B, u128_32 N, uint ninv) {
    uint a[NLIMBS], b[NLIMBS], n[NLIMBS];
    ld(a, A); ld(b, B); ld(n, N);
    uint t[NLIMBS + 2];
    for (int i = 0; i < NLIMBS + 2; i++) t[i] = 0;

    for (int i = 0; i < NLIMBS; i++) {
        ulong C = 0;
        uint bi = b[i];
        for (int j = 0; j < NLIMBS; j++) {
            ulong p = (ulong)a[j] * bi + t[j] + C;
            t[j] = (uint)p;
            C = p >> 32;
        }
        ulong s = (ulong)t[NLIMBS] + C;
        t[NLIMBS] = (uint)s;
        t[NLIMBS + 1] = (uint)(s >> 32);

        uint m = t[0] * ninv;
        C = ((ulong)m * n[0] + t[0]) >> 32;   /* low limb becomes 0 */
        for (int j = 1; j < NLIMBS; j++) {
            ulong p = (ulong)m * n[j] + t[j] + C;
            t[j - 1] = (uint)p;
            C = p >> 32;
        }
        s = (ulong)t[NLIMBS] + C;
        t[NLIMBS - 1] = (uint)s;
        t[NLIMBS] = t[NLIMBS + 1] + (uint)(s >> 32);
    }

    uint r[NLIMBS];
    for (int i = 0; i < NLIMBS; i++) r[i] = t[i];
    if (t[NLIMBS] || geq(r, n)) sub_in(r, n);
    return st(r);
}

inline u128_32 m32_sqr(u128_32 a, u128_32 n, uint ninv) { return m32_mul(a, a, n, ninv); }

/* binary gcd of two 128-bit integers (not in Montgomery form) */
inline int is_zero(const uint *a) { return (a[0]|a[1]|a[2]|a[3]) == 0; }
inline void shr1(uint *a) {
    for (int i = 0; i < NLIMBS - 1; i++) a[i] = (a[i] >> 1) | (a[i+1] << 31);
    a[NLIMBS-1] >>= 1;
}
inline void shl1(uint *a) {
    for (int i = NLIMBS - 1; i > 0; i--) a[i] = (a[i] << 1) | (a[i-1] >> 31);
    a[0] <<= 1;
}
inline u128_32 m32_gcd(u128_32 A, u128_32 B) {
    uint a[NLIMBS], b[NLIMBS];
    ld(a, A); ld(b, B);
    if (is_zero(a)) return B;
    if (is_zero(b)) return A;
    int shift = 0;
    while (((a[0] | b[0]) & 1) == 0) { shr1(a); shr1(b); shift++; }
    while ((a[0] & 1) == 0) shr1(a);
    do {
        while ((b[0] & 1) == 0) shr1(b);
        if (geq(a, b)) { uint tmp[NLIMBS]; for(int i=0;i<NLIMBS;i++){tmp[i]=a[i];a[i]=b[i];b[i]=tmp[i];} }
        sub_in(b, a);   /* b = b - a, b was >= a */
    } while (!is_zero(b));
    for (int i = 0; i < shift; i++) shl1(a);
    return st(a);
}
