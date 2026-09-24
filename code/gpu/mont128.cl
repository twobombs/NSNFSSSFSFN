/*
 * 128-bit fixed-width Montgomery arithmetic for OpenCL, matching CADO's
 * utils/arith/modredc_2ul2 representation: modulus n < 2^127 (odd), two
 * 64-bit limbs little-endian (lo, hi), R = 2^128, invm = -n^{-1} mod 2^64.
 * Residues are kept in Montgomery form a*R mod n.
 *
 * This header is shared by the ECM kernels and by the arithmetic unit test.
 * mul_hi(ulong,ulong) supplies the high 64 bits of a 64x64 product.
 */

typedef ulong2 u128;   /* .x = low limb, .y = high limb */

inline uint addc64(ulong a, ulong b, uint cin, ulong *out) {
    ulong s = a + cin;
    uint c1 = (s < a);
    s += b;
    *out = s;
    return c1 + (s < b);
}

inline int geq128(u128 a, u128 b) {
    return (a.y > b.y) || (a.y == b.y && a.x >= b.x);
}

inline u128 sub128(u128 a, u128 b) {   /* assumes a >= b */
    u128 r;
    r.x = a.x - b.x;
    r.y = a.y - b.y - (a.x < b.x);
    return r;
}

inline u128 add128(u128 a, u128 b) {   /* mod 2^128 */
    u128 r;
    r.x = a.x + b.x;
    r.y = a.y + b.y + (r.x < a.x);
    return r;
}

/* Montgomery add: (a+b) mod n */
inline u128 mont_add(u128 a, u128 b, u128 n) {
    u128 r = add128(a, b);
    /* a,b < n < 2^127 so a+b < 2^128: a carry out or r>=n both mean subtract */
    int carry = (r.y < a.y) || (r.x < a.x && r.y == a.y);
    if (carry || geq128(r, n))
        r = sub128(r, n);
    return r;
}

inline u128 mont_sub(u128 a, u128 b, u128 n) {
    if (geq128(a, b))
        return sub128(a, b);
    return sub128(add128(a, n), b);   /* a + n - b, a<b<n so a+n < 2^128 */
}

/*
 * CIOS Montgomery multiplication for L=2 limbs.
 * Returns a*b*R^{-1} mod n, with 0 <= result < n.
 */
inline u128 mont_mul(u128 a, u128 b, u128 n, ulong invm) {
    ulong an[2] = {a.x, a.y};
    ulong bn[2] = {b.x, b.y};
    ulong nn[2] = {n.x, n.y};
    ulong t0 = 0, t1 = 0, t2 = 0;

    for (int i = 0; i < 2; i++) {
        ulong bi = bn[i];
        /* t += a * bi */
        ulong lo0 = an[0] * bi, hi0 = mul_hi(an[0], bi);
        ulong lo1 = an[1] * bi, hi1 = mul_hi(an[1], bi);
        uint c;
        c = addc64(t0, lo0, 0, &t0);          /* t0 += lo0 */
        ulong carry_hi0 = hi0 + c;             /* hi0 < 2^64-1 so no overflow */
        c = addc64(t1, lo1, 0, &t1);          /* t1 += lo1 */
        ulong carry = hi1 + c;
        c = addc64(t1, carry_hi0, 0, &t1);    /* t1 += hi0 + c0 */
        carry += c;
        t2 += carry;                           /* t2 += hi1 + carries */

        /* m = t0 * invm mod 2^64; t += m * n; then shift right one limb */
        ulong m = t0 * invm;
        ulong ml0 = m * nn[0], mh0 = mul_hi(m, nn[0]);
        ulong ml1 = m * nn[1], mh1 = mul_hi(m, nn[1]);
        ulong tmp;
        c = addc64(t0, ml0, 0, &tmp);          /* t0 + ml0, low limb -> 0, keep carry */
        ulong carry_mh0 = mh0 + c;
        c = addc64(t1, ml1, 0, &t1);
        ulong carry2 = mh1 + c;
        c = addc64(t1, carry_mh0, 0, &t1);
        carry2 += c;
        /* shift: new t0=t1, t1=t2+carry2, t2=carry out */
        ulong newt1;
        c = addc64(t2, carry2, 0, &newt1);
        t0 = t1;
        t1 = newt1;
        t2 = c;
    }

    u128 r; r.x = t0; r.y = t1;
    /* t2 holds a possible carry limb; reduce so that 0 <= r < n */
    if (t2 || geq128(r, n))
        r = sub128(r, n);
    return r;
}

inline u128 mont_sqr(u128 a, u128 n, ulong invm) { return mont_mul(a, a, n, invm); }

/* binary GCD of two 128-bit integers (not in Montgomery form) */
inline u128 gcd128(u128 a, u128 b) {
    if (a.x == 0 && a.y == 0) return b;
    if (b.x == 0 && b.y == 0) return a;
    int shift = 0;
    while (((a.x | b.x) & 1) == 0) {
        a.x = (a.x >> 1) | (a.y << 63); a.y >>= 1;
        b.x = (b.x >> 1) | (b.y << 63); b.y >>= 1;
        shift++;
    }
    while ((a.x & 1) == 0) { a.x = (a.x >> 1) | (a.y << 63); a.y >>= 1; }
    do {
        while ((b.x & 1) == 0) { b.x = (b.x >> 1) | (b.y << 63); b.y >>= 1; }
        if (geq128(a, b)) { u128 t = a; a = b; b = t; }   /* a = min, b = max */
        b = sub128(b, a);
    } while (!(b.x == 0 && b.y == 0));
    /* a <<= shift */
    for (int i = 0; i < shift; i++) { a.y = (a.y << 1) | (a.x >> 63); a.x <<= 1; }
    return a;
}
