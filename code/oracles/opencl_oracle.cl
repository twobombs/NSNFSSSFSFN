/*
 * Batch RSA "raw" signing oracle: compute a^d mod N for many bases a, with N
 * and d fixed across the batch. One work-item per base.
 *
 * Fixed-width big integers of LIMBS 32-bit limbs, little-endian. Modular
 * multiplication is CIOS Montgomery multiplication (Koc et al.), so the host
 * supplies n0inv = -N^{-1} mod 2^32 and R2 = R^2 mod N with R = 2^(32*LIMBS).
 *
 * LIMBS and DBITS are supplied by the host as -D defines when building.
 */

#ifndef LIMBS
#error "define LIMBS"
#endif
#ifndef DBITS
#error "define DBITS"
#endif

typedef uint limb;

/* t = a*b*R^{-1} mod n, all LIMBS-limb; needs n0inv. CIOS, t scratch LIMBS+2. */
static void mulmont(const limb *a, const limb *b, const limb *n,
                    limb n0inv, limb *out)
{
    limb t[LIMBS + 2];
    for (int i = 0; i < LIMBS + 2; i++)
        t[i] = 0;

    for (int i = 0; i < LIMBS; i++) {
        ulong C = 0;
        for (int j = 0; j < LIMBS; j++) {
            ulong s = (ulong)t[j] + (ulong)a[j] * (ulong)b[i] + C;
            t[j] = (limb)s;
            C = s >> 32;
        }
        ulong s = (ulong)t[LIMBS] + C;
        t[LIMBS] = (limb)s;
        t[LIMBS + 1] = (limb)(s >> 32);

        limb m = (limb)((ulong)t[0] * (ulong)n0inv);
        ulong cs = (ulong)t[0] + (ulong)m * (ulong)n[0];
        C = cs >> 32;
        for (int j = 1; j < LIMBS; j++) {
            cs = (ulong)t[j] + (ulong)m * (ulong)n[j] + C;
            t[j - 1] = (limb)cs;
            C = cs >> 32;
        }
        cs = (ulong)t[LIMBS] + C;
        t[LIMBS - 1] = (limb)cs;
        t[LIMBS] = t[LIMBS + 1] + (limb)(cs >> 32);
    }

    /* conditional final subtraction: if t >= n then t -= n */
    limb borrow = 0, diff[LIMBS];
    for (int j = 0; j < LIMBS; j++) {
        long d = (long)t[j] - (long)n[j] - (long)borrow;
        diff[j] = (limb)d;
        borrow = (d < 0) ? 1 : 0;
    }
    /* subtract when there was no borrow out of the top, or t had a carry limb */
    limb ge = (t[LIMBS] != 0) | (borrow == 0);
    for (int j = 0; j < LIMBS; j++)
        out[j] = ge ? diff[j] : t[j];
}

__kernel void powmod(__global const limb *bases,   /* count * LIMBS */
                     __global const limb *n_g,     /* LIMBS */
                     __global const limb *r2_g,    /* LIMBS, = R^2 mod n */
                     __global const limb *d_g,     /* LIMBS, exponent */
                     const limb n0inv,
                     const uint count,
                     __global limb *out)           /* count * LIMBS */
{
    uint gid = get_global_id(0);
    if (gid >= count)
        return;

    limb n[LIMBS], r2[LIMBS], d[LIMBS], base[LIMBS];
    for (int j = 0; j < LIMBS; j++) {
        n[j] = n_g[j];
        r2[j] = r2_g[j];
        d[j] = d_g[j];
        base[j] = bases[(size_t)gid * LIMBS + j];
    }

    /* aM = base * R mod n  (= mulmont(base, R^2)) */
    limb aM[LIMBS];
    mulmont(base, r2, n, n0inv, aM);

    /* result = R mod n = Montgomery form of 1 = mulmont(1, R^2) */
    limb one[LIMBS], result[LIMBS], tmp[LIMBS];
    for (int j = 0; j < LIMBS; j++)
        one[j] = 0;
    one[0] = 1;
    mulmont(one, r2, n, n0inv, result);

    /* left-to-right square-and-multiply over the bits of d, MSB first */
    for (int bit = DBITS - 1; bit >= 0; bit--) {
        mulmont(result, result, n, n0inv, tmp);
        for (int j = 0; j < LIMBS; j++)
            result[j] = tmp[j];
        if ((d[bit >> 5] >> (bit & 31)) & 1u) {
            mulmont(result, aM, n, n0inv, tmp);
            for (int j = 0; j < LIMBS; j++)
                result[j] = tmp[j];
        }
    }

    /* convert out of Montgomery form: mulmont(result, 1) */
    mulmont(result, one, n, n0inv, tmp);
    for (int j = 0; j < LIMBS; j++)
        out[(size_t)gid * LIMBS + j] = tmp[j];
}
