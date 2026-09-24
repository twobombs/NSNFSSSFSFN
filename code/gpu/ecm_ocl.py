#!/usr/bin/env python3
"""
OpenCL ECM stage-1 driver + validation.

Runs the Montgomery-curve ECM stage-1 ladder (ecm.cl) on a batch of 128-bit
cofactors, one work item per (cofactor, sigma), using the Brent-Suyama
parameterization computed exactly on the host so curves match CADO's `-ecm`
(BRENT12). Validates the kernel against a pure-Python reference ECM using the
same curve and scalar, and checks every emitted factor divides its cofactor.
"""
import os
import random
import numpy as np
import pyopencl as cl
from sympy import primerange

HERE = os.path.dirname(os.path.abspath(__file__))
MASK64 = (1 << 64) - 1


def u128_np(x):
    return np.array([x & MASK64, (x >> 64) & MASK64], dtype=np.uint64)


def np_u128(a):
    return int(a[0]) | (int(a[1]) << 64)


def stage1_E(B1):
    """E = prod_{p<=B1} p^floor(log_p B1), the standard ECM stage-1 multiplier."""
    E = 1
    for p in primerange(2, B1 + 1):
        pk = p
        while pk * p <= B1:
            pk *= p
        E *= pk
    return E


def brent_suyama(n, sigma):
    """CADO's Brent-Suyama parameterization. Returns (x0, z0, b) mod n, or a
    ('factor', g) if a non-trivial gcd appears during the inversion."""
    s = sigma % n
    v = (4 * s) % n
    u = (s * s - 5) % n
    x0 = pow(u, 3, n)
    z0 = pow(v, 3, n)
    t1 = (16 * x0 * v) % n           # 16 * u^3 * v
    bnum = ((v - u) ** 3 * (3 * u + v)) % n
    g = np.gcd(t1, n) if False else _gcd(t1, n)
    if g != 1:
        return ("factor", g)
    b = (bnum * pow(t1, -1, n)) % n   # (A+2)/4
    return (x0 % n, z0 % n, b % n)


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


# ---- pure-Python reference ECM (same Montgomery formulas as the kernel) ----
def ref_dbl(P, n, b):
    x, z = P
    u = (x + z) ** 2 % n
    v = (x - z) ** 2 % n
    x2 = u * v % n
    w = (u - v) % n
    uu = (w * b + v) % n
    return (x2, w * uu % n)


def ref_add(P, Q, D, n):
    u = (P[0] - P[1]) * (Q[0] + Q[1]) % n
    v = (P[0] + P[1]) * (Q[0] - Q[1]) % n
    w = (u + v) ** 2 % n
    vv = (u - v) ** 2 % n
    return (w * D[1] % n, vv * D[0] % n)


def ref_ladder(P0, ebits, n, b):
    R0, R1 = P0, ref_dbl(P0, n, b)
    for bit in ebits[1:]:
        if bit:
            R0 = ref_add(R0, R1, P0, n)
            R1 = ref_dbl(R1, n, b)
        else:
            R1 = ref_add(R0, R1, P0, n)
            R0 = ref_dbl(R0, n, b)
    return R0


def ebits_msb(E):
    return [int(c) for c in bin(E)[2:]]


def run(cofactors, sigmas, B1, ctx=None):
    """Return list of (cofactor_index, sigma, factor) the kernel found."""
    E = stage1_E(B1)
    ebits = ebits_msb(E)
    items = []          # (n, invm, x0M, z0M, bM, cof_idx, sigma)
    trivial = []
    for ci, n in enumerate(cofactors):
        R = 1 << 128
        invm = (-pow(n, -1, 1 << 64)) & MASK64
        for sg in sigmas:
            bs = brent_suyama(n, sg)
            if bs[0] == "factor":
                trivial.append((ci, sg, bs[1]))
                continue
            x0, z0, b = bs
            items.append((n, invm, x0 * R % n, z0 * R % n, b * R % n, ci, sg))
    if not items:
        return trivial

    cnt = len(items)
    n_np = np.stack([u128_np(it[0]) for it in items]).astype(np.uint64)
    invm_np = np.array([it[1] for it in items], dtype=np.uint64)
    x0_np = np.stack([u128_np(it[2]) for it in items]).astype(np.uint64)
    z0_np = np.stack([u128_np(it[3]) for it in items]).astype(np.uint64)
    b_np = np.stack([u128_np(it[4]) for it in items]).astype(np.uint64)
    eb_np = np.array(ebits, dtype=np.uint8)

    ctx = ctx or cl.create_some_context()
    q = cl.CommandQueue(ctx)
    src = open(os.path.join(HERE, "mont128.cl")).read() + "\n" + \
        open(os.path.join(HERE, "ecm.cl")).read()
    prog = cl.Program(ctx, src).build()
    mf = cl.mem_flags

    def buf(a):
        return cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a)
    d_n, d_invm, d_x0, d_z0, d_b, d_eb = map(
        buf, (n_np, invm_np, x0_np, z0_np, b_np, eb_np))
    xz = np.empty((cnt * 2, 2), dtype=np.uint64)
    g = np.empty((cnt, 2), dtype=np.uint64)
    d_xz = cl.Buffer(ctx, mf.WRITE_ONLY, xz.nbytes)
    d_g = cl.Buffer(ctx, mf.WRITE_ONLY, g.nbytes)
    k = cl.Kernel(prog, "ecm_stage1")
    k.set_args(d_n, d_invm, d_x0, d_z0, d_b, d_eb, np.uint32(len(ebits)), d_xz, d_g)
    cl.enqueue_nd_range_kernel(q, k, (cnt,), None)
    cl.enqueue_copy(q, xz, d_xz)
    cl.enqueue_copy(q, g, d_g)
    q.finish()

    found = list(trivial)
    for idx, it in enumerate(items):
        n = it[0]
        gg = np_u128(g[idx])
        if 1 < gg < n and n % gg == 0:
            found.append((it[5], it[6], gg))
    return found, items, xz, ebits


def _self_test():
    rng = random.Random(2024)

    def randprime(bits):
        from sympy import nextprime
        return int(nextprime(rng.getrandbits(bits)))

    # cofactors: semiprimes with a ~25-35 bit factor, product < 2^127
    cof = []
    for _ in range(300):
        p = randprime(rng.randint(24, 34))
        q = randprime(rng.randint(60, 80))
        if (p * q).bit_length() < 127:
            cof.append(p * q)
    sigmas = list(range(6, 26))
    B1 = 600

    res = run(cof, sigmas, B1)
    found, items, xz, ebits = res
    E = stage1_E(B1)

    # 1) every emitted factor divides its cofactor
    bad = [f for f in found if not (1 < f[2] < cof[f[0]] and cof[f[0]] % f[2] == 0)]
    print("emitted factors all valid divisors:", not bad, "(%d found)" % len(found))

    # 2) kernel stage-1 point matches the pure-Python reference for every item
    R = 1 << 128
    mism = 0
    for idx, it in enumerate(items):
        n, invm = it[0], it[1]
        x0, z0, b = it[2] * pow(R, -1, n) % n, it[3] * pow(R, -1, n) % n, it[4] * pow(R, -1, n) % n
        rx, rz = ref_ladder((x0, z0), ebits, n, b)
        kx = np_u128(xz[2 * idx]) * pow(R, -1, n) % n
        kz = np_u128(xz[2 * idx + 1]) * pow(R, -1, n) % n
        if (kx, kz) != (rx, rz):
            mism += 1
            if mism <= 3:
                print("  ladder mismatch item", idx)
    print("kernel ladder == python reference for all %d items: %s" % (len(items), mism == 0))
    return not bad and mism == 0


if __name__ == "__main__":
    ok = _self_test()
    print("SELF-TEST OK" if ok else "SELF-TEST FAILED")
    raise SystemExit(0 if ok else 1)
