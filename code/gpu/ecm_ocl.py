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


def stage2_pairs(B1, B2, D):
    """For each prime in (B1,B2], (j,k) with p = 2D*j +/- k, 1<=k<=D."""
    pj, pk = [], []
    twoD = 2 * D
    for p in primerange(B1 + 1, B2 + 1):
        j = (p + D) // twoD           # nearest giant index
        k = abs(p - twoD * j)
        if k == 0 or k > D:
            continue                  # (shouldn't happen for D>=1, prime>2D)
        pj.append(j)
        pk.append(k)
    n_giant = max(pj) if pj else 1
    return pj, pk, n_giant


def _ref_stage2(Qxz, n, b, D, pj, pk):
    """Pure-Python reference matching ecm_stage2.cl."""
    def dbl(P):
        return ref_dbl(P, n, b)

    def add(P, Q, Dp):
        return ref_add(P, Q, Dp, n)
    Q = Qxz
    B = {1: Q}
    if D >= 2:
        B[2] = dbl(Q)
    for k in range(3, D + 1):
        B[k] = add(B[k - 1], Q, B[k - 2])
    step = dbl(B[D])
    G = {1: step}
    ng = max(pj) if pj else 1
    if ng >= 2:
        G[2] = dbl(step)
    for j in range(3, ng + 1):
        G[j] = add(G[j - 1], step, G[j - 2])
    acc = None
    for j, k in zip(pj, pk):
        diff = (G[j][0] * B[k][1] - B[k][0] * G[j][1]) % n
        acc = diff if acc is None else acc * diff % n
    if acc is None:
        return 1
    return _gcd(acc, n)


def run_full(cofactors, sigmas, B1, B2, D=32, ctx=None):
    """Stage 1 then stage 2 on the device. Returns (found, details)."""
    ctx = ctx or cl.create_some_context()
    q = cl.CommandQueue(ctx)
    src = (open(os.path.join(HERE, "mont128.cl")).read() + "\n" +
           open(os.path.join(HERE, "ecm.cl")).read() + "\n" +
           open(os.path.join(HERE, "ecm_stage2.cl")).read())

    E = stage1_E(B1)
    ebits = ebits_msb(E)
    pj, pk, n_giant = stage2_pairs(B1, B2, D)
    # size the per-work-item tables to the actual need (private memory is scarce)
    prog = cl.Program(ctx, src).build(
        options=["-DMAXD=%d" % (D + 1), "-DMAXG=%d" % (n_giant + 1)])
    mf = cl.mem_flags

    items = []
    for ci, n in enumerate(cofactors):
        R = 1 << 128
        invm = (-pow(n, -1, 1 << 64)) & MASK64
        for sg in sigmas:
            bs = brent_suyama(n, sg)
            if bs[0] == "factor":
                continue
            x0, z0, b = bs
            items.append((n, invm, x0 * R % n, z0 * R % n, b * R % n, ci, sg))
    cnt = len(items)

    def col(f, dt):
        return np.array([f(it) for it in items], dtype=dt)
    n_np = np.stack([u128_np(it[0]) for it in items]).astype(np.uint64)
    invm_np = col(lambda it: it[1], np.uint64)
    x0_np = np.stack([u128_np(it[2]) for it in items]).astype(np.uint64)
    z0_np = np.stack([u128_np(it[3]) for it in items]).astype(np.uint64)
    b_np = np.stack([u128_np(it[4]) for it in items]).astype(np.uint64)
    eb_np = np.array(ebits, dtype=np.uint8)

    pj_np = np.array(pj, dtype=np.uint32)
    pk_np = np.array(pk, dtype=np.uint32)
    d_pj = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pj_np)
    d_pk = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pk_np)
    eb_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=eb_np)
    k1 = cl.Kernel(prog, "ecm_stage1")
    k2 = cl.Kernel(prog, "ecm_stage2")

    xz = np.empty((cnt * 2, 2), dtype=np.uint64)
    g1 = np.empty((cnt, 2), dtype=np.uint64)
    g2 = np.empty((cnt, 2), dtype=np.uint64)

    # Dispatch in chunks: PoCL is unstable with very large private-memory
    # kernels over huge global sizes, and a real batch would be chunked too.
    CHUNK = 1024
    for s in range(0, cnt, CHUNK):
        e_ = min(s + CHUNK, cnt)
        m = e_ - s

        def cbuf(a):
            return cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                             hostbuf=np.ascontiguousarray(a))
        d_n = cbuf(n_np[s:e_]); d_invm = cbuf(invm_np[s:e_])
        d_x0 = cbuf(x0_np[s:e_]); d_z0 = cbuf(z0_np[s:e_]); d_b = cbuf(b_np[s:e_])
        d_xz = cl.Buffer(ctx, mf.READ_WRITE, m * 2 * 16)
        d_g1 = cl.Buffer(ctx, mf.WRITE_ONLY, m * 16)
        d_g2 = cl.Buffer(ctx, mf.WRITE_ONLY, m * 16)
        k1.set_args(d_n, d_invm, d_x0, d_z0, d_b, eb_buf,
                    np.uint32(len(ebits)), d_xz, d_g1)
        cl.enqueue_nd_range_kernel(q, k1, (m,), None)
        k2.set_args(d_n, d_invm, d_b, d_xz, np.uint32(D), np.uint32(n_giant),
                    d_pj, d_pk, np.uint32(len(pj)), d_g2)
        cl.enqueue_nd_range_kernel(q, k2, (m,), None)
        cl.enqueue_copy(q, xz[s * 2:e_ * 2], d_xz)
        cl.enqueue_copy(q, g1[s:e_], d_g1)
        cl.enqueue_copy(q, g2[s:e_], d_g2)
        q.finish()

    found = []
    for idx, it in enumerate(items):
        n = it[0]
        for garr in (g1, g2):
            gg = np_u128(garr[idx])
            if 1 < gg < n and n % gg == 0:
                found.append((it[5], it[6], gg))
                break
    return found, items, xz, g1, g2, (ebits, pj, pk, D)


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

    # ---- stage 2 ----
    B2 = 5000
    D = 32
    full = run_full(cof, sigmas, B1, B2, D=D)
    f_found, f_items, f_xz, g1, g2, (ebits2, pj, pk, Dd) = full
    R = 1 << 128
    # kernel stage-2 gcd == python reference stage-2 gcd, per item
    s2mis = 0
    for idx, it in enumerate(f_items):
        n, invm = it[0], it[1]
        b = it[4] * pow(R, -1, n) % n
        Qx = np_u128(f_xz[2 * idx]) * pow(R, -1, n) % n
        Qz = np_u128(f_xz[2 * idx + 1]) * pow(R, -1, n) % n
        ref = _ref_stage2((Qx, Qz), n, b, D, pj, pk)
        ker = np_u128(g2[idx])
        # both reduce to the same divisor structure (ref is a gcd, ker a gcd)
        if _gcd(ref, n) != ker and ref != ker:
            s2mis += 1
            if s2mis <= 3:
                print("  stage2 mismatch item", idx, ref, ker)
    print("kernel stage2 == python reference for all %d items: %s" % (len(f_items), s2mis == 0))

    # stage 1+2 finds strictly more (cofactor,sigma) than stage 1 alone
    s1_hits = {(it[5], it[6]) for idx, it in enumerate(f_items)
               if 1 < np_u128(g1[idx]) < it[0] and it[0] % np_u128(g1[idx]) == 0}
    s12_hits = {(ci, sg) for ci, sg, f in f_found}
    print("stage1 hits=%d  stage1+2 hits=%d  (stage2 added %d)"
          % (len(s1_hits), len(s12_hits), len(s12_hits - s1_hits)))

    return not bad and mism == 0 and s2mis == 0 and len(s12_hits) > len(s1_hits)


if __name__ == "__main__":
    ok = _self_test()
    print("SELF-TEST OK" if ok else "SELF-TEST FAILED")
    raise SystemExit(0 if ok else 1)
