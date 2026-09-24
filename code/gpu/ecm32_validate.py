#!/usr/bin/env python3
"""
Validate the 32-bit-limb ECM stage-1 kernel (ecm32.cl) against the pure-Python
reference ECM in ecm_ocl.py, bit-for-bit, and check the two kernels (64-bit
mont128 and 32-bit mont32) agree on the factors they find.
"""
import os
import random
import numpy as np
import pyopencl as cl

import ecm_ocl as e

HERE = os.path.dirname(os.path.abspath(__file__))
MASK32 = (1 << 32) - 1


def to4(x):
    return np.array([(x >> (32 * i)) & MASK32 for i in range(4)], dtype=np.uint32)


def from4(a):
    return sum(int(v) << (32 * i) for i, v in enumerate(a))


def run32(cofactors, sigmas, B1, ctx):
    q = cl.CommandQueue(ctx)
    src = open(os.path.join(HERE, "mont32.cl")).read() + "\n" + \
        open(os.path.join(HERE, "ecm32.cl")).read()
    prog = cl.Program(ctx, src).build()
    mf = cl.mem_flags
    E = e.stage1_E(B1)
    ebits = e.ebits_msb(E)
    R = 1 << 128
    items = []
    for ci, n in enumerate(cofactors):
        ninv = (-pow(n, -1, 1 << 32)) & MASK32
        for sg in sigmas:
            bs = e.brent_suyama(n, sg)
            if bs[0] == "factor":
                continue
            x0, z0, b = bs
            items.append((n, ninv, x0 * R % n, z0 * R % n, b * R % n, ci, sg))
    cnt = len(items)
    n_np = np.stack([to4(it[0]) for it in items])
    ninv_np = np.array([it[1] for it in items], dtype=np.uint32)
    x0_np = np.stack([to4(it[2]) for it in items])
    z0_np = np.stack([to4(it[3]) for it in items])
    b_np = np.stack([to4(it[4]) for it in items])
    eb_np = np.array(ebits, dtype=np.uint8)

    def buf(a):
        return cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.ascontiguousarray(a))
    d_n, d_ninv, d_x0, d_z0, d_b, d_eb = map(buf, (n_np, ninv_np, x0_np, z0_np, b_np, eb_np))
    xz = np.empty((cnt * 2, 4), dtype=np.uint32)
    g = np.empty((cnt, 4), dtype=np.uint32)
    d_xz = cl.Buffer(ctx, mf.WRITE_ONLY, xz.nbytes)
    d_g = cl.Buffer(ctx, mf.WRITE_ONLY, g.nbytes)
    k = cl.Kernel(prog, "ecm32_stage1")
    k.set_args(d_n, d_ninv, d_x0, d_z0, d_b, d_eb, np.uint32(len(ebits)), d_xz, d_g)
    cl.enqueue_nd_range_kernel(q, k, (cnt,), None)
    cl.enqueue_copy(q, xz, d_xz)
    cl.enqueue_copy(q, g, d_g)
    q.finish()
    return items, xz, g, ebits


def run32_full(cofactors, sigmas, B1, B2, D, ctx):
    """32-bit stage 1 -> stage 2 on the device. Returns (found, items, g1, g2, plan)."""
    q = cl.CommandQueue(ctx)
    src = (open(os.path.join(HERE, "mont32.cl")).read() + "\n" +
           open(os.path.join(HERE, "ecm32.cl")).read() + "\n" +
           open(os.path.join(HERE, "ecm32_stage2.cl")).read())
    E = e.stage1_E(B1)
    ebits = e.ebits_msb(E)
    pj, pk, n_giant = e.stage2_pairs(B1, B2, D)
    prog = cl.Program(ctx, src).build(
        options=["-DMAXD=%d" % (D + 1), "-DMAXG=%d" % (n_giant + 1)])
    mf = cl.mem_flags
    R = 1 << 128
    items = []
    for ci, n in enumerate(cofactors):
        ninv = (-pow(n, -1, 1 << 32)) & MASK32
        for sg in sigmas:
            bs = e.brent_suyama(n, sg)
            if bs[0] == "factor":
                continue
            x0, z0, b = bs
            items.append((n, ninv, x0 * R % n, z0 * R % n, b * R % n, ci, sg))
    cnt = len(items)
    n_np = np.stack([to4(it[0]) for it in items])
    ninv_np = np.array([it[1] for it in items], dtype=np.uint32)
    x0_np = np.stack([to4(it[2]) for it in items])
    z0_np = np.stack([to4(it[3]) for it in items])
    b_np = np.stack([to4(it[4]) for it in items])
    eb_np = np.array(ebits, dtype=np.uint8)

    def buf(a):
        return cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.ascontiguousarray(a))
    d_n, d_ninv, d_x0, d_z0, d_b, d_eb = map(buf, (n_np, ninv_np, x0_np, z0_np, b_np, eb_np))
    d_pj = buf(np.array(pj, dtype=np.uint32))
    d_pk = buf(np.array(pk, dtype=np.uint32))
    k1 = cl.Kernel(prog, "ecm32_stage1")
    k2 = cl.Kernel(prog, "ecm32_stage2")
    xz = np.empty((cnt * 2, 4), dtype=np.uint32)
    g1 = np.empty((cnt, 4), dtype=np.uint32)
    g2 = np.empty((cnt, 4), dtype=np.uint32)
    CHUNK = 1024
    for s in range(0, cnt, CHUNK):
        m = min(CHUNK, cnt - s)

        def cb(a):
            return cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                             hostbuf=np.ascontiguousarray(a[s:s + m]))
        dn, dv = cb(n_np), cb(ninv_np)
        dx, dz, db = cb(x0_np), cb(z0_np), cb(b_np)
        dxz = cl.Buffer(ctx, mf.READ_WRITE, m * 2 * 16)
        dg1 = cl.Buffer(ctx, mf.WRITE_ONLY, m * 16)
        dg2 = cl.Buffer(ctx, mf.WRITE_ONLY, m * 16)
        k1.set_args(dn, dv, dx, dz, db, d_eb, np.uint32(len(ebits)), dxz, dg1)
        cl.enqueue_nd_range_kernel(q, k1, (m,), None)
        k2.set_args(dn, dv, db, dxz, np.uint32(D), np.uint32(n_giant),
                    d_pj, d_pk, np.uint32(len(pj)), dg2)
        cl.enqueue_nd_range_kernel(q, k2, (m,), None)
        cl.enqueue_copy(q, xz[s * 2:(s + m) * 2], dxz)
        cl.enqueue_copy(q, g1[s:s + m], dg1)
        cl.enqueue_copy(q, g2[s:s + m], dg2)
        q.finish()
    found = []
    for idx, it in enumerate(items):
        n = it[0]
        for garr in (g1, g2):
            gg = from4(garr[idx])
            if 1 < gg < n and n % gg == 0:
                found.append((it[5], it[6], gg))
                break
    return found, items, xz, g1, g2, (ebits, pj, pk, D)


def _validate_stage2(cof, sigmas, B1, ctx):
    B2, D = 5000, 32
    found, items, xz, g1, g2, (ebits, pj, pk, Dd) = run32_full(cof, sigmas, B1, B2, D, ctx)
    R = 1 << 128
    # 32-bit stage-2 gcd == python reference stage-2 gcd
    mis = 0
    for idx, it in enumerate(items):
        n = it[0]
        b = it[4] * pow(R, -1, n) % n
        Qx = from4(xz[2 * idx]) * pow(R, -1, n) % n
        Qz = from4(xz[2 * idx + 1]) * pow(R, -1, n) % n
        ref = e._ref_stage2((Qx, Qz), n, b, D, pj, pk)
        ker = from4(g2[idx])
        if e._gcd(ref, n) != ker and ref != ker:
            mis += 1
    print("mont32 stage2 == python reference for all %d items: %s" % (len(items), mis == 0))
    # same combined finds as the 64-bit full pipeline
    r64 = e.run_full(cof, sigmas, B1, B2, D=D, ctx=ctx)
    f64 = {(ci, sg) for ci, sg, f in r64[0]}
    f32 = {(ci, sg) for ci, sg, f in found}
    print("mont32 stage1+2 finds %d, mont128 finds %d, identical set: %s"
          % (len(f32), len(f64), f32 == f64))
    return mis == 0 and f32 == f64


def main():
    rng = random.Random(77)
    from sympy import nextprime
    cof = []
    for _ in range(300):
        p = int(nextprime(rng.getrandbits(rng.randint(24, 34))))
        qq = int(nextprime(rng.getrandbits(rng.randint(60, 80))))
        if (p * qq).bit_length() < 127:
            cof.append(p * qq)
    sigmas = list(range(6, 26))
    B1 = 600
    ctx = cl.create_some_context()

    items, xz, g, ebits = run32(cof, sigmas, B1, ctx)
    R = 1 << 128
    # 1) 32-bit kernel ladder == python reference, bit-for-bit
    mism = 0
    for idx, it in enumerate(items):
        n = it[0]
        x0 = it[2] * pow(R, -1, n) % n
        z0 = it[3] * pow(R, -1, n) % n
        b = it[4] * pow(R, -1, n) % n
        rx, rz = e.ref_ladder((x0, z0), ebits, n, b)
        kx = from4(xz[2 * idx]) * pow(R, -1, n) % n
        kz = from4(xz[2 * idx + 1]) * pow(R, -1, n) % n
        if (kx, kz) != (rx, rz):
            mism += 1
    print("mont32 kernel ladder == python reference for all %d items: %s"
          % (len(items), mism == 0))

    # 2) factors valid, and same finds as the 64-bit kernel
    found32 = {(it[5], it[6]) for idx, it in enumerate(items)
               if 1 < from4(g[idx]) < it[0] and it[0] % from4(g[idx]) == 0}
    r64 = e.run(cof, sigmas, B1, ctx=ctx)
    found64 = {(ci, sg) for ci, sg, f in r64[0]}
    print("mont32 finds %d, mont128 finds %d, identical set: %s"
          % (len(found32), len(found64), found32 == found64))
    s2ok = _validate_stage2(cof, sigmas, B1, ctx)
    ok = (mism == 0 and found32 == found64 and s2ok)
    print("32-BIT VALIDATION OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
