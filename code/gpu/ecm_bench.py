#!/usr/bin/env python3
"""
Throughput benchmark for the OpenCL ECM cofactorization kernels, to run on a
real target device (e.g. a Radeon Pro V340 GPU) and replace the estimates in
docs/las_gpu_notes.md with a measured curves/sec number.

It reports two numbers, because they scale differently:

  * kernel-only curves/sec  -- pure device dispatch time (enqueue -> finish),
    everything already resident on the device. This is what scales with GPU
    lanes and across a cluster.
  * end-to-end curves/sec   -- including the host-side Brent-Suyama
    parameterization (currently pure-Python bignum) and the host<->device
    copies. The parameterization is negligible next to the ladder and can be
    moved on-device or pipelined; this number is the pessimistic bound today.

Select the device with PYOPENCL_CTX (e.g. PYOPENCL_CTX=0). With no GPU it runs
on any OpenCL platform including PoCL (CPU), which is how it was developed.

Examples:
  PYOPENCL_CTX=0 python3 ecm_bench.py --list-devices
  PYOPENCL_CTX=0 python3 ecm_bench.py --curves 200000 --b1 600
  PYOPENCL_CTX=0 python3 ecm_bench.py --curves 200000 --b1 600 --b2 5000
"""
import argparse
import os
import random
import time

import numpy as np
import pyopencl as cl

import ecm_ocl as e   # same directory


def list_devices():
    for pi, p in enumerate(cl.get_platforms()):
        for di, d in enumerate(p.get_devices()):
            gmhz = d.max_clock_frequency
            print("platform %d / device %d: %s | %s | %d CU @ %d MHz | "
                  "%.1f GB | OpenCL %s"
                  % (pi, di, d.name.strip(), cl.device_type.to_string(d.type),
                     d.max_compute_units, gmhz,
                     d.global_mem_size / 2**30, d.opencl_c_version))


def _rand_cofactors(count, rng):
    """count semiprimes < 2^127 with a ~26-34 bit factor (realistic leftovers)."""
    from sympy import nextprime
    out = []
    while len(out) < count:
        p = int(nextprime(rng.getrandbits(rng.randint(26, 34))))
        q = int(nextprime(rng.getrandbits(70)))
        n = p * q
        if n.bit_length() < 127:
            out.append(n)
    return out


def bench(curves, B1, B2, D, reps, ctx):
    """Time the device kernels over `curves` (cofactor,sigma) work items."""
    q = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
    src = (open(os.path.join(e.HERE, "mont128.cl")).read() + "\n" +
           open(os.path.join(e.HERE, "ecm.cl")).read() + "\n" +
           open(os.path.join(e.HERE, "ecm_stage2.cl")).read())

    rng = random.Random(1234)
    # one sigma per cofactor keeps host prep proportional to `curves`
    cofs = _rand_cofactors(curves, rng)
    sigmas_per = [11]

    do_stage2 = B2 > B1
    E = e.stage1_E(B1)
    ebits = e.ebits_msb(E)
    pj, pk, n_giant = e.stage2_pairs(B1, B2, D) if do_stage2 else ([], [], 1)

    opts = ["-DMAXD=%d" % (D + 1), "-DMAXG=%d" % (n_giant + 1)]
    prog = cl.Program(ctx, src).build(options=opts)
    mf = cl.mem_flags

    # ---- host-side parameterization (timed as "host prep") ----
    t0 = time.perf_counter()
    R = 1 << 128
    items = []
    for n in cofs:
        invm = (-pow(n, -1, 1 << 64)) & e.MASK64
        for sg in sigmas_per:
            bs = e.brent_suyama(n, sg)
            if bs[0] == "factor":
                continue
            x0, z0, b = bs
            items.append((n, invm, x0 * R % n, z0 * R % n, b * R % n))
    cnt = len(items)
    n_np = np.stack([e.u128_np(it[0]) for it in items]).astype(np.uint64)
    invm_np = np.array([it[1] for it in items], dtype=np.uint64)
    x0_np = np.stack([e.u128_np(it[2]) for it in items]).astype(np.uint64)
    z0_np = np.stack([e.u128_np(it[3]) for it in items]).astype(np.uint64)
    b_np = np.stack([e.u128_np(it[4]) for it in items]).astype(np.uint64)
    eb_np = np.array(ebits, dtype=np.uint8)
    host_prep = time.perf_counter() - t0

    def buf(a):
        return cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                         hostbuf=np.ascontiguousarray(a))
    d_n = buf(n_np); d_invm = buf(invm_np)
    d_x0 = buf(x0_np); d_z0 = buf(z0_np); d_b = buf(b_np); d_eb = buf(eb_np)
    d_xz = cl.Buffer(ctx, mf.READ_WRITE, cnt * 2 * 16)
    d_g1 = cl.Buffer(ctx, mf.WRITE_ONLY, cnt * 16)
    k1 = cl.Kernel(prog, "ecm_stage1")
    k1.set_args(d_n, d_invm, d_x0, d_z0, d_b, d_eb, np.uint32(len(ebits)), d_xz, d_g1)

    k2 = d_g2 = d_pj = d_pk = None
    if do_stage2:
        d_pj = buf(np.array(pj, dtype=np.uint32))
        d_pk = buf(np.array(pk, dtype=np.uint32))
        d_g2 = cl.Buffer(ctx, mf.WRITE_ONLY, cnt * 16)
        k2 = cl.Kernel(prog, "ecm_stage2")
        k2.set_args(d_n, d_invm, d_b, d_xz, np.uint32(D), np.uint32(n_giant),
                    d_pj, d_pk, np.uint32(len(pj)), d_g2)

    # Dispatch in chunks via a global-work offset. A real GPU could take the
    # whole batch at once, but PoCL is unstable launching a huge global size
    # over the stage-2 kernel's private-memory tables; chunking is also how a
    # streaming batch would feed the device. All buffers stay resident, so this
    # still measures device throughput, not host copies.
    CHUNK = 4096

    def dispatch():
        for s in range(0, cnt, CHUNK):
            m = min(CHUNK, cnt - s)
            cl.enqueue_nd_range_kernel(q, k1, (m,), None, global_work_offset=(s,))
            if do_stage2:
                cl.enqueue_nd_range_kernel(q, k2, (m,), None, global_work_offset=(s,))
        q.finish()

    dispatch()   # warm-up

    best = float("inf")
    total = 0.0
    for _ in range(reps):
        t = time.perf_counter()
        dispatch()
        dt = time.perf_counter() - t
        best = min(best, dt)
        total += dt
    kernel_avg = total / reps

    return {
        "curves": cnt, "B1": B1, "B2": B2, "stage2": do_stage2,
        "kernel_best_s": best, "kernel_avg_s": kernel_avg,
        "host_prep_s": host_prep,
        "kernel_cps": cnt / best,
        "endtoend_cps": cnt / (best + host_prep),
    }


def bench32(curves, B1, reps, ctx):
    """Stage-1 throughput of the 32-bit-limb kernel (ecm32.cl)."""
    MASK32 = (1 << 32) - 1
    q = cl.CommandQueue(ctx)
    src = (open(os.path.join(e.HERE, "mont32.cl")).read() + "\n" +
           open(os.path.join(e.HERE, "ecm32.cl")).read())
    prog = cl.Program(ctx, src).build()
    mf = cl.mem_flags
    rng = random.Random(1234)
    cofs = _rand_cofactors(curves, rng)
    E = e.stage1_E(B1)
    ebits = e.ebits_msb(E)
    R = 1 << 128

    def to4(x):
        return np.array([(x >> (32 * i)) & MASK32 for i in range(4)], dtype=np.uint32)

    t0 = time.perf_counter()
    items = []
    for n in cofs:
        ninv = (-pow(n, -1, 1 << 32)) & MASK32
        bs = e.brent_suyama(n, 11)
        if bs[0] == "factor":
            continue
        x0, z0, b = bs
        items.append((to4(n), ninv, to4(x0 * R % n), to4(z0 * R % n), to4(b * R % n)))
    cnt = len(items)
    n_np = np.stack([it[0] for it in items])
    ninv_np = np.array([it[1] for it in items], dtype=np.uint32)
    x0_np = np.stack([it[2] for it in items])
    z0_np = np.stack([it[3] for it in items])
    b_np = np.stack([it[4] for it in items])
    eb_np = np.array(ebits, dtype=np.uint8)
    host_prep = time.perf_counter() - t0

    def buf(a):
        return cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.ascontiguousarray(a))
    d_n, d_ninv, d_x0, d_z0, d_b, d_eb = map(buf, (n_np, ninv_np, x0_np, z0_np, b_np, eb_np))
    d_xz = cl.Buffer(ctx, mf.WRITE_ONLY, cnt * 2 * 16)
    d_g = cl.Buffer(ctx, mf.WRITE_ONLY, cnt * 16)
    k = cl.Kernel(prog, "ecm32_stage1")
    k.set_args(d_n, d_ninv, d_x0, d_z0, d_b, d_eb, np.uint32(len(ebits)), d_xz, d_g)
    CHUNK = 4096

    def dispatch():
        for s in range(0, cnt, CHUNK):
            m = min(CHUNK, cnt - s)
            cl.enqueue_nd_range_kernel(q, k, (m,), None, global_work_offset=(s,))
        q.finish()
    dispatch()
    best = min((_timed(dispatch) for _ in range(reps)))
    return {"curves": cnt, "B1": B1, "stage2": False, "kernel_best_s": best,
            "host_prep_s": host_prep, "kernel_cps": cnt / best,
            "endtoend_cps": cnt / (best + host_prep)}


def _timed(fn):
    t = time.perf_counter(); fn(); return time.perf_counter() - t


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--list-devices", action="store_true")
    ap.add_argument("--curves", type=int, default=100000,
                    help="number of (cofactor,sigma) work items")
    ap.add_argument("--b1", type=int, default=600)
    ap.add_argument("--b2", type=int, default=0,
                    help="stage-2 bound; >b1 enables stage 2 (0 = stage 1 only)")
    ap.add_argument("--d", type=int, default=32, help="stage-2 giant/baby size")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--limb", choices=["64", "32"], default="64",
                    help="limb width: 64 (mont128) or 32 (mont32, GPU-friendly). "
                         "32 is stage-1 only for now.")
    args = ap.parse_args()

    ctx = cl.create_some_context()
    dev = ctx.devices[0]
    print("device: %s (%s, %d CU @ %d MHz)"
          % (dev.name.strip(), cl.device_type.to_string(dev.type),
             dev.max_compute_units, dev.max_clock_frequency))
    if args.list_devices:
        list_devices()
        return 0

    if args.limb == "32":
        r = bench32(args.curves, args.b1, args.reps, ctx)
        r["B2"] = 0
    else:
        r = bench(args.curves, args.b1, args.b2, args.d, args.reps, ctx)
    stage = "stage 1+2" if r["stage2"] else "stage 1"
    print("\n%s  (%s-bit limbs)  B1=%d%s  work items=%d  reps=%d"
          % (stage, args.limb, r["B1"], ("  B2=%d" % r["B2"]) if r["stage2"] else "",
             r["curves"], args.reps))
    print("  kernel-only : best %.4f s  -> %s curves/sec"
          % (r["kernel_best_s"], f"{r['kernel_cps']:,.0f}"))
    print("  host prep   : %.4f s (pure-Python Brent-Suyama, movable on-device)"
          % r["host_prep_s"])
    print("  end-to-end  : %s curves/sec (kernel + host prep)"
          % f"{r['endtoend_cps']:,.0f}")
    print("\nCompare to the CPU baseline in docs/las_gpu_notes.md")
    print("(CADO ECM ~92k curves/sec/core at B1=600 on the dev Xeon).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
