#!/usr/bin/env python3
"""Unit-test the 128-bit Montgomery ops in mont128.cl against Python."""
import os
import random
import numpy as np
import pyopencl as cl

HERE = os.path.dirname(os.path.abspath(__file__))
MASK64 = (1 << 64) - 1


def u128_to_np(x):
    return np.array([x & MASK64, (x >> 64) & MASK64], dtype=np.uint64)


def np_to_u128(a):
    return int(a[0]) | (int(a[1]) << 64)


KERNEL = open(os.path.join(HERE, "mont128.cl")).read() + r"""
__kernel void t_mul(__global const ulong2* a, __global const ulong2* b,
                    ulong2 n, ulong invm, __global ulong2* out) {
    int i = get_global_id(0); out[i] = mont_mul(a[i], b[i], n, invm);
}
__kernel void t_addsub(__global const ulong2* a, __global const ulong2* b,
                       ulong2 n, __global ulong2* add, __global ulong2* sub) {
    int i = get_global_id(0);
    add[i] = mont_add(a[i], b[i], n); sub[i] = mont_sub(a[i], b[i], n);
}
__kernel void t_gcd(__global const ulong2* a, __global const ulong2* b,
                    __global ulong2* out) {
    int i = get_global_id(0); out[i] = gcd128(a[i], b[i]);
}
"""


def main():
    ctx = cl.create_some_context()
    q = cl.CommandQueue(ctx)
    prog = cl.Program(ctx, KERNEL).build()
    mf = cl.mem_flags
    rng = random.Random(12345)
    N = 20000

    # random odd modulus < 2^127
    def rand_mod():
        return (rng.getrandbits(127) | 1)

    ok = True
    for trial in range(6):
        n = rand_mod()
        R = 1 << 128
        invm = (-pow(n, -1, 1 << 64)) & MASK64
        Rmodn = R % n
        a = [rng.randrange(n) for _ in range(N)]
        b = [rng.randrange(n) for _ in range(N)]
        # Montgomery form
        aM = [(x * R) % n for x in a]
        bM = [(x * R) % n for x in b]
        a_np = np.stack([u128_to_np(x) for x in aM]).astype(np.uint64)
        b_np = np.stack([u128_to_np(x) for x in bM]).astype(np.uint64)
        d_a = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a_np)
        d_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_np)
        n_vec = np.array([n & MASK64, n >> 64], dtype=np.uint64)

        out = np.empty_like(a_np)
        d_out = cl.Buffer(ctx, mf.WRITE_ONLY, out.nbytes)
        prog.t_mul(q, (N,), None, d_a, d_b, n_vec, np.uint64(invm), d_out)
        cl.enqueue_copy(q, out, d_out); q.finish()
        # mont_mul(aM,bM) = aM*bM/R = a*b*R mod n  (Montgomery product of aM,bM)
        bad = 0
        for i in range(N):
            exp = (aM[i] * bM[i] * pow(R, -1, n)) % n
            if np_to_u128(out[i]) != exp:
                bad += 1
                if bad <= 3:
                    print("  MUL mismatch", i, hex(np_to_u128(out[i])), hex(exp))
        ok &= bad == 0

        add = np.empty_like(a_np); sub = np.empty_like(a_np)
        d_add = cl.Buffer(ctx, mf.WRITE_ONLY, add.nbytes)
        d_sub = cl.Buffer(ctx, mf.WRITE_ONLY, sub.nbytes)
        prog.t_addsub(q, (N,), None, d_a, d_b, n_vec, d_add, d_sub)
        cl.enqueue_copy(q, add, d_add); cl.enqueue_copy(q, sub, d_sub); q.finish()
        badas = 0
        for i in range(N):
            if np_to_u128(add[i]) != (aM[i] + bM[i]) % n: badas += 1
            if np_to_u128(sub[i]) != (aM[i] - bM[i]) % n: badas += 1
        ok &= badas == 0

        # gcd on plain (non-Montgomery) integers
        ga = np.stack([u128_to_np(x) for x in a]).astype(np.uint64)
        gb = np.stack([u128_to_np(x) for x in b]).astype(np.uint64)
        import math
        d_ga = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=ga)
        d_gb = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gb)
        g = np.empty_like(ga); d_g = cl.Buffer(ctx, mf.WRITE_ONLY, g.nbytes)
        prog.t_gcd(q, (N,), None, d_ga, d_gb, d_g)
        cl.enqueue_copy(q, g, d_g); q.finish()
        badg = sum(1 for i in range(N) if np_to_u128(g[i]) != math.gcd(a[i], b[i]))
        ok &= badg == 0
        print(f"trial {trial}: n={n.bit_length()}b  mul_bad={bad} addsub_bad={badas} gcd_bad={badg}")

    print("ALL OK" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
