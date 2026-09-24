#!/usr/bin/env python3
"""Unit-test the 32-bit-limb 128-bit Montgomery ops in mont32.cl vs Python."""
import os
import random
import numpy as np
import pyopencl as cl

HERE = os.path.dirname(os.path.abspath(__file__))
MASK32 = (1 << 32) - 1


def to4(x):
    return np.array([(x >> (32 * i)) & MASK32 for i in range(4)], dtype=np.uint32)


def from4(a):
    return sum(int(v) << (32 * i) for i, v in enumerate(a))


KERNEL = open(os.path.join(HERE, "mont32.cl")).read() + r"""
__kernel void t_mul(__global const uint4* a, __global const uint4* b,
                    uint4 n, uint ninv, __global uint4* out) {
    int i = get_global_id(0); out[i] = m32_mul(a[i], b[i], n, ninv);
}
__kernel void t_as(__global const uint4* a, __global const uint4* b,
                   uint4 n, __global uint4* add, __global uint4* sub) {
    int i = get_global_id(0); add[i]=m32_add(a[i],b[i],n); sub[i]=m32_sub(a[i],b[i],n);
}
__kernel void t_gcd(__global const uint4* a, __global const uint4* b, __global uint4* o){
    int i=get_global_id(0); o[i]=m32_gcd(a[i],b[i]);
}
"""


def main():
    ctx = cl.create_some_context()
    q = cl.CommandQueue(ctx)
    prog = cl.Program(ctx, KERNEL).build()
    mf = cl.mem_flags
    rng = random.Random(999)
    N = 20000
    import math
    ok = True
    for trial in range(6):
        n = rng.getrandbits(rng.choice([64, 96, 127])) | 1
        R = 1 << 128
        ninv = (-pow(n, -1, 1 << 32)) & MASK32
        a = [rng.randrange(n) for _ in range(N)]
        b = [rng.randrange(n) for _ in range(N)]
        aM = [(x * R) % n for x in a]
        bM = [(x * R) % n for x in b]
        A = np.stack([to4(x) for x in aM]); B = np.stack([to4(x) for x in bM])
        dA = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=A)
        dB = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=B)
        nvec = to4(n)
        out = np.empty_like(A); dOut = cl.Buffer(ctx, mf.WRITE_ONLY, out.nbytes)
        cl.Kernel(prog, "t_mul")(q, (N,), None, dA, dB, nvec, np.uint32(ninv), dOut)
        cl.enqueue_copy(q, out, dOut); q.finish()
        Rinv = pow(R, -1, n)
        bad = sum(1 for i in range(N) if from4(out[i]) != (aM[i]*bM[i]*Rinv) % n)
        add = np.empty_like(A); sub = np.empty_like(A)
        dAdd = cl.Buffer(ctx, mf.WRITE_ONLY, add.nbytes); dSub = cl.Buffer(ctx, mf.WRITE_ONLY, sub.nbytes)
        cl.Kernel(prog, "t_as")(q, (N,), None, dA, dB, nvec, dAdd, dSub)
        cl.enqueue_copy(q, add, dAdd); cl.enqueue_copy(q, sub, dSub); q.finish()
        bas = sum(1 for i in range(N)
                  if from4(add[i]) != (aM[i]+bM[i]) % n or from4(sub[i]) != (aM[i]-bM[i]) % n)
        gA = np.stack([to4(x) for x in a]); gB = np.stack([to4(x) for x in b])
        dGA = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gA)
        dGB = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gB)
        g = np.empty_like(gA); dG = cl.Buffer(ctx, mf.WRITE_ONLY, g.nbytes)
        cl.Kernel(prog, "t_gcd")(q, (N,), None, dGA, dGB, dG)
        cl.enqueue_copy(q, g, dG); q.finish()
        bg = sum(1 for i in range(N) if from4(g[i]) != math.gcd(a[i], b[i]))
        ok &= (bad == 0 and bas == 0 and bg == 0)
        print(f"trial {trial}: n={n.bit_length()}b mul_bad={bad} addsub_bad={bas} gcd_bad={bg}")
    print("ALL OK" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
