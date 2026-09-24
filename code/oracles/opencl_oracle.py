#!/usr/bin/env python3
"""
OpenCL batch signing oracle: compute a^d mod N for many bases a on an OpenCL
device (GPU, or any CPU device such as PoCL). Drop-in for sage_oracle.py --
same -d/-N/--in/--out interface and the same JSON output {str(a): str(a^d)}.

The heavy per-query work (a fixed-N, fixed-d modular exponentiation) is done
on the device with Montgomery (CIOS) big-integer arithmetic; one work-item
per base. N and d are the SOFTWARE-SIMULATED private key, exactly as in
sage_oracle.py -- this only speeds up the local simulation used for
experiments, not any real oracle.

Needs pyopencl and numpy. Select the device with PYOPENCL_CTX, e.g.
PYOPENCL_CTX=0 (or run interactively to be prompted). Falls back to a pure
Python pow() if OpenCL is unavailable, so it always produces correct output.
"""

import argparse
import os
import sys
from datetime import datetime
from math import ceil
from time import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from misc_tools import fast_json_dump  # noqa: E402

LIMB_BITS = 32
LIMB_MASK = (1 << LIMB_BITS) - 1
_KERNEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "opencl_oracle.cl")


def timestamp(ts=None):
    return datetime.fromtimestamp(ts or time()).strftime("%Y-%m-%d %H:%M:%S")


def timeprint(*args):
    print("%s:" % timestamp(), *args)
    sys.stdout.flush()


def n_limbs(x):
    return max(1, ceil(x.bit_length() / LIMB_BITS))


def to_limbs(x, limbs):
    return [(x >> (LIMB_BITS * i)) & LIMB_MASK for i in range(limbs)]


def from_limbs(arr):
    return sum(int(v) << (LIMB_BITS * i) for i, v in enumerate(arr))


def mont_n0inv(n):
    """-n^{-1} mod 2^32 (n odd)."""
    inv = pow(n % (1 << LIMB_BITS), -1, 1 << LIMB_BITS)
    return (-inv) & LIMB_MASK


def run_python(d, N, bases):
    """Reference/fallback path."""
    return {str(a): str(pow(a, d, N)) for a in bases}


def run_opencl(d, N, bases):
    import numpy as np
    import pyopencl as cl

    if N % 2 == 0:
        raise ValueError("N must be odd for Montgomery arithmetic")

    limbs = n_limbs(N)
    dbits = max(1, d.bit_length())
    R = 1 << (LIMB_BITS * limbs)
    r2 = (R * R) % N
    n0inv = mont_n0inv(N)

    count = len(bases)
    base_flat = np.zeros(count * limbs, dtype=np.uint32)
    for i, a in enumerate(bases):
        base_flat[i * limbs:(i + 1) * limbs] = to_limbs(a % N, limbs)
    n_arr = np.array(to_limbs(N, limbs), dtype=np.uint32)
    r2_arr = np.array(to_limbs(r2, limbs), dtype=np.uint32)
    d_arr = np.array(to_limbs(d, limbs), dtype=np.uint32)

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    with open(_KERNEL) as fh:
        src = fh.read()
    prog = cl.Program(ctx, src).build(
        options=["-DLIMBS=%d" % limbs, "-DDBITS=%d" % dbits])

    mf = cl.mem_flags
    d_base = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=base_flat)
    d_n = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=n_arr)
    d_r2 = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=r2_arr)
    d_d = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=d_arr)
    out = np.empty(count * limbs, dtype=np.uint32)
    d_out = cl.Buffer(ctx, mf.WRITE_ONLY, out.nbytes)

    prog.powmod(queue, (count,), None, d_base, d_n, d_r2, d_d,
                np.uint32(n0inv), np.uint32(count), d_out)
    cl.enqueue_copy(queue, out, d_out)
    queue.finish()

    rdict = {}
    for i, a in enumerate(bases):
        rdict[str(a)] = str(from_limbs(out[i * limbs:(i + 1) * limbs]))
    return rdict


def main():
    parser = argparse.ArgumentParser(
        prog="opencl_oracle.py",
        description="OpenCL batch a^d mod N signing oracle (sage_oracle.py drop-in)")
    parser.add_argument("-d", "--d", dest="d", required=True)
    parser.add_argument("-N", "--N", dest="N", required=True)
    parser.add_argument("--in", dest="infile", required=True)
    parser.add_argument("--out", dest="outfile", required=True)
    parser.add_argument("--force-cpu", action="store_true",
                        help="skip OpenCL, use the pure-Python reference path")
    args = parser.parse_args()

    d = int(args.d)
    N = int(args.N)
    with open(args.infile) as fh:
        bases = [int(line) for line in fh if line.strip()]

    use_opencl = not args.force_cpu
    if use_opencl:
        try:
            timeprint("Running OpenCL oracle on %d queries" % len(bases))
            rdict = run_opencl(d, N, bases)
        except Exception as ex:  # noqa: BLE001 -- correctness must not depend on the device
            timeprint("OpenCL unavailable (%s); falling back to Python pow()" % ex)
            use_opencl = False
    if not use_opencl:
        timeprint("Running Python oracle on %d queries" % len(bases))
        rdict = run_python(d, N, bases)

    timeprint("Finished oracle queries")
    fast_json_dump(rdict, args.outfile)


if __name__ == "__main__":
    main()
