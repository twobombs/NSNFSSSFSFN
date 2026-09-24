# GPU ECM cofactorization for CADO `las`

Toward GPU-accelerating the dominant cost of the computation (CADO `las`,
~98.5% of the work; see `../../docs/las_gpu_notes.md`). The GPU-amenable part
of `las` is **ECM cofactorization** of sieve survivors: many small, fixed-width
factoring attempts, run over `(survivors × ncurves)` independent work items.

This directory holds a standalone, validated OpenCL implementation:

| file | what |
|------|------|
| `mont128.cl` | 128-bit Montgomery arithmetic (matches CADO `modredc_2ul2`), `mul`/`add`/`sub`/`sqr`/`gcd` |
| `ecm.cl` | ECM **stage 1**: Montgomery-curve differential add/double + ladder `[E]P0`, then `gcd(z, n)` |
| `ecm_stage2.cl` | ECM **stage 2**: baby-step/giant-step "product of point differences" over primes in `(B1,B2]`, one `gcd` |
| `test_mont128.py` | unit tests for the 128-bit ops vs Python (6 moduli × 20k samples) |
| `ecm_ocl.py` | host driver + Brent-Suyama parameterization (exact, matches CADO `-ecm`), stage-2 plan, and a self-test |
| `ecm_bench.py` | throughput benchmark (curves/sec) for any OpenCL device (run it on your GPU) |
| `las_split.py` + `las_profiling.md` | measure the sieving-vs-cofactorization split in `las` (the Amdahl ceiling) |
| `mont32.cl` | 128-bit Montgomery arithmetic with **32-bit limbs** (GPU-friendly: no 64-bit `mul_hi`) |
| `ecm32.cl` | ECM stage 1 on the 32-bit-limb field (same formulas as `ecm.cl`) |
| `ecm32_stage2.cl` | ECM stage 2 on the 32-bit-limb field (counterpart of `ecm_stage2.cl`) |
| `test_mont32.py`, `ecm32_validate.py` | unit tests + bit-for-bit validation of the 32-bit kernel |

Curves use the Brent-Suyama parameterization (CADO's `BRENT12`), computed
exactly on the host so a given sigma yields the same curve CADO uses. Stage 1
multiplies `P0` by `E = prod p^k ≤ B1`.

## 64-bit vs 32-bit limbs

Two arithmetic backends compute the same results: `mont128.cl` uses 2x64-bit
limbs (`mul_hi(ulong,ulong)`), `mont32.cl` uses 4x32-bit limbs (only
`uint*uint` products). On a **CPU** (incl. PoCL) the 64-bit backend is faster
because the CPU multiplies 64-bit natively — measured here ~13k vs ~9k
curves/sec. On **GCN/Vega and most GPUs**, 64-bit integer multiply is
synthesized from 32-bit ops, so the **32-bit backend is expected to be
faster**; that is the whole reason it exists. Benchmark both on the target
(`ecm_bench.py --limb 32` vs `--limb 64`, with `--b2` for stage 2) to see
which wins there. The 32-bit backend now covers stage 1 **and** stage 2.

## Validation (PoCL CPU OpenCL — no GPU needed to develop)

- `test_mont128.py`: all 128-bit ops match Python across 6 random 127-bit
  moduli, 20k samples each.
- `ecm_ocl.py`: the OpenCL stage-1 ladder **and** stage-2 accumulator each
  match a pure-Python reference **bit-for-bit** across 6000 `(cofactor, sigma)`
  items, and every factor the kernel emits is a true non-trivial divisor.
  Stage 2 more than doubles the yield (e.g. 790 → 1935 finds on one batch).
- Cross-check against CADO's own `sieve/ecm/testbench`:
  - **stage 1 only** (`-ecm B1 0`): 99.3% agreement, and the kernel never
    over-reports (kernel-only = 0). The ~0.7% `testbench` catches and stage 1
    alone does not are curves whose point order has a prime square near `B1`
    (a CADO stage-1 setup nuance).
  - **stage 1+2** (`-ecm B1 B2`, B2=5000): **99.9% agreement**, and the kernel
    catches everything `testbench` does (testbench-only = 0); the only
    differences are a couple of extra valid factors the kernel finds.

## Scope / honesty

- PoCL runs OpenCL on the CPU; it is for **correctness and development**, not
  to beat CADO's hand-tuned C on the same CPU. The same kernel is what would
  run on a GPU with thousands of concurrent work items.
- The Brent-Suyama parameterization is done on the host in Python here; it is
  negligible next to the ladder and can move on-device later.
- Integration into `las`'s batch cofactorization path (`ecm/batch.cpp`) is a
  separate, larger effort (see the notes doc).

## Command-line usage

Install the dependencies once (PoCL provides a CPU OpenCL device, so no GPU is
required to develop):

```bash
pip install pyopencl numpy sympy pocl-binary-distribution
```

Every tool selects its OpenCL device from the `PYOPENCL_CTX` environment
variable (e.g. `PYOPENCL_CTX=0`); omit it to be prompted interactively.

### `ecm_bench.py` — throughput benchmark

The benchmark is this script itself (there is no `--bench` flag). Each run
prints **kernel-only** curves/sec (pure device time — this is what scales
across a GPU or cluster) and **end-to-end** curves/sec (including the
host-side Brent-Suyama parameterization).

| option | default | meaning |
|--------|---------|---------|
| `--list-devices` | off | list OpenCL platforms/devices (name, type, compute units, clock, memory), then exit |
| `--curves N` | 100000 | number of `(cofactor, sigma)` work items |
| `--b1 N` | 600 | stage-1 smoothness bound B1 |
| `--b2 N` | 0 | stage-2 bound; `> --b1` enables stage 2, `0` = stage 1 only |
| `--d N` | 32 | stage-2 giant/baby-step size |
| `--reps N` | 5 | timed repetitions (reports the best) |
| `--limb {64,32}` | 64 | arithmetic backend: `64` = `mont128` (2×64-bit), `32` = `mont32` (4×32-bit, GPU-friendly) |

```bash
# see what OpenCL devices exist and pick an index for PYOPENCL_CTX
PYOPENCL_CTX=0 python3 ecm_bench.py --list-devices

# default stage-1 benchmark (100k curves, B1=600, 64-bit limbs)
PYOPENCL_CTX=0 python3 ecm_bench.py

# the comparison that matters on a GPU (e.g. Radeon Pro V340): 64- vs 32-bit
PYOPENCL_CTX=0 python3 ecm_bench.py --curves 200000 --b1 600 --limb 64
PYOPENCL_CTX=0 python3 ecm_bench.py --curves 200000 --b1 600 --limb 32

# full stage 1+2 (enable stage 2 with --b2 > --b1), both backends
PYOPENCL_CTX=0 python3 ecm_bench.py --curves 100000 --b1 600 --b2 5000 --limb 64
PYOPENCL_CTX=0 python3 ecm_bench.py --curves 100000 --b1 600 --b2 5000 --limb 32

# sweep B1 to see throughput drop as stage-1 work grows
for b1 in 315 600 3000; do PYOPENCL_CTX=0 python3 ecm_bench.py --b1 $b1; done

# a quick, low-variance smoke run
PYOPENCL_CTX=0 python3 ecm_bench.py --curves 20000 --reps 3
```

### Validation / self-tests (no options)

These take no CLI options; each runs a fixed check and prints a pass/fail
summary. Device is still `PYOPENCL_CTX`.

```bash
PYOPENCL_CTX=0 python3 ecm_ocl.py         # 64-bit stage 1+2 self-test (bit-exact vs reference)
PYOPENCL_CTX=0 python3 ecm32_validate.py  # 32-bit stage 1+2 validation vs reference & 64-bit
PYOPENCL_CTX=0 python3 test_mont128.py    # 64-bit Montgomery arithmetic unit tests
PYOPENCL_CTX=0 python3 test_mont32.py     # 32-bit Montgomery arithmetic unit tests
```

### `las_split.py` — sieving-vs-cofactorization split

Takes a `las -v -v` log on **stdin** or as a **file argument** (no other
options); prints the split and the Amdahl ceiling. See `las_profiling.md` for
producing the input.

```bash
<las> ... -v -v | python3 las_split.py     # pipe a live run
python3 las_split.py las_output.log        # or parse a saved log
```
