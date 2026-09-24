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

Curves use the Brent-Suyama parameterization (CADO's `BRENT12`), computed
exactly on the host so a given sigma yields the same curve CADO uses. Stage 1
multiplies `P0` by `E = prod p^k ≤ B1`.

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

## Run

```bash
pip install pyopencl numpy sympy pocl-binary-distribution   # PoCL = CPU OpenCL
PYOPENCL_CTX=0 python3 test_mont128.py
PYOPENCL_CTX=0 python3 ecm_ocl.py   # stage 1+2 self-test
```
