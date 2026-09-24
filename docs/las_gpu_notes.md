# Moving `las` toward the GPU: target, build recipe, and baseline

Working notes for accelerating the dominant cost of the computation. The
container these were produced in is ephemeral, so everything needed to
reproduce the measurements is written down here.

## Where the time goes

From `timings_clean/timing.n1024` (~12.16M CPU-hours total):

| Stage | CPU-hours | Share | Program |
|---|---:|---:|---|
| Extension sieving | 9,240,758 | 76.0% | CADO `sieve/las` |
| Descent (mostly `descent_bottom`) | 1,421,195 | 11.7% | CADO `las` (descent) + ECM |
| Algebraic query sieving | 1,311,975 | 10.8% | CADO `las` |
| Linear algebra (`bwc`) | 146,264 | 1.2% | CADO block Wiedemann |
| Polynomial selection | 39,992 | 0.3% | CADO `polyselect` |
| e-th root, filter | ~2,235 | <0.02% | Python in this repo |

~98.5% of the work is CADO `las`. So "GPU work on this repo" means GPU work
inside the `code/cado` submodule, not the Python here.

`las` itself splits into two parts per special-q:
1. **Sieving** — marking the sieve region. Bucket sieving is scatter-heavy
   and memory-bound; a poor GPU fit, and a large rewrite.
2. **Cofactorization** — for each sieve survivor, finish factoring the
   leftover norms with P-1 / P+1 / **ECM**. This is many small, independent,
   fixed-width factoring attempts: an excellent GPU fit.

The logs don't record the sieving-vs-cofactorization split (las ran under the
Slurm helpers and its per-job stdout wasn't captured), so it must be measured
with `las -v` or with `testbench` (below).

## The GPU target: ECM cofactorization (`facul`)

Entry points (all in `code/cado/sieve/`):
- `las-cofactor.hpp` / `.cpp`: `factor_leftover_norms(...)` → `facul(...)`.
- `ecm/batch.{hpp,cpp}`: the **batch cofactorization** path. Survivors are
  collected as `cofac_candidate` and factored in bulk at the end of `las`
  (`las.cpp` "batch cofactorization"). This is the ideal GPU seam — the work
  is already gathered into a batch.
- `ecm/facul.cpp`, `facul_doit*.cpp`: run an ordered strategy of methods
  (P-1, P+1, then several ECM curves) against each cofactor.

Parallelism available: (number of survivors) × (`ncurves` per side). The
n1024 descent uses up to `ncurves = 600` (`config/n1024.config`), and query
sieving produces very many survivors — easily enough work items to fill a GPU.

Code a kernel must reimplement (all header-only C, templated over width):
- Fixed-width Montgomery modular arithmetic:
  `utils/arith/modredc_ul.h` (64-bit modulus),
  `modredc_15ul.h` (96-bit), `modredc_2ul2.h` (128-bit).
  Cofactors after sieving fit these widths (leftover norm ≤ `mfb`, split into
  ≤ `lpb`-bit pieces), so no bignum is needed on the device.
- Curve arithmetic: `sieve/ecm/ec_arith_Edwards.h` (twisted Edwards, torsion
  12 — the first ECM curve CADO tries), `ec_arith_Montgomery.h`,
  `ec_arith_common.h`; parameterization in `ec_parameterization.h`.
- Stage 1 is a **bytecode** chain (double-base / precomputed) interpreted per
  curve: `ecm.c` (`bytecode_dbchain_interpret_edwards*`). The bytecode depends
  only on B1, so it is computed once on the host and the same stream drives
  every work item.
- Stage 2 uses a baby-step/giant-step table (`ecm.c` stage 2, `batch.cpp`).

## Build recipe (validated in this container)

CADO needs dev packages the base image lacks. What was required to build the
cofactorization harness:

```bash
apt-get update
apt-get install -y libgmp-dev python3-flask python3-requests
# (full las additionally needs libhwloc-dev)
cd code/cado
make testbench          # cmake-configures, builds utils + facul + testbench
# -> build/<host>/sieve/ecm/testbench
```

`sieve/ecm/testbench` is CADO's standalone cofactorization driver and is the
**ground-truth oracle** for any GPU kernel: give it the same cofactors,
curve, B1, B2 and sigma, and compare the factors found. Key options:
`-ecmem12 <B1> <B2> <s>` (Edwards T12), `-ecm/-ecmm12/-ecmm16`,
`-cof <n>` (multiply each input into a realistic composite),
`-inp <file>`, `-vf` (print factors), `-p` (primes only), `-q` (quiet).

## CPU baseline (single core, this box: Xeon @ 2.10GHz)

ECM, Edwards torsion-12 curve, on ~100-bit composite cofactors, 200k curves:

| B1 | curves/sec/core |
|---:|---:|
| 315 | ~137,000 |
| 600 | ~92,000 |
| 3000 | ~20,000 |

This is the number a GPU kernel must beat (per core; ×cores for the socket).
A real B1 for this computation is small (the descent/cofactorization B1s are
in the hundreds–low thousands), so the per-curve work is short and the win
comes from throughput across many work items, exactly what a GPU provides.

## Plan for the OpenCL kernel

Mirror the approach already used for `oracles/opencl_oracle.py`: implement,
validate against a CADO ground truth (here `testbench`), then benchmark on
PoCL (CPU OpenCL, so no GPU is needed to develop and test).

1. **Width first.** Port `modredc_2ul2` (128-bit) Montgomery arithmetic to an
   OpenCL header. Unit-test each op (add/sub/mul/redc/to/from-Montgomery)
   against the C backend over random inputs.
2. **Curve + stage 1.** Port `ec_arith_Edwards` (T12) and a host-side
   bytecode generator for a fixed B1. Kernel = one work item per
   (cofactor, sigma): build the curve, run the stage-1 bytecode chain,
   take gcd(result, N). Validate factors found against
   `testbench -ecmem12 B1 0 s` on the same inputs (stage 1 only, B2=0).
3. **Stage 2.** Add the baby-step/giant-step stage 2; validate against
   `testbench -ecmem12 B1 B2 s`.
4. **Batch driver.** A standalone tool that reads cofactors (the `-inp`
   format testbench uses), runs N curves each on the device, and prints
   factors — diffable against `testbench`. This is the throughput benchmark.
5. **Integration (separate, larger effort).** Feed the batch
   cofactorization path (`ecm/batch.cpp` `cofac_candidate` list) to the
   device instead of the CPU `facul` loop. Keeps `las`'s sieving untouched.

### Risks / open questions
- Matching CADO's exact curve parameterization (sigma → curve) so factors
  are bit-identical to `testbench`; the parameterization headers make this
  tractable but fiddly.
- Divergence: cofactors have varying sizes/strategies; a fixed-width,
  fixed-B1 kernel handles one bucket, so the host must bin candidates by
  width and strategy before dispatch.
- Whether the sieving/cofactorization split makes this worth it — measure
  with `las -v` on a small config (e.g. `config/n192.config`) before
  committing to integration.
