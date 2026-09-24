# Measuring the sieving vs cofactorization split in `las`

The GPU ECM kernels here accelerate **cofactorization** only. How much that
helps end to end is capped by cofactorization's share of `las` time (Amdahl).
That share is **not fixed** — it grows sharply with `mfb` (the large-prime
cofactor bound), and the real configs use large `mfb`. Measure it on your own
config rather than assuming; `las_split.py` does the arithmetic.

## Build `las`

On top of the `testbench` recipe in `../../docs/las_gpu_notes.md`, `las` also
needs hwloc:

```bash
apt-get install -y libgmp-dev libhwloc-dev python3-flask python3-requests
cd code/cado
make las makefb          # -> build/<host>/sieve/{las,makefb}
```

`las` prints the split only when built **without** `-production` (the default
`make` build is fine) and run at **verbose level 2** (`-v -v`).

## Run and parse

Quickest self-contained check, using a bundled test polynomial:

```bash
B=code/cado/build/vm/sieve
P=code/cado/tests/misc/c60.poly
$B/makefb -poly $P -lim 111342 -maxbits 10 -out /tmp/c60.fb1
$B/las -poly $P -fb1 /tmp/c60.fb1 -lim0 78682 -lim1 111342 \
       -lpb0 18 -lpb1 19 -mfb0 17 -mfb1 80 -I 10 \
       -q0 200000 -q1 202000 -v -v -t 1 \
  | python3 code/gpu/las_split.py
```

`las_split.py` reads the `# Total cpu time … sieving S … factor F …` line and
prints the cofactorization fraction plus the Amdahl ceiling for a
cofactorization-only GPU speedup.

To profile a **real** config from this repo, run `las` with that config's
`lpb`/`mfb`/`A`/`lim` (see `code/config/nNNN.config` and the `las` calls in
`code/helpers.py` `do_algebraic_query_sieving` / `do_fb_extension_sieving`) on
its generated poly + factor base, adding `-v -v`, and pipe to `las_split.py`.

## Measured on this box (c60 poly, one q-block) — the split is mfb-driven

| `mfb1` | sieving | cofactor | cofactor share |
|------:|--------:|---------:|---------------:|
| 38 | 1.1 s | 0.1 s | **8%** |
| 60 | 1.2 s | 1.0 s | 45% |
| 80 | 1.2 s | 4.6 s | **79%** |

At `mfb1=80` cofactorization is ~79% of sieve+cofactor, so moving only
cofactorization to the GPU has an Amdahl ceiling of `1/(1-0.79) ≈ 4.8×`
end to end. The production configs use larger `mfb` still (e.g. n1024:
`sieve.mfb1=120`, `desc.mfb1=150`), so cofactorization is expected to be the
**dominant** `las` cost there — which is exactly where a fast GPU cofactorizer
pays off. Measure your real config to get the exact ceiling.
