#!/usr/bin/env python3
"""
Parse a CADO `las` verbose run and report the sieving-vs-cofactorization
split -- the Amdahl fraction that caps how much GPU cofactorization can help
end to end.

`las` (built WITHOUT -production, run at verbose level 2, e.g. `-v -v`) prints:

  # Total cpu time T s, useful U s [norm a+b, sieving S (...), factor F (...),
    rest R], wasted+waited W s, rest ...

where `sieving S` is total sieve time and `factor F` is cofactorization time.
The cofactorization fraction of the *useful* CPU time is F / (norm + S + F + R);
the fraction of sieve+cofac is F / (S + F).

Usage:
  <las> -v -v <las-args...>  2>&1 | python3 las_split.py
  python3 las_split.py las_output.log

See las_profiling.md for how to produce the input on a real config.
"""
import re
import sys

LINE = re.compile(
    r"# Total cpu time ([\d.]+)s, useful ([\d.]+)s \[norm ([\d.]+)\+([\d.]+), "
    r"sieving ([\d.]+).*?factor ([\d.]+)")


def parse(text):
    rows = []
    for m in LINE.finditer(text):
        total, useful, n0, n1, sieving, factor = map(float, m.groups())
        rows.append(dict(total=total, useful=useful, norm=n0 + n1,
                         sieving=sieving, factor=factor))
    return rows


def main():
    text = open(sys.argv[1]).read() if len(sys.argv) > 1 else sys.stdin.read()
    rows = parse(text)
    if not rows:
        sys.stderr.write(
            "no '# Total cpu time ... sieving ... factor ...' line found.\n"
            "Build las WITHOUT -production and run it at verbose level 2 "
            "(-v -v); see las_profiling.md.\n")
        return 2
    # aggregate (a run may print several special-q blocks)
    norm = sum(r["norm"] for r in rows)
    siev = sum(r["sieving"] for r in rows)
    fact = sum(r["factor"] for r in rows)
    useful = sum(r["useful"] for r in rows)
    sf = siev + fact
    print("las timing over %d block(s):" % len(rows))
    print("  norm       : %10.2f s" % norm)
    print("  sieving    : %10.2f s" % siev)
    print("  cofactor   : %10.2f s" % fact)
    print("  useful cpu : %10.2f s" % useful)
    print()
    if sf > 0:
        print("  cofactorization / (sieving+cofactor) = %.1f%%" % (100 * fact / sf))
    if useful > 0:
        f = fact / useful
        print("  cofactorization / useful cpu         = %.1f%%" % (100 * f))
        print()
        print("  Amdahl ceiling if ONLY cofactorization is moved to the GPU:")
        for sp in (5, 10, 50, 1e9):
            cap = 1.0 / ((1 - f) + f / sp)
            tag = "inf" if sp > 1e8 else "%gx" % sp
            print("    kernel %-4s faster -> %.2fx end-to-end" % (tag, cap))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
