from sage.all import *

from helpers import LinalgOutput
from cado_sage import CadoPolyFile
from misc_tools import fast_persistent_save, fast_persistent_load
from timing import timeprint

def find_inert_primes(f, TTplus=[], TTminus=[], avoid=set()):
    # Let's find 100 inert primes, then filter out any that fail the TTplus and TTminus checks
    timeprint("Finding candidate inert primes")
    candidates = []
    for i, p in enumerate(primes(2, infinity)):
        if i > 5000:
            raise ValueError("Cannot find an inert prime."
                             " Maybe this polynomial has a peculiar Galois group, or maybe need to search longer")
        fp = f.change_ring(GF(p))
        if fp.degree() < f.degree():
            continue
        if p in avoid:
            print(f"not using p={p} because"
                  " we listed it among the primes to avoid in the lift")
            continue
        if not fp.is_irreducible():
            continue
        candidates.append(p)
        if len(candidates) >= 100:
            break
    timeprint("Candidate inert primes:", candidates)

    out = []
    for p in candidates:
        timeprint(f"Checking {p}...")
        fp = f.change_ring(GF(p))
        Fpn = GF(p**(f.degree()), 'alpha_p', modulus=f)
        if any(a%p == 0 and b%p == 0 for (a,b),k in TTplus):
            continue
        if any(a%p == 0 and b%p == 0 for (a,b),k in TTminus):
            continue
        out.append(p)
    timeprint("Primes to use:", out)
    return out

if __name__ == "__main__":
    from sys import argv
    if len(argv) < 5:
        print("Usage: hybrid_root_find_primes.py {f.poly} {TTplus.sobj} {TTminus.sobj} {primelist.out}")
        exit(1)
    polyfile = argv[1]
    ttplusfile = argv[2]
    ttminusfile = argv[3]
    outfile = argv[4]
    poly = CadoPolyFile(polyfile)
    timeprint("Reading f.poly...")
    poly.read()
    timeprint("Done reading f.poly")
    timeprint("Loading TTplus...")
    TTplus = fast_persistent_load(ttplusfile)
    timeprint("Loading TTminus...")
    TTminus = fast_persistent_load(ttminusfile)
    timeprint("Done loading everything")
    primes_to_use = find_inert_primes(poly.f[1], TTplus, TTminus, avoid={65537})
    with open(outfile, "w") as f:
        print(primes_to_use, file=f)
