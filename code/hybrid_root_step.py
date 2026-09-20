from sage.all import *

from os.path import exists
from os import makedirs, mkdir
from subprocess import run
import json

from multiprocessing import Pool

from cado_sage import CadoPolyFile
import hybrid_root_crt as crt
from hybrid_root_find_primes import find_inert_primes

from misc_tools import fast_persistent_load
from timing import timeprint

SAGE="/usr/local/sagemath/10.7/bin/sage"

# Assumes the following files already exist:
# TTplus.sobj and TTminus.sobj (created by "sage run.py ... serialize_ttplus_ttminus")
# gamma-fac.sobj

# STEPS
# 1) Create crt_primes if it doesn't exist by calling hybrid_root_find_primes.py
# 2) Create one workdir per prime, in {datadir}/padic_jobs/{jobno}
#    Populate each jobdir with p
# 3) Compute necessary value of lg_ell
# 4) Launch each padic job (hybrid_padic_job.py). Wait for them all to finish. Now in each jobdir there are residue_0 through residue_5
# 5) Run CRT precomp to produce crt basis elts and M (= prod p_i^ell), if it hasn't already been run
# 6) Run CRT reconstruction on each coefficient separately, giving a polynomial mod M (and mod f(x))
# 7) Do rational reconstruction on the polynomial as follows:
# Return a(x)/b(x)

def run_hybrid_root_step(bits_to_reconstruct, datadir, limit_jobs=0):
    with open(f"{datadir}/params.json","r") as f:
        params = json.load(f)
    polyfile = params['files']['POLYFILE']
    poly = CadoPolyFile(polyfile)
    poly.read()
    f = poly.f[1]

    ps = create_crt_primes(datadir, f, params['parameters']['e'])

    create_workdirs(datadir, ps)

    ell = Integer(ceil(bits_to_reconstruct / prod(ps).nbits()))
    lg_ell = Integer(ceil(log(ell,2)))
    ell = 2**lg_ell
    assert sum(ell*log(p,2) for p in ps) > bits_to_reconstruct
    timeprint(f"Will run p-adic steps up to p^(2^{lg_ell}),"
          f" with {len(ps)} primes in parallel,"
          f" for a total of {sum(ell*log(p,2) for p in ps).n():.1f} bits")
    
    if exists(f"{datadir}/padic_jobs/M"):
        timeprint("M file exists, assuming CRT precomp has already been done")
    else:
        timeprint("Running CRT precomp...")
        crt.precomp(ps, lg_ell, datadir=f"{datadir}/padic_jobs")
        timeprint("Done with CRT precomp")
    M = crt.read_bigint(f"{datadir}/padic_jobs/M")
    
    timeprint("Launching parallel p-adic jobs")
    num_processes = limit_jobs if limit_jobs > 0 else None
    with Pool(num_processes) as pool:
        def print_from_job(i,p,msg,skip_x=False):
            return lambda x : timeprint(f"Job {i} (p={p}) {msg}", x if not skip_x else "")
        for (i,p) in enumerate(ps):
            pool.apply_async(run_padic_job,
                (i, p, lg_ell, datadir),
                callback=print_from_job(i,p,"Done!", True),
                error_callback=print_from_job(i,p,"Error:"),
            )
        pool.close()
        pool.join()
    timeprint("All parallel p-adic jobs done")

    timeprint("Start CRT reconstruction")
    cx = [None] * f.degree()
    for i in range(f.degree()):
        timeprint(f"Reconstructing x^{i} coefficient")
        cx[i] = crt.reconstruct(len(ps), datadir=f"{datadir}/padic_jobs", suffix=f"_{i}")
        timeprint(f"Done reconstructing x^{i} coefficient")

    timeprint("Finished CRT reconstruction")
    timeprint("Starting rational reconstruction")
    return rational_reconstruction(cx, f, M)

def rational_reconstruction(cx, f, M):
    """
    cx: list of coefficients of a polynomial mod M, f(x)
    f: polynomial
    M: modulus
    """
    # We have c(x) = a(x)/b(x) mod (M, f(x)) for small unknown a(x),b(x).
    # Represent the polynomials as square coefficient matrices (where row i = x^i c(x) mod f(x))
    # Do LLL on the matrix
    # [c(x) 1]
    # [  M   ]
    # to recover
    # [a(x) b(x)]
    R = Zmod(M).extension(f)
    cx = R(cx)
    alpha = R.gen()
    n = f.degree()
    L = block_matrix(2, 2,
        [
            matrix(ZZ,[(cx * alpha**i).list() for i in range(n)]), 1,
            M, 0
        ]).LLL()
    coeffs = L[0].list()
    ZP = ZZ['x']
    rn = ZP(coeffs[:n])
    rd = ZP(coeffs[n:])
    assert rn(alpha) - cx * rd(alpha) == 0
    return rn / rd

def create_crt_primes(datadir, f, e):
    if not exists(f"{datadir}/crt_primes"):
        ps = find_inert_primes(f,
            fast_persistent_load(f"{datadir}/TTplus.sobj"),
            fast_persistent_load(f"{datadir}/TTminus.sobj"),
            avoid={e}
        )
        with open(f"{datadir}/crt_primes","w") as fd:
            print(ps, file=fd)
        return ps
    else:
        timeprint("using existing crt_primes")
        with open(f"{datadir}/crt_primes","r") as fd:
            return [Integer(x) for x in fd.read().strip().lstrip("[").rstrip("]").split(", ")]

def create_workdirs(datadir, ps):
    makedirs(f"{datadir}/padic_jobs", exist_ok=True)
    for (i, p) in enumerate(ps):
        makedirs(f"{datadir}/padic_jobs/{i}", exist_ok=True)
        with open(f"{datadir}/padic_jobs/{i}/p","w") as f:
            print(p, file=f)
    

def run_padic_job(i, p, lg_ell, datadir):
    with open(f"{datadir}/padic_jobs/{i}/padic.stdout",'a') as stdout:
        with open(f"{datadir}/padic_jobs/{i}/padic.stderr",'a') as stderr:
            timeprint("----- run_padic_job start -----", file=stdout, flush=True)
            timeprint("----- run_padic_job start -----", file=stderr, flush=True)
            run([SAGE,
                "hybrid_root_padic_job.py",
                str(p), str(lg_ell),
                f"{datadir}/padic_jobs/{i}",
                f"{datadir}/params.json",
                f"{datadir}/TTplus.txt",
                f"{datadir}/TTminus.txt",
                f"{datadir}/gamma-fac.sobj"],
                stdout=stdout, stderr=stderr
            )

