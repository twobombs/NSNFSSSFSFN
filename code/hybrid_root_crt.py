from sage.all import *
import multiprocessing as mp
from os import makedirs, urandom, listdir
import time

timeprint = lambda *args : print(f"{time.strftime("%Y-%M-%D")}:", *args)

take_first_n = lambda iterator, n : list(zip(*zip(iterator, range(n))))[0]

CRT_DIR = "/tmp/crt_demo"
PRECOMP_DIR = f"{CRT_DIR}/precomp"
RECONSTRUCT_DIR = f"{CRT_DIR}/reconstruct"

def lift(a, ainv, n):
    """
    Given a and its inverse ainv mod n,
    compute a^2 and its inverse mod n^2.
    """
    y = 3 - 2*a*ainv
    n2 = n**2
    return (a**2, (y*ainv**2)%n2, n2)

def write_bigint(bigint, file):
    with open(file, "wb") as f:
        f.write(bigint.to_bytes(bigint.nbits() // 8 + 1, byteorder="little"))

def read_bigint(file):
    with open(file, "rb") as f:
        return Integer(int.from_bytes(f.read(), byteorder="little"))

def write_precomp(p, basiselt, jobdir):
    makedirs(jobdir, exist_ok=True)
    with open(f"{jobdir}/p", "w") as f:
        print(p, file=f)
    write_bigint(basiselt, f"{jobdir}/basiselt")

def read_precomp(jobdir):
    with open(f"{jobdir}/p", "r") as f:
        p = int(f.read().strip())
    return (p, read_bigint(f"{jobdir}/basiselt"))

def precomputation_job(p, lg_ell, m, jobnum, datadir=PRECOMP_DIR):
    """
    Compute the CRT basis coefficient Q * (q^{-ell} mod p^ell),
    which is 1 mod p_i^ell (for this i) and 0 mod p_j^ell (for all j != i)
    lg_ell is the base-2 log of ell (which we assume is a power of 2)
    m is the product of all the p_i
    The output basiselt gets written (as raw bytes, little-endian) to {datadir}/{jobnum}/basiselt,
    and p gets written (in decimal) to {datadir}/{jobnum}/p.
    If demo == True, also write a random value somewhat less than p_i^ell to {datadir}/{jobnum}/residue,
    to be used in reconstruction
    """
    q = Integer(m // p)
    q_inv = pow(q, -1, p)
    Q = q**(2**lg_ell)
    a, ainv, modulus = (q, q_inv, p)
    for z in range(lg_ell):
        (a, ainv, modulus) = lift(a, ainv, modulus)
    # ainv is now q^-ell mod p^ell
    del a, modulus
    basiselt = Q * ainv
    del Q, ainv
    # TODO: might be better for basiselts to be in [-M/2, M/2) instead of [0,M).
    write_precomp(p, basiselt, f"{datadir}/{jobnum}")

def precomputation_M_job(m, ell, datadir=PRECOMP_DIR):
    """
    Precomputes the modulus M = m^ell, where m is the product of all the primes.
    Writes output (as raw bytes, little-endian) to {datadir}/M
    """
    makedirs(datadir, exist_ok=True)
    M = Integer(m)**ell
    with open(f"{datadir}/M", "wb") as f:
        f.write(M.to_bytes(M.nbits()//8+1, byteorder="little"))

def reconstruction_mult_job(jobnum, value=None, workdir=RECONSTRUCT_DIR, precompdir=PRECOMP_DIR, suffix="residue"):
    makedirs(f"{workdir}/{jobnum}", exist_ok=True)
    p, basiselt = read_precomp(f"{precompdir}/{jobnum}")
    if value is None:
        value = read_bigint(f"{precompdir}/{jobnum}/{suffix}")
    write_bigint(basiselt * value, f"{workdir}/{jobnum}/basiselt_times_{suffix}")

#####

def precomp(ps, lg_ell, datadir):
    ell = 2**lg_ell
    m = prod(ps)
    # Eventually we may use slurm and/or MPI, but for now just python multiprocessing
    pool = mp.Pool()
    try:
        pool.apply_async(precomputation_M_job, args=(m, ell, datadir))
        for i, p in enumerate(ps):
            pool.apply_async(precomputation_job, args=(p, lg_ell, m, i, datadir))
    finally:
        pool.close()
        pool.join()


def reconstruct(num_ps, datadir, suffix=""):
    pool = mp.Pool()
    try:
        for i in range(num_ps):
            pool.apply_async(reconstruction_mult_job, args=(i, None, datadir, datadir, f"residue{suffix}"))
    finally:
        pool.close()
        pool.join()
    timeprint("multiplications all done, starting singlethreaded addition")
    out = 0
    M = read_bigint(f"{datadir}/M")
    for i in range(num_ps):
        out += read_bigint(f"{datadir}/{i}/basiselt_times_residue{suffix}")
    # For our sizes, this is actually fine
    out %= M
    print("Result:", out)
    write_bigint(out, f"{datadir}/result{suffix}")
    with open(f"{datadir}/result{suffix}.txt","w") as f:
        print(out, file=f)
    return out

#####

def demo_precomputation(n=50_000, k=50, lg_ell=8, primebits=10):
    """
    n: bitsize of value to reconstruct
    k: number of primes (and number of jobs - 1)
    lg_ell: base 2 log of exponent ell
    primebits: bitsize of primes to use
    """
    ell = 2**lg_ell # for simplicity, a power of 2
    ps = take_first_n(primes(start=2**primebits, stop=2**(primebits+1)), k)
    m = prod(ps)
    assert prod(ps).nbits() * ell >= n, f"CRT will reconstruct {prod(ps).nbits() * ell} bits but need {n}"
    precomp(ps, lg_ell, PRECOMP_DIR)

def demo_make_residues(n=50_000, k=50, lg_ell=8):
    digits = ceil(n / log(10,2)) - 2
    value = Integer("7" * digits)
    assert n - 10 < value.nbits() < n-1
    for i in range(k):
        p, _ = read_precomp(f"{PRECOMP_DIR}/{i}")
        pn = p**(2**lg_ell)
        residue = value % pn
        write_bigint(residue, f"{PRECOMP_DIR}/{i}/residue")

def demo_reconstruction():
    """
    For now we're using the filesystem to communicate between jobs,
    which hopefully can be improved with openmp
    """
    num_ps = len(listdir(PRECOMP_DIR)) - 1
    pool = mp.Pool()
    try:
        for i in range(num_ps):
            pool.apply_async(reconstruction_mult_job, args=(i, None))
    finally:
        pool.close()
        pool.join()
    # TODO: the addition tree part, which should be pretty fast compared to the multiplications
    # For now we'll just do it singlethreaded, because the additions shouldn't be *that* large
    timeprint("multiplications all done, starting singlethreaded addition")
    out = 0
    M = read_bigint(f"{PRECOMP_DIR}/M")
    for i in range(num_ps):
        out += read_bigint(f"{RECONSTRUCT_DIR}/{i}/basiselt_times_residue")
    # For our sizes, this is actually fine
    out %= M
    print("Result:", out)
    write_bigint(out, f"{RECONSTRUCT_DIR}/result")
    return out
    

if __name__ == "__main__":
    from sys import argv
    if len(argv) >= 2 and argv[1] == "precomp":
        demo_precomputation()
    elif len(argv) >= 2 and argv[1] == "demo_padic":
        demo_make_residues()
    elif len(argv) >= 2 and argv[1] == "reconstruct":
        demo_reconstruction()
    else:
        print("usage: crt.py {precomp|demo_padic|reconstruct}")
        exit(1)
