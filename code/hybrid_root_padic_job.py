from sage.all import *
from timing import timeprint as timeprint_

timeprint = lambda *args, **kwargs : timeprint_(*args, flush=True, **kwargs)

def load_TTplusorminus(TTplusorminus_txt):
    """ Returns an iterator that yields ((a,b),k) tuples"""
    with open(TTplusorminus_txt,"r") as f:
        for line in f:
            a,b,k = [Integer(x) for x in line.strip().split()]
            yield ((a,b),k)

def padic_root_job(TTplus, TTminus, f, p, e, gamma_fac, lg_ell):
    timeprint("padic_root_job start")
    # will be run as its own process, with its own workdir
    # inputs:
    #   e, f
    #   TTplus, TTminus, gamma_fac
    #   p, lg_ell
    # output:
    #  list of the coefficients mod p^(2^lg_ell) of the eth root of (prod(TTplus)/prod(TTminus))/prod(gamma_fac)**e mod f(x)
    
    ZP = ZZ['x']

    d = f.degree()
    # Ri is (Z/p^(2^i))[x]/f(x)
    R0 = GF(p**(f.degree()), 'alpha_p', modulus=f)
    alpha0 = R0.gen()
    # We know how far we need to go, so let's just compute y and z at that precision from the beginning
    Rl = Zmod(p**(2**lg_ell)).extension(f)
    alphal = Rl.gen()
    yl = prod((a-b*alphal)**k for (a,b),k in TTplus)
    zl = prod((a-b*alphal)**k for (a,b),k in TTminus)
    y0 = R0(yl.list())
    z0 = R0(zl.list())
    u0 = 1/z0
    # Our iteration will compute the *inverse* e-th root of y/z
    rf0 = 1
    if gamma_fac is not None:
        for fac in gamma_fac:
            nf_elt = fac[0]
            exp = fac[1]
            poly_alpha0 = nf_elt.polynomial().change_ring(R0)(alpha0)
            rf0 *= (1/(poly_alpha0**exp))       # whole product is 1/gamma

    r0 = (z0 / y0 / rf0**e).nth_root(e)
    # r0 = gamma * (z / y)^(1/e)

    y,z,u,r,R,pk,i,rf = y0,z0,u0,r0,R0,p,0,rf0
    
    while i < lg_ell:
        # invariants:
        assert pk == p**(2**i)
        # R = (Z/p^(2^i))[x]/f(x)
        assert y in R
        assert z in R
        assert u in R
        assert r in R
        assert u * z == 1
        assert r**e * y * rf**e * u == 1

        i += 1
        pk = pk * pk
        R = Integers(pk).extension(f)

        timeprint(f"Lifting to {p}^(2^{i})")

        r = ZP(r.list())(R.gen())
        # lift the preinverse of z, too.
        u = ZP(u.list())(R.gen())

        # Newton step
        y = R(yl.list())
        z = R(zl.list())

        # if u*z = 1+pb, then the higher order inverse of z is u*(1-pb)
        # IOW, a Newton step on u: u becomes u*(1-(u*z-1))
        u = u * (1 - (u * z - 1))

        rf = 1

        if gamma_fac is not None:
            for fac in gamma_fac:
                nf_elt = fac[0]
                exp = fac[1]
                poly_R_gen = nf_elt.polynomial().change_ring(R)(R.gen())
                rf *= (1/(poly_R_gen**exp))     # whole product is 1/gamma

        # This is the Newton iteration on the function f(x) = y/z - x^-e
        r += r * (1 - r**e * rf**e * y * u) / e

    # remember that we have computed the _inverse_ of the e-th
    # root.
    root = 1/r

    # there are d coefficients, we'll separately crt-reconstruct each one
    out = root.list()
    return out


if __name__ == "__main__":
    from sys import argv
    if len(argv) < 8:
        print(f"Usage: {argv[0]} p lg_ell workdir params.json TTplus.txt TTminus.txt gamma_fac.sobj")
        exit(1)

    import json
    from cado_sage import CadoPolyFile
    from misc_tools import fast_persistent_load
    from hybrid_root_crt import write_bigint
    from os import makedirs

    timeprint("START!")

    p = Integer(argv[1])
    assert p.is_prime()
    lg_ell = Integer(argv[2])
    workdir = argv[3]
    params_file = argv[4]
    TTplus_file = argv[5]
    TTminus_file = argv[6]
    gammafac_file = argv[7]
    
    makedirs(workdir, exist_ok=True)
        
    with open(params_file, 'r') as f:
        params = json.load(f)
    e = params['parameters']['e']
    polyfile = params['files']['POLYFILE']

    timeprint("Loading f.poly...")
    poly = CadoPolyFile(polyfile)
    poly.read()
    f = poly.f[1]

    timeprint("Opening TTplus...")
    TTplus = load_TTplusorminus(TTplus_file)
    timeprint("Opening TTminus...")
    TTminus = load_TTplusorminus(TTminus_file)
    timeprint("Loading gamma_fac...")
    gamma_fac = fast_persistent_load(gammafac_file)
    timeprint("Everything is loaded!")

    out = padic_root_job(TTplus, TTminus, f, p, e, gamma_fac, lg_ell)

    timeprint("Finished padic_root_job!")

    for (i, residue) in enumerate(out):
        write_bigint(residue.lift(), f"{workdir}/residue_{i}")

    timeprint("FINISHED!")
