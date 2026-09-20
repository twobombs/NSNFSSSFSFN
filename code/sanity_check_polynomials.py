from run import *
from sage.all import *
import os
import glob
import time
import argparse
from cado_sage import CadoPolyFile
import padic_eth_root


ABSOLUTELY_NOT = 0
SM_APPEND_WILL_BE_WEIRD = 1
SEEMS_OK = 2


def check_polyfile(polyfile):
    poly = CadoPolyFile(polyfile)
    poly.read()

    e = 65537
    K = poly.K[1]
    alpha = K.gen()
    Kw = poly.nt[1]
    sm_maps = Kw.schirokauer_maps(e)

    OK = K.maximal_order()
    f = K.defining_polynomial()
    lc = f.leading_coefficient()
    alpha_hat = K.gen() * lc
    f_hat = alpha_hat.minpoly()
    disc = f.discriminant()

    if gcd(lc, e) != 1:
        return ABSOLUTELY_NOT

    if valuation(disc, e) >= 2:
        return ABSOLUTELY_NOT

    try:
       p = padic_eth_root.find_inert_prime(poly.f[1])
       assert p > 2
    except ValueError:
        return ABSOLUTELY_NOT

    ideal_fac = OK.fractional_ideal(e).factor()

    print(str(list(ideal_fac)))

    if len(list(ideal_fac)) > 1:
        return SM_APPEND_WILL_BE_WEIRD

    assert gcd(lc, e) == 1
    assert valuation(disc, e) < 2

    return SEEMS_OK


def check_sm_append(cands):
    for c in cands:
        command = [
            "build-457bd1173/build-no-mpi-large-lpb/filter/sm_append",
            "-ell", "65537",
            "-in", "n1024/select/good/smallsmtest.forsm",
            "-out",
            f"n1024/select/good/smallsmtest.forsm.{c}.withsm",
            "-poly",
            f"n1024/select/good/polyselect.out.poly.option.{c}",
            "-b", "256",
            "-t", "88",
            "-nsm", "0,6"
        ]
        subprocess.run(command)


def make_fbs(cands):
    for c in cands:
        command = [
            "build-457bd1173/build-no-mpi-large-lpb/sieve/makefb",
            "-lim", str(2**31),
            "-t", "88",
            "-poly",
            f"n1024/select/good/polyselect.out.poly.option.{c}",
            "-out",
            f"n1024/select/good/fb.option.{c}"
        ]
        subprocess.run(command)


if __name__=='__main__':

    candidates = [
        64, 26, 88, 57, 48, 90, 46, 76, 50, 42, 54, 11, 81, 36, 94, 0, 1, 2
    ]

    # check_sm_append(candidates)

    make_fbs(candidates)

    #ngood = 0
    #candidates = glob.glob('n1024/select/polyselect.out.poly.option.[0-9][0-9]')
    #for cand in candidates:
    #    print("\n")
    #    res = check_polyfile(cand)
    #    print(str(cand))
    #    print(str(res))

    #    if res == SEEMS_OK:
    #        ngood += 1
    #        subprocess.run(
    #            ["cp", cand, "n1024/select/good/"]
    #        )
    #        subprocess.run(
    #            ["cp", cand + ".only-side1", "n1024/select/good/"]
    #        )

    #print(f"found {ngood} reasonable polynomials")


    # We should run tests for:
    # las at reasonable parameters
    # las_descent at large parameters, to see how skewed of a lattice we can have
    # sm_append, to make sure things are nonzero, and cado_sage while we're at it
    # Montgomery / padic to the extent that it's doable
