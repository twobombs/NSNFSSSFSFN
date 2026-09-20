#!/usr/bin/env sage

from sage.all import *
from cado.scripts import descent
from helpers import timeprint
import multiprocessing
import argparse
import json
import json_custom
import re
import itertools
import os
import random
import glob
from helpers import silent_remove, masked_target, construct_S, CadoExplainRenumberFile, generate_or_load_target
from helpers import truncate_S, sanity_check_descent_outfile, handle_bottom_special_q, handle_bottom_special_q_composites
from contextlib import redirect_stdout,redirect_stderr
from candy import print_command_line,warning_message,major_message,error_message
from cado_sage import CadoPolyFile
from cado_nfs_binaries import CadoNFS,CadoNFSBinaries
import time
import descent_ecm_utils
import functools


# XXX This is not unified with Params in run.py, which is unfortunate.
# This type definition is quite like argparse.Namespace, by the way.
class Params(dict):
    __getattr__ = dict.get

if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='descent_rb_helper.py')
    topparser.add_argument('--params',dest='params',required=True)
    topparser.add_argument('--outfile',dest='outfile',required=True)
    topparser.add_argument('--side',dest='side',required=True)
    topparser.add_argument('--q',dest='q',required=True)
    topparser.add_argument('--rho',dest='rho')
    topparser.add_argument('--seed',dest='seed')
    topparser.add_argument('--strategy',dest='strategy')
    topparser.add_argument('--overwrite-lpb0',dest='overwrite_lpb0')
    topparser.add_argument('--overwrite-lpb1',dest='overwrite_lpb1')
    topargs = topparser.parse_args()

    assert topargs.strategy in ['P','C','C1','C2']
    # either Polynomials or Composites

    params = Params(json.loads(open(topargs.params,'r').read()))

    e = params.parameters['e']
    N = params.parameters['N']
    ZN = Integers(N)
    LPB0 = params.parameters['LPB0']
    LPB1 = params.parameters['LPB1']
    BOUNDR = params.BOUNDR
    BOUNDA = params.BOUNDA

    side = int(topargs.side)
    assert side in [0,1]

    q = Integer(topargs.q)
    #if side == 0:
        #assert q > BOUNDR
    #else:
        #assert q > BOUNDA

    if topargs.rho is not None:
        rho = Integer(topargs.rho)
    else:
        rho = None

    if topargs.seed is not None:
        working_pfx = params.dirs['DESC'] + 'seed' + topargs.seed[:10] + "/rbjobs/"
    else:
        working_pfx = params.dirs['DESC']

    CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])
    og_poly = CadoPolyFile(params.files['POLYFILE']); og_poly.read()

    f = og_poly.f[1]
    g = og_poly.f[0]

    if topargs.overwrite_lpb0:
        used_lpb0 = int(topargs.overwrite_lpb0)
        BOUNDR = 2**used_lpb0
    else:
        used_lpb0 = None

    if topargs.overwrite_lpb1:
        used_lpb1 = int(topargs.overwrite_lpb1)
        BOUNDA = 2**used_lpb1
    else:
        used_lpb1 = None

    # Polynomials strategy
    if topargs.strategy == 'P':

        a, b, fac0, fac1, winner, newg, newf, rr, new_shared_root = handle_bottom_special_q(
            params, q, side, working_pfx, rho, None, used_lpb0, used_lpb1
        )

        a = Integer(a)
        b = Integer(b)

        if side == 0:
            assert prod(fac0)*q == ZZ(g(a/b)*b**g.degree()).abs()
            assert prod(fac1) == ZZ(f(a/b)*b**f.degree()).abs()
            timeprint("Yay! The relation works under the original polynomials.")
            fac0.append(q)

        else:
            assert prod(fac0) == ZZ(g(a/b)*b**g.degree()).abs()
            assert prod(fac1)*q == ZZ(f(a/b)*b**f.degree()).abs()
            timeprint("Yay! The relation works under the original polynomials.")
            fac1.append(q)

    # Composites strategy
    elif topargs.strategy in ['C','C1','C2']:

        a, b, fac0, fac1 = handle_bottom_special_q_composites(
            params, q, side, working_pfx, rho, topargs.strategy, used_lpb0, used_lpb1
        )

        print("a", str(a))
        print("b", str(b))
        print("fac0", str(fac0))
        print("fac1", str(fac1))

        # In fact casting to Integer is necessary to pass the assertions. Blegh!
        a = Integer(a)
        b = Integer(b)

        fac0_int = [ Integer(y,16) for y in fac0 ]
        fac1_int = [ Integer(y,16) for y in fac1 ]

        assert prod(fac0_int) == ZZ(g(a/b)*b**g.degree()).abs()
        assert prod(fac1_int) == ZZ(f(a/b)*b**f.degree()).abs()
        timeprint("Yay! The relation works.")

        fac0 = fac0_int
        fac1 = fac1_int

    assert len(fac0) > 0 and len(fac1) > 0

    for rat_fac in fac0:
        assert (rat_fac == q and side == 0) or (rat_fac < BOUNDR)

    for alg_fac in fac1:
        assert (alg_fac == q and side == 1) or (alg_fac < BOUNDA)

    # At this point, we have: a, b, fac0, fac1
    # We write it in the way that the las_descent parser will recognize
    with open(topargs.outfile, "w") as out:

        if side == 0:
            rho = g.roots(GF(q))[0][0]

        out.write("# Taking decision on ")
        out.write(str(q.nbits()))
        out.write(f"@{side} side-{side} q=")
        out.write(str(q))
        out.write("; rho=")
        out.write(str(rho))
        out.write("\n")

        out.write("Taken: " + str(a) + "," + str(b))
        out.write(":")
        out.write(",".join([ hex(x)[2:] for x in fac0 ]))
        out.write(":")
        out.write(",".join([ hex(x)[2:] for x in fac1 ]))
        out.write("\n")



# -------------------------
# We have a few strategies.
# 1. For a given q, sieve over qq' for small primes q'.
#    Factor out q' at the end.
#    The benefit here is we can try some fresh sieving lattices.
# 2. Do the transform_polys_by_q trick.
#    This should work for either side.
# In general we could also bias the las_descent hintfile to descend side 1
# and deal with more side 0 qs in this special step.
