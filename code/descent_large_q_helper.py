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
from helpers import truncate_S, sanity_check_descent_outfile, handle_very_large_special_q
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
    topparser = argparse.ArgumentParser(prog='descent_large_q_helper.py')
    topparser.add_argument('--params',dest='params')
    topparser.add_argument('--existing-init-data',dest='existing_init_data',required=True)
    topparser.add_argument('--lq',dest='lq',required=True)
    topargs = topparser.parse_args()

    params = Params(json.loads(open(topargs.params,'r').read()))

    e = params.parameters['e']
    N = params.parameters['N']
    ZN = Integers(N)
    target = generate_or_load_target(params)

    init_starttime = time.time()
    timeprint("Starting large q: " + str(topargs.lq))

    CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])
    og_poly = CadoPolyFile(params.files['POLYFILE']); og_poly.read()

    try:
        with open(topargs.existing_init_data, "r") as fp:
            init_dict = json.load(fp)

        todofile = init_dict['todofilename']
        u = Integer(init_dict['u'])
        v = Integer(init_dict['v'])
        u_fac = init_dict['u_fac']
        v_fac = init_dict['v_fac']
        mask = Integer(init_dict['mask'])
        seed = Integer(init_dict['seed'])
        firstrelsfile = None    # Note this is always None anyway
        largeq_rels_file = init_dict['largeq_rels_file']
    except Exception as ex:
        print(f"Uh oh! There was a problem parsing the init_data in {topargs.existing_init_data}.")
        raise ex
        sys.exit(1)
        os._exit(1)

    pfx = params.dirs['DESC'] + 'largeq.' + str(seed)[:10] + "."
    f = og_poly.f[1]
    g = og_poly.f[0]

    lq = Integer(topargs.lq)
    assert lq.nbits() > params.parameters.get('LARGEQ', 90)
    assert lq in (u_fac + v_fac)

    # rho = ZZ(g.roots(GF(lq))[0][0])   # done in handle_very_large_special_q
    timeprint(f"Handling the large q: {lq} of {lq.nbits()} bits")
    a, b, fac0, fac1, winning_rel, newg, newf, rr, new_shared_root = handle_very_large_special_q(params, lq, 0, pfx, rho=None)
    print("winning_rel: ", str(winning_rel))
    a = Integer(a)
    b = Integer(b)

    assert prod(fac0)*lq == ZZ(g(a/b)*b**g.degree()).abs()
    assert prod(fac1) == ZZ(f(a/b)*b**f.degree()).abs()

    todo_qs = []
    for fac in fac0:
        if fac > params.BOUNDR:
            todo_qs.append((0, fac, 0))
    for fac in fac1:
        if fac > params.BOUNDA:
            bb = GF(fac)(b)
            if bb == 0:
                continue
            r = ZZ(a/bb)
            todo_qs.append((1, fac, r))

    with open(todofile, "a") as f:
        for (side, q, r) in todo_qs:
            if side == 0:
                f.write(f"{side} {q}")
            elif side == 1:
                f.write(f"{side} {q} {r}")
            f.write("\n")

    with open(largeq_rels_file, "a") as f:
        f.write("# Taking decision on ")
        f.write(str(lq.nbits()))
        f.write("@0 side-0 q=")
        f.write(str(lq))
        f.write("; rho=")
        rho = g.roots(GF(lq))[0][0]
        f.write(str(rho))
        f.write("\n")

        f.write("Taken: " + str(a) + "," + str(b))
        f.write(":")
        f.write(",".join([ hex(x)[2:] for x in fac0 ]))
        f.write("," + hex(lq)[2:] + ":")
        f.write(",".join([ hex(x)[2:] for x in fac1 ]))
        f.write("\n")

    init_endtime = time.time()
    descent_init_time = round(init_endtime-init_starttime)

    timeprint("Success! We found a relation for large q: " + str(lq))
    timeprint("Data has been appended to: " + str(todofile) + " and " + str(largeq_rels_file))
