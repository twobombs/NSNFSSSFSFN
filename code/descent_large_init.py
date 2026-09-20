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
import sys
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
    topparser = argparse.ArgumentParser(prog='descent_large_init.py')
    topparser.add_argument('--params',dest='params')
    topparser.add_argument('--seed',dest='seedval')
    topargs = topparser.parse_args()

    params = Params(json.loads(open(topargs.params,'r').read()))
    # assert(params.descent_ecm_init)

    e = params.parameters['e']
    N = params.parameters['N']
    ZN = Integers(N)
    target = generate_or_load_target(params)

    init_starttime = time.time()
    timeprint("Starting large descent initialization!")

    q0 = params.parameters.get('desc.ecm.q0', 2147483648)
    max_ecm_trials = params.parameters.get('desc.ecm.max_ecm_trials', 50)
    nq = params.parameters.get('desc.ecm.nq', 10)
    qdiff = params.parameters.get('desc.ecm.qdiff', 10000000)
    ecm_nthreads = params.parameters.get('desc.ecm.ecm_nthreads', 24)
    init_tkewness = params.parameters.get('desc.ecm.tkewness', 2000000000)

    ecminit_lim0 = params.parameters.get('desc.ecm.lim0', 2147483648)
    ecminit_lim1 = params.parameters.get('desc.ecm.lim1', 2147483648)
    ecminit_lpb0 = params.parameters.get('desc.ecm.lpb0', 110)
    ecminit_lpb1 = params.parameters.get('desc.ecm.lpb1', 110)
    ecminit_mfb0 = params.parameters.get('desc.ecm.mfb0', 300)
    ecminit_mfb1 = params.parameters.get('desc.ecm.mfb1', 300)
    ecminit_I_sieving = params.parameters.get('desc.ecm.I_sieving', 12)
    ecminit_ncurves0 = params.parameters.get('desc.ecm.ncurves0', 1)
    ecminit_ncurves1 = params.parameters.get('desc.ecm.ncurves1', 1)
    ecm_smoothB = 2 ** params.parameters.get('desc.ecm.smoothB_log', 135)
    ecm_cofacB = 2 ** params.parameters.get('desc.ecm.cofacB_log', 260)
    ecm_B1 = params.parameters.get('desc.ecm.B1', 500000)
    ecm_ncurves = params.parameters.get('desc.ecm.ncurves', 600)

    found = False
    out = None
    ntrial = 1

    CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])

    init_poly_bound = N.bit_length() // 2 + 20

    random.seed(topargs.seedval)

    init_polyfile = params.dirs['DESC'] + 'ecminit.' + str(topargs.seedval)[:10] + '.poly'
    init_fbfile = params.dirs['DESC'] + 'ecminit.' + str(topargs.seedval)[:10] + '.fb'
    mask = 1
    zz = 1
    gg = None

    while True:
        mask = random.randrange(N)
        zz = (pow(mask, e, N) * target) % N
        gg = descent_ecm_utils.myxgcd(int(zz), int(N), int(init_tkewness))
        if (gg[0][0].bit_length() < init_poly_bound and
                gg[1][0].bit_length() < init_poly_bound and
                gg[0][1].bit_length() < init_poly_bound and
                gg[1][1].bit_length() < init_poly_bound):
            # we're happy
            break
        print("Skewed reconstruction. Let's randomize the input.")

    with open(init_polyfile, 'w') as f:
        f.write("n: %d\n" % N)
        f.write("skew: 1\n")
        f.write("c1: %d\n" % gg[0][0])
        f.write("c0: %d\n" % gg[1][0])
        f.write("Y1: %d\n" % gg[0][1])
        f.write("Y0: %d\n" % gg[1][1])

    init_poly = CadoPolyFile(init_polyfile); init_poly.read()
    og_poly = CadoPolyFile(params.files['POLYFILE']); og_poly.read()
    timeprint("Successfully chose an initialization polynomial.")
    timeprint("Beginning las and ecm filtering...")

    while not found:
        las_outfile = params.dirs['DESC'] + 'ecminit.' + str(topargs.seedval)[:10] + '.las.trial' + str(ntrial)
        print("**** Trial number " + str(ntrial) + " ****")
        q = q0 + ZZ.random_element(qdiff)
        print(f"Trying with {nq} q's after {q}")

        CadoNFS("sieve/las",
            "-poly", 'POLY',
            '-lim0', str(ecminit_lim0),
            '-lim1', str(ecminit_lim1),
            '-lpb0', str(ecminit_lpb0),
            '-lpb1', str(ecminit_lpb1),
            '-mfb0', str(ecminit_mfb0),
            '-mfb1', str(ecminit_mfb1),
            '-ncurves0', str(ecminit_ncurves0),
            '-ncurves1', str(ecminit_ncurves1),
            '-I', str(ecminit_I_sieving),
            '-nq', str(nq),
            '-q0', str(q),
            '-t', str(ecm_nthreads),
            "-v",
            #"-bkthresh1", str(int(min(ecminit_lim0, ecminit_lim1)/2)),
            '-batch-print-survivors', las_outfile + '.survivors',
            "-out", 'OUTFILE',
            outputs={
                'OUTFILE': las_outfile,},
            inputs={
                'POLY': init_polyfile,
            }
        )

        # May produce a bunch of .int output files
        wildcard = las_outfile + '.survivors*'
        for filename in glob.glob(wildcard):
            print("TRYING FILE", filename)
            found, out, survivor_line = descent_ecm_utils.try_file(filename,
                                                                    ecm_smoothB, ecm_cofacB,
                                                                    ecm_B1, ecm_ncurves,
                                                                    44)

            good_large_qs = False

            if found:
                n1 = out[0]
                n2 = out[1]

                cofac0 = Integer(survivor_line.split()[2].strip())
                cofac1 = Integer(survivor_line.split()[3].strip())

                abcd = survivor_line.strip().split()
                a = int(abcd[0], 10)
                b = int(abcd[1], 10)

                Num = a * gg[0][0] + b * gg[1][0]
                Den = a * gg[0][1] + b * gg[1][1]

                factNum = []
                for ff in n2:
                    if ff.is_prime():
                        factNum.append(ff)
                    else:
                        # weirdly, this does happen
                        for fff in ff.factor():
                            for mult in range(0, fff[1]):
                                factNum.append(fff[0])

                factDen = []
                for ff in n1:
                    if ff.is_prime():
                        factDen.append(ff)
                    else:
                        # weirdly, this does happen
                        for fff in ff.factor():
                            for mult in range(0, fff[1]):
                                factDen.append(fff[0])

                for ff in (Num/cofac1).factor():
                    for mult in range(0, ff[1]):
                        factNum.append(ff[0])

                for ff in (Den/cofac0).factor():
                    for mult in range(0, ff[1]):
                        factDen.append(ff[0])

                large_q = []
                for fact in factNum + factDen:
                    if fact.nbits() > params.parameters.get('LARGEQ', 90):
                        large_q.append(fact)

                og_polyfile = params.files['POLYFILE']
                og_poly = CadoPolyFile(og_polyfile); og_poly.read()
                m = og_poly.m
                f = og_poly.f[1]
                g = og_poly.f[0]

                good_large_qs = True

                for q in large_q:
                    rho = ZZ(g.roots(GF(q))[0][0])
                    newg, newf, coeff = descent_ecm_utils.transform_polys_by_q(g, f, q, rho, side=0)
                    rr = Integer(newg.resultant(newf))
                    c0 = Integer(newg.list()[0])
                    c1 = Integer(newg.list()[1])
                    if gcd(c1, rr) > 1:
                        # sad
                        good_large_qs = False

            if found and good_large_qs:
                break

        ntrial = ntrial+1
        if ntrial >= max_ecm_trials:
            timeprint(f"Uh oh! Failed over {ntrial-1} attempts at descent ecm init.")
            sys.exit(0)

    time_ecm_done = time.time()
    timeprint(f"Finished ecm filtering! So far we've taken time {time_ecm_done-init_starttime}")

    print("found", str(found))
    print("out", str(out))
    print("survivor_line", survivor_line)

    n1 = out[0]
    n2 = out[1]

    cofac0 = Integer(survivor_line.split()[2].strip())
    cofac1 = Integer(survivor_line.split()[3].strip())

    assert(prod(n1) == cofac0)
    assert(prod(n2) == cofac1)

    abcd = survivor_line.strip().split()
    a = int(abcd[0], 10)
    b = int(abcd[1], 10)

    timeprint(f"The winning (a,b) pair is: {a},{b}")

    #Num = ZZ(init_poly.f[0](a/b)*b)
    #Den = ZZ(init_poly.f[1](a/b)*b)
    #assert ((mask**e * target)*Den-Num) % N == 0

    Num = a * gg[0][0] + b * gg[1][0]
    Den = a * gg[0][1] + b * gg[1][1]
    # zz = Num/Den mod N

    assert (Num % cofac1) == 0
    assert (Den % cofac0) == 0
    assert (zz * Den - Num) % N == 0

    factNum = []
    for ff in n2:
        if ff.is_prime():
            factNum.append(ff)
        else:
            # weirdly, this does happen
            for fff in ff.factor():
                for mult in range(0, fff[1]):
                    factNum.append(fff[0])

    assert prod(factNum) == prod(n2)

    factDen = []
    for ff in n1:
        if ff.is_prime():
            factDen.append(ff)
        else:
            # weirdly, this does happen
            for fff in ff.factor():
                for mult in range(0, fff[1]):
                    factDen.append(fff[0])

    assert prod(factDen) == prod(n1)

    for ff in (Num/cofac1).factor():
        for mult in range(0, ff[1]):
            factNum.append(ff[0])

    for ff in (Den/cofac0).factor():
        for mult in range(0, ff[1]):
            factDen.append(ff[0])

    print("factNum", str(factNum))
    print("factDen", str(factDen))

    assert(prod(factNum) == abs(Num))
    assert(prod(factDen) == abs(Den))

    factNum.sort()
    factDen.sort()
    print("bits in factors")
    print(str([ x.nbits() for x in factNum ]))
    print(str([ x.nbits() for x in factDen ]))

    large_q = []
    small_q = []
    for fact in factNum + factDen:
        if fact.nbits() > params.parameters.get('LARGEQ', 90):
            large_q.append(fact)
        else:
            small_q.append(fact)

    todofile = params.dirs['DESC'] + 'ecminit.' + str(topargs.seedval)[:10] + '.todo'
    init_dict = dict()
    init_dict['todofilename'] = todofile
    init_dict['u'] = Num
    init_dict['v'] = Den
    init_dict['u_fac'] = factNum
    init_dict['v_fac'] = factDen
    init_dict['mask'] = mask
    init_dict['seed'] = topargs.seedval
    # the below file will exist AFTER descent_large_las
    init_dict['DRELS_FILE'] = params.dirs['DESC'] + 'seed' + str(topargs.seedval)[:10] + '/desc.total.rels'
    # Note that desc.total.rels should contain all Taken lines from
    # 1. All the largeq descents
    # 2. All the rounds of las_descent

    largeq_rels_file = params.dirs['DESC'] + 'ecminit.' + str(topargs.seedval)[:10] + '.largeq.rels'
    init_dict['largeq_rels_file'] = largeq_rels_file

    init_dict_file = params.dirs['DESC'] + 'ecminit.' + str(topargs.seedval)[:10] + '.initdata'
    with open(init_dict_file, "w") as fp:
        json_custom.dump(init_dict, fp, indent=True)

    with open(todofile, "w") as f:
        for sq in small_q:
            if sq > params.BOUNDR:
                f.write(f"{0} {sq}")
                f.write("\n")

    timeprint("Saved initial split and data in: " + str(init_dict_file))
    timeprint("Next large qs need to be handled by descent_large_q.")
