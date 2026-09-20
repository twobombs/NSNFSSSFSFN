from sage.all import *
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from crt_ethroot import *
from timing import *
import tqdm
import mmap

from sage.misc.persist import SagePickler, SageUnpickler
from sage.libs.libecm import ecmfactor

import abc
import argparse
import copy
import functools
import glob
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tempfile
import concurrent.futures
import signal
import heapq
import bisect
import logging
import math

from collections import defaultdict, namedtuple
from pathlib import Path
import time
from contextlib import redirect_stdout,redirect_stderr
from misc_tools import fast_persistent_save, fast_persistent_load

from cado.scripts import descent
from cado_sage import CadoPolyFile
from misc_tools import cat_or_zcat, fast_json_load, fast_json_dump, find_factors_close_to_square_root
from wait_for_file import wait_for_file, wait_for_file_content

from cado_nfs_binaries import CadoNFS, CadoNFSBinaries
from candy import print_command_line,warning_message,major_message
from relations import las_relation, las_relations_from_file
from relations import indexed_relations_from_file
from relations import swap_parts_of_relations
from relations import keep_only_one_relation_per_q, parse_fb_extension_relations
from relations import convert_to_indexed_relation, big_convert_to_indexed_relation

import tocfile
from descent_ecm_utils import transform_polys_by_q, filter_with_ecm
from bwc_helpers import write_ascii_vector
from partial_R import PartialRenumber1024
from concurrent.futures import ProcessPoolExecutor

LinalgOutput = namedtuple('LinalgOutput', ['sol', 'ST_list', 'ST_alg_vector', 'T_list', 'row_to_aquery', 'S_rat_vector', 'indexed_relations_file'])

x = polygen(QQ, 'x')

def silent_remove(filename):
    if os.path.exists(filename):
        timeprint("Deleting",filename)
        os.remove(filename)

def glob_remove(pattern):
    for filename in glob.glob(pattern):
        timeprint("Deleting",filename)
        os.remove(filename)

def make_and_clean(dirname):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    Path(dirname).mkdir(exist_ok=True)

def get_uv_fac(u, v, target_info):
    try:
        u_fac = target_info['u_fac']
        v_fac = target_info['v_fac']

        uv_fac_list = []
        for fac in u_fac:
            uv_fac_list.append( (ZZ(fac),1) )
        for fac in v_fac:
            uv_fac_list.append( (ZZ(fac),-1) )

        uv_fac = Factorization(uv_fac_list, unit=sign(u)*sign(v))
    except KeyError:
        uv_fac = Factorization([(p, sign(e)) for p, e in factor(ZZ(u) / ZZ(v)) for _ in range(abs(e))], unit=sign(u)*sign(v))

    print("u/v", str(u/v))
    print("uv_fac", str(uv_fac.prod()))
    assert uv_fac.prod() == u/v
    return uv_fac

def masked_target(target, e):
    """
    given target as an element of Z/NZ, return the quadruple
    (mask, mask^e*target, u, v) such that mask^e*target = u/v
    """
    ZN = target.parent()
    mask = ZN.random_element()
    h = target * mask**e
    u, v = choose_uv(h)
    return mask, h, u, v

def parse_polyfile(POLYFILE):
    poly = CadoPolyFile(POLYFILE)
    poly.read()
    N = poly.N
    m = poly.m
    f = poly.f[1]
    return (N, m, f)

def check_N_d_and_e(params):
    parameters = params.parameters
    N = parameters['N']
    d = Integer(parameters['d'])
    e = Integer(parameters['e'])

    if d == 0:
        if 'p' in parameters and 'q' in parameters:
            p = Integer(parameters['p'])
            q = Integer(parameters['q'])
            phi = (p-1) * (q-1)
        else:
            timeprint("d=0 in config file,"
                      " factoring N in order to find correct d")
            phi = euler_phi(N)
        assert gcd(e, phi) == 1
        d = ZZ(1/Integers(phi)(e))
        params.parameters['d'] = d
        timeprint(f"Computed d={d}")

    ZN = Integers(N)
    # This is to make sure that we clean up all files if we want to
    # experiment with other exponents.
    assert ZN(2)**(d*e) == 2

    return N, d, e

def generate_or_load_target(params):
    # Should be run ONCE for a computation.
    # I'm not sure what exactly reads/writes to params.target
    if os.path.isfile(params.files['THE_TARGET']):
        timeprint("Using existing target in: " + params.files['THE_TARGET'])
        with open(params.files['THE_TARGET'], "r") as f:
            return Integer(f.readline().strip())

    N = params.poly.N
    ZN = Integers(N)
    e = params.parameters['e']
    nbits = Integer(N).nbits()

    tgt_d = 1
    target = 1

    while True:
        tgt_d = ZN.random_element()
        target = ZN(tgt_d ** e)

        # We want to make sure we get a full-size target
        if Integer(tgt_d).nbits() == nbits and Integer(target).nbits() == nbits:
            break

    with open(params.files['THE_TARGET'], "w") as f:
        f.write(str(target))
        f.write("\n")
        f.write("only for debugging:")
        f.write("\n" + str(tgt_d) + "\n")

    return target

def write_hintfile(infilename,outfilename, BOUNDR, BOUNDA, I):
    # https://sympa.inria.fr/sympa/arc/cado-nfs/2018-07/msg00044.html
    # This file needs to exist and seems to be crucial for the efficiency of the descent
    #sizes = range(12, 24)
    with open(outfilename, "w") as outfile, open(infilename, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            foo = re.match(r"(^.*I=)(\d+)\s+(\d+)(,[\d.,]+)"
                           r"\s+(\d+)(,[\d.,]+)$",
                           line)
            if not foo:
                timeprint("Warning, parse error in input hint file",
                          "on line:\n" + line)
                continue
            prolog,hint_I,lim0,params0,lim1,params1 = foo.groups()

            # If the I or lim values are bad, we should give a warning
            if int(hint_I) > I:
                s = (
                    "Attempting to use hintfile: ",
                    str(infilename),
                    " but encountered a line with I=",
                    str(hint_I),
                    " which is larger than the allowed I=",
                    str(I)
                )
                raise RuntimeError(s)
            if int(lim0) != BOUNDR:
                s = (
                    "Attempting to use hintfile: ",
                    str(infilename),
                    " but encountered a line with lim0=",
                    str(lim0),
                    " which is different than the expected BOUNDR=",
                    str(BOUNDR)
                )
                raise RuntimeError(s)
            if int(lim1) != BOUNDA:
                s = (
                    "Attempting to use hintfile: ",
                    str(infilename),
                    " but encountered a line with lim1=",
                    str(lim1),
                    " which is different than the expected BOUNDA=",
                    str(BOUNDA)
                )
                raise RuntimeError(s)

            outfile.write(prolog + str(hint_I) + " " + str(lim0) + params0 + " " + str(lim1) + params1 + "\n")
    return

def write_hintfile2(params, **kwargs):
    """
    very much WIP at this point. Does not do anything interesting yet. Do
    not use.
    """

    ratio = params.parameters.get('initial_descent_smoothness_ratio', 2.25)
    ratio = float(ratio)
    initial_smoothness_maxbits = params.parameters['MODULUS_BITS'] / 2 / ratio
    initial_smoothness_maxbits = int(initial_smoothness_maxbits)

    I = kwargs.get('I', params.I_sieving)

    # We have prepared factor bases up to this size.
    BOUNDR = params.BOUNDR
    BOUNDA = params.BOUNDA
    LPB0 = params.parameters['LPB0']
    LPB1 = params.parameters['LPB1']

    s = params.poly.skewness

    # We have an affine relation that, for a given value of I and a
    # special-q of bitsize b, gives the approximate bitsize of the two
    # norms.
    multipliers = []
    for side in range(2):
        f = params.poly.f[side]
        d = f.degree()
        f_offset = max([log(abs(f[i])*s**(i-d/2),2) for i in range(d+1)])
        multipliers.append((d, f_offset))

    for b in range(min(LPB0, LPB1), initial_smoothness_maxbits+1):
        raw_norms = [(I+b/2)*m1 + m0 for m1, m0 in multipliers]
        for side in range(2):
            n = [x.round() for x in raw_norms]
            n[side] -= b
            print(f"With I={I}, a {b}@{side} special-q"
                  f" will have norms of {n[0]} and {n[1]} bits")


def polysel_lattice(params,N,e,deg,force_monic=True):
    x = polygen(ZZ)
    ZP = x.parent()
    # Return (root mod N, [coeffs of the algebraic-side polynomial])
    if force_monic:
        timeprint("Generating monic polynomial with f irreducible mod e")
        m = Integer(int(floor(Integer(N)**(Integer(1)/deg))))
        coeffs = Integer(N).digits(m)
        f = sum(ci * x**i for (i,ci) in enumerate(coeffs))
        assert f(m) % N == 0 and f.degree() == deg
        assert f.is_irreducible() # if not, we've already factored N.
        assert f.is_monic()
        # Make sure f(m) has no roots mod e
        done = f.change_ring(GF(e)).is_irreducible()
        tries = 0
        if not done:
            m = m-1
        while not done:
            m += 1
            timeprint("trying a new f,m")
            N_diff = N - m**deg
            B = matrix(deg+1,deg+2)
            w = random_prime(floor(sqrt(N)))
            wi = random_prime(1000)
            for i in range(deg):
                B[i,i] = wi**i
                B[i,deg+1] = w*(m**i)
            B[deg,deg] =  next_prime(m)
            B[deg,deg+1] = -w*N_diff
            B_red = B.LLL()
            for v in B_red:
                if v[-1] == 0:
                    f = sum(Integer(ci/(wi**i)) * x**i for (i,ci) in enumerate(v) if i < deg)
                    if f(m) == 0:
                        timeprint("f(m)!=0")
                        done = False
                        continue
                    f += x**deg
                    print(f)
                    if not f.degree() == deg:
                        done = False
                        continue
                    if not (f(m) % N == 0):
                        timeprint("f(m) % N != 0")
                        done = False
                        continue
                    if not f.is_irreducible():
                        done = False
                        continue
                    if f.change_ring(GF(e)).is_irreducible():
                        done = True
                        break
                    timeprint("f(x) = ",f, "not irreducible mod e")
                    tries += 1
                    if tries > 200:
                        timeprint("Giving up, using f that is not irreducible mod e")
                        done = True
                        break
        timeprint("exited loop with",f,m)
    else:
        m = randint(int(round(N**(1/(deg+1)))), int(round(N**(1/deg))))
        coeffs = []
        curr_deg = deg
        curr_val = N
        while curr_deg >= 0:
            next_coeff = floor(curr_val/(m**curr_deg))
            coeffs.append(next_coeff)
            curr_val -= next_coeff*(m**curr_deg)
            curr_deg -= 1
        coeffs.reverse()
        f = ZP(coeffs)
    return f,m


@timing
def polysel_higheffort_select(params):
    assert 'slurm.numjobs' in params.parameters
    numjobs = int(params.parameters['slurm.numjobs'])
    admin_min = float(params.parameters['poly.admin'])
    admax_max = float(params.parameters['poly.admax'])
    ad_diff = floor((admax_max - admin_min)/numjobs)

    params.save_to_file()

    processes = []
    for jobnum in range(numjobs):
        admin = admin_min + jobnum*ad_diff
        admax = admin + ad_diff
        outfile = params.files['POLYSEL_SELECT_PFX'] + "." + str(jobnum)
        command_list = [
            params.files['SAGE'],
            "polyselect_helper.py",
            "--params", params.files['PARAMS'],
            "--admin", str(admin),
            "--admax", str(admax),
            "--outfile", outfile,
        ]
        processes.append(slurmit(params, " ".join(command_list), params.prefix[:-1]+"-polyselect", jobnum))

    finished_processes, cputime_slurm = slurm_wait(processes)
    overall_cputime.add(cputime_slurm)
    print(f"Finished running {numjobs} polyselects!")


@timing
def polysel_higheffort_ropt(params):
    assert 'slurm.numjobs' in params.parameters
    numjobs = int(params.parameters['slurm.numjobs'])

    params.save_to_file()

    processes = []
    for jobnum in range(numjobs):
        infile = params.files['POLYSEL_SELECT_PFX'] + "." + str(jobnum)
        outfile = params.files['POLYSEL_SELECT_PFX'] + "." + str(jobnum) + ".ropt"
        command_list = [
            params.files['SAGE'],
            "poly_ropt_helper.py",
            "--params", params.files['PARAMS'],
            "--infile", infile,
            "--outfile", outfile,
        ]
        processes.append(slurmit(params, " ".join(command_list), params.prefix[:-1]+"-polyropt", jobnum))

    finished_processes, cputime_slurm = slurm_wait(processes)
    overall_cputime.add(cputime_slurm)
    print(f"Finished running {numjobs} poly_ropts!")


@timing
def polysel_higheffort_candidates(params, ncandidates=10):
    assert 'slurm.numjobs' in params.parameters
    numjobs = int(params.parameters['slurm.numjobs'])
    output = ""

    deg = int(params.parameters['POLY_DEG'])
    N = Integer(params.parameters['N'])

    best = [(None,None,None,None) for i in range(ncandidates)]
    current = None
    x = polygen(ZZ, 'x')
    ZN = Integers(N)

    for jobnum in range(numjobs):
        outfilename = params.files['POLYSEL_SELECT_PFX'] + "." + str(jobnum) + ".ropt"
        with open(outfilename, "r") as outfile:
            output += outfile.read()
            output += "\n"

    for line in output.split("\n"):
        if m := re.match(r"^#* root-optimized polynomial (\d+) #*", line):
            current = {}
            current['index'] = int(m.group(1))
        elif m := re.match(r"^(\w+): (-?\d+(\.\d*)?)", line):
            current[m.group(1)] = m.group(2)
        elif m := re.match(r"^# side 1 MurphyE\(.*\)=([\d\.e-]+)", line):
            score = float(m.group(1))
            f = sum([int(current[f'c{i}'])*x**i for i in range(deg+1)])
            g = sum([int(current[f'Y{i}'])*x**i for i in range(2)])
            m = Integer(-ZN(g[0])/ZN(g[1]))
            best_to_replace = -1
            for i in range(len(best)):
                if best[i][0] is None or best[i][0] < score:
                    best_to_replace = i
                    break
            if best_to_replace < 0:
                continue
            else:
                best[best_to_replace] = (score, g, f, m)

    for i in range(len(best)):
        score = best[i][0]
        g = best[i][1]
        f = best[i][2]
        m = best[i][3]

        if score is None:
            continue

        filename = params.files['POLYSEL_SELECT_PFX'] + ".poly.option." + str(i)

        with open(filename, "w") as file:
            print(f"n: {N}", file=file)
            print(f"m: {m}", file=file)
            print("poly0: ", ", ".join([str(x) for x in g.list()]), file=file)
            print("poly1: ", ", ".join([str(x) for x in f.list()]), file=file)
        with open(filename + ".only-side1", "w") as file:
            print(f"n: {N}", file=file)
            print(f"m: {m}", file=file)
            print("poly0: ", ", ".join([str(x) for x in f.list()]), file=file)

        try:
            skew = float(CadoNFS("polyselect/skewness", "poly",
                         inputs={"poly": filename},
                         capture=True))
        except CadoNFSBinaries.BinaryNotFound:
            skew = 0.7
            warning_message(f"binary polyselect/skewness not found, using phony skewness value {skew} instead")

        with open(filename, "a+") as file:
            print(f"skew: {skew}", file=file)

        with open(filename + ".Me", "w") as efile:
            print(f"\nscore: {score}\n", file=efile)

    print("Finished selecting best " + str(ncandidates) + " candidate polynomials!")


@timing
def polysel_cado(params,N,e,deg):
    output = CadoNFS("polyselect/polyselect",
                     "-P",params.parameters['poly.P'],
                     "-N",N,
                     "-degree",deg,
                     "-t", params.polyselect_nthreads_or_auto,
                     "-admax",params.parameters['poly.admax'],
                     "-incr",params.parameters['poly.incr'],
                     "-nq",625,
                     capture=True)
    surviving_polynomials = []
    coeffs = {}
    polylines = []

    lines = iter(output.decode('utf-8').split("\n"))
    candidate_polys = []
    for line in lines:
        if line.startswith("#"):
            continue
        if line.startswith("n:"):
            polylines = [line]
            while not line.startswith("#"):
                line = next(lines)
                polylines.append(line)
            matches = re.search(r"lognorm ([\d\.]+),", line)
            lognorm = matches.groups(0)[0]
            heapq.heappush(candidate_polys, (float(lognorm),polylines))

    params.save_to_file()

    if not params.slurm:
        timeprint("Running polyselect_ropt without slurm, setting numjobs = 1")
        numjobs = 1
    else:
        if 'slurm.numjobs' in params.parameters:
            numjobs = int(params.parameters['slurm.numjobs'])
        else:
            numjobs = 1
        timeprint("Running polyselect_ropt with",numjobs,"jobs")
    jobsize = min(floor(len(candidate_polys)/numjobs),1000)
    filenames = []
    processes = []
    for jobnum in range(numjobs):
        outfilename = params.files['POLYSELECT']+"."+str(jobnum)
        filenames.append(outfilename)
        with open(outfilename,"w") as outfile:
            #outfile.write(output.decode('utf-8'))
            for i in range(jobsize):
                if not candidate_polys:
                    break
                lognorm, polylines = heapq.heappop(candidate_polys)
                for line in polylines:
                    outfile.write(line)
                    outfile.write("\n")
                outfile.write("\n")

        command_list = [
            params.files['SAGE'],
            "poly_ropt_helper.py",
            "--params", params.files['PARAMS'],
            "--infile", outfilename,
            "--outfile", outfilename+".out"
        ]

        if params.slurm:
            processes.append(slurmit(params,
                " ".join(command_list),
                params.prefix[:-1]+"-polysel",
                jobnum))
        else:
            processes.append(subprocess.Popen(["time", "-p"] + command_list, stderr=subprocess.PIPE, text=True))

    expected_filenames = [f"{filename}.out" for filename in filenames]
    if params.slurm:
        finished_processes, cputime_slurm = slurm_wait(processes)
        overall_cputime.add(cputime_slurm)
        if not wait_for_file(expected_filenames):
            logging.error("Time out waiting for files: ", ', '.join(expected_filenames))
            exit(1)
    else:
        for process in processes:
            process.wait()
            stderr = process.communicate()[1]
            overall_cputime.add(extract_time(stderr))

    # output = CadoNFS("polyselect/polyselect_ropt",
    #                  "-t",params.polyselect_nthreads_or_auto,
    #                  "-inputpolys",params.files['POLYSELECT'],
    #                  "-area", params.parameters['poly.area'],
    #                  "-Bf",params.parameters['poly.Bf'],
    #                  "-Bg",params.parameters['poly.Bg'],
    #                  "-ropteffort",params.parameters.get('poly.ropteffort',5),
    #                  capture=True)

    best = (None,)
    current = None
    x = polygen(ZZ, 'x')
    ZN = Integers(N)
    output = ""
    for outfilename in expected_filenames:
        with open(outfilename,"r") as outfile:
            output += outfile.read()
        output += "\n"

    for line in output.split("\n"):
        if m := re.match(r"^#* root-optimized polynomial (\d+) #*", line):
            current = {}
            current['index'] = int(m.group(1))
        elif m := re.match(r"^(\w+): (-?\d+(\.\d*)?)", line):
            current[m.group(1)] = m.group(2)
        elif m := re.match(r"^# side 1 MurphyE\(.*\)=([\d\.e-]+)", line):
            score = float(m.group(1))
            if best[0] is not None and score < best[0]:
                continue
            f = sum([int(current[f'c{i}'])*x**i for i in range(deg+1)])
            if not f.change_ring(GF(e)).is_irreducible():
                continue
            g = sum([int(current[f'Y{i}'])*x**i for i in range(2)])
            m = Integer(-ZN(g[0])/ZN(g[1]))
            best = (score, g, f, m)
    timeprint("Best poly has score:",score)
    return best[1:]

@timing
def write_polyfile(params,algorithm,N, e, deg, filename, force_monic=True, hardcoded_m=None, hardcoded_polys=None):
    x = polygen(ZZ, 'x')
    if algorithm == 'custom':
        f,m = polysel_lattice(params,N,e,deg,force_monic=True)
        g = x - m
    if algorithm == 'cado':
        g,f,m = polysel_cado(params,N,e,deg)
    if algorithm == 'snfs538':
        m = 12**30
        g = x - m
        f = 12*x**5 - 1
    if algorithm == 'hardcoded':
        assert hardcoded_m is not None
        assert hardcoded_polys is not None
        m = hardcoded_m
        g = hardcoded_polys[0]
        f = hardcoded_polys[1]

    with open(filename, "w") as file:

        print(f"n: {N}", file=file)
        print(f"m: {m}", file=file)
        # It's a bit of a pity: we could at least theoretically write the
        # polynomials in polynomial form, but cado's descent.py won't
        # understand it.
        #print(f"poly0: -{m}, 1", file=file)
        print("poly0: ", ", ".join([str(x) for x in g.list()]), file=file)
        print("poly1: ", ", ".join([str(x) for x in f.list()]), file=file)

    with open(filename + ".only-side1", "w") as file:
        #x = polygen(ZZ, 'x')
        print(f"n: {N}", file=file)
        print(f"m: {m}", file=file)
        print("poly0: ", ", ".join([str(x) for x in f.list()]), file=file)

    # We have an external tool that can compute the skewness
    try:
        skew = float(CadoNFS("polyselect/skewness", "poly",
                             inputs={"poly": filename},
                             capture=True))

        if skew < 0.0000000001:
            skew = 1.0  # fake, but 2e-155 won't be parsed correctly

    except CadoNFSBinaries.BinaryNotFound:
        skew = 0.7
        warning_message(f"binary polyselect/skewness not found, using phony skewness value {skew} instead")

    with open(filename, "a+") as file:
        print(f"skew: {skew}", file=file)

@timing
def run_oracle(params, todofilename, jsonfilename):
    if params.oracle == "sage":
        N, d, e = check_N_d_and_e(params)

        command_line = [params.files['SAGE'],
            f"oracles/sage_oracle.py",
            "-d", str(d),
            "-N", str(N),
            "--in", todofilename,
            "--out", jsonfilename
        ]

        if params.mpi:
            command_line = [
                params.mpi_mpirun_bin,
                *(params.mpi_extra_args.split(" ")),
                "-n", "1",
                *command_line,
                "--mpi",
                "--thr", str(params.parameters['mpi.thr'])
            ]
    elif params.oracle == "luna_k6":
        bits = params.parameters['MODULUS_BITS']
        N = params.parameters['N']
        e = params.parameters['e']
        command_line = [params.files['PYTHON'],
            f"oracles/luna_k6_hsm_oracle_client.py",
            "--bits", str(bits),
            "--label", f"id_oracle_rsa_{bits}_exp_{e}",
            "-N", str(N),
            "-e", str(e),
            "--logs_dir", params.dirs['TEMP_OUTPUT_DIR'].rstrip("/"),
            "--in", todofilename,
            "--out", jsonfilename
        ]
    elif params.oracle == "remote_luna_S750":
        bits = params.parameters['MODULUS_BITS']
        N = params.parameters['N']
        e = params.parameters['e']
        command_line = [params.files['PYTHON'],
            "oracles/remote_queries/local_query_handler.py",
            "--oracle", "luna_S750",
            "--label", f"id_oracle_rsa_{bits}_exp_{e}",
            "-N", str(N),
            "--work-dir", params.dirs['TEMP_OUTPUT_DIR'].rstrip("/"),
            "--in", todofilename,
            "--out", jsonfilename
        ]
    else:
        raise Exception(f"Unsupported oracle {params.oracle}")

    print_command_line(*command_line)
    p = subprocess.run(["time", "-p"] + command_line, check=True, stderr=subprocess.PIPE, text=True)
    overall_cputime.add(extract_time(p.stderr))


@timing
def do_rational_queries(params):
    boundr = params.BOUNDR
    todofilename = params.files['RQUERIES_TODO']
    jsonfilename = params.files['RQUERIES_FILE']

    nqueries = 0
    do_generate = True
    if os.path.isfile(todofilename):
        with open(todofilename, "r") as todofile:
            do_generate = len(todofile.readlines()) < prime_pi(boundr)

    if do_generate:
        with open(todofilename, "w") as todofile:
            primesieve_bin = "/usr/bin/primesieve"
            if os.path.exists(primesieve_bin):
                command_line = [
                    "/usr/bin/primesieve",
                    str(boundr),
                    "--print"
                ]
                print_command_line(*command_line)
                p = subprocess.run(
                    ["time", "-p"] + command_line,
                    check=True,
                    stdout=todofile,
                    stderr=subprocess.PIPE,
                    text=True
                )
                overall_cputime.add(extract_time(p.stderr))
            else:
                for p in prime_range(boundr):
                    print(p, file=todofile)
                    nqueries += 1

    with open(todofilename, "r+") as todofile:
        with open(todofilename, "r") as todofile:
            queries = todofile.readlines()
            nqueries = len(queries)
        last_primes = frozenset([int(p.strip()) for p in queries])

        # I'm not sure it's useful on the algebraic side. Quite possibly
        # not. It doesn't hurt, though.
        for side in range(2):
            F = params.poly.f[side]
            for p, e in F.leading_coefficient().factor(limit=boundr):
                if p < boundr or p in last_primes:
                    continue
                timeprint(f"including factor {p} of the leading coefficient of"
                          f" f{side} = {F} among the rational queries")
                print(p, file=todofile)
                nqueries += 1

    if not do_generate and os.path.isfile(jsonfilename):
        timeprint(f"Skipping calls to oracle for rational queries and instead reusing the ones already in {jsonfilename}.")
    else:
        run_oracle(params, todofilename, jsonfilename)

    return nqueries

def slurmit(params, command, jobname, jobnum, exclusive=True, setup_commands=[], partition=None, nested_mpi=False):
    jobfile = f"{params.dirs['SLURMJOBS']}{jobname}.{jobnum}.job"
    outfile = f"{params.dirs['SLURMJOBS']}{jobname}.{jobnum}.%j.out"
    errfile = f"{params.dirs['SLURMJOBS']}{jobname}.{jobnum}.%j.err"

    if nested_mpi:
        m = re.search(r'--mpi[\s=]+(\d+)x(\d+)', command)
        if m:
            groups = m.groups()
            nnodes = int(groups[0]) * int(groups[1])
        else:
            nnodes = 1

    with open(jobfile,"w") as f:
        f.writelines("#!/bin/bash\n")
        if nested_mpi:
            f.writelines(f"#SBATCH -N {nnodes}\n")
        else:
            f.writelines(f"#SBATCH -n 1\n")
        f.writelines(f"#SBATCH --job-name {jobname}\n")
        f.writelines("#SBATCH --requeue\n")
        if len(params.slurm_exclude) > 0:
            f.writelines(f"#SBATCH --exclude={params.slurm_exclude}\n")
        if exclusive:
            f.writelines("#SBATCH --exclusive\n")
        if not partition:
            partition = params.slurm_job_partition
        f.writelines(f"#SBATCH --partition={partition}\n")
        f.writelines(f"#SBATCH --output={outfile}\n")
        f.writelines(f"#SBATCH -e {errfile}\n")
        f.writelines("\n")
        f.writelines("export DOT_SAGE=/tmp/$USER.sage/\n")
        f.writelines("set -e\n")
        f.writelines(setup_commands)

        if nested_mpi:
            # Remark: Possible ways to call MPI are:
            #   - "sbatch/script/mpirun" (this file),
            #   - "sbatch/mpirun script" (OpenMPI recommendation but can't have further mpirun calles inside script)
            #   - srun --mpi=pmix script (can't have further mpirun calles inside script)
            f.writelines("""
# srun cannot be used here since the script below has internal calls to MPI.
# Using srun loses the node/tasks information from sbatch and thus tries to
# run the below script from all tasks at the same time, causing utter mayhem.
""")
            f.writelines(f"time -p " + command)
        else:
            f.writelines("\n# srun is necessary so that SIGINT signals are forwarded to the executed job, time is here to measure cputime.\n")
            f.writelines(f"srun time -p " + command)

        f.writelines("\n")

    result = subprocess.run(["sbatch", jobfile],stdout=subprocess.PIPE)
    timeprint("Submitted slurm job:")
    major_message(command)
    pattern = rb'Submitted batch job (\d+)'
    match = re.search(pattern,result.stdout)
    if match:
        jobnum = int(match.group(1))
        timeprint("job number is",jobnum)
    return jobnum, errfile.replace("%j", str(jobnum))

def slurm_wait(process_list, throw_error_on_failed_job=False, params=None):
    finished_processes = []
    launched_processes = []
    cputime_slurm = 0
    while process_list:
        for jobnum, errfile in process_list:
            result = subprocess.run(["sacct","-j",str(jobnum),"--format=state","--noheader"],stdout=subprocess.PIPE)

            # Remark: the order of the below matters. Some batch jobs have multiple
            # subjobs that are listed in the sacct output.

            # If any (sub)jobs failed we want to fail the whole job.
            if re.search(b'FAILED',result.stdout):
                result2 = subprocess.run(["sacct","-j",str(jobnum),"--format=nodelist,JobName%50","--noheader","-X"],stdout=subprocess.PIPE)
                nodename, jobname = [entry.decode('ascii') for entry in result2.stdout.split()]
                exception_msg = f"job {jobname} ({jobnum}) failed on {nodename} :("
                process_list.remove((jobnum, errfile))
                cputime_slurm += extract_time_from_file(errfile)
                timeprint(exception_msg)

                if throw_error_on_failed_job:
                    if os.path.isfile(errfile):
                        with open(errfile, "r") as fp:
                            exception_msg += f"\nHere is the job's error output from file {errfile}:\n" + fp.read()
                    raise Exception(exception_msg)

            # If any (sub)jobs is still running, we want to wait.
            elif re.search(b'RUNNING',result.stdout):
                if jobnum in launched_processes:
                    continue
                result2 = subprocess.run(["sacct","-j",str(jobnum),"--format=nodelist","--noheader"],stdout=subprocess.PIPE)
                nodename = result2.stdout.split()[0].decode('ascii')
                launched_processes.append(jobnum)
                timeprint(f"job {jobnum} launched on node {nodename}")

            # If all (sub)jobs are either complete or cancelled, we want to mark the whole job as cancelled.
            elif re.search(b'CANCELLED',result.stdout):
                result2 = subprocess.run(["sacct","-j",str(jobnum),"--format=elapsed","--noheader"],stdout=subprocess.PIPE)
                elapsed = result2.stdout.split()[0].decode('ascii')
                finished_processes.append(jobnum)
                process_list.remove((jobnum, errfile))
                job_cputime_slurm = extract_time_from_file(errfile)
                cputime_slurm += job_cputime_slurm
                timeprint(f"job {jobnum} cancelled :/ elapsed: {elapsed} (cputime: {job_cputime_slurm:2.4f}s)")

            # If all (sub)jobs completed, we finally are done!
            elif re.search(b'COMPLETED',result.stdout):
                result2 = subprocess.run(["sacct","-j",str(jobnum),"--format=elapsed","--noheader"],stdout=subprocess.PIPE)
                elapsed = result2.stdout.split()[0].decode('ascii')
                finished_processes.append(jobnum)
                process_list.remove((jobnum, errfile))
                job_cputime_slurm = extract_time_from_file(errfile)
                cputime_slurm += job_cputime_slurm
                timeprint(f"job {jobnum} succeeded :) elapsed: {elapsed}s (cputime: {job_cputime_slurm:2.4f}s)")
    return finished_processes, cputime_slurm

def get_slurm_job_status(slurm_job):
    result = subprocess.run(["sacct","-j",str(slurm_job),"--format=state","--noheader"],stdout=subprocess.PIPE)
    if re.search(b'RUNNING',result.stdout):
        return 'RUNNING'
    if re.search(b'COMPLETED',result.stdout):
        return 'COMPLETED'
    if re.search(b'FAILED',result.stdout):
        return 'FAILED'
    if re.search(b'PENDING',result.stdout):
        return 'PENDING'
    return None

@timing
def call_todo_sieving(params, is_alg_or_ext=True):
    # True for alg, False for extension

    R = params.R

    if is_alg_or_ext:
        shortname = "alg"
        basefile = params.files['AQRELS_FILE']
        base_outfile = params.dirs['TEMP_OUTPUT_DIR'] + "algrels/alg"
        outstanding_renumber_indices = set(
                [R.column_to_renumber(i)
                 for i in range(R.number_of_algebraic_columns())
                 if i == 0
                 or R.side_and_index_to_ideal(1, i)[1] < params.BOUNDA_queries
                 ])
        A_start = int(params.parameters.get('algebraic_todo_sieving.A_min', params.parameters['A_sieving']))
        A_max = int(params.parameters.get('algebraic_todo_sieving.A_max', A_start+1))
    else:
        shortname = "ext"
        basefile = params.files['EXTRELS_FILE']
        base_outfile = basefile
        outstanding_renumber_indices = set(
                [R.column_to_renumber(i)
                 for i in range(R.number_of_algebraic_columns())
                 if i > 0
                 and R.side_and_index_to_ideal(1, i)[1] >= params.BOUNDA_queries
                 ])
        A_start = int(params.parameters.get('extension_todo_sieving.A_min', params.parameters['A_sieving']))
        A_max = int(params.parameters.get('extension_todo_sieving.A_max', A_start+1))

    num_total_ideals = len(outstanding_renumber_indices)

    for r in convert_to_indexed_relation(las_relations_from_file(basefile), params):
        for ii in r.indices:
            outstanding_renumber_indices -= {ii}

    extra_flags = []
    round_info = dict(num=0,
                      A = A_start,
                      mfb1 = int(params.parameters['sieve.mfb1']))

    #if not params.slurm:
    #    raise RuntimeError("Haven't implemented non-slurm todo sieving.")
    #    exit(0)

    if 'extension_sieving.numjobs' in params.parameters:
        numjobs = int(params.parameters['extension_sieving.numjobs'])
    elif 'sieve.numjobs' in params.parameters:
        numjobs = int(params.parameters['sieve.numjobs'])
    elif 'slurm.numjobs' in params.parameters:
        numjobs = int(params.parameters['slurm.numjobs'])
    else:
        timeprint("slurm.numjobs not set; defaulting to numjobs = 24")
        numjobs = 24

    if not params.slurm:
        numjobs = 1

    while len(outstanding_renumber_indices) and round_info['A'] <= A_max:
        outstanding_qs = [R.side_and_index_to_ideal(1,
                                                    R.renumber_to_column(i))
                      for i in outstanding_renumber_indices]
        timeprint(f"Missing relations for {len(outstanding_qs)} ideals")
        if len(outstanding_qs) < 20:
            timeprint("Complete list:", outstanding_qs)

        pct_missing = round(100.0 * len(outstanding_qs) / num_total_ideals, 10)
        timeprint(f"That means the percentage of missing {shortname} ideals is: {pct_missing}")

        sieve_processes = []
        nonlinear_todo_qs = []
        qs_per_job = ceil(len(outstanding_qs)/numjobs)
        last_cmd = []

        for job in range(numjobs):
            outfile_suffix = f".round{round_info['num']}.job{job}"
            todofilename = base_outfile + outfile_suffix + ".todo"
            start = job*qs_per_job
            stop = min((job+1)*qs_per_job, len(outstanding_qs))

            with open(todofilename, "w") as f:
                for side, q, rho in outstanding_qs[start:stop]:
                    assert side == 1
                    if "alpha" in str(rho) and (q not in nonlinear_todo_qs):
                        # This means it's a nonlinear ideal, so rho is not an integer,
                        # but the todofile only accepts integers. Ugh.
                        nonlinear_todo_qs.append(q)
                    else:
                        print(f"0 {q} {rho}", file=f)

            command_list = [params.files['SAGE'],
                            "todo_sieving_helper.py",
                            "--params",params.files['PARAMS'],
                            "--roundA", str(round_info['A']),
                            "--roundMfb", str(round_info['mfb1']),
                            "--outfile", base_outfile + outfile_suffix,
                            "--todofile", todofilename,
                            "--q0",str(-1),
                            "--q1",str(-1)]

            jobname = f"{params.prefix[:-1]}-todosieve-{round_info['num']}-{job}"
            if params.slurm:
                sieve_processes.append(
                    slurmit(
                        params,
                        " ".join(command_list),
                        jobname, job))
            last_cmd = command_list

        todo_start = time.time()
        if params.slurm:
            finished_processes, cputime_slurm = slurm_wait(sieve_processes)
            overall_cputime.add(cputime_slurm)
        else:
            assert numjobs == 1
            sp = subprocess.run(last_cmd, stdout=subprocess.PIPE)
            todo_fin = time.time()
            overall_cputime.add(int(round(todo_fin-todo_start)))

        timeprint(f"This round of todo jobs all finished!")

        # TODO: Actually use extra_flags. Add to timing dictionary.

        outfile_suffix = f".round{round_info['num']}"
        with open(base_outfile + outfile_suffix, "w") as round_out:
            for job in range(numjobs):
                jobfile = f"{base_outfile}.round{round_info['num']}.job{job}"
                if not wait_for_file([jobfile]):
                    continue

                try:
                    infile = open(jobfile, "r")
                    for line in infile.readlines():
                        _ = round_out.write(line)
                    infile.close()
                except FileNotFoundError:
                    timeprint(f"{jobfile} does not exist, skipping")

        if not is_alg_or_ext:
            keep_only_one_relation_per_q(base_outfile + outfile_suffix,
                                        (1, params.BOUNDA_queries, params.BOUNDA),
                                        keep_file=True)

        round_info['num'] += 1
        for r in convert_to_indexed_relation(las_relations_from_file(base_outfile + outfile_suffix), params):
            for ii in r.indices:
                outstanding_renumber_indices -= {ii}

        if len(outstanding_renumber_indices):
            timeprint("Increasing A, mfb1 and trying todo sieving again!")
            round_info['A'] += 1
            round_info['mfb1'] += 2
            extra_flags = [ "--adjust-strategy", 2 ]

    if len(outstanding_renumber_indices) == 0:
        timeprint("No outstanding ideals! All done with todo sieving.")
    else:
        timeprint("Reached the maximum A value of " + str(round_info['A']))
        timeprint("Giving up on todo sieving. The number of missing ideals is " + str(len(outstanding_renumber_indices)))

    with open(basefile, "a") as outfile:
        for r in range(round_info['num']):
            timeprint(f"Adding relations from round {r} todo sieving")
            with open(f"{base_outfile}.round{r}", "r") as infile:
                for line in infile.readlines():
                    if not line.startswith("#"):
                        _ = outfile.write(line)

    if is_alg_or_ext:
        ulfile = params.files['UNLINKED_IDEALS']
    else:
        ulfile = params.files['EXT_UNLINKED_IDEALS']

    # This would be the first encounter with unlinked ideals
    with open(ulfile, "w") as f:
        for ii in outstanding_renumber_indices:
            f.write(str(ii) + "\n")

    num_remaining = len(outstanding_renumber_indices)
    pct_missing = round(100.0 * num_remaining / num_total_ideals, 10)
    timeprint(f"After todo sieving, the percentage of missing {shortname} ideals is: {pct_missing}")

@timing
def call_extra_sieving_granular(params, prefix, given_A, given_mfb1, is_alg_or_ext=True, is_filter=False,
                                mult_by_nprimes=0):
    # True for alg, False for extension
    # mult_by_nprimes: For a given q, add todos for q*q' for q' the first n primes.
    # mult_by_nprimes can be helpful when A becomes too large. It will require --allow_compsq.
    # If used for extension sieving, the de-duplication by q certainly should be done.

    #R = params.R
    extra_flags = []

    if is_alg_or_ext:
        basefile = params.files['AQRELS_FILE']
        base_outfile = params.dirs['TEMP_OUTPUT_DIR'] + "algrels/alg"
        if is_filter:
            unlinkedfile = params.files['FILT_UNLINKED_IDEALS']
            unlinked_todofile = params.files['FILT_UNLINKED_TODOS']
        else:
            unlinkedfile = params.files['UNLINKED_IDEALS']
            unlinked_todofile = params.files['UNLINKED_IDEALS_TODOS']
    else:
        basefile = params.files['EXTRELS_FILE']
        base_outfile = basefile
        unlinkedfile = params.files['EXT_UNLINKED_IDEALS']
        unlinked_todofile = params.files['EXT_UNLINKED_TODOS']

    if not params.slurm:
        raise RuntimeError("Haven't implemented non-slurm todo sieving.")
        exit(0)

    if 'extension_sieving.numjobs' in params.parameters and (not is_alg_or_ext):
        numjobs = int(params.parameters['extension_sieving.numjobs'])
    elif 'sieve.numjobs' in params.parameters:
        numjobs = int(params.parameters['sieve.numjobs'])
    elif 'slurm.numjobs' in params.parameters:
        numjobs = int(params.parameters['slurm.numjobs'])
    else:
        timeprint("slurm.numjobs not set; defaulting to numjobs = 72")
        numjobs = 72

    #if is_filter:
        # more reasonable
        # should be quite a small list
    #    numjobs = 4

    #outstanding_qs = []
    #with open(unlinkedfile, "r") as f:
    #    for line in f.readlines():
    #        renumber_index = int(line.strip())
    #        ideal = R.side_and_index_to_ideal(1, R.renumber_to_column(renumber_index))
    #        outstanding_qs.append(ideal)

    outstanding_todo_lines = []
    with open(unlinked_todofile, 'r') as f:
        for line in f:
            outstanding_todo_lines.append(line)

    timeprint(f"Missing relations for {len(outstanding_todo_lines)} ideals")

    sieve_processes = []
    nonlinear_todo_qs = []
    qs_per_job = ceil(len(outstanding_todo_lines)/numjobs)
    #if mult_by_nprimes > 0:
    #    primeset = primes_first_n(mult_by_nprimes)  #[4:]  # skip first few primes
    #else:
    #    primeset = []

    primeset = primes_first_n(100)[60:]

    og_poly = CadoPolyFile(params.files['POLYFILE']); og_poly.read()
    alg_poly = og_poly.f[1]

    for job in range(numjobs):
        outfile = f"{base_outfile}.{prefix}.job{job}"
        todofilename = outfile + ".todo"
        start = job*qs_per_job
        stop = min((job+1)*qs_per_job, len(outstanding_todo_lines))

        with open(todofilename, "w") as f:
            for line in outstanding_todo_lines[start:stop]:
                #f.write(line)

                lineinfo = line.strip().split()
                assert lineinfo[0] == '0'
                q = int(lineinfo[1])
                rho = int(lineinfo[2])

                for q_prime in primeset:
                    roots = alg_poly.roots(GF(q_prime))
                    if len(roots) > 0:
                        rho_prime = roots[0][0]
                        try:
                            rho_q_qprime = crt(Integer(rho), Integer(rho_prime), Integer(q), Integer(q_prime))
                            assert rho_q_qprime % q == rho
                            assert rho_q_qprime % q_prime == rho_prime
                            f.write(f"0 {q*q_prime} {rho_q_qprime}\n")
                        except ValueError:
                            continue

            #for side, q, rho in outstanding_qs[start:stop]:
            #    assert side == 1
            #    if "alpha" in str(rho) and (q not in nonlinear_todo_qs):
                    # This means it's a nonlinear ideal, so rho is not an integer,
                    # but the todofile only accepts integers. Ugh.
            #        nonlinear_todo_qs.append(q)
            #    else:
                    # side 0 since it's 1-sided sieving
            #        print(f"0 {q} {rho}", file=f)

            #        if mult_by_nprimes > 0:
            #            for q_prime in primeset:
            #                roots = alg_poly.roots(GF(q_prime))
            #                if len(roots) > 0:
            #                    rho_prime = roots[0][0]
            #                    rho_q_qprime = crt(Integer(rho), Integer(rho_prime), Integer(q), Integer(q_prime))
            #                    assert rho_q_qprime % q == rho
            #                    assert rho_q_qprime % q_prime == rho_prime
            #                    print(f"0 {q*q_prime} {rho_q_qprime}", file=f)

        #if len(nonlinear_todo_qs) > 0:
        #    major_message("Warning: some nonlinear todo qs.")
        #    if len(nonlinear_todo_qs) < 25:
        #        print("Here they are:")
        #        print(str(nonlinear_todo_qs))

        command_list = [params.files['SAGE'],
                        "todo_sieving_helper.py",
                        "--params",params.files['PARAMS'],
                        "--roundA", str(given_A),
                        "--roundMfb", str(given_mfb1),
                        "--outfile", outfile,
                        "--todofile", todofilename,
                        "--q0",str(-1),
                        "--q1",str(-1)]

        jobname = f"{params.prefix[:-1]}-extrasieve-{prefix}-{job}"
        sieve_processes.append(
            slurmit(
                params,
                " ".join(command_list),
                jobname, job))

    todo_start = time.time()
    finished_processes, cputime_slurm = slurm_wait(sieve_processes)
    overall_cputime.add(cputime_slurm)
    todo_fin = time.time()
    timeprint("Done doing extra sieving!")

    # TODO: Actually use extra_flags. Add to timing dictionary.
    # Also handle the nonlinear ideals, using q0 and q1.

    #with open(basefile, "a") as ff:
    #    for job in range(numjobs):
    #        jobfile = f"{base_outfile}.{prefix}.job{job}"
    #        with open(jobfile, "r") as infile:
    #            for line in infile.readlines():
    #                if not line.startswith("#"):
    #                    _ = ff.write(line)

    #timeprint("Done adding extra relations to " + basefile + "!")


@timing
def call_algebraic_query_sieving(params, is_extra=False):
    params.save_to_file()

    c0 = QQ(params.parameters.get("algebraic_query_sieving.q0_ratio", 1/4))
    c1 = 1 #4
    q0 = floor(c0*params.BOUNDA_queries)
    q1 = c1*params.BOUNDA_queries

    if is_extra and 'extra_alg.q0' in params.parameters and 'extra_alg.q1' in params.parameters:
        q0 = int( params.parameters['extra_alg.q0'] )
        q1 = int( params.parameters['extra_alg.q1'] )
        assert q1 > q0

    #R = params.R

    sieve_processes = []

    if not params.slurm:
        timeprint("Running sieving without slurm, setting numjobs = 1")
        numjobs = 1
        if is_extra:
            jobstr = params.extra_prefix + "0"
        else:
            jobstr = "0"
        command_list = ["time", "-p",
                        params.files['SAGE'],
                        "algebraic_query_sieving_helper.py",
                        "--params",params.files['PARAMS'],
                        "--jobnum",jobstr,
                        "--q0",str(q0),
                        "--q1",str(q1)]
        process = subprocess.run(command_list, stderr=subprocess.PIPE, text=True)
        overall_cputime.add(extract_time(process.stderr))
        jobrange = range(numjobs)

    else:
        if 'sieve.numjobs' in params.parameters:
            numjobs = int(params.parameters['sieve.numjobs'])
        elif 'slurm.numjobs' in params.parameters:
            numjobs = int(params.parameters['slurm.numjobs'])
        else:
            timeprint("slurm.numjobs not set; defaulting to numjobs = 24")
            numjobs = 24
        step = floor((q1-q0)/numjobs)
        intervals = [q0+i*step for i in range(numjobs)]+[q1]
        timeprint(intervals)

        if 'alg_sieving.jobnum.min' in params.parameters and 'alg_sieving.jobnum.max' in params.parameters:
            jobrange = range(
                int(params.parameters['alg_sieving.jobnum.min']),
                int(params.parameters['alg_sieving.jobnum.max']),
            )
        else:
            jobrange = range(numjobs)

        for jobnum in jobrange:

            if is_extra:
                jobstr = params.extra_prefix + str(jobnum)
            else:
                jobstr = str(jobnum)

            if os.path.exists(f'alg.{jobstr}.tmp'):
                continue

            command_list = [params.files['SAGE'],
                            "algebraic_query_sieving_helper.py",
                            "--params",params.files['PARAMS'],
                            "--jobnum",jobstr,
                            "--q0",str(intervals[jobnum]),
                            "--q1",str(intervals[jobnum+1])]
            if is_extra:
                jobname = f"{params.prefix[:-1]}-{params.extra_prefix}-algsieve-{q0}-{q1}"
            else:
                jobname = f"{params.prefix[:-1]}-algsieve-{q0}-{q1}"
            sieve_processes.append(
                    slurmit(
                        params,
                        " ".join(command_list),
                        jobname, jobnum))

        aqrels_start = time.time()
        finished_processes, cputime_slurm = slurm_wait(sieve_processes)
        overall_cputime.add(cputime_slurm)
        aqrels_fin = time.time()
        if is_extra:
            params.timing["aqrels_extra_slurm_sieving"] = aqrels_fin-aqrels_start
        else:
            params.timing["aqrels_slurm_sieving"] = aqrels_fin-aqrels_start
        timeprint("Aqrels slurm jobs all finished!")

    AQRELS_FILE = params.files['AQRELS_FILE']

    if is_extra:
        mode = "a"
    else:
        mode = "w"
    with open(AQRELS_FILE, mode) as outfile:
        for jobnum in jobrange:
            if is_extra:
                jobstr = params.extra_prefix + str(jobnum)
            else:
                jobstr = str(jobnum)

            #jobfile = f"{AQRELS_FILE}.{jobstr}"
            jobfile = params.dirs['TEMP_OUTPUT_DIR'] + f"algrels/alg.{jobstr}"

            if not wait_for_file([jobfile]):
                timeprint(f"{jobfile} does not exist, skipping")
                continue

            with open(jobfile, "r") as infile:
                for line in infile.readlines():
                    _ = outfile.write(line)

    timeprint("Finished writing to AQRELS_FILE!")
    if not is_extra:
        call_todo_sieving(params, is_alg_or_ext=True)


# No timing decorator, fails with fake outsourced params
def do_algebraic_query_sieving(params,jobnum=0,q0=None,q1=None):
    BOUNDA_queries = params.BOUNDA_queries
    AQRELS_FILE = params.files['AQRELS_FILE']
    FBFILE = params.files['CAPPED_FBGZ']
    POLYFILE = params.files['POLYFILE']
    LPB0 = params.parameters['LPB0']
    LPB1_queries = params.parameters['LPB1_queries']

    # AQRELS_FILE_tmp = f"{AQRELS_FILE}.{jobnum}.tmp"
    AQRELS_FILE_tmp = params.dirs['TEMP_OUTPUT_DIR'] + f"algrels/alg.{jobnum}.tmp"

    c0 = 1/4 # 1/2 # Magic constants
    c1 = 1 #4

    if not q0:
        q0 = floor(c0*params.BOUNDA_queries)
    if not q1:
        q1 = c1*params.BOUNDA_queries

    A_used = params.parameters.get('algebraic_query_sieving.A', params.parameters['A_sieving'])
    if int(A_used) > 32 and 'las.bigA.hwloc_job_binding_policy' in params.parameters:
        t_used = params.parameters['las.bigA.hwloc_job_binding_policy']
    else:
        t_used = params.las_job_binding_policy

    lim0_t = params.parameters.get('algebraic_query_sieving.lim', params.BOUNDA_queries)
    lim0 = min(lim0_t, params.BOUNDA_queries, 2**31)

    # We're doing sieving on one side only, so we only use lim0, lpb0, et
    # caetera. The special-q side becomes 0 as well.
    # NOTE: If you update this las call, please also correspondingly update the two las calls in estimate.py so that the estimator remains accurate
    CadoNFS("sieve/las",
            "-sqside", 0,
            "-B", params.parameters.get('algebraic_query_sieving.B', 16),
            "-A", str(A_used),
            "--adjust-strategy", params.parameters.get('algebraic_query_sieving.adjust_strategy', 0),
            "-q0", q0,
            "-q1", q1,
            "-skew", params.poly.skewness,
            "-lpb0", LPB1_queries,
            "-mfb0", params.parameters['sieve.mfb1'],
            "-powlim", params.parameters['sieve.powlim'],
            "-bkmult", params.parameters.get('sieve.bkmult', "1s:1.1"),
            "-bkthresh1", params.parameters.get('sieve.bkthresh1', params.BOUNDA_queries),
            "-poly", 'POLY',
            "-fb0", 'FB1',
            "-out", 'AQRELS',
            "-lim0", str(lim0),
            "-t", t_used,
            "--memory-margin", params.parameters.get('las.memory_margin', 20),
            outputs={'AQRELS': AQRELS_FILE_tmp},
            inputs={
                'FB1': FBFILE,
                'POLY': POLYFILE + ".only-side1",
            }
            )

    # Will need to do later -- on machines with enough memory
    swap_parts_of_relations(AQRELS_FILE_tmp, remove_tmp_suffix=True, do_ext_dedup=None,
                           extra_relation_check=(True, params.files['DEBUG_RENUMBER_FILE']))
    return

def do_todo_query_sieving(params, round_A, round_mfb1, outfile, qmin, qmax, todofilename, extra_flags=[]):
    # Do sieving for this particular A and mfb1.
    # If qmin and qmax are given, go through all special q's in the range.
    # Else, if a list of todo_qs is given, go through the specified special q's.
    # This function works for either algebraic or extension sieving.
    FBFILE = params.files['CAPPED_FBGZ']
    POLYFILE = params.files['POLYFILE']

    outfile_tmp = f"{outfile}.tmp"

    if "extrels" in outfile:
        do_ext_dedup = (True, params.BOUNDA_queries, params.BOUNDA)
    else:
        do_ext_dedup = None

    if int(qmin) < 0:
        qrange = ["-todo", todofilename]
    else:
        qrange = ["-q0", qmin, "-q1", qmax]

    if int(round_A) > 32 and 'las.bigA.hwloc_job_binding_policy' in params.parameters:
        t_used = params.parameters['las.bigA.hwloc_job_binding_policy']
    else:
        t_used = params.las_job_binding_policy

    lim0_t = params.parameters.get('algebraic_query_sieving.lim', params.BOUNDA_queries)
    lim0 = min(lim0_t, params.BOUNDA_queries, 2**31)

    CadoNFS("sieve/las",
        "--memory-margin", params.parameters.get('las.memory_margin', 20),
        "-sqside", 0,
        "-B", params.parameters.get('algebraic_query_sieving.B', 16),
        "-A", round_A,
        "--adjust-strategy", params.parameters.get('algebraic_query_sieving.adjust_strategy', 0),
        "-skew", params.poly.skewness,
        "-lpb0", params.parameters['LPB1_queries'],
        "-mfb0", round_mfb1,
        "-powlim", params.parameters['sieve.powlim'],
        "-bkmult", params.parameters.get('sieve.bkmult', "1s:1.1"),
        "-bkthresh1", params.parameters.get('sieve.bkthresh1', params.BOUNDA_queries),
        "-poly", 'POLY',
        "-fb0", 'FB1',
        "-out", 'OUT',
        "-lim0", str(lim0),
        "--allow-largesq",
        "--allow_compsq",
        *qrange,
        "--never-discard",
        "-t", t_used,
        "-sync",
        "-exit-early", 1,
        *extra_flags,
        outputs={'OUT': outfile_tmp},
        inputs={
            'FB1': FBFILE,
            'POLY': POLYFILE + ".only-side1",
        }
    )
    #swap_parts_of_relations(outfile_tmp, remove_tmp_suffix=True, do_ext_dedup=do_ext_dedup,
    #                        extra_relation_check=(True, params.files['DEBUG_RENUMBER_FILE']))
    return

#alg_or_ext = "algebraic" or "extension"
@timing
def do_algebraic_queries(params, alg_or_ext):
    if alg_or_ext == "algebraic":
        jsonfilename = params.files['AQUERIES_FILE']
        input_rels = params.files['AQRELS_INDEXED']
        todofilename = params.files['AQUERIES_TODO']
    elif alg_or_ext == "extension":
        jsonfilename = params.files['EXT_QUERIES_FILE']
        input_rels = params.files['EXTRELS_INDEXED']
        todofilename = params.files['EXT_QUERIES_TODO']
    else:
        raise ValueError(f"do_algebraic_queries got invalid keyword {alg_or_ext}")

    # N = params.poly.N
    # ZN = Integers(N)

    nqueries = 0

    with open(todofilename, "w") as todofile:
        for rel in indexed_relations_from_file(input_rels):

            # Note -- we don't have check() defined for indexed relations.
            #if not params.cadopoly:
                # TODO: fix this call if it's useful?
                #rel.check(params.poly)
            # ideal = K.ideal(a-b*alpha)
            # ideal_fac = ideal.prime_factors()
            # is_smooth = True
            # for fac in list(ideal_fac):
            #     if fac.norm() > BOUNDA_queries:
            #         is_smooth = False
            #         break

            # if (not is_smooth):
            #     pass

            # print(ZN(rel.norm(params.poly, 0)), file=todofile)
            print(rel.a-rel.b*params.poly.m, file=todofile)

            nqueries += 1

    run_oracle(params, todofilename, jsonfilename)

    return nqueries

def do_fb_extension_sieving(params,jobnum=0,q0=None,q1=None):
    EXTRELS_FILE = params.files['EXTRELS_FILE']
    FBFILE = params.files['CAPPED_FBGZ']
    POLYFILE = params.files['POLYFILE']
    LPB0 = params.parameters['LPB0']
    LPB1_queries = params.parameters['LPB1_queries']

    EXTRELS_FILE_tmp = f"{EXTRELS_FILE}.{jobnum}.tmp"

    if not q0: q0 = str(params.BOUNDA_queries)
    if not q1: q1 = str(params.BOUNDA)

    switch_never_discard = ""
    if params.parameters.get('extension_sieving.never_discard', 0) == 1:
        switch_never_discard = "--never-discard"

    A_used = params.parameters.get('extension_sieving.A', params.parameters['A_sieving'])
    if int(A_used) > 32 and 'las.bigA.hwloc_job_binding_policy' in params.parameters:
        t_used = params.parameters['las.bigA.hwloc_job_binding_policy']
    else:
        t_used = params.las_job_binding_policy

    CadoNFS("sieve/las",
            "-sqside", 0,
            "--memory-margin", params.parameters.get('las.memory_margin', 20),
            "-B", params.parameters.get('extension_sieving.B', 16),
            "-A", str(A_used),
            "-skew", params.poly.skewness,
            "-lpb0", LPB1_queries,
            "-mfb0",
                params.parameters.get('extension_sieving.mfb',
                    params.parameters['sieve.mfb1']),
            "-powlim", params.parameters['sieve.powlim'],
            "-poly", 'POLY',
            "-fb0", 'FB1',
            "-out", 'OUT',
            "-lim0", params.parameters.get('extension_sieving.lim',
                                           params.BOUNDA_queries),
            "--allow-largesq",
            "--adjust-strategy", params.parameters.get('extension_sieving.adjust_strategy', 0),
            "-bkmult", params.parameters.get('sieve.bkmult', "1s:1.1"),
            "-bkthresh1", params.parameters.get('sieve.bkthresh1', params.BOUNDA_queries),
            "-q0", q0,
            "-q1", q1,
            switch_never_discard,
            "-t", t_used,
            "-sync",
            "-exit-early", 1,
            outputs={
                'OUT': EXTRELS_FILE_tmp,
            },
            inputs={
                'FB1': FBFILE,
                'POLY': POLYFILE + ".only-side1",
            }
            )
    swap_parts_of_relations(EXTRELS_FILE_tmp, remove_tmp_suffix=True,
                            do_ext_dedup = (True, params.BOUNDA_queries, params.BOUNDA),
                            extra_relation_check=(True, params.files['DEBUG_RENUMBER_FILE']))
    return

@timing
def call_fb_extension_sieving(params):
    # By the end of this function, all of the extension relations should be in EXTRELS_FILE.
    # That file should contain ONE relation for EACH extension q. Because sieving sometimes
    # doesn't find a relation for some of the q's, this may mean rerunning the one-round
    # function a few times. (Each round produces its own output file, but everything should
    # be coalesced into EXTRELS_FILE, which is used elsewhere.)

    EXTRELS_FILE = params.files['EXTRELS_FILE']

    q0 = params.BOUNDA_queries
    q1 = params.BOUNDA
    #R = params.R

    params.save_to_file()
    if not params.slurm:
        timeprint("Running extension sieving without slurm, setting numjobs = 1")
        numjobs = 1
        command_list = ["time", "-p",
                        params.files['SAGE'],
                        "fb_extension_sieving_helper.py",
                        "--params",params.files['PARAMS'],
                        "--jobnum","0",
                        "--q0",str(q0),
                        "--q1",str(q1)]
        print_command_line(*command_list)
        process = subprocess.run(command_list, stderr=subprocess.PIPE, text=True)
        overall_cputime.add(extract_time(process.stderr))
    else:
        par = params.parameters
        numjobs = par.get('extension_sieving.numjobs',
                          par.get('sieve.numjobs',
                                  par.get('slurm.numjobs')))
        if numjobs is None:
            timeprint("slurm.numjobs not set; defaulting to numjobs = 24")
            numjobs = 24

        sieve_processes = []
        step = floor((q1-q0)/numjobs)
        intervals = [q0+i*step for i in range(numjobs)]+[q1]
        timeprint(intervals)

        for jobnum in range(numjobs):
            command_list = [params.files['SAGE'],
                            "fb_extension_sieving_helper.py",
                            "--params",params.files['PARAMS'],
                            "--jobnum",str(jobnum),
                            "--q0",str(intervals[jobnum]),
                            "--q1",str(intervals[jobnum+1])]
            jobname = f"{params.prefix[:-1]}-extsieve-{q0}-{q1}"
            sieve_processes.append(
                    slurmit(params," ".join(command_list),
                            jobname, jobnum))

        extslurm_start = time.time()
        _, cputime_slurm = slurm_wait(sieve_processes)
        overall_cputime.add(cputime_slurm)
        extslurm_fin = time.time()
        params.timing["extension_slurm_sieving"] = extslurm_fin-extslurm_start
        timeprint(f"extension_slurm_sieving finished in {params.timing['extension_slurm_sieving']}s!")

    with open(EXTRELS_FILE, "w") as outfile:
        for jobnum in range(numjobs):
            jobfile = f"{EXTRELS_FILE}.{jobnum}"

            if not wait_for_file([jobfile]):
                timeprint(f"{jobfile} does not exist, skipping")
                continue

            with open(jobfile, "r") as infile:
                for line in infile.readlines():
                    _ = outfile.write(line)

    keep_only_one_relation_per_q(EXTRELS_FILE, (1, params.BOUNDA_queries, params.BOUNDA), keep_file=True)
    call_todo_sieving(params, is_alg_or_ext=False)


def choose_uv(h):
    # Find small u,v such that u/v = h mod N
    # For this we use lattice basis (0, N) (1, -h)
    # in spirit, it should be just
    # h.rational_reconstruction().as_integer_ratio(), but unfortunately
    # the latter is specified with strict bounds and sometimes fails.
    M = Matrix(ZZ, 2, 2)
    M[0,0] = 0
    M[0,1] = h.parent().characteristic()
    M[1,0] = 1
    M[1,1] = (-1)*ZZ(h)
    Mr = M.LLL()

    u = Mr[0][1]
    v = Mr[0][0]*(-1)
    assert h.parent()(u/v) == h
    return (u, v)

def do_descent_sieving(params, DESCENT_PREFIX):

    parser = argparse.ArgumentParser(description="Descent sieving")
    parser.add_argument("--target",
                        help="Element whose DL is wanted",
                        type=str,
                        required=True)
    descent.GeneralClass.declare_args(parser)
    descent.DescentMiddleClass.declare_args(parser)

    output_name = params.dirs['DESC'] + DESCENT_PREFIX + ".descent.tgt.middle.rels"

    timeprint("output_name:", output_name)

    silent_remove(output_name)
    silent_remove(output_name+".cond")

    inputs = [
        "--poly", params.files['POLYFILE'],
        "--fb1", params.files['CAPPED_FBGZ'],
        "--lim1", params.BOUNDA,
        "--lim0", params.BOUNDR,
        "--I", params.parameters['I_sieving'],
        "--B", 16,
        "--lpb0", params.parameters['LPB0'],
        "--lpb1", params.parameters['LPB1'],
        "--mfb0", params.parameters['desc.mfb0'],
        "--mfb1", params.parameters['desc.mfb1'],
        "--descent-hint", params.files['HINTFILE'],
        "--threads", multiprocessing.cpu_count(),
        "--cadobindir", params.dirs['CADO_BUILD_DIR'],
        "--prefix", DESCENT_PREFIX,
        "--datadir", params.dirs['DESC'],
        "--ell", params.poly.N,
        "--renumber", params.files['RENUMBERFILE'],
        "--no-logs",
        "--target", "tgt"
    ]
    outfile = params.dirs['DESC']+DESCENT_PREFIX
    with open(outfile+".out","w") as fout, open(outfile+".err","w") as ferr:
        fout.reconfigure(line_buffering=True)
        ferr.reconfigure(line_buffering=True)
        with redirect_stdout(fout):
            with redirect_stderr(ferr):
                args = parser.parse_args([str(c) for c in inputs])
                general = descent.GeneralClass(args)
                middle = descent.DescentMiddleClass(general, args)
                relsfile = middle.do_descent(params.files['SPECQ_TODO_FILE']+"."+DESCENT_PREFIX)
    # sanity_check_descent_outfile(relsfile)
    # DEAD code
    timeprint("relsfile:",relsfile)
    return relsfile


class descent_processes_pool():
    def __init__(self, params, seedval=None):
        self.params = params
        self.n_max_descent_processes = params.parameters.get("desc.nprocesses_max", 1)
        self.descent_processes = []
        self.descent_count = 0
        if seedval is not None:
            self.seedval = seedval
        else:
            # we'll do seedval+i, so let's avoid the case where we try
            # the same seeds when we rerun the same job a few seconds
            # later.
            self.seedval = int(time.time()) * 10**6
        self.max_runs = params.parameters['MAX_DESCENT_TRIES']

    def schedule_one_process(self):
        cmd = ["time", "-p",
               self.params.files['SAGE'],
               "descent_helper.py",
               "--params", self.params.files['PARAMS'],
               "--seed", str(self.seedval + self.descent_count)]
        major_message(' '.join(cmd))
        self.descent_processes.append(
                subprocess.Popen(cmd,
                                 shell=False,
                                 preexec_fn=os.setsid,
                                 stderr=subprocess.PIPE, text=True)
                )
        self.descent_count += 1

    def __enter__(self):
        while len(self.descent_processes) < self.n_max_descent_processes:
            self.schedule_one_process()
        return self

    def __exit__(self, *args):
        for P in self.descent_processes:
            if P.poll() is None:
                timeprint("killing",P.pid)
                yy = os.getpgid(P.pid)
                os.kill(P.pid,signal.SIGTERM)
                os.killpg(yy,signal.SIGTERM)
                P.terminate()
                P.wait()
                stderr = P.communicate()[1]
                overall_cputime.add(extract_time(stderr))

    def __len__(self):
        return len(self.descent_processes)

    def poll(self):
        assert len(self)
        finished = []
        pending = []
        for i, P in enumerate(self.descent_processes):
            if P.poll() is None:
                pending.append(P)
            else:
                finished.append(P)
                stderr = P.communicate()[1]
                overall_cputime.add(extract_time(stderr))

        self.descent_processes = pending
        while self.descent_count < self.max_runs and len(self) < self.n_max_descent_processes:
               self.schedule_one_process()
        if not finished:
            time.sleep(1)
        return finished


class descent_slurm_pool():
    def __init__(self, params, existing_init_data=None, seedval=None):
        self.params = params
        self.descent_processes = []
        self.descent_count = 0
        self.existing_init_data = existing_init_data
        if 'slurm.numjobs' in self.params.parameters:
            self.numslurm = self.params.parameters['slurm.numjobs']
        else:
            self.numslurm = 4
        if seedval is not None:
            self.seedval = seedval
        else:
            self.seedval = int(time.time()) * 10**6
        self.max_runs = params.parameters['MAX_DESCENT_TRIES']
        self.cputime_slurm_pool = 0

    def schedule_one_process(self):
        cmd = [self.params.files['SAGE'],
               "descent_helper.py",
               "--params", self.params.files['PARAMS'],
               "--seed", str(self.seedval + self.descent_count)]
        if self.existing_init_data is not None:
            cmd.append("--existing-init-data")
            cmd.append(self.existing_init_data)
        major_message(' '.join(cmd))
        jobnum = self.descent_count
        jobname = f"{self.params.prefix[:-1]}-desc-{jobnum}"
        self.descent_processes.append(
            slurmit(
                self.params,
                ' '.join(cmd),
                jobname,
                jobnum
            )
        )
        self.descent_count += 1

    def __enter__(self):
        while len(self.descent_processes) < self.numslurm:
            self.schedule_one_process()
        return self

    def __exit__(self, *args):
        cancelled_jobs_errfiles = []
        for slurmjob, errfile in self.descent_processes:
            status = get_slurm_job_status(slurmjob)
            if status in ['RUNNING', 'PREEMPTED', 'SUSPENDED']:
                timeprint("Sending SIGINT to slurm job", str(slurmjob))
                p = subprocess.Popen([
                    "scancel",
                    "--signal=SIGINT",
                    "--hurry",
                    str(slurmjob)
                ],stdout=subprocess.PIPE)
                cancelled_jobs_errfiles.append((p, slurmjob, errfile))
            else:
                # Kill forcefully, don't extract any timing since this job was never run in the first place.
                timeprint("Sending SIGKILL to slurm job", str(slurmjob))
                subprocess.Popen([
                    "scancel",
                    str(slurmjob)
                ])

        for p, slurmjob, errfile in cancelled_jobs_errfiles:
            p.wait()
            wait_for_file_content(errfile, ["real ", "user ", "sys "])
            cputime_slurm = extract_time_from_file(errfile)
            timeprint(f"Cancelled slurm job {slurmjob} spent {cputime_slurm}s cputime")
            self.cputime_slurm_pool += cputime_slurm

        overall_cputime.add(self.cputime_slurm_pool)

    def __len__(self):
        return len(self.descent_processes)

    def poll(self):
        assert len(self.descent_processes)
        finished = []
        pending = []
        for slurmjob, errfile in self.descent_processes:
            status = get_slurm_job_status(slurmjob)
            if status in ['RUNNING', 'PENDING', None]:
                pending.append((slurmjob, errfile))
            elif status == 'COMPLETED':
                finished.append(slurmjob)
                self.cputime_slurm_pool += extract_time_from_file(errfile)
            elif status == 'FAILED':
                self.cputime_slurm_pool += extract_time_from_file(errfile)
        self.descent_processes = pending
        while self.descent_count < self.max_runs and len(self.descent_processes) < self.numslurm:
            self.schedule_one_process()
        if not finished:
            time.sleep(1)
        return finished

@timing
def call_descents(target, params):
    params.save_to_file()

    # just _any_ run that succeeds is good. We'll pick the result from
    # files named like this. So we start by removing potential traces of
    # the older ones.
    wildcard = params.dirs['TEMP_OUTPUT_DIR']+"tgt.json*"

    glob_remove(wildcard)
    with descent_processes_pool(params) as D:
        while len(D):
            for P in D.poll():
                if P.returncode != 0:
                    warning_message("Descent process finished"
                                    f" with return code {P.returncode}:",
                                    *P.args)
                    continue
                for filename in glob.glob(wildcard):
                    shutil.copyfile(filename, params.files['TGT_INFO'])
                    with open(filename, "r") as fp:
                        D = json.load(fp)

                    def cast(D, field, parent):
                        D[field] = parent(D[field])
                    ZN = Integers(params.poly.N)
                    cast(D, 'tgt', ZN)
                    cast(D, 'mask', ZN)
                    cast(D, 'h', ZN)
                    cast(D, 'u', ZZ)
                    cast(D, 'v', ZZ)
                    cast(D, 'init_time', ZZ)
                    cast(D, 'middle_time', ZZ)
                    params.target_info = D
                    params.timing['descent_init'] = D['init_time']
                    params.timing['descent_middle'] = D['middle_time']
                    return params.target_info
    return None

@timing
def call_descents_slurm(target, params, existing_init_data=None):
    params.save_to_file()
    wildcard = params.dirs['TEMP_OUTPUT_DIR']+"tgt.json*"
    glob_remove(wildcard)

    with descent_slurm_pool(params, existing_init_data=existing_init_data) as Pool:
        while len(Pool):
            for completed_jobs in Pool.poll():
                for filename in glob.glob(wildcard):
                    shutil.copyfile(filename, params.files['TGT_INFO'])
                    with open(filename, "r") as fp:
                        tgtfile = json.load(fp)

                    def cast(tgtfile, field, parent):
                        tgtfile[field] = parent(tgtfile[field])
                    ZN = Integers(params.poly.N)
                    cast(tgtfile, 'tgt', ZN)
                    cast(tgtfile, 'mask', ZN)
                    cast(tgtfile, 'h', ZN)
                    cast(tgtfile, 'u', ZZ)
                    cast(tgtfile, 'v', ZZ)
                    cast(tgtfile, 'init_time', ZZ)
                    cast(tgtfile, 'middle_time', ZZ)
                    params.target_info = tgtfile
                    params.timing['descent_init'] = tgtfile['init_time']
                    params.timing['descent_middle'] = tgtfile['middle_time']
                    return params.target_info
    return None

def call_descent_helper(target,params,seedval):
    seedval = int(seedval)
    cmdline = [
        "time", "-p",
        params.files['SAGE'],
        "descent_helper.py",
        "--params", params.files['PARAMS'],
        "--seed", str(seedval)
    ]
    print(f"Run descent helper for target {target} with command:\n{' '.join(cmdline)}")
    completed_descent = subprocess.run(cmdline, stderr=subprocess.PIPE, text=True)
    overall_cputime.add(extract_time(completed_descent.stderr))

    print(f"yay finished with seedval {seedval}")
    print("returncode",completed_descent.returncode)

    if completed_descent.returncode == 0:
        timeprint("Descent success")
        for filename in glob.glob(params.dirs['TEMP_OUTPUT_DIR']+"tgt.json*"):
            with open(filename, "r") as fp:
                return json.load(fp)
    else:
        timeprint("Descent failed")
    return None

class TakenLineMissing(Exception):
    def __init__(self, missed, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.missed = missed
    def __str__(self):
        m = "; ".join([f"{side},{q},{rho}" for side,q,rho in self.missed])
        return f"\"Taken\" line missing for special-qs: [{m}]"

def sanity_check_descent_outfile(g, descent_file, todofile=None):
    # We want to make sure the descent output file has a Taken line for each special-q.
    # Technically, we only care about this for special-q's that go into constructing S,
    # but that would essentially recreate the logic of construct_S. For now, just check
    # that each special-q job has a Taken line.
    timeprint("Starting sanity_check_descent_outfile")

    all_qs = set()
    taken_qs = set()

    if todofile is not None:
        with open(todofile, "r") as f:
            for line in f.readlines():
                # Note: most todofiles have only side-0 primes
                # but there could be side-1 primes from large descent
                ss = line.strip().split()
                ssi = tuple([Integer(i) for i in ss])
                assert len(ssi) in [2,3]
                assert ssi[0] in [0,1]

                if ssi[0] == 1:
                    all_qs.add(ssi)
                else:
                    qq = ssi[1]
                    rho = g.roots(GF(qq))[0][0]
                    all_qs.add((0, qq, rho))

    for line in open(descent_file).readlines():
        line = line.strip()

        if m := re.match(r"^# Now sieving side-(\d+) q=(\d+); rho=(\d+)", line):
            side, sq, rho = (int(c) for c in m.groups())
            all_qs.add((side,sq,rho))
            continue

        if m := re.match(r"^# [descent] pushing side-(\d+) q=(\d+); rho=(\d+) .* to todo list .*", line):
            side, sq, rho = (int(c) for c in m.groups())
            all_qs.add((side,sq,rho))
            continue

        if m := re.match(r"# Taking decision on .* side-(\d+) q=(\d+); rho=(\d+)", line):
            # To be stored in taken_qs if next line matches if condition below
            side, sq, rho = (int(c) for c in m.groups())

        if (m := re.match("^Taken: (.*)", line)):
            taken_qs.add((side, sq, rho))

    if len(all_qs) > len(taken_qs):
        raise TakenLineMissing(all_qs - taken_qs)
    return

def fast_ab_parsing(rel_lines):
    l = []
    for rel_line in rel_lines:
        if not rel_line.startswith(b'#'):
            a, b = rel_line.split(b':')[0].split(b',')
            l.append((int(a, 16), int(b, 16)))
    return l

class LinearAlgebraMatrix(abc.ABC):
    def __init__(self, params, indexed_relations_file):

        timeprint("Inside init of LinearAlgebraMatrix")

        # Only store relative path persistently to allow copying
        self.indexed_relations_file_relative = self.convert_absolute_to_relative_path(params, indexed_relations_file)

        timeprint("Inside init of LinearAlgebraMatrix; starting to load row_to_aquery")

        # pre-read this.
        self.row_to_aquery = []
        # works for linux and macos
        mem_total_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        if params.nthreads > 1 and 2 * os.path.getsize(indexed_relations_file) < mem_total_bytes:
            # This is less memory efficient but faster than the one thread version.
            # Since indexed_relations_file was ~17GB for n1024, it's fine to hold both
            # the full file and the parsed a-b-pairs in memory at the same time.
            timeprint(f"loading {indexed_relations_file} into memory")
            with open(indexed_relations_file, mode="r", encoding="ascii") as indexed_rels_fp:
                with mmap.mmap(indexed_rels_fp.fileno(), length=0, access=mmap.ACCESS_READ) as indexed_rels_mmap:
                    timeprint(f"Starting to parse {indexed_relations_file} with {params.nthreads} processes")
                    indexed_rels_lines = []
                    l = 0
                    while True:
                        line = indexed_rels_mmap.readline()
                        l += 1
                        if line == b"":
                            break
                        indexed_rels_lines.append(line)
                    timeprint(f"Finished splitting {indexed_relations_file} into lines.")

                with ProcessPoolExecutor(max_workers=params.nthreads) as executor:
                    batch_size = math.ceil(l / params.nthreads)
                    batch_args = [indexed_rels_lines[i * batch_size : (i+1) * batch_size] for i in range(params.nthreads)]
                    for row_to_aquery_batch in tqdm.tqdm(executor.map(fast_ab_parsing, batch_args)):
                        self.row_to_aquery += row_to_aquery_batch
                timeprint(f"Done parsing")
        else:
            for irel in tqdm.tqdm(indexed_relations_from_file(indexed_relations_file)):
                self.row_to_aquery.append((irel.a, irel.b))

        timeprint("Finished loading row_to_aquery")

        self.Ze = IntegerModRing(params.parameters['e'])

    def convert_absolute_to_relative_path(self, params, abs_path):
        return abs_path.replace(params.dirs['TEMP_OUTPUT_DIR'], "")

    def convert_relative_to_absolute_path(self, params, rel_path):
        return os.path.join(params.dirs['TEMP_OUTPUT_DIR'], rel_path)

    def base_ring(self):
        return self.Ze

    @abc.abstractmethod
    def ab_per_row(self, params):
        pass

    def indexed_relations_filename(self, params):
        # XXX: for backwards compatibility with matrices persistently stored on disk
        # before making them path independent.
        if not hasattr(self, 'indexed_relations_file_relative'):
            prefix_idx = self.indexed_relations_file.index(params.prefix)
            self.indexed_relations_file_relative = self.indexed_relations_file[prefix_idx + len(params.prefix):]

        return self.convert_relative_to_absolute_path(params, self.indexed_relations_file_relative)

    def all_ab(self):
        """
        this iterator returns all the (a,b) pairs that were used to
        construct the matrix.
        """
        for ab in self.row_to_aquery:
            yield ab


def construct_sage_matrix(params, indexed_relations_file):
    e = params.parameters['e']
    R = params.R
    Ze = Integers(e)

    num_cols = R.number_of_fb_valuations(params.BOUNDA_queries)

    M_entries = defaultdict(Ze)

    row_to_aquery = []
    for r, rel in enumerate(indexed_relations_from_file(indexed_relations_file)):
        for idx in rel.indices:
            col = R.renumber_to_column(idx)
            # FIXME J_valuation_inconsistency
            # to be activated someday
            if False and col == 0 and not params.poly.f[1].is_monic():
                M_entries[r, col] -= 1
            elif col >= 0:
                M_entries[r, col] += 1
        row_to_aquery.append((rel.a,rel.b))

    M0 = matrix(Ze, len(row_to_aquery), num_cols, dict(M_entries), sparse=True)
    Z = matrix(Ze, len(row_to_aquery), 1)

    # XXX oddly enough, the existing code is padding M with two zero
    # columns...
    # M = block_matrix(1,3,[Z,M0,Z])
    M = M0

    return M, row_to_aquery

class LinearAlgebraMatrix_IndexedFileOnly(LinearAlgebraMatrix):
    def __init__(self, params, indexed_relations_file):
        super().__init__(params, indexed_relations_file)
        self.M, self.row_to_aquery = construct_sage_matrix(params,
                                                           indexed_relations_file)
        self.nrows = self.M.nrows()
        self.ncols = self.M.ncols()

    def matrix(self):
        return self.M

    def ab_per_row(self, params):
        """
        For each row of the matrix, return the combination of all (a,b)
        pairs that led to it, with exponents.

        Of course in the case of the current class, it's a fairly trivial
        thing.
        """
        return [{x:1} for x in self.row_to_aquery]

    def __str__(self):
        return "matrix built from the full set of " \
               f"{len(self.row_to_aquery)} indexed relations"



class LinearAlgebraMatrix_Filtered(LinearAlgebraMatrix):
    def __init__(self,
                 params,
                 indexed_relations_file,
                 purged_file,
                 ideals_file, index_file,
                 matrix_file):

        timeprint("Inside init of LinearAlgebraMatrix_Filtered")

        super().__init__(params, indexed_relations_file)

        self.params = params
        self.matrix_file = self.convert_absolute_to_relative_path(params, matrix_file)
        self.ideals_file = self.convert_absolute_to_relative_path(params, ideals_file)
        self.index_file = self.convert_absolute_to_relative_path(params, index_file)
        self.purged_file = self.convert_absolute_to_relative_path(params, purged_file)

        rwfile = re.sub(r'\.bin$', '.rw.bin', matrix_file)
        cwfile = re.sub(r'\.bin$', '.cw.bin', matrix_file)
        self.nrows = os.stat(rwfile).st_size // 4
        self.ncols = os.stat(cwfile).st_size // 4
        self.ncoeffs = ((os.stat(matrix_file).st_size // 4) - self.nrows) // 2

        from cado.tests.sagemath.cado_sage import bwc

        timeprint("Inside init of LinearAlgebraMatrix_Filtered; starting bwc setup")

        self.bwcparams = bwc.BwcParameters(m=6, n=6,
                                           p=ZZ(self.params.parameters['e']))

        sage_or_scipy = True        # sage
        if params.scipy_matrix:
            sage_or_scipy = False   # scipy

        self.sage_or_scipy = sage_or_scipy
        self.M = bwc.BwcMatrix(self.bwcparams, matrix_file, nthreads=self.params.nthreads, sage_or_scipy=sage_or_scipy)

        timeprint("Inside init of LinearAlgebraMatrix_Filtered; starting M.read()")

        self.M.read()
        timeprint(self.M)

        self._expand_map=[]
        with open(ideals_file) as f:
            line = next(f)
            m = re.match(r'# (\d+)$', line)
            assert m
            assert self.ncols == int(m.group(1))
            for i, line in enumerate(f.readlines()):
                ii, j = line.strip().split()
                assert int(ii) == i
                self._expand_map.append(int(j, 16))
        self._shrink_map={xj: j for j,xj in enumerate(self._expand_map)}

    def get_matrix_file(self, params):
        return self.convert_relative_to_absolute_path(params, self.matrix_file)

    def get_ideals_file(self, params):
        return self.convert_relative_to_absolute_path(params, self.ideals_file)

    def get_index_file(self, params):
        return self.convert_relative_to_absolute_path(params, self.index_file)

    def get_purged_file(self, params):
        return self.convert_relative_to_absolute_path(params, self.purged_file)

    def matrix(self):
        return self.M.M

    def __str__(self):
        return   f"cado-nfs matrix of size {self.nrows} x {self.ncols}" \
               + f" with {self.ncoeffs} non-zero coefficients"

    def ab_per_row(self, params):
        """
        For each row of the matrix, return the combination of all (a,b)
        pairs that led to it, with exponents.
        """
        # Note that row_to_aquery is mapping indices/rows in the *full* list of relations.
        # The index file is mapping indices/rows in the *purged* file.

        timeprint("Starting ab_per_row")

        purge_row_to_ab = []
        with open(self.get_purged_file(params)) as pf:
            nlines = int(next(pf).split()[1])
            for line in pf.readlines():
                # ed,58:11,1a,ad,224,42a,482
                a = int(line.split(":")[0].split(",")[0], 16)
                b = int(line.split(":")[0].split(",")[1], 16)
                purge_row_to_ab.append((a,b))

        assert(nlines == len(purge_row_to_ab))

        timeprint("Finished reading purged file. Starting to read index file")

        with open(self.get_index_file(params)) as f:
            nrows = int(next(f))
            assert nrows == self.nrows
            for line in f.readlines():
                length, *coeffs = line.strip().split()
                assert int(length) == len(coeffs)
                yield {purge_row_to_ab[int(j, 16)]: int(c)
                       for j, c in [x.split(':') for x in coeffs]}


    def column_shrink_map(self):
        """
        this maps the column range of the indexed relations to the column
        range of the matrix in matrix_file
        """
        return self._shrink_map

    def column_expand_map(self):
        return self._expand_map




@timing
def construct_S(params, DRELS_FILE, DRELS_INDEXED,
                u, v, already_indexed=False, partial_R=False, uv_fac=None):
    """
    Return the dictionary {(a,b):exp)} terms as well as the S_alg_vector and S_rat_vector
    """
    BOUNDR = params.BOUNDR
    BOUNDA = params.BOUNDA
    #u = params.target_info['u']
    #v = params.target_info['v']

    if partial_R:
        R = PartialRenumber1024(params.poly)
        partial_renum_filename = params.dirs['TEMP_OUTPUT_DIR'] + 'partial.renumber.map'
        R.load(partial_renum_filename)
        R.do_sanity_checks()
    else:
        R = params.R

    S_rat_vector = vector(ZZ, R.number_of_rational_queries(), sparse=True)
    timeprint(f"rat vector in dimension {len(S_rat_vector)}")
    S_alg_vector = vector(ZZ, R.number_of_fb_valuations(params.BOUNDA), sparse=True)
    timeprint(f"alg vector in dimension {len(S_alg_vector)} (bound={params.BOUNDA})")

    # (side, q) --> (a,b)
    special_qs = dict()

    # (a,b) --> [rough factors]
    rough_alg_factors = dict()
    rough_rat_factors = dict()

    taken = list(extract_taken_relations(DRELS_FILE))

    ext_unlinked_qs = set()

    if os.path.exists(params.files['EXT_UNLINKED_TODOS']):

        with open(params.files['EXT_UNLINKED_TODOS'], "r") as extfile:
            for line in extfile.readlines():
                if line.strip() == "":
                    continue

                # it's a todo file
                lineinfo = line.strip().split(" ")
                assert lineinfo[0] == '0'
                q = Integer(lineinfo[1])
                r = Integer(lineinfo[2])
                ext_unlinked_qs.add((1,q,r))

    # Here we use pull_large_factors to __modify__ the taken relations,
    # and remove the factors above BOUNDR and BOUNDA. This means in
    # particular that the indexed relations that we produce are no longer
    # valid relations: they only make sense if we consider them together
    # with these rough factors.
    for rel, (side, sq, rho) in taken:
        Rf, Af = rel.pull_large_factors((BOUNDR, BOUNDA), (side, sq), ext_unlinked_qs)
        rough_rat_factors[int(side), int(sq), rel.a, rel.b] = Rf
        rough_alg_factors[int(side), int(sq), rel.a, rel.b] = Af
        special_qs[side, sq] = (rel.a, rel.b)

    rels = [rel for rel, q in taken]

    if already_indexed:
        # DRELS_INDEXED is a file of indexed relations,
        # by some other means.
        indexed_relations_file = indexed_relations_from_file(DRELS_INDEXED)
    else:
        # The usual case
        indexed_relations_file = big_convert_to_indexed_relation(rels,
                                               params,
                                               DRELS_INDEXED)

    indexed = { (irel.a,irel.b): irel for irel in indexed_relations_file }

    from collections import deque

    # rough, and outstanding (same items, for the moment): the set of
    # special-q's that need to be descended.
    rough = deque()
    outstanding = defaultdict(int)

    # indices: what we find in the relations. Some of them are rational
    # primes, but they only get converted in the end.
    # As regards the rest, because we work with relations from which we
    # pulled all factors above the linear algebra bounds, these really
    # all go in the matrix.
    indices = defaultdict(int)

    # rat_primes: the rational primes in the factorization. Initially, we
    # only have those from the factorization of u/v
    rat_primes = defaultdict(int)

    # rel_combination: how we arrange relations together.
    rel_combination = defaultdict(int)

    if uv_fac is not None:
        fac_target = uv_fac
    else:
        # silly to do, since we have this information saved in a file,
        # but it is fast since u and v are very smooth.
        fac_target = factor(ZZ(u)/ZZ(v))

    timeprint("Factorization of target is", fac_target)

    for p, k in fac_target:
        if p >= BOUNDR:
            outstanding[0, p] += k
            rough.append((0, p))
        else:
            S_rat_vector[R.rational_prime_to_prime_index(p)] += k

    while rough:
        side, q = rough.popleft()
        exponent = outstanding[side, q]

        S_exponent = -exponent

        a, b = special_qs[side, q]

        print(f"killer relation for {(side,q)} is", indexed[a, b])
        if rough_rat_factors:
            print("attached rough factors (rat)", rough_rat_factors[int(side), int(q), a,b])
        if rough_alg_factors:
            print("attached rough factors (alg)", rough_alg_factors[int(side), int(q), a,b])

        rel_combination[a,b] += S_exponent

        for fac in indexed[a,b].indices:
            # FIXME J_valuation_inconsistency
            # to be activated someday
            if False and fac == 0 and not params.poly.f[1].is_monic():
                indices[fac] -= S_exponent
            else:
                indices[fac] += S_exponent

        for aa in rough_alg_factors[(int(side), int(q), a,b)]:
            outstanding[1, aa] += S_exponent
            rough.append((1, aa))

        for rr in rough_rat_factors[(int(side), int(q), a,b)]:
            outstanding[0, rr] += S_exponent
            rough.append((0, rr))

        outstanding[side, q] += S_exponent
        assert outstanding[side, q] == 0

    for fac, v in indices.items():
        col = R.renumber_to_column(fac)
        if col >= 0:
            S_alg_vector[col] += v
        else:
            S_rat_vector[R.renumber_to_rational_prime_index(fac)] += v

    # The S vector is over ideals in the extended factor base.
    # We look for the BOUNDA_queries unlinked ideals in 'FILT_UNLINKED_IDEALS'
    # which is the total set after filtering.
    # We look for the BOUNDA unlinked ideals in 'EXT_UNLINKED_IDEALS'.
    ul_ideal_cols = set()

    if os.path.exists(params.files['FILT_UNLINKED_IDEALS']):
        with open(params.files['FILT_UNLINKED_IDEALS'], "r") as ul_file:
            for renum_index in ul_file.readlines():
                ul_col = R.renumber_to_column(int(renum_index))
                if ul_col >= 0:
                    ul_ideal_cols.add(ul_col)

    if os.path.exists(params.files['EXT_UNLINKED_IDEALS']):
        with open(params.files['EXT_UNLINKED_IDEALS'], "r") as ul_file:
            for line in ul_file.readlines():
                renumber_index = int(line.strip())
                ul_col = R.renumber_to_column(renumber_index)
                if ul_col >= 0:
                    ul_ideal_cols.add(ul_col)

    for ul in ul_ideal_cols:
        if S_alg_vector[ul] != 0:
            ul_renum = R.column_to_renumber(ul)
            ul_renum_hex = hex(ul_renum)
            raise RuntimeError("construct_S:"
                                f" The relation involves ideal number {ul_renum}"
                                f" or {ul_renum_hex}"
                                " (in renumber coordinates),"
                                f" which is {ul}"
                                " (in column coordinates),"
                                " however it is not linked to the others."
                                " Too bad, really.")
            exit(0)

    return rel_combination, S_alg_vector, S_rat_vector


def extract_taken_relations(DRELS_FILE):
    """
    This is a bit like condense_descent_relations, except that we don't
    rendezvous on a file, and we just yield the relations (with the
    special-q attached) instead
    """
    for line in open(DRELS_FILE).readlines():
        line = line.strip()

        if m := re.match(r"# Taking decision on .* side-(\d+) q=(\d+); rho=(\d+)", line):
            # To be stored in taken_qs if next line matches if condition below
            side, sq, rho = (int(c) for c in m.groups())

        if (m := re.match("^Taken: (.*)", line)):
            yield (las_relation(m.group(1)), (side, sq, rho))

@timing
def filter_relation_file(params):
    #filename = params.files['AQRELS_FILE']
    #TEMP_OUTPUT_DIR = dirs['TEMP_OUTPUT_DIR']
    POLYFILE = params.files['POLYFILE']
    RENUMBERFILE = params.files['RENUMBERFILE']
    # Return the filename of the filtered relations

    # For 1024: we can use the 31-version of the renumber file,
    # which should be faster to load.
    # ^ removed for now

    basename = os.path.basename(params.files['AQRELS_FILE'])
    indexed_relations_file = params.files['AQRELS_FILE'] + ".indexed"
    purged_file            =  indexed_relations_file + ".purged"
    relsdel_file           =  indexed_relations_file + ".relsdel"
    history_file           =  indexed_relations_file + ".history"
    index_file             =  indexed_relations_file + ".index"
    ideals_file            =  indexed_relations_file + ".ideals"
    matrix_file            =  indexed_relations_file + ".matrix.bin"

    make_and_clean(params.dirs['DUP'])

    # we don't need to bother with dup1 doing anything fun. So really,
    # it's just about making a copy (and stripping out comments)

    CadoNFS("filter/dup1",
            "-out", 'TMPDIR',
            "-prefix", basename + ".dup",
            "-n", 0, 'DATA',
            inputs={'TMPDIR': params.dirs['TEMP_OUTPUT_DIR'][:-1],
                    'DATA': params.files['AQRELS_FILE']})

    # count the number of relations in all files from L
    wc_l = lambda L: sum([len(open(f).readlines()) for f in L])

    files_to_dup2 = glob.glob(params.dirs['DUP']+'*.dup.*')
    nrels = wc_l(files_to_dup2)

    CadoNFS("filter/dup2",
            "-dl",
            "-nrels", nrels,
            "-poly", 'POLY',
            "-renumber", 'RENUMBER',
            *files_to_dup2,
            inputs={'POLY': POLYFILE, 'RENUMBER': RENUMBERFILE})


    # We no longer call this function to renumber the relations from
    # extrels or from the descent.
    assert "desc" not in params.files['AQRELS_FILE']
    assert "extrels" not in params.files['AQRELS_FILE']

    # so really, it's only about algebraic query relations.
    assert "aqrels" in params.files['AQRELS_FILE']


    # record the collection of renumbered relations in a single file, for
    # potential future use.
    with open(indexed_relations_file, "w") as out:
        for f in files_to_dup2:
            with open(f) as slice:
                for line in slice.readlines():
                    out.write(line)


    if not params.cado_nfs_filter:
        # do *not* using the cado-nfs filtering tools.
        warning_message("***** skipping purge, copying rels to",
                        indexed_relations_file)

        return LinearAlgebraMatrix_IndexedFileOnly(params,
                                                indexed_relations_file)

    nrels = wc_l(files_to_dup2)
    major_message(f"Number of algebraic relations after dup2: {nrels}")

    CadoNFS("filter/purge",
            "-col-max-index", params.BOUNDA_queries,
            "-col-min-index", 0,
            "-nrels", nrels,
            "-out", 'PURGED',
            "-outdel", 'RELSDEL',
            # It's important that we "keep" at least as many relations as the
            # number of Schirokauer maps.
            "-keep", #10,
            params.parameters.get('purge.keep', 10),
            *files_to_dup2,
            outputs={'PURGED': purged_file, 'RELSDEL': relsdel_file})


    CadoNFS("filter/merge-dl",
            "-mat", 'PURGED',
            "-out", 'HIST',
            # 10 is good for the n60, but otherwise it's really too small!
            "-target_density",
            params.parameters.get('merge.target_density', 10),
            "-t", str(params.parameters.get('nthreads', 8)),
            inputs={'PURGED': purged_file},
            outputs={'HIST': history_file})


    CadoNFS("filter/replay-dl",
            "-purged", 'PURGED',
            "-his", 'HIST',
            "-index", 'INDEX',
            "-ideals", 'IDEALS',
            "-out", 'MATRIX',
            inputs={'PURGED': purged_file,
                    'HIST': history_file},
            outputs={'INDEX': index_file,
                    'IDEALS': ideals_file,
                    'MATRIX': matrix_file})

    M = LinearAlgebraMatrix_Filtered(params,
                                     indexed_relations_file,
                                     purged_file,
                                     ideals_file, index_file,
                                     matrix_file)

    return M


class CadoExplainRenumberFile:
    """
    This class takes inspiration (and copies much code!) from
    CadoIdealsDebugFile, except that we don't want to waste time
    computing the ideals.

    (now that we cheat on the max order computation, there's probably not
    much point in duplicating code)
    """
    def __init__(self, poly, filename):
        """
        poly is a CadoPolyFile object.
        filename is expected to be the file produced by
        explain_indexed_relation -all -raw
        (yes, explain_indexed_relation is really poorly named, since in
        effect it's used as a renumber table dump tool)
        """
        self.poly = poly
        self.filename = filename
        self.__clear_fields_for_read()
        self.has_merged_J = False
        self.index_of_J = []

    def __clear_fields_for_read(self):
        self._ideals = None
        self._indices_per_side = None

    def __repr__(self):
        return ("CadoExplainRenumberFile("
                + f"CadoPolyFile(\"{self.poly.filename}\")"
                + f", \"{self.filename}\")")

    def __str__(self):
        rep = "cado-nfs proxy to the renumber table"
        if self.filename:
            rep += f" (path={self.filename})"
        else:
            rep += " (transient)"
        if not self.relsets:
            rep += ", no data read yet"
            return rep
        else:
            rep += f", {len(self.ideals)} ideals"
            return rep

    def __len__(self):
        return len(self._ideals)

    # def __iter__(self):
    #     return iter(self._ideals)
    #
    # def __getitem__(self, i):
    #     return self._ideals[i]

    def index(self, i):
        return self._ideals.index(i)

    def read(self):
        self.__clear_fields_for_read()

        if get_verbose():
            timeprint(f"Reading {self.filename}")

        K = self.poly.K

        self._ideals = []
        #self._ideals_to_index = dict()
        self._indices_per_side = [[] for f in K]

        with cat_or_zcat(self.filename) as fp:
            for t in fp:
                if t.startswith('#'):
                    continue

                parser, *data = t.split()

                i = len(self._ideals)

                if parser == 'J':
                    # the "data" field is actually a bit of a lie, b
                    self.has_merged_J = len(data) > 1
                    if self.has_merged_J:
                        # We're going to cheat, and not return an ideal
                        side = tuple(int(s) for s in data)
                        I = tuple([f"J{s}" for s in data])
                    else:
                        side = int(data[0])
                        I = (f"J{side}",)
                    self.index_of_J.append(len(self._ideals))
                elif parser == 'rat':
                    side, p = data
                    side = int(side)
                    p = ZZ(p)
                    I = (p,)
                    #self._ideals_to_index[side,p] = i
                elif parser == 'proj':
                    side, p = data
                    side = int(side)
                    p = ZZ(p)
                    I = (p,)
                    # we're going to have problems here if we have projective
                    # primes.
                    #self._ideals_to_index[side,p] = i
                elif parser == 'easy':
                    side, p, r = data
                    side = int(side)
                    p = ZZ(p)
                    r = ZZ(r)
                    I = (p, r)
                    #self._ideals_to_index[side,p,r] = i
                elif parser == 'generic':
                    side, p, denom, *coeffs = data
                    side = int(side)
                    p = ZZ(p)
                    denom = ZZ(denom)
                    try:
                        theta = K[side]([ZZ(c) for c in coeffs]) / denom
                    except Exception as e:
                        print(side, p, denom, coeffs)
                        raise e
                    I = (p, theta)
                    #self._ideals_to_index[side,p,theta] = i

                # _indices_per_side[side] is always a list of indices in the
                # renumber table (that is, in self._ideals) of all ideals
                # that have something to do with the given {side}.

                # reciprocally, self._ideals[i] has info on the index of this
                # among ideals on the same side. Watch out for the special
                # treatment of the J ideals, though!

                if type(side) is tuple:
                    # This is a special case for J, which is sometimes
                    # attached to two sides.
                    side_restricted_index = tuple(len(self._indices_per_side[s]) for s in side)
                    for s in side:
                        self._indices_per_side[s].append(i)
                else:
                    side_restricted_index = len(self._indices_per_side[side])
                    self._indices_per_side[side].append(i)

                self._ideals.append((parser, side, I, side_restricted_index))

                END_READING_EARLY = True
                if END_READING_EARLY and parser != 'J' and int(I[0]) > 2**31:
                    break
                    # When handling only LPB1_queries ideals,
                    # we don't need to read beyond 2**31,
                    # and it would be expensive to do so.
                    # Comment this out when we do extension things.

            assert not self.has_merged_J or len(self.index_of_J) == 1

    def ideal_to_index(self, q):
        """
        q is a tuple (side, *things)
        """
        raise NotImplementedError("AFAIK ideal_to_index is never called.")
        #return self._ideals_to_index[q]

    def index_to_ideal(self, i):
        parser, side, I, col_index = self._ideals[i]
        # This is really only something we can understand if we're away
        # from the special cases such as J, bad ideals, and so on.
        return side, *I

    def side_and_index_to_ideal(self, side, i):
        rside, *I = self.index_to_ideal(self._indices_per_side[side][i])
        assert rside == side
        return side, *I

    def ideal_to_side_and_index(self, q):
        """
        q is a tuple (side, *things)
        """
        raise NotImplementedError("AFAIK ideal_to_side_and_index is never called.")
        #i = self._ideals_to_index[q]
        #parser, side, I, col_index = self._ideals[i]
        #assert side == q[0]
        #return side, col_index


    def renumber_to_column(self, idx):
        """
        given an index to the renumber table, return the column index in
        the valuation matrix.
        This returns -1 if the index corresponds to a rational prime.
        """

        parser, side, I, col_index = self._ideals[idx]

        if type(side) is tuple:
            assert len(side) == 2
            # it had better be 0 on both sides, otherwise we don't really
            # know what to return...
            assert col_index[0] == col_index[1]
            return col_index[0]

        return -1 if side == 0 else col_index

    def column_to_renumber(self, col_index):
        return self._indices_per_side[1][col_index]

    def column_to_sage_ideal(self, col_index, side_hint=None):
        # XXX: Note that there's a parallelizable version of this function in
        # montgomery_ethroot.py that has close to identical code. Any updates
        # to this function should likely also be applied there.

        parser, side, Idata, _col_index = self._ideals[self.column_to_renumber(col_index)]

        if parser == 'J' and self.has_merged_J:
            # This case is special, really.
            assert side_hint is not None    # what can we do?
            assert len(side) == len(_col_index)
            assert _col_index == (_col_index[0],)*len(side)
            side = side[side_hint]
            _col_index = _col_index[side_hint]

        assert _col_index == col_index

        # XXX this assert is a bit excessive in full generality.
        # Perhaps we'd like to assert side == side_hint, at most.
        assert side == 1

        K = self.poly.K
        J = self.poly.nt.J()
        OK = self.poly.nt.maximal_orders()

        if parser == 'J':
            return J[side]
        elif parser == 'rat':
            (p,) = Idata
            I = OK[side].fractional_ideal(p)
        elif parser == 'proj':
            (p,) = Idata
            I = OK[side].fractional_ideal(p) + J[side]
        elif parser == 'easy':
            (p, r) = Idata
            I = OK[side].fractional_ideal(p, K[side].gen() - r) * J[side]
        elif parser == 'generic':
            (p, theta) = Idata
            I = OK[side].fractional_ideal(p, theta)
        else:
            raise AssertionError("Unknown parser:", parser)
        return I

    def renumber_to_rational_prime(self, idx):
        parser, side, I, col_index = self._ideals[idx]
        if parser == 'J':
            # can't we? Sure it's not a prime, we still would have
            # something to return, wouldn't we?
            raise RuntimeError("We can really not convert J to a rational prime")
        if side == 1:
            raise KeyError
        assert parser == 'rat'
        return I[0]

    def renumber_to_rational_prime_index(self, idx):
        parser, side, I, col_index = self._ideals[idx]
        if parser == 'J':
            raise RuntimeError("J does not have an index as a rational prime")
        if side == 1:
            raise KeyError
        assert parser == 'rat'
        return col_index - int(self.has_merged_J)

    def _rational_prime_to_prime_index_raw(self, p):
        """
        This returns the prime index, possibly offset by 1 if we have a
        merged J
        """
        lo = int(self.has_merged_J)
        b = bisect.bisect_left(self._indices_per_side[0], p, lo=lo,
                               key=lambda x: self._ideals[x][2][0])
        if b == len(self._indices_per_side[0]):
            raise KeyError(f"rational prime p={p} is not in the renumber table")

        j = self._indices_per_side[0][b]
        if self._ideals[j] == ('rat', 0, (p,), b):
            pass
        elif not p.is_prime(proof=False):
            raise ValueError(f"p={p} is not prime!")
        else:
            raise RuntimeError(f"rational prime p={p} is not"
                               " in the renumber table, which looks like"
                               " a data consistency bug."
                               " First entry above:"
                               f"{self._ideals[j]}")
        return b, j


    def rational_prime_to_prime_index(self, p):
        b, r = self._rational_prime_to_prime_index_raw(p)
        return b - int(self.has_merged_J)

    def rational_prime_to_renumber(self, p):
        b, r = self._rational_prime_to_prime_index_raw(p)
        return r

    def rational_prime_index_to_prime(self, i):
        lo = int(self.has_merged_J)
        b = i + lo
        j = self._indices_per_side[0][b]
        p = self._ideals[j][2][0]
        assert self._ideals[j] == ('rat', 0, (p,), b)
        return p

    def rational_ideal_index_to_prime(self, b):
        """
        This is _almost_ the same as rational_prime_index_to_prime, with
        the difference that the rational J ideal is taken into account
        (when there is one, of course). So:
        self.rational_ideal_index_to_prime(0) is 1/leading coeff of f[0]
        self.rational_ideal_index_to_prime(1) is 2
        etc.
        """
        if b == 0 and self.has_merged_J:
            return 1/self.poly.f[0].leading_coefficient()
        else:
            p = self.rational_prime_index_to_prime(b - int(self.has_merged_J))
            return p


    def _number_of_valuations(self, bound=None):
        # This takes the
        # corresponding bound as a parameter.
        # if the given bound is None, then we return the full number of
        # valuations (i.e., the total number of columns in the renumber
        # table)
        if bound is None:
            return len(self._indices_per_side[1])
        index = None
        for parser, side, I, col_index in self._ideals:
            if side != 1:
                continue
            assert type(I) is tuple
            if I[0] >= bound:
                break
            index = col_index
        try:
            assert index is not None
            return index + 1
        except AssertionError as e:
            raise RuntimeError(f"_number_of_valuations({bound}) did not find _any_ ideal below {bound}, which is nuts, really")

    def number_of_fb_valuations(self, BOUNDA_queries):
        """
        returns the number of prime ideals that end up in the
        linear algebra matrix (before filtering).
        """
        return self._number_of_valuations(BOUNDA_queries)

    def number_of_algebraic_columns(self):
        return self._number_of_valuations()

    def number_of_rational_queries(self):
        """
        return the number of rational queries that we need to make in
        order to find all the unknowns on the rational side. When we have
        a merged J ideal, it does not count since this one already
        appears as an unknown on the algebraic side, and thus we have a
        coordinate for it with the linear algebra solution
        """
        return len(self._indices_per_side[0]) - self.has_merged_J


@timing
def run_S_sanity_checks(params, S_list, S_alg_vector, S_rat_vector, u, v, partial_R=False):
    N = params.poly.N
    f0 = params.poly.f[0]
    x = f0.parent().gen()
    m = f0.roots(QQ,multiplicities=False)[0]
    assert m == params.poly.K[0].gen()
    ZN = Integers(N)
    assert ZN(m) == ZN(params.poly.m)

    # note that (a-b*x).resultant(f0) is a-b*m if f0 is monic. Otherwise
    # it's lc(f0)*(a-b*m)

    # The following checks are only valid if we kept track of the
    # rational side, which in fact we'd like to avoid!

    if partial_R:
        R = PartialRenumber1024(params.poly)
        partial_renum_filename = params.dirs['TEMP_OUTPUT_DIR'] + 'partial.renumber.map'
        R.load(partial_renum_filename)
        R.do_sanity_checks()
    else:
        R = params.R

    pi2p = R.rational_prime_index_to_prime

    rat_factored = Factorization([(ZN(pi2p(i)),k)
                                  for i,k in S_rat_vector.dict().items()])

    # Check u/v * S(m)
    rational_smooth_product = rat_factored.prod()

    S_m = ZN(u) / ZN(v)

    # see remark above. S_m isn't exactly S(m)
    # S_m *= prod([ZN((a-b*x).resultant(f0))**k for (a, b), k in S_list.items()])
    S_m *= prod([ZN((a-b*m))**k for (a, b), k in S_list.items()])
    lc_fix = sum(S_list.values())
    S_m *= ZN(f0.leading_coefficient())**lc_fix

    print("S_m", str(S_m))
    print("rational_smooth_product", str(rational_smooth_product))
    print("lc_fix", str(lc_fix))

    assert(rational_smooth_product == S_m or
           rational_smooth_product == -S_m)


    # another way to put it. Maybe this one is a bit expensive because it
    # computes over the integers and not in ZN (but on the other hand
    # everything remains in sparse form so it isn't that bad)
    fac_elems = u.factor()/v.factor()*prod([(a-b*x).resultant(f0).factor()**k
                                            for (a, b), k in S_list.items()])
    fac_primes = Factorization([(pi2p(i),k)
                                for i,k in S_rat_vector.dict().items()])
    ff = fac_elems / fac_primes
    assert ff.prod() == ff.unit()


    # also check what we can check on the algebraic side.

    f1 = params.poly.f[1]
    K = params.poly.K[1]
    alpha = K.gen()
    OK = K.maximal_order()

    # note that we can only check the ideal factorizations!
    S_alpha = OK.fractional_ideal(prod([(a-b*alpha)**k
                                        for (a, b), k in S_list.items()]))

    # FIXME J_valuation_inconsistency. This is the same as the situation
    # we encounted in the Montgomery e-th root computation (see "FIXME
    # J_valuation_inconsistency" there)
    if True and not f1.is_monic():
        J1 = R.column_to_sage_ideal(0, side_hint=1)
        S_alpha *= J1**(lc_fix*2)

        print("J1 norm", str(J1.norm()))
        print("lc_fix", str(lc_fix))

    c2a = lambda c: R.column_to_sage_ideal(c, side_hint=1)
    alg_factored = Factorization([(c2a(i),k)
                                  for i,k in S_alg_vector.dict().items()])

    diff = S_alpha.norm() -  alg_factored.prod().norm()
    ratio = S_alpha.norm() / alg_factored.prod().norm()
    ratio2 = alg_factored.prod().norm() / S_alpha.norm()

    print("diff", diff)
    print("ratio", ratio)
    print("ratio2", ratio2)

    assert S_alpha.norm() == alg_factored.prod().norm()
    assert S_alpha == alg_factored.prod()

@timing
def run_ST_sanity_checks(params, S_list, T_list, ST_alg_vector, S_rat_vector, u, v, partial_R=False):
    N = params.poly.N
    f0 = params.poly.f[0]
    x = f0.parent().gen()
    m = f0.roots(QQ,multiplicities=False)[0]
    assert m == params.poly.K[0].gen()
    ZN = Integers(N)
    assert ZN(m) == ZN(params.poly.m)

    if partial_R:
        desc_R = PartialRenumber1024(params.poly)
        partial_renum_filename = params.dirs['TEMP_OUTPUT_DIR'] + 'partial.renumber.map'
        desc_R.load(partial_renum_filename)
        desc_R.do_sanity_checks()
        params_R = params.R     # has primes up to 2**31
    else:
        desc_R = params.R
        params_R = params.R

    # note that (a-b*x).resultant(f0) is a-b*m if f0 is monic. Otherwise
    # it's lc(f0)*(a-b*m)

    # The following checks are only valid if we kept track of the
    # rational side, which in fact we'd like to avoid!
    pi2p = desc_R.rational_prime_index_to_prime
    rat_factored = Factorization([(ZN(pi2p(i)),k)
                                  for i,k in S_rat_vector.dict().items()])
    rational_smooth_product = rat_factored.prod()

    # Check u/v * S(m)
    S_m = ZN(u) / ZN(v)

    # see remark above. S_m isn't exactly S(m)
    # S_m *= prod([ZN((a-b*x).resultant(f0))**k for (a, b), k in S_list.items()])
    S_m *= prod([ZN((a-b*m))**k for (a, b), k in S_list.items()])
    lc_fix = sum(S_list.values())
    S_m *= ZN(f0.leading_coefficient())**lc_fix

    assert(rational_smooth_product == S_m or
           rational_smooth_product == -S_m)

    # NOTE: We have ignored T_list for now: these elements are not
    # smooth, and we only expect to find their roots via additional
    # queries:
    # [(a-b*x).resultant(f0).factor()**k for (a, b), k in T_list.items()]


    # also check what we can check on the algebraic side. Here, we use
    # T_list as well!

    f1 = params.poly.f[1]
    K = params.poly.K[1]
    alpha = K.gen()
    OK = K.maximal_order()

    # note that we can only check the ideal factorizations!
    S_alpha = OK.fractional_ideal(prod([(a-b*alpha)**k
                                        for (a, b), k in S_list.items()]))
    T_alpha = OK.fractional_ideal(prod([(a-b*alpha)**k
                                        for (a, b), k in T_list.items()]))
    if True and not f1.is_monic():
        # FIXME: This is becoming very annoying because here our
        # correcting factor is different from the one on the rational
        # side...
        J1 = params_R.column_to_sage_ideal(0, side_hint=1)
        S_alpha *= J1**(2*sum(S_list.values()))
        T_alpha *= J1**(2*sum(T_list.values()))

    c2a = lambda c: params_R.column_to_sage_ideal(c, side_hint=1)
    alg_factored = Factorization([(c2a(i),k)
                                  for i,k in ST_alg_vector.dict().items()])

    S_times_T = S_alpha*T_alpha

    diff = S_times_T.norm() -  alg_factored.prod().norm()
    ratio = S_times_T.norm() / alg_factored.prod().norm()
    ratio2 = alg_factored.prod().norm() / S_times_T.norm()

    timeprint("diff", str(diff))
    timeprint("ratio", str(ratio))
    timeprint("ratio2", str(ratio2))

    assert S_times_T == alg_factored.prod()


@timing
def run_sol_sanity_checks(params, sol, ST_alg_vector, indexed_relations_file):
    """
    This re-reads all the algebraic query relations, and checks that the
    valuations at the algebraic prime ideals of the combination match the
    valuations found in ST_alg_vector (which depends on the target)
    """
    e = params.parameters['e']
    Ze = Integers(e)
    sol_prime_ideal_vals = vector(Ze, len(ST_alg_vector))
    rels = indexed_relations_from_file(indexed_relations_file)
    current_row = 0
    while (current_row < len(sol)):
        sol_r = sol[current_row]
        irel = next(rels)
        for ii in irel.indices:
            col = params.R.renumber_to_column(ii)
            # FIXME J_valuation_inconsistency
            if False and col == 0 and not params.poly.f[1].is_monic():
                sol_prime_ideal_vals[col] -= sol_r
            elif col >= 0:
                sol_prime_ideal_vals[col] += sol_r
        current_row += 1

    # print(str(sol_prime_ideal_vals))
    # print(str(ST_alg_vector))

    for i in range(len(ST_alg_vector)):
        assert(-1*sol_prime_ideal_vals[i] % e == ST_alg_vector[i] % e)

@timing
def algebraic_consistency_check(params, S_list, S_alg_vector):
    alpha1 = params.poly.K[1].gen()
    alg_thing = ZZ(1)
    for (a,b),e in S_list.items():
        alg_thing *= (a-b*alpha1).norm()**e

    # print(alg_thing)
    # print(alg_thing.factor())

    biggie = ZZ(1)
    for i,v in enumerate(S_alg_vector):
        if not v:
            continue
        side,q,rho = params.R.side_and_index_to_ideal(1, i)
        biggie *= q**v

    # print(biggie)
    # print(biggie.factor())

    return biggie / alg_thing in [1, -1]

@timing
def truncate_S(params, S_list, S_alg_vector, partial_R=False):

    EXTRELS_FILE = params.files['EXTRELS_FILE']
    EXTRELS_INDEXED = params.files['EXTRELS_INDEXED']

    if partial_R:
        desc_R = PartialRenumber1024(params.poly)
        partial_renum_filename = params.dirs['TEMP_OUTPUT_DIR'] + 'partial.renumber.map'
        desc_R.load(partial_renum_filename)
        desc_R.do_sanity_checks()
        params_R = params.R     # has primes up to 2**31
    else:
        desc_R = params.R
        params_R = params.R

    # num_cols = R.number_of_fb_valuations(params.BOUNDA)
    # for irel in indexed_relations_from_file(EXTRELS_INDEXED):
    #     for idx in irel.indices:

    #print("S_alg_vector", ', '.join([ f"{i}:{S_alg_vector[i]}" for i in S_alg_vector.nonzero_positions()]))
    if params.debug:
        assert algebraic_consistency_check(params, S_list, S_alg_vector)

    starttime = time.time()
    # this is messy, really. We should be able to avoid some of these
    # detours.
    per_q = parse_fb_extension_relations(EXTRELS_FILE, (1, params.BOUNDA_queries, params.BOUNDA))
    endtime = time.time()
    params.timing['parse_fb_extension_relations'] = endtime-starttime
    timeprint('parse_fb_extension_relations took', endtime-starttime)
    starttime = time.time()
    indexed_relations_file = big_convert_to_indexed_relation(per_q.values(),
                                               params,
                                               EXTRELS_INDEXED)
    endtime = time.time()
    params.timing['big_convert_to_indexed_relation'] = endtime-starttime
    timeprint('big_convert_to_indexed_relation', endtime-starttime)

    indexed = { (irel.a,irel.b): irel for irel in indexed_relations_file }

    T_list = defaultdict(int)
    ST_list = copy.copy(S_list)
    ST_alg_vector = copy.copy(S_alg_vector)

    num_cols = params_R.number_of_fb_valuations(params.BOUNDA_queries)
    timeprint("num_cols", str(num_cols))
    timeprint("num_cols_2", str(desc_R.number_of_fb_valuations(params.BOUNDA_queries)))
    #assert num_cols == params_R.number_of_fb_valuations(params.BOUNDA_queries)

    num_missing_ideals = 0

    for i in S_alg_vector.nonzero_positions():
        if i < num_cols:
            continue

        exponent = S_alg_vector[i]

        side, *I = desc_R.side_and_index_to_ideal(1, i)
        if side != 1 or len(I) != 2:
            raise RuntimeError(f"ideal {i} in the renumber table points to {(side, *I)}")
        q, rho = I

        if not partial_R:
            # _indices_per_side should be defined as usual
            # Not sure why we use that instead of column_to_renumber.
            timeprint(f"now killing index {i}->0x{params.R._indices_per_side[1][i]:x} == {(side, q, rho)}")
        else:
            timeprint(f"now killing index {i}->0x{desc_R.column_to_renumber(i)} == {(side, q, rho)}")

        if (side, q, rho) not in per_q:

            num_missing_ideals += 1
            print(f"Missing prime ideal {side} {q} {rho}")
            print(f"So far missing: {num_missing_ideals}")

            continue

            #raise RuntimeError(f"Uh oh! The prime ideal of index 0x{i:x}"
            #                   " in the extension factor base,"
            #                   f" a.k.a.  {(side,q,rho)}"
            #                   " didn't show up in any relations to"
            #                   " the query factor base."
            #                   " This bug is supposedly gone."
            #                   f" See {EXTRELS_INDEXED}")

        las_rel = per_q[side, q, rho]
        irel = indexed[las_rel.a, las_rel.b]

        # We divide both sides by (a-b*alpha) = q * smooth_prod
        ST_list[las_rel.a, las_rel.b] -= exponent
        T_list[las_rel.a, las_rel.b] -= exponent

        found_me = False
        for f in irel.indices:
            col_desc = desc_R.renumber_to_column(f)
            try:
                col_params = params_R.renumber_to_column(f)
            except IndexError:
                col_params = -1

            if col_desc < 0 and col_params < 0:
                # rational
                continue
            elif col_desc == i:
                found_me = True
                # our extension prime
                assert (col_params < 0) or (col_desc == col_params)
                assert col_desc >= num_cols
                ST_alg_vector[col_desc] -= exponent
            elif col_params >= 0:
                assert (col_desc == col_params) or (col_desc < 0)
                # prime below LPB1_queries
                ST_alg_vector[col_params] -= exponent
        assert found_me

    timeprint(f"At least made it through main loop of truncate_S")
    timeprint(f"Beyond this we will need all the extension relations.")

    assert num_missing_ideals == 0
    assert ST_alg_vector[num_cols:].nonzero_positions() == []

    ST_alg_vector = vector(ST_alg_vector[:num_cols])

    if params.debug:
        assert algebraic_consistency_check(params, ST_list, ST_alg_vector)

    # After truncation, we should only be within BOUNDA_queries, so only
    # need to look for unlinked ideals among the FILT_UNLINKED_IDEALS.
    ul_ideal_cols = set()
    if os.path.exists(params.files['FILT_UNLINKED_IDEALS']):
        with open(params.files['FILT_UNLINKED_IDEALS'], "r") as ul_file:
            for renum_index in ul_file.readlines():
                ul_col = params_R.renumber_to_column(int(renum_index))
                if ul_col >= 0:
                    ul_ideal_cols.add(ul_col)

    for ul in ul_ideal_cols:
        if ST_alg_vector[ul] != 0:
            ul_renum = params_R.column_to_renumber(ul)
            ul_renum_hex = hex(ul_renum)
            raise RuntimeError("truncate_S:"
                                f" The relation involves ideal number {ul_renum}"
                                f" or {ul_renum_hex}"
                                " (in renumber coordinates),"
                                f" which is {ul}"
                                " (in column coordinates),"
                                " however it is not linked to the others."
                                " Too bad, really.")
            exit(0)

    return T_list, ST_list, ST_alg_vector

@timing
def get_RSU_m(S_list, sol, row_to_aquery, f, e, N, m):
    # R**e = S*U
    ZN = Integers(N)
    abm_list = []
    extra_abks = []
    S_m = ZN(1)
    U_m = ZN(1)

    # Note that the crt_ethroot script uses the (a+b*m) convention, not (a-b*m).

    for (a, b), exponent in S_list.items():
        S_m = ZN(S_m * power_mod(a-b*m, exponent, N))
        if exponent > 0:
            abm_list.append((a, -1*b, exponent))
        elif exponent < 0:
            # Find k such that 0 <= m + k*e
            k = 0
            while(0 > exponent + k*e):
                k += 1
            extra_abks.append((a, b, k))
            abm_list.append((a, -1*b, k*e + exponent))

    for i in range(len(sol)):
        a,b = row_to_aquery[i]
        exponent = Integer(sol[i])
        U_m = ZN(U_m * power_mod(a-b*m, exponent, N))
        if exponent > 0:
            abm_list.append((a, -1*b, exponent))
        elif exponent < 0:
            k = 0
            while(0 > exponent + k*e):
                k += 1
            extra_abks.append((a, b, k))
            abm_list.append((a, -1*b, k*e + exponent))

    for a,b,m in abm_list:
        assert(m >= 0)

    ethroot_gen = ethroot(abm_list, f, e, N)
    ### For debugging, uncomment the following lines to write the eth-root step's input to a file
    #with open('/tmp/ethroot_data','w') as file:
    #    file.write(f"f: {f}\n")
    #    file.write(f"e: {e}\n")
    #    file.write(f"N: {N}\n")
    #    file.write(f"len(abm_list): {len(abm_list)}\n")
    #    file.write(f"abm_list: {abm_list}\n")


    if False:
        assert U_m == prod([ZN(a-b*m)^k
                            for (a,b),k in [
                                (linalg_output.row_to_aquery[i],s)
                                for i,s in enumerate(linalg_output.sol)]])

        assert S_m == prod([ZN(a-b*m)^k
                            for (a,b),k in linalg_output.ST_list.items()])

        TT = [(linalg_output.row_to_aquery[i],s) for i,s in enumerate(linalg_output.sol)]
        TT += list(linalg_output.ST_list.items())
        assert S_m*U_m == prod([ZN(a-b*m)^k for (a,b),k in TT])
        assert S_m*U_m == prod([ZN(a+minus_b*m)^k for a,minus_b,k in abm_list])/prod([ZN(a-b*m)^k for a,b,k in extra_abks])^e



    return S_m, U_m, ethroot_gen, extra_abks, abm_list


def extract_decimal_ab_from_file(outfilename, infilename):
    with open(outfilename, "w") as fw:
        timeprint("Reading", infilename)
        for rel in indexed_relations_from_file(infilename):
            print(f"{rel.a},{rel.b}", file=fw)


def fast_union(s1, s2):
    s1.update(s2)
    return s1

def get_dict_keyset(d):
    return set(d.keys())

class SchirokauerMapsAppender(abc.ABC):
    def __init__(self, params):
        self.params = params
        self.N = params.poly.N
        self.e = int(params.parameters['e'])
        self.Ze = IntegerModRing(self.e)

    @abc.abstractmethod
    def maps_from_more_ab(self, ab_pairs):
        pass

    @abc.abstractmethod
    def compute_sm_block_for_matrix(self, mat):
        pass


def cado_sage_maps_from_ab_parallel(inputs):
    (a, b, e, alpha_pkl, sm_maps_pkl) = inputs
    alpha = SageUnpickler.loads(alpha_pkl)
    sm_maps = SageUnpickler.loads(sm_maps_pkl)
    return vector(Integers(e),
                  sum([s(a-b*alpha).list() for s in sm_maps], [])
                  )


class SchirokauerMapsAppender_cado_sage(SchirokauerMapsAppender):
    def _maps_from_ab(self, a, b):
        return vector(self.Ze,
                      sum([s(a-b*self.alpha).list() for s in self.sm_maps], [])
                      )

    def __init__(self, params):
        super().__init__(params)
        self.K = params.poly.K[1]
        self.alpha = self.K.gen()
        self.Kw = params.poly.nt[1]
        self.sm_maps = self.Kw.schirokauer_maps(self.e)
        self.nthr = params.parameters.get('sm_append.thr', multiprocessing.cpu_count())
        self.e = params.parameters['e']

    @timing
    def compute_sm_block_for_matrix(self, mat):
        """
        Form the block of Schirokauer map values that we're going to paste
        into the matrix.
        """

        # This uses the mat.ab_per_row() interface, which can actually
        # return stuff that is more complicated than the 1-1 mapping that
        # we have in a LinearAlgebraMatrix_IndexedFileOnly object.

        # precompute sm maps for all a,b, in case some a,b appear several
        # times.

        timeprint("Starting compute_sm_block_for_matrix")

        try:
            import mr4mp
            all_ab = mr4mp.pool().mapreduce(get_dict_keyset, fast_union, list(mat.ab_per_row(self.params)))
        except ModuleNotFoundError:
            timeprint("WARNING: Install Python package 'mr4mp' for faster mapreduce operation")
            all_ab = functools.reduce(set.union, [R.keys() for R in mat.ab_per_row(self.params)], set())

        #all_ab = functools.reduce(set.union, [R.keys() for R in mat.ab_per_row(self.params)], set())

        timeprint("Finished making all_ab set")

        pool = multiprocessing.Pool(processes=self.nthr)
        all_inputs = []

        DO_PARALLEL_CADO_SAGE = False
        DO_SLURM_CADO_SAGE = True

        if DO_PARALLEL_CADO_SAGE:
            for (a,b) in all_ab:
                all_inputs.append(
                    (a, b, self.e, SagePickler.dumps(self.alpha), SagePickler.dumps(self.sm_maps))
                )
            res = pool.map(cado_sage_maps_from_ab_parallel, all_inputs)

            D = { (all_inputs[i][0], all_inputs[i][1]) : res[i] for i in range(len(all_inputs)) }

        elif DO_SLURM_CADO_SAGE:
            timeprint("starting slurm allocations")

            all_ab_list = list(all_ab)
            job_line_to_ab = dict()

            numjobs = int(self.params.parameters['slurm.numjobs'])
            ab_per_job = ceil( len(all_ab) / numjobs)

            timeprint("ab_per_job", str(ab_per_job))

            processes = []

            for j in range(numjobs):
                infile = f"{self.params.dirs['TEMP_OUTPUT_DIR']}sm/abs.job{j}.in"
                outfile = f"{self.params.dirs['TEMP_OUTPUT_DIR']}sm/abs.job{j}.out"
                jobname = f"sm-768-{j}"

                linenum = 0

                with open(infile, "w") as abfile:
                    for i in range( j*ab_per_job, min( (j+1)*ab_per_job, len(all_ab) ) ):
                        (a,b) = all_ab_list[i]
                        job_line_to_ab[(j,linenum)] = (a,b)
                        abfile.write(f"{a},{b}\n")
                        linenum += 1

                command_list = [
                    self.params.files['SAGE'],
                    "sm_cado_sage_helper.py",
                    "--infile", infile,
                    "--outfile", outfile,
                    "--e", str(self.params.parameters['e']),
                    "--poly", self.params.files['POLYFILE']
                ]
                processes.append(slurmit(self.params, " ".join(command_list), jobname, j))

            finished_processes, cputime_slurm = slurm_wait(processes)
            overall_cputime.add(cputime_slurm)

            D = dict()

            for j in range(numjobs):
                outfile = f"{self.params.dirs['TEMP_OUTPUT_DIR']}sm/abs.job{j}.out"

                with open(outfile, "r") as smfile:
                    i = 0
                    for line in smfile.readlines():
                        line = line.strip()
                        a,b = job_line_to_ab[(j,i)]
                        D[(a,b)] = vector(self.Ze, line.split(","))
                        i += 1

        else:
            D = { (a,b):self._maps_from_ab(a, b) for (a,b) in all_ab }

        timeprint("Finished constructing cado_sage dictionary D. Now making the matrix to return.")

        return matrix(self.Ze, [
            sum([D[a, b] * k for (a,b),k in R.items()])
            for R in mat.ab_per_row(self.params)])

    @timing
    def maps_from_more_ab(self, ab_pairs):
        """
        We typically use this function with the (a,b) pairs that we
        encountered in the descent: those are typically not found in the
        file AQRELS_INDEXED
        """
        return { x: self._maps_from_ab(*x) for x in ab_pairs }


class SchirokauerMapsAppender_sm_simple(SchirokauerMapsAppender):
    def __init__(self, params):
        super().__init__(params)


    def _call_sm_simple(self, inputname):

        outputname = inputname + ".sm"

        CadoNFS("filter/sm_simple",
                "-ell", self.e,
                "-inp", 'IN',
                "-out", 'OUT',
                "-poly", 'POLY',
                outputs={'OUT': outputname},
                inputs={
                    'IN': inputname,
                    'POLY': self.params.files['POLYFILE'],
                    }
                )
        return [vector(self.Ze, line.strip().split()) for line in
                open(outputname).readlines()]

    @timing
    def compute_sm_block_for_matrix(self, mat):
        """
        Form the block of Schirokauer map values that we're going to paste
        into the matrix.
        """

        # Each row of the matrix is for a linear combination of (a,b) pairs (because filtering)
        # but filter/sm_simple wants as input one (a,b) pair per line.

        # Make a list of all (a,b) used in the filtered matrix, removing duplicates.
        # Note that mat.ab_per_row() returns an iterator of {(a,b): coeff} dicts, one per row of mat.
        all_ab = functools.reduce(set.union, [R.keys() for R in mat.ab_per_row(self.params)], set())
        # Write to a file, and remember which (a,b) are on which line
        base = mat.indexed_relations_filename(self.params)
        aqrels_decimal_filename = base + ".filtered.deci"
        ab_to_line = dict()
        with open(aqrels_decimal_filename, "w") as fw:
            for (lineno, (a,b)) in enumerate(all_ab):
                ab_to_line[(a,b)] = lineno
                print(f"{a},{b}", file=fw)
        block_unfiltered = matrix(self.Ze, self._call_sm_simple(aqrels_decimal_filename))
        # Now the characters for (a,b) are row (ab_to_line[(a,b)]) of block_unfiltered

        return matrix(self.Ze, [
            sum([block_unfiltered[ab_to_line[(a, b)]] * k for (a,b),k in R.items()])
            for R in mat.ab_per_row(self.params)
        ])

    @timing
    def maps_from_more_ab(self, ab_pairs):
        tmp = self.params.dirs['TEMP_OUTPUT_DIR']
        with tempfile.NamedTemporaryFile(dir=tmp, mode="w") as f:
            for a,b in ab_pairs:
                print(f"{a},{b}", file=f.file)

            f.flush()
            sms = self._call_sm_simple(f.name)
            return { (a,b):s for (a,b),s in zip(ab_pairs, sms) }

def parse_sm_append_lines_to_vector(args):
    Ze, lines = args
    vectors = []
    for line_bytes in lines:
        line = line_bytes.decode()
        if not line.startswith('#'):
            vectors.append(vector(Ze, line.strip().split(':')[1].split(',')))
    return vectors

class SchirokauerMapsAppender_sm_append(SchirokauerMapsAppender):
    def __init__(self, params):
        super().__init__(params)


    def _call_sm_append(self, inputname):
        outputname = inputname + ".withsm"

        nsm = []
        if nc := self.params.parameters.get('num_character_columns'):
            nsm = ["-nsm", f"0,{nc}"]

        sm_append_mpi_exec = []
        sm_append_add_args = []
        if self.params.mpi:
            timeprint(f"Run sm_append with MPI enabled.")
            sm_append_mpi_exec = [self.params.mpi_mpirun_bin] + self.params.mpi_extra_args.split(" ")
        else:
            sm_append_thr = self.params.parameters.get('sm_append.thr',
                    multiprocessing.cpu_count())
            timeprint(f"Run sm_append with {sm_append_thr} threads.")
            sm_append_add_args = [
                "-t", sm_append_thr
            ]

        CadoNFS("filter/sm_append",
                "-ell", self.e,
                "-in", 'IN',
                "-out", 'OUT',
                "-poly", 'POLY',
                "-b", "256",
                *sm_append_add_args,
                *nsm,
                mpi_exec=sm_append_mpi_exec,
                outputs={ 'OUT': outputname },
                inputs={
                    'IN': inputname,
                    'POLY': self.params.files['POLYFILE'],
                    })

        if self.params.nthreads > 1:
            timeprint(f"loading {outputname} into memory")
            with open(outputname, mode="r", encoding="ascii") as output_fp:
                with mmap.mmap(output_fp.fileno(), length=0, access=mmap.ACCESS_READ) as output_mmap:
                    timeprint(f"Starting to parse {outputname} with {self.params.nthreads} processes")
                    output_lines = []
                    l = 0
                    while True:
                        line = output_mmap.readline()
                        l += 1
                        if line == b"":
                            break
                        output_lines.append(line)
                    timeprint(f"Finished splitting {outputname} into lines.")

                vectors = []
                with ProcessPoolExecutor(max_workers=self.params.nthreads) as executor:
                    batch_size = math.ceil(l / self.params.nthreads)
                    batch_args = [(self.Ze, output_lines[i * batch_size : (i+1) * batch_size]) for i in range(self.params.nthreads)]
                    for vectors_batch in tqdm.tqdm(executor.map(parse_sm_append_lines_to_vector, batch_args)):
                        vectors += vectors_batch
                timeprint(f"Done parsing")
                return vectors
        else:
            return [vector(self.Ze, line.strip().split(':')[1].split(','))
                    for line in open(outputname).readlines()
                    if not line.startswith('#')]

    @timing
    def compute_sm_block_for_matrix(self, mat):
        """
        Form the block of Schirokauer map values that we're going to paste
        into the matrix.
        """

        # Each row of the matrix is for a linear combination of (a,b) pairs (because filtering)
        # but filter/sm_append wants as input one (a,b) pair per line.

        timeprint("Making list of all (a,b)'s used")

        # Make a list of all (a,b) used in the filtered matrix, removing duplicates.
        # Note that mat.ab_per_row() returns an iterator of {(a,b): coeff} dicts, one per row of mat.
        try:
            import mr4mp
            all_ab = mr4mp.pool().mapreduce(get_dict_keyset, fast_union, list(mat.ab_per_row(self.params)))
        except ModuleNotFoundError:
            timeprint("WARNING: Install Python package 'mr4mp' for faster mapreduce operation")
            all_ab = functools.reduce(set.union, [R.keys() for R in mat.ab_per_row(self.params)], set())

        # Write to a file, and remember which (a,b) are on which line
        base = mat.indexed_relations_filename(self.params)
        abpairs_for_sm_filename = base + ".forsm"
        timeprint(f"Reading {base} into memory and"
                  f" copying to {abpairs_for_sm_filename}")
        ab_to_line = dict()
        with open(abpairs_for_sm_filename, "w") as fw:
            for (lineno, (a,b)) in enumerate(all_ab):
                ab_to_line[(a,b)] = lineno

                #if Integer(a).nbits() >= 64 or Integer(b).nbits() >= 64:
                #    raise RuntimeError(f"Large a,b passed to matrix sm_append: {a}, {b}")

                print(f"{a:x},{b:x}", file=fw)

        timeprint("Calling sm_append")
        block_unfiltered = matrix(self.Ze, self._call_sm_append(abpairs_for_sm_filename))
        # Now the characters for (a,b) are row (ab_to_line[(a,b)]) of block_unfiltered

        timeprint("Recovering sm block with a matrix multiplication")
        return matrix(self.Ze, [
            sum([block_unfiltered[ab_to_line[(a, b)]] * k for (a,b),k in R.items()])
            for R in mat.ab_per_row(self.params)
        ])

    @timing
    def maps_from_more_ab(self, ab_pairs):
        """
        this one is really hacky. sm_append is not meant for that, even
        though we promised to use it here. Anyway, here goes.
        """
        tmp = self.params.dirs['TEMP_OUTPUT_DIR']
        with tempfile.NamedTemporaryFile(dir=tmp, mode="w") as f:
            for a,b in ab_pairs:
                #if Integer(a).nbits() >= 64 or Integer(b).nbits() >= 64:
                #    raise RuntimeError(f"Large a,b passed to sm_append: {a}, {b}")

                print(f"{a:x},{b:x}", file=f)

            f.file.flush()
            sms = self._call_sm_append(f.name)
            return { (a,b):s for (a,b),s in zip(ab_pairs, sms) }

@timing
def make_linalg_system(params, matrix, ST_list, ST_alg_vector):
    # We have several options for the computation of the Schirokauer maps
    # block that goes in the matrix. Here, the cado_sage method has
    # higher priority because it works more generally.

    timeprint("Starting make_linalg_system")

    compute_sm_method = {
        'sm_simple': SchirokauerMapsAppender_sm_simple,
        'sm_append': SchirokauerMapsAppender_sm_append,
        'cado_sage': SchirokauerMapsAppender_cado_sage,
    }

    # starttime = time.time()
    # sm_simple_init = SchirokauerMapsAppender_sm_simple(params)
    # endtime = time.time()
    # print("sm_simple_init took",endtime-starttime)
    # params.timing["sm_simple_init"] = endtime-starttime
    # starttime = time.time()
    # sm_append_init = SchirokauerMapsAppender_sm_append(params)
    # endtime = time.time()
    # print("sm_append_init took",endtime-starttime)
    # params.timing["sm_append_init"] = endtime-starttime
    # starttime = time.time()
    # cado_sage_init = SchirokauerMapsAppender_cado_sage(params)
    # endtime = time.time()
    # print("cado_sage_init took",endtime-starttime)
    # params.timing["cado_sage_init"] = endtime-starttime

    # compute_sm_flags = {
    #     'sm_simple': sm_simple_init,
    #     'sm_append': sm_append_init,
    #     'cado_sage': cado_sage_init
    # }

    timeprint("Starting sm_init")
    starttime = time.time()
    compute_sm_flags = {params.sm_alg: compute_sm_method[params.sm_alg](params)}
    endtime = time.time()
    timeprint("Finished sm_init")
    params.timing["sm_init"] = endtime-starttime

    MC_saved_file = params.dirs['TEMP_OUTPUT_DIR']+"MC.sobj"
    S_block_saved_file = params.dirs['TEMP_OUTPUT_DIR']+"S_block.sobj"
    if os.path.exists(MC_saved_file) and os.path.exists(S_block_saved_file) and (not params.overwrite_MC):
        # The matrix already exists,
        # and we do not want to overwrite it.
        # The target vector must still be made.
        write_matrix = False
    else:
        write_matrix = True

    variants = dict()

    #ST_list_ab_small = [
    #    (a,b) for (a,b) in ST_list.keys() if Integer(a).nbits() < 64 and Integer(b).nbits() < 64
    #]
    #ST_list_ab_big = [
    #    (a,b) for (a,b) in ST_list.keys() #if Integer(a).nbits() >= 64 or Integer(b).nbits() >= 64
    #]

    timeprint("Starting compute_sm_block")
    for method, APP in compute_sm_flags.items():
        if write_matrix:
            S_block = APP.compute_sm_block_for_matrix(matrix)
        else:
            major_message("Loading existing S_block from file.")
            S_block = fast_persistent_load(S_block_saved_file)

        D = APP.maps_from_more_ab(ST_list.keys())

        # OLD DEBUGGING:
        #D = APP.maps_from_more_ab(ST_list_ab_small)
        #D = dict()
        #timeprint("Finished the ST_list_ab_small")
        #timeprint("Starting the ST_list_ab_big")

        #K = params.poly.K[1]
        #alpha = K.gen()
        #Kw = params.poly.nt[1]
        #sm_maps = Kw.schirokauer_maps(params.parameters['e'])

        #print("Len of sm_maps", str(len(sm_maps)))

        #for (big_a, big_b) in ST_list_ab_big:
        #    print("Handling a big a,b...")
        #    vv = vector(Integers(params.parameters['e']),
        #                sum([s(big_a-big_b*alpha).list() for s in sm_maps], [])
        #                )
        #    D[(big_a,big_b)] = vv

        starttime = time.time()
        C_block = sum([m * D[a,b] for (a,b),m in ST_list.items()])
        endtime = time.time()
        params.timing['C_block'] = endtime-starttime
        timeprint("C_block took",endtime-starttime)
        variants[method] = (S_block, C_block)
    timeprint("Finished compute_sm_block")

    for method, (S_block, C_block) in variants.items():
        timeprint(f"S block computed by {method}:",
                  f"{S_block.nrows()}x{S_block.ncols()}")
        timeprint(f"C block computed by {method}:",
                  f"{len(C_block)}")

    starttime = time.time()
    all_variants = list(variants.items())
    for i in range(len(all_variants)):
        k0, (S0, C0) = all_variants[i]
        for j in range(i + 1,len(all_variants)):
            k, (S, C) = all_variants[j]
            if S != S0:
                tail = ""
                c = min(S.ncols(), S0.ncols())
                if S[:,:c] == S0[:,:c]:
                    tail = f' (but they agree on the first {c} coordinates)'
                warning_message(f"matrix S_block computed by {k} differs from",
                                f"the one computed by {k0}" + tail)
            if C != C0:
                tail = ""
                c = min(len(C), len(C0))
                if C[:c] == C0[:c]:
                    tail = f' (but they agree on the first {c} coordinates)'
                warning_message(f"vector C_block computed by {k} differs from",
                                f"the one computed by {k0}" + tail)

#    S_block, C_block =  variants['cado_sage']
    S_block, C_block = variants[params.sm_alg]
#    endtime = time.time()
#    params.timing['variants'] = endtime-starttime
#    print("variants took",endtime-starttime)

    timeprint(f"M: {matrix.nrows} x {matrix.ncols}")
    timeprint(f"S_block: {S_block.nrows()} x {S_block.ncols()}")

    timeprint(type(matrix.matrix()))
    timeprint(type(S_block))
    timeprint(type(C_block))
    timeprint(repr(matrix.matrix()))
    timeprint(repr(S_block))
    timeprint(repr(C_block))

    if write_matrix:
        starttime = time.time()

        S_block_ZZ = S_block.change_ring(ZZ)
        timeprint("type of S_block_ZZ", str(type(S_block_ZZ)))
        timeprint(repr(S_block_ZZ))

        if params.scipy_matrix:
            # No need to make MC at this moment.
            # We will get it later as M || S_block_ZZ
            from matrix_helpers import MC_wrapper
            MC = MC_wrapper(matrix.matrix().nrows(), matrix.matrix().ncols() + S_block_ZZ.ncols())
        else:
            MC = block_matrix(1, 2, [ matrix.matrix(), S_block_ZZ ], sparse=True)

        endtime = time.time()
        params.timing['block_matrix'] = endtime-starttime
        timeprint("block_matrix took", endtime-starttime)
        timeprint("type of MC", str(type(MC)))
    else:
        major_message("Loading existing MC from file.")
        MC = fast_persistent_load(MC_saved_file)

    timeprint(f"ST_alg_vector: {len(ST_alg_vector)}")

    starttime = time.time()
    assert ST_alg_vector.is_sparse()
    SC_vector = vector(matrix.base_ring(),
                       len(ST_alg_vector) + len(C_block),
                       ST_alg_vector.dict())
    for j,c in enumerate(C_block):
        SC_vector[len(ST_alg_vector) + j] = c
    endtime = time.time()
    params.timing['SC_vector'] = endtime-starttime
    timeprint("SC_vector took", endtime-starttime)

    timeprint(f"SC_vector: {len(SC_vector)}")

    return MC, SC_vector, S_block, C_block

@timing
def build_killer_rels_dict(params, MM):
    num_og_prime_ideals = params.R.number_of_fb_valuations(params.BOUNDA_queries)

    killers = { i:None for i in MM.column_shrink_map().keys() }
    indexed_relations_file = params.files['AQRELS_FILE'] + ".indexed"
    purged_file            = indexed_relations_file + ".purged"
    relsdel_file           = indexed_relations_file + ".relsdel"
    unlinked = {i:None for i in range(num_og_prime_ideals) if (not params.R.column_to_renumber(i) in killers)}

    rels  = {(i.a, i.b): ([ c for c in i.indices if (not c in killers) and (params.R.renumber_to_column(c)>0) ],i)
             for i in indexed_relations_from_file(relsdel_file) }
    rels |= {(i.a, i.b): ([ c for c in i.indices if (not c in killers) and (params.R.renumber_to_column(c)>0) ],i)
             for i in indexed_relations_from_file(purged_file) }

    still_cancelling = True
    num_iters = 0
    order_of_cancellation = []

    while still_cancelling:
        timeprint(f"Still cancelling (iteration {num_iters})...")
        timeprint(f"unknowns:  {num_og_prime_ideals - len(killers)}")
        timeprint(f"num rels:  {len(rels)}")
        still_cancelling = False
        num_iters += 1

        drops = []
        for ab,(r,rr) in rels.items():
            # get the indices in o for which we don't have a killer rel yet.
            o = [c for c in r if c not in killers]
            if len(o) == 1:
                i = o[0]
                # keep the killers[] dict in renumber coordinates for the
                # time being, since it's what we use in the loop.
                killers[i] = rr
                still_cancelling = True
                o = []
                order_of_cancellation.append(params.R.renumber_to_column(i))
                del unlinked[params.R.renumber_to_column(i)]
            if not o:
                drops.append(ab)
                still_cancelling = True
            elif len(o) < len(r):
                # we don't remove this relation just yet, but its
                # collection of outstanding ideals got shortened.
                rels[ab] = (o,rr)
                still_cancelling = True
        for ab in drops:
            del rels[ab]

        timeprint("Made progress?", still_cancelling)

    killers = {str(params.R.renumber_to_column(i)):str(k)
               for i,k in killers.items()
               if k is not None}
    order_of_cancellation = list(reversed(order_of_cancellation))

    fast_json_dump(killers, params.files['KILLER_RELS_DICT'])

    fast_json_dump(order_of_cancellation, params.files['KILLER_RELS_ORDER'])

    num_unlinked = 0
    with open(params.files['FILT_UNLINKED_IDEALS'], "w") as ul_file:
        for col in unlinked.keys():
            renum_index = params.R.column_to_renumber(col)
            ul_file.write(str(renum_index) + "\n")
            num_unlinked += 1

    with open(params.files['FILT_UNLINKED_TODOS'], "w") as todo_file:
        for col in unlinked.keys():
            ideal = params.R.side_and_index_to_ideal(1, col)
            side, q, rho = ideal
            if "alpha" in str(rho):
                major_message("Warning: a nonlinear todo q.")
            else:
                print(f"0 {q} {rho}", file=todo_file)

    assert(num_unlinked == num_og_prime_ideals - len(killers) - len(MM.column_shrink_map()))
    print(f"Unknowns at the end of build_killer_rels: {num_unlinked}")

    pct_unknown = round(100.0 * num_unlinked / num_og_prime_ideals, 10)
    timeprint(f"That means, of ideals up to BOUNDA_queries, the percent unlinked is: {pct_unknown}")

    return killers, order_of_cancellation


@timing
def update_collection_of_unlinked_ideals(params, which_ones):
    # Make sure the UNLINKED_IDEALS file accurately reflects which ideals are
    # NOT covered in relations.
    # Run this if we are doing extra algebraic sieving. It's also useful for coming
    # up with a sieving todo list.
    # For our purposes in this function unlinked ideals come from call_todo_sieving.
    # If we are running this, build_killer_rels and filtering will happen again anyway.
    # Yeah, this could be slow, so we shouldn't run it all the time.

    assert (which_ones in ['alg', 'ext'])

    R = params.R

    # TEMPORARY, just testing something.
    print("R.rational_prime_to_prime_index(2)", R.rational_prime_to_prime_index(2))
    print("R.rational_prime_to_prime_index(3)", R.rational_prime_to_prime_index(3))
    print("R.rational_prime_to_prime_index(5)", R.rational_prime_to_prime_index(5))
    print("R.rational_prime_to_prime_index(7)", R.rational_prime_to_prime_index(7))

    print("R.renumber_to_rational_prime_index(9)", R.renumber_to_rational_prime_index(9))
    print("R.renumber_to_rational_prime_index(11)", R.renumber_to_rational_prime_index(11))
    print("R.renumber_to_rational_prime_index(14)", R.renumber_to_rational_prime_index(14))

    for jj in range(0, 17):
        print(f"R.renumber_to_column({jj})", R.renumber_to_column(jj))

    if which_ones == 'alg':
        basefile = params.files['AQRELS_FILE']
        outstanding_renumber_indices = set(
            [R.column_to_renumber(i)
                for i in range(R.number_of_algebraic_columns())
                if i == 0
                or R.side_and_index_to_ideal(1, i)[1] < params.BOUNDA_queries
            ])
        ul_file = params.files['UNLINKED_IDEALS']
        todo_file = params.files['UNLINKED_IDEALS_TODOS']
    else:
        basefile = params.files['EXTRELS_FILE']
        outstanding_renumber_indices = set(
            [R.column_to_renumber(i)
                for i in range(R.number_of_algebraic_columns())
                if i > 0
                and R.side_and_index_to_ideal(1, i)[1] >= params.BOUNDA_queries
            ])
        ul_file = params.files['EXT_UNLINKED_IDEALS']
        todo_file = params.files['EXT_UNLINKED_TODOS']

    num_total = len(outstanding_renumber_indices)

    for r in convert_to_indexed_relation(las_relations_from_file(basefile), params):
        for ii in r.indices:
            outstanding_renumber_indices -= {ii}

    n = 0

    # we do overwrite an existing UNLINKED_IDEALS file
    with open(ul_file, "w") as f:
        for ii in outstanding_renumber_indices:
            f.write(str(ii) + "\n")
            n += 1

    # also create a TODO file, in the {0 q rho} format
    with open(todo_file, "w") as f:
        for ii in outstanding_renumber_indices:
            ideal = R.side_and_index_to_ideal(1, R.renumber_to_column(ii))
            side, q, rho = ideal
            if "alpha" in str(rho):
                major_message("Warning: a nonlinear todo q.")
            else:
                print(f"0 {q} {rho}", file=f)

    timeprint(f"Wrote {n} unlinked ideals to {ul_file}.")
    pct_unlinked = round(100.0 * n / num_total, 10)
    timeprint(f"The percent of unlinked {which_ones} ideals is: {pct_unlinked}")


@timing
def apply_filtering_to_target(params, M, MM, MC, SC_vector, killer_relations, order_of_cancellation):
    num_character_cols = MC.ncols()-M.ncols()
    len_tgt_full = len(SC_vector) - num_character_cols
    tgt_full = [SC_vector[i] for i in range(len_tgt_full)]

    e = params.parameters['e']
    #Kw = params.poly.nt[1]
    #sm_maps = Kw.schirokauer_maps(e)
    #K = params.poly.K[1]
    #alpha = K.gen()

    sm_computer = {
        'sm_simple': SchirokauerMapsAppender_sm_simple,
        'sm_append': SchirokauerMapsAppender_sm_append,
        'cado_sage': SchirokauerMapsAppender_cado_sage,
    }[params.sm_alg](params)

    COMPUTE_SMS_ONE_AT_A_TIME = False
    # If we compute the characters for killer relations all at once at the beginning,
    # then we only invoke the (e.g.) sm_simple binary once, rather than once per killer rel.
    # This is likely a little faster, and it logs way less unnecessary output.
    if not COMPUTE_SMS_ONE_AT_A_TIME:
        killer_rel_a = (
            int(killer_relations[str(i)].split(":")[0].split(",")[0],16)
            for i in order_of_cancellation
        )
        killer_rel_b = (
            int(killer_relations[str(i)].split(":")[0].split(",")[1],16)
            for i in order_of_cancellation
        )
        killer_rel_ab_pairs = list(zip(killer_rel_a, killer_rel_b))
        assert len(killer_rel_ab_pairs) == len(order_of_cancellation)
        _a, _b = killer_rel_ab_pairs[0]

        maps_ab_start = time.time()
        killer_rel_ab_sm_dict = sm_computer.maps_from_more_ab(killer_rel_ab_pairs)
        params.timing["killer relations maps_from_more_ab"] = time.time() - maps_ab_start

    updated_characters = vector(Integers(e), num_character_cols)
    # Start with the target characters, then iteratively update in agreement
    # with the killer relations.
    for i in range(num_character_cols):
        updated_characters[i] = SC_vector[len_tgt_full + i]
    assert(len_tgt_full + num_character_cols == len(SC_vector))

    saved_kr_rels = []

    for col_i in order_of_cancellation:
        if tgt_full[col_i] == 0:
            # No need to do anything here
            continue
        else:
            coeff = tgt_full[col_i]
            killer_rel = killer_relations[str(col_i)]
            saved_kr_rels.append((coeff, killer_rel))
            # Add (-coeff)*killer_rel to tgt_full
            # We only chose killer_rel if there was one factor of ideal_i
            ideals = killer_rel.split(":")[1].split(",")
            for j in ideals:
                col_j = params.R.renumber_to_column(int(j, 16))
                tgt_full[col_j] = tgt_full[col_j] - coeff
            # Also need to update the characters in the SC vector
            a = int(killer_rel.split(":")[0].split(",")[0], 16)
            b = int(killer_rel.split(":")[0].split(",")[1], 16)
            #ab_sm_vector = vector(Integers(e), sum([s(a-b*alpha).list() for s in sm_maps], []))
            if COMPUTE_SMS_ONE_AT_A_TIME:
                ab_sm_vector = sm_computer.maps_from_more_ab([(a,b)])[(a,b)]
            else:
                ab_sm_vector = killer_rel_ab_sm_dict[(a,b)]
            assert(len(ab_sm_vector) == num_character_cols)
            for jj in range(num_character_cols):
                updated_characters[jj] = updated_characters[jj] - coeff*ab_sm_vector[jj]

    SC_vector_filtered = vector(MM.base_ring(), MC.ncols())

    for i in range(len(tgt_full)):
        if tgt_full[i] == 0:
            continue
        ideal_i = params.R.column_to_renumber(i)
        if ideal_i in MM.column_shrink_map().keys():
            # The ideal has survived filtering
            filtered_col = MM.column_shrink_map()[ideal_i]
            SC_vector_filtered[filtered_col] = tgt_full[i]
        elif ideal_i not in MM.column_shrink_map().keys():
            # Bad! tgt_full[i] is nonzero but this should be a
            # cancelled ideal.
            hex_ideal_i = hex(ideal_i)
            col_ideal_i = params.R.renumber_to_column(ideal_i)
            raise RuntimeError("tgt_full was not completely cancelled:"
                               f" The relation involves ideal number {ideal_i}"
                               f" or {hex_ideal_i}"
                               " (in renumber coordinates),"
                               f" which is {col_ideal_i}"
                               " (in column coordinates),"
                               " however it is not linked to the others."
                               " Too bad, really.")
            exit(0)

    num_filtered_ideals = len(SC_vector_filtered) - num_character_cols
    for j in range(num_character_cols):
        SC_vector_filtered[num_filtered_ideals + j] = updated_characters[j]

    return SC_vector_filtered, saved_kr_rels


@timing
def solve_filtering_plus_sage(params, M, MM, MC, SC_vector):
    '''
    Here, we are finding the sol vector (of coefficients) given that we have
    already done all filtering, and want to use sage solve_left.
    This returns a vector sol, of the full length, not the filtered length.
    sol should be a solution for the input target (SC_vector).
    '''
    # SC_vector is the target over all prime ideals. We are going to
    # iteratively update it so that SC_vector[i] is only nonzero if
    # ideal_i has survived filtering (has a column in the matrix).
    num_character_cols = MC.ncols()-M.ncols()
    len_tgt_full = len(SC_vector) - num_character_cols
    e = params.parameters['e']

    #killer_relations, order_of_cancellation = build_killer_rels_dict(params, MM, len_tgt_full)

    killer_relations = fast_json_load(params.files['KILLER_RELS_DICT'])
    order_of_cancellation = fast_json_load(params.files['KILLER_RELS_ORDER'])

    SC_vector_filtered, saved_kr_rels = apply_filtering_to_target(
        params, M, MM, MC, SC_vector, killer_relations, order_of_cancellation
    )

    # At this point, killer_relations should have all our necessary info.
    # We go in the order of order_of_cancellation, where the first item in the list
    # was the most recently cancelled. That is, we can write the ideal at index i
    # in terms of ideals at indices from i+1 on.
    # The surviving ideals should never show up in order_of_cancellation.

    timeprint("Updated target! Created SC_vector_filtered.")
    # print("Number of ideals that survived filtering is " + str(num_filtered_ideals))
    timeprint("The SC filtered vector is of length " + str(len(SC_vector_filtered)))
    sol_filtered = (MC.dense_matrix()).solve_left(SC_vector_filtered)
    assert(M.nrows()==len(sol_filtered))
    assert sol_filtered * MC == SC_vector_filtered
    timeprint("Solved the filtered linalg system using sage.")


    row_to_aquery = MM.row_to_aquery
    ab_per_row = MM.ab_per_row(params)

    aquery_to_row = {q:r for r,q in enumerate(row_to_aquery)}

    sol = vector(Integers(e), len(row_to_aquery))
    for row_in_filtered_M, (coeff, abs) in enumerate(zip(sol_filtered, ab_per_row)):
        for ab in abs.keys():
            c = abs[ab]
            a = ab[0]
            b = ab[1]
            row_in_og_M = aquery_to_row[(a,b)]
            sol[row_in_og_M] += (c * coeff)

    for (coeff, killer_rel) in saved_kr_rels:
        a = int(killer_rel.split(":")[0].split(",")[0], 16)
        b = int(killer_rel.split(":")[0].split(",")[1], 16)
        r = aquery_to_row[(a,b)]
        sol[r] += coeff

    return sol


@timing
def run_linalg_nobwc(params, use_target_info_file):

    timeprint("Loading LinearAlgebraMatrix_Filtered...")
    indexed_relations_file = params.files['AQRELS_FILE'] + ".indexed"
    purged_file = indexed_relations_file + ".purged"
    ideals_file = indexed_relations_file + ".ideals"
    index_file = indexed_relations_file + ".index"
    matrix_file = indexed_relations_file + ".matrix.bin"
    MM = LinearAlgebraMatrix_Filtered(params,
                                      indexed_relations_file,
                                      purged_file,
                                      ideals_file, index_file,
                                      matrix_file)
    timeprint("Finished loading LinearAlgebraMatrix_Filtered.")

    target_info = json.load(open(use_target_info_file))
    seed_ten = str(target_info['seed'])[:10]

    M = MM.matrix()

    u = ZZ(target_info['u'])
    v = ZZ(target_info['v'])
    uv_fac = get_uv_fac(u, v, target_info)

    DRELS_FILE = target_info['DRELS_FILE']
    DRELS_INDEXED = DRELS_FILE + ".cond.indexed"

    if 'n1024' in params.files['RENUMBERFILE']:
        partial_R=True
    else:
        partial_R=False

    timeprint("Start make_linalg_system")
    S_list, S_alg_vector, S_rat_vector = construct_S(
        params, DRELS_FILE, DRELS_INDEXED, u, v, already_indexed=False, partial_R=partial_R, uv_fac=uv_fac
    )
    T_list, ST_list, ST_alg_vector = truncate_S(
        params, S_list, S_alg_vector, partial_R=partial_R
    )
    MC, SC_vector, S_block, C_block = make_linalg_system(
        params, MM, ST_list, ST_alg_vector
    )
    timeprint("Finished make_linalg_system")

    # Now we have all of M, MM, MC, SC_vector, S_block, C_block

    num_character_cols = MC.ncols()-M.ncols()
    len_tgt_full = len(SC_vector) - num_character_cols
    e = params.parameters['e']

    killer_relations = fast_json_load(params.files['KILLER_RELS_DICT'])
    order_of_cancellation = fast_json_load(params.files['KILLER_RELS_ORDER'])

    SC_vector_filtered, saved_kr_rels = apply_filtering_to_target(
        params, M, MM, MC, -SC_vector, killer_relations, order_of_cancellation
    )

    SC_filtered_chars = SC_vector_filtered[-num_character_cols:]
    SC_filtered_nochars = SC_vector_filtered[:-num_character_cols]
    timeprint("Updated target! Created SC_vector_filtered.")
    timeprint("The rhs vector is of length " + str(len(SC_filtered_nochars)))

    matrix_file = MM.get_matrix_file(params)
    timeprint(f"Matrix file MM ({type(MM)}) is {matrix_file} of size {MM.nrows}x{MM.ncols}")
    timeprint(f"Matrix MC ({type(MC)}) of size {MC.nrows()}x{MC.ncols()}")
    timeprint(f"Matrix M ({type(M)}) of size {M.nrows()}x{M.ncols()}")

    assert(M.ncols()==len(SC_filtered_nochars))
    FILES_DIR = params.dirs['TEMP_OUTPUT_DIR']
    WORKDIR = os.path.join(FILES_DIR,"bwc/")

    bwc_mn = min(params.parameters['POLY_DEG']+1, num_character_cols+1)
    if 'bwc.numsols' in params.parameters and int(params.parameters['bwc.numsols']) > bwc_mn:
        bwc_mn = int(params.parameters['bwc.numsols'])

    timeprint("bwc completed! Starting some sage matrix multiplications.")
    Ze = Integers(e)
    start_makesolmatrix = time.time()
    sol_matrix_rows = []
    for j in range(bwc_mn):
        with open(f"{WORKDIR}/K.sols{j}-{j+1}.0.txt",'r') as f:
            sol_matrix_rows.append([int(line) for line in f])

    sol_matrix = matrix(Ze, sol_matrix_rows, sparse=True)
    end_makesolmatrix = time.time()

    # The optimization only works if all occurring integers can be represented as float64
    # (i.e., without losing precision).
    # This is the case for all integers up to 2^53.
    do_parallelized_mat_mult = (max(M.ncols(), S_block.ncols()) * e * e).bit_length() <= 53

    timeprint("do_parallelized_mat_mult", str(do_parallelized_mat_mult))
    check_value = (max(M.ncols(), S_block.ncols()) * e * e).bit_length()
    timeprint("check_value", str(check_value))

    # assert do_parallelized_mat_mult   # Fails but it may still be ok
    assert params.scipy_matrix

    from sparse_dot_mkl import dot_product_mkl
    import numpy as np
    from scipy.sparse import csr_array, bmat

    timeprint("Use scipy and Intel MKL for parallelized sparse matrix multiplication (interpreting entries in Z(e) as float64 and reducing mod e later; the code ensures this precision is not an issue, i.e., all integers are < 2^53)")

    def convert_sparse_matrix_scipy_to_sage(R, m):
        """
        Convert scipy.sparse array to a sparse sage matrix *without*
        constructing the full dense matrix during conversion.
        Forces the scipy matrix entries to be integers, then converts
        them to ring R (for us, Ze = Integers(e)).
        """
        m_coo = m.astype(int).tocoo()
        row_col_val_dict = {(r, c): v for r, c, v in zip(m_coo.row, m_coo.col, m_coo.data)}
        return matrix(R, *(m.shape), row_col_val_dict, sparse=True)

    def convert_sparse_matrix_sage_to_scipy(m, dtype):
        """
        Convert sparse sage matrix to scipy.sparse.csr_array *without*
        constructing the full dense matrix during conversion.
        Sparse sage vectors of dimension D are converted to a csr_array
        of form 1xD (because that's how the blocks work out with
        scipy.sparse.bmat later on).
        """
        m_dict = m.dict()
        l = len(m_dict)

        rows = [None] * l
        cols = [None] * l
        data = [None] * l
        i = 0

        if hasattr(m, "is_vector"):
            shape = (1, m.degree())
            r = 0
        else:
            shape = m.dimensions()

        for k, d in m_dict.items():
            if hasattr(m, "is_vector"):
                c = k
            else:
                r, c = k
            rows[i] = r
            cols[i] = c
            data[i] = d
            i += 1

        return csr_array((data, (rows, cols)), shape=shape, dtype=dtype)

    ### Time assert
    start_assert = time.time()

    timeprint("Convert sparse sage matrices to scipy (sol_matrix, block_matrix = [[M], [SC_filtered_nochars]])")
    timeprint("type of M", type(M))
    sol_matrix_sc = convert_sparse_matrix_sage_to_scipy(sol_matrix, dtype=np.float64)
    SC_filtered_nochars_sc = convert_sparse_matrix_sage_to_scipy(SC_filtered_nochars, dtype=np.float64)
    if params.scipy_matrix:
        # M is already stored as a scipy matrix
        M_sc = M.scipy_M()
    else:
        assert False    # shouldn't hit this case
    block_matrix_filtered_nochars_sc = bmat([
        [M_sc],
        [SC_filtered_nochars_sc]
    ])

    timeprint("Multiply sol_matrix * block_matrix = sol_times_filtered_nochars")
    sol_times_filtered_nochars_sc = dot_product_mkl(sol_matrix_sc, block_matrix_filtered_nochars_sc)
    sol_times_filtered_nochars = convert_sparse_matrix_scipy_to_sage(Ze, sol_times_filtered_nochars_sc)
    assert sol_times_filtered_nochars.is_zero()
    end_assert = time.time()

    ### Time matonright
    start_matonright = time.time()
    timeprint("Convert sparse sage matrices to scipy (block_matrix_2 = [[S_block], [SC_filtered_chars]])")
    S_block_sc = convert_sparse_matrix_sage_to_scipy(S_block, dtype=np.float64)
    SC_filtered_chars_sc = convert_sparse_matrix_sage_to_scipy(SC_filtered_chars, dtype=np.float64)
    block_matrix_filtered_chars_sc = bmat([
        [S_block_sc],
        [SC_filtered_chars_sc]
    ])

    timeprint("Multiply sol_matrix * block_matrix_2 = small_matrix_on_the_right")
    small_matrix_on_the_right_sc = dot_product_mkl(sol_matrix_sc, block_matrix_filtered_chars_sc)
    small_matrix_on_the_right = convert_sparse_matrix_scipy_to_sage(Ze, small_matrix_on_the_right_sc)
    end_matonright = time.time()

    params.timing["making sol_matrix"] = end_makesolmatrix - start_makesolmatrix
    params.timing["assert sol_matrix"] = end_assert - start_assert
    params.timing["small_matrix_on_the_right"] = end_matonright-start_matonright

    timeprint(f"found sol_matrix, {sol_matrix.nrows()}x{sol_matrix.ncols()}")

    # At this point, SC_filtered_chars is the target characters.  We need
    # to compute the characters for each of the rows in the sol_matrix.
    # Then, some combination of sol_matrix should give us
    # SC_filtered_chars, and we use this as our overall solution
    # (something in the nullspace of M|t).  Recall also that S_block is
    # the character portion of the matrix.  C_block we totally ignore,
    # since that is from the pre-filtered target vector.

    # First get the left nullspace of small_matrix_on_the_right
    second_solve = small_matrix_on_the_right
    start_secondsolve = time.time()
    ns2 = second_solve.left_kernel().basis_matrix()
    end_secondsolve = time.time()
    params.timing["second_solve.left_kernel"] = end_secondsolve-start_secondsolve
    timeprint("Started the second solve!")

    timeprint("ns2 nrows = " + str(ns2.nrows()))
    timeprint("ns2 ncols = " + str(ns2.ncols()))

    if ns2.nrows() == 0:
        raise RuntimeError("Cannot solve the small linear system after bwc.")

    start_lastsolve = time.time()
    # find a row combination that reaches -1 on the last coordinate
    assert not (ns2 * sol_matrix[:,-1:]).is_zero(), "The last column of ns2*sol_matrix is zero, which (I think) means that no combination of the solutions bwc found to the non-character part of the matrix will be able to cancel out the target characters"
    combine = (ns2 * sol_matrix[:,-1:]).solve_left(vector([-1]))
    end_lastsolve = time.time()
    params.timing["final solve_left"] = end_lastsolve-start_lastsolve
    timeprint("Completed the final solve_left!")

    start_solf = time.time()
    sol_f = combine * ns2 * sol_matrix[:,:-1]

    # This is the same assert that we have in solve_filtering_plus_sage.
    # It must hold!
    timeprint("type of sol_f:", str(type(sol_f)))
    timeprint("type of MC:", str(type(MC)))
    timeprint("type of SC_vector_filtered:", str(type(SC_vector_filtered)))

    def convert_sparse_vector_scipy_to_sage_dense(R, v):
        v_coo = v.astype(int).tocoo()
        row_col_val_dict = {(r, c): v for r, c, v in zip(v_coo.row, v_coo.col, v_coo.data)}
        assert v.shape[0] == 1
        v_len = v.shape[1]
        v_s = vector(R, v_len)  # not sparse

        for k in row_col_val_dict.keys():
            (r,c) = k
            v = row_col_val_dict[k]
            assert r == 0
            v_s[c] = R(v)

        return v_s

    def convert_sage_dense_vector_to_scipy(v, dtype):
        # vector is not sparse!
        l = len(v)
        rows = [None] * l
        cols = [None] * l
        data = [None] * l
        i = 0

        for c in range(len(v)):
            r = 0
            d = int(v[c])
            rows[i] = r
            cols[i] = c
            data[i] = d
            i += 1

        return csr_array((data, (rows, cols)), shape=(1,len(v)), dtype=dtype)

    sol_f_scipy = convert_sage_dense_vector_to_scipy(sol_f, dtype=np.float64)
    M_scipy = M.scipy_M()
    S_block_scipy = convert_sparse_matrix_sage_to_scipy(S_block, dtype=np.float64)
    assert M_scipy.shape[0] == S_block_scipy.shape[0]
    MC_scipy = bmat([[
        M_scipy,
        S_block_scipy
    ]])
    leftside_scipy = dot_product_mkl(sol_f_scipy, MC_scipy)
    leftside_sage = convert_sparse_vector_scipy_to_sage_dense(Ze, leftside_scipy)

    timeprint("type of leftside_sage", str(type(leftside_sage)))
    timeprint("type of SC_vector_filtered", str(type(SC_vector_filtered)))
    assert leftside_sage == SC_vector_filtered

    end_solf = time.time()
    params.timing["sol_f computations"] = end_solf-start_solf

    start_filtering = time.time()
    row_to_aquery = MM.row_to_aquery
    ab_per_row = MM.ab_per_row(params)
    aquery_to_row = {q:r for r,q in enumerate(row_to_aquery)}

    sol = vector(Integers(e), len(row_to_aquery))
    for row_in_filtered_M, (coeff, abs) in enumerate(zip(sol_f, ab_per_row)):
        for ab in abs.keys():
            c = abs[ab]
            a = ab[0]
            b = ab[1]
            row_in_og_M = aquery_to_row[(a,b)]
            sol[row_in_og_M] += (c * coeff)

    for (coeff, killer_rel) in saved_kr_rels:
        a = int(killer_rel.split(":")[0].split(",")[0], 16)
        b = int(killer_rel.split(":")[0].split(",")[1], 16)
        r = aquery_to_row[(a,b)]
        sol[r] += coeff

    end_filtering = time.time()
    params.timing["recover sol from sol_f"] = end_filtering - start_filtering

    # save stuff to files, not sure if these are all small enough...
    sol_sparse = sol.sparse_vector()
    fast_persistent_save(sol_sparse, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-sol_sparse.sobj")
    fast_persistent_save(sol, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-sol.sobj")
    fast_persistent_save(ST_list, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-ST_list.sobj")
    fast_persistent_save(T_list, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-T_list.sobj")
    fast_persistent_save(ST_alg_vector, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-ST_alg_vector.sobj")
    fast_persistent_save(S_rat_vector, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-S_rat_vector.sobj")

    toreturn = LinalgOutput(sol=sol_sparse, ST_list=ST_list, ST_alg_vector=ST_alg_vector, T_list=T_list, row_to_aquery=row_to_aquery, S_rat_vector=S_rat_vector, indexed_relations_file=indexed_relations_file)

    fast_persistent_save(toreturn, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-linalgoutput.sobj")

    run_sol_sanity_checks(params, sol, ST_alg_vector, indexed_relations_file)
    timeprint("Found a solution! We have sol*M = ST mod e.")

    timeprint("Making indexed relations table-of-contents file...")
    tocfile.maketoc(indexed_relations_file, indexed_relations_file + ".toc")
    timeprint("Done making indexed relations table-of-contents file")

    return toreturn


@timing
def run_linalg_nobwc_noairbags(params, use_target_info_file):

    timeprint("Loading LinearAlgebraMatrix_Filtered...")
    indexed_relations_file = params.files['AQRELS_FILE'] + ".indexed"
    purged_file = indexed_relations_file + ".purged"
    ideals_file = indexed_relations_file + ".ideals"
    index_file = indexed_relations_file + ".index"
    matrix_file = indexed_relations_file + ".matrix.bin"
    MM = LinearAlgebraMatrix_Filtered(params,
                                      indexed_relations_file,
                                      purged_file,
                                      ideals_file, index_file,
                                      matrix_file)
    timeprint("Finished loading LinearAlgebraMatrix_Filtered.")

    target_info = json.load(open(use_target_info_file))
    seed_ten = str(target_info['seed'])[:10]

    M = MM.matrix()

    u = ZZ(target_info['u'])
    v = ZZ(target_info['v'])
    uv_fac = get_uv_fac(u, v, target_info)

    DRELS_FILE = target_info['DRELS_FILE']
    DRELS_INDEXED = DRELS_FILE + ".cond.indexed"

    if 'n1024' in params.files['RENUMBERFILE']:
        partial_R=True
    else:
        partial_R=False

    timeprint("Start make_linalg_system")
    S_list, S_alg_vector, S_rat_vector = construct_S(
        params, DRELS_FILE, DRELS_INDEXED, u, v, already_indexed=False, partial_R=partial_R, uv_fac=uv_fac
    )
    T_list, ST_list, ST_alg_vector = truncate_S(
        params, S_list, S_alg_vector, partial_R=partial_R
    )
    MC, SC_vector, S_block, C_block = make_linalg_system(
        params, MM, ST_list, ST_alg_vector
    )
    timeprint("Finished make_linalg_system")

    # Now we have all of M, MM, MC, SC_vector, S_block, C_block

    num_character_cols = MC.ncols()-M.ncols()
    len_tgt_full = len(SC_vector) - num_character_cols
    e = params.parameters['e']

    killer_relations = fast_json_load(params.files['KILLER_RELS_DICT'])
    order_of_cancellation = fast_json_load(params.files['KILLER_RELS_ORDER'])

    SC_vector_filtered, saved_kr_rels = apply_filtering_to_target(
        params, M, MM, MC, -SC_vector, killer_relations, order_of_cancellation
    )

    SC_filtered_chars = SC_vector_filtered[-num_character_cols:]
    SC_filtered_nochars = SC_vector_filtered[:-num_character_cols]
    timeprint("Updated target! Created SC_vector_filtered.")
    timeprint("The rhs vector is of length " + str(len(SC_filtered_nochars)))

    matrix_file = MM.get_matrix_file(params)
    timeprint(f"Matrix file MM ({type(MM)}) is {matrix_file} of size {MM.nrows}x{MM.ncols}")
    timeprint(f"Matrix MC ({type(MC)}) of size {MC.nrows()}x{MC.ncols()}")
    timeprint(f"Matrix M ({type(M)}) of size {M.nrows()}x{M.ncols()}")

    assert(M.ncols()==len(SC_filtered_nochars))
    FILES_DIR = params.dirs['TEMP_OUTPUT_DIR']
    WORKDIR = os.path.join(FILES_DIR,"bwc/")

    bwc_mn = min(params.parameters['POLY_DEG']+1, num_character_cols+1)
    if 'bwc.numsols' in params.parameters and int(params.parameters['bwc.numsols']) > bwc_mn:
        bwc_mn = int(params.parameters['bwc.numsols'])

    timeprint("bwc completed! Starting some sage matrix multiplications.")
    Ze = Integers(e)
    start_makesolmatrix = time.time()
    sol_matrix_rows = []
    for j in range(bwc_mn):
        with open(f"{WORKDIR}/K.sols{j}-{j+1}.0.txt",'r') as f:
            sol_matrix_rows.append([int(line) for line in f])

    sol_matrix = matrix(Ze, sol_matrix_rows, sparse=True)
    end_makesolmatrix = time.time()

    start_matonright = time.time()
    small_matrix_on_the_right = sol_matrix * block_matrix(2,1,[S_block,matrix([SC_filtered_chars])], sparse=True)
    end_matonright = time.time()

    params.timing["making sol_matrix"] = end_makesolmatrix - start_makesolmatrix
    params.timing["small_matrix_on_the_right"] = end_matonright-start_matonright

    timeprint(f"found sol_matrix, {sol_matrix.nrows()}x{sol_matrix.ncols()}")

    # At this point, SC_filtered_chars is the target characters.  We need
    # to compute the characters for each of the rows in the sol_matrix.
    # Then, some combination of sol_matrix should give us
    # SC_filtered_chars, and we use this as our overall solution
    # (something in the nullspace of M|t).  Recall also that S_block is
    # the character portion of the matrix.  C_block we totally ignore,
    # since that is from the pre-filtered target vector.

    # First get the left nullspace of small_matrix_on_the_right
    second_solve = small_matrix_on_the_right
    start_secondsolve = time.time()
    ns2 = second_solve.left_kernel().basis_matrix()
    end_secondsolve = time.time()
    params.timing["second_solve.left_kernel"] = end_secondsolve-start_secondsolve
    timeprint("Started the second solve!")

    timeprint("ns2 nrows = " + str(ns2.nrows()))
    timeprint("ns2 ncols = " + str(ns2.ncols()))

    if ns2.nrows() == 0:
        raise RuntimeError("Cannot solve the small linear system after bwc.")

    start_lastsolve = time.time()
    # find a row combination that reaches -1 on the last coordinate
    assert not (ns2 * sol_matrix[:,-1:]).is_zero(), "The last column of ns2*sol_matrix is zero, which (I think) means that no combination of the solutions bwc found to the non-character part of the matrix will be able to cancel out the target characters"
    combine = (ns2 * sol_matrix[:,-1:]).solve_left(vector([-1]))
    end_lastsolve = time.time()
    params.timing["final solve_left"] = end_lastsolve-start_lastsolve
    timeprint("Completed the final solve_left!")

    start_solf = time.time()
    sol_f = combine * ns2 * sol_matrix[:,:-1]
    end_solf = time.time()
    params.timing["sol_f computations"] = end_solf-start_solf

    start_filtering = time.time()
    row_to_aquery = MM.row_to_aquery
    ab_per_row = MM.ab_per_row(params)
    aquery_to_row = {q:r for r,q in enumerate(row_to_aquery)}

    sol = vector(Integers(e), len(row_to_aquery))
    for row_in_filtered_M, (coeff, abs) in enumerate(zip(sol_f, ab_per_row)):
        for ab in abs.keys():
            c = abs[ab]
            a = ab[0]
            b = ab[1]
            row_in_og_M = aquery_to_row[(a,b)]
            sol[row_in_og_M] += (c * coeff)

    for (coeff, killer_rel) in saved_kr_rels:
        a = int(killer_rel.split(":")[0].split(",")[0], 16)
        b = int(killer_rel.split(":")[0].split(",")[1], 16)
        r = aquery_to_row[(a,b)]
        sol[r] += coeff

    end_filtering = time.time()
    params.timing["recover sol from sol_f"] = end_filtering - start_filtering

    # save stuff to files, not sure if these are all small enough...
    sol_sparse = sol.sparse_vector()
    fast_persistent_save(sol_sparse, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-sol_sparse.sobj")
    fast_persistent_save(sol, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-sol.sobj")
    fast_persistent_save(ST_list, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-ST_list.sobj")
    fast_persistent_save(T_list, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-T_list.sobj")
    fast_persistent_save(ST_alg_vector, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-ST_alg_vector.sobj")
    fast_persistent_save(S_rat_vector, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-S_rat_vector.sobj")

    run_sol_sanity_checks(params, sol, ST_alg_vector, indexed_relations_file)
    timeprint("Found a solution! We have sol*M = ST mod e.")

    timeprint("Making indexed relations table-of-contents file...")
    tocfile.maketoc(indexed_relations_file, indexed_relations_file + ".toc")
    timeprint("Done making indexed relations table-of-contents file")

    toreturn = LinalgOutput(sol=sol_sparse, ST_list=ST_list, ST_alg_vector=ST_alg_vector, T_list=T_list, row_to_aquery=row_to_aquery, S_rat_vector=S_rat_vector, indexed_relations_file=indexed_relations_file)

    fast_persistent_save(toreturn, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-linalgoutput.sobj")

    return toreturn


@timing
def solve_filtering_plus_bwc(params, M, MM, MC, SC_vector, S_block, C_block, seed_ten=None):
    # Notes: S_block is the character portion of the matrix. C_block is the
    # character portion of the target vector, BEFORE filtering. The BWC-output
    # matrix (the ...matrix.bin file) does NOT include character columns.

    num_character_cols = MC.ncols()-M.ncols()
    len_tgt_full = len(SC_vector) - num_character_cols
    e = params.parameters['e']

    #killer_relations, order_of_cancellation = build_killer_rels_dict(params, M, MM, MC, len_tgt_full)

    killer_relations = fast_json_load(params.files['KILLER_RELS_DICT'])
    order_of_cancellation = fast_json_load(params.files['KILLER_RELS_ORDER'])

    if seed_ten is None:
        seed_pfx = ""
    else:
        seed_pfx = seed_ten + "-"

    SC_vector_filtered_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"SC_vector_filtered.sobj"
    saved_kr_rels_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"saved_kr_rels.sobj"

    if os.path.exists(SC_vector_filtered_file) and os.path.exists(saved_kr_rels_file):
        SC_vector_filtered = fast_persistent_load(SC_vector_filtered_file)
        saved_kr_rels = fast_persistent_load(saved_kr_rels_file)
    else:
        SC_vector_filtered, saved_kr_rels = apply_filtering_to_target(
            params, M, MM, MC, -SC_vector, killer_relations, order_of_cancellation
        )

        if 'n1024' not in params.files['RENUMBERFILE']:
            fast_persistent_save(SC_vector_filtered, params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"SC_vector_filtered.sobj")
            fast_persistent_save(saved_kr_rels, params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"saved_kr_rels.sobj")

    # separate out character columns
    SC_filtered_chars = SC_vector_filtered[-num_character_cols:]
    SC_filtered_nochars = SC_vector_filtered[:-num_character_cols]

    timeprint("Updated target! Created SC_vector_filtered.")
    timeprint("The rhs vector is of length " + str(len(SC_filtered_nochars)))

    matrix_file = MM.get_matrix_file(params)
    timeprint(f"Matrix file MM ({type(MM)}) is {matrix_file} of size {MM.nrows}x{MM.ncols}")
    timeprint(f"Matrix MC ({type(MC)}) of size {MC.nrows()}x{MC.ncols()}")
    timeprint(f"Matrix M ({type(M)}) of size {M.nrows()}x{M.ncols()}")
    vector_file = params.dirs['TEMP_OUTPUT_DIR'] + "tgt.ascii"

    start_writetgt = time.time()
    assert(M.ncols()==len(SC_filtered_nochars))
    write_ascii_vector(vector_file, SC_filtered_nochars, e)
    end_writetgt = time.time()
    params.timing["write_tgt_vector"] = end_writetgt-start_writetgt

    FILES_DIR = params.dirs['TEMP_OUTPUT_DIR']
    WORKDIR = os.path.join(FILES_DIR,"bwc/")
    e = int(params.parameters['e'])
    Path(WORKDIR).mkdir(exist_ok=True)
    make_and_clean(WORKDIR)

    bwc_mn = min(params.parameters['POLY_DEG']+1, num_character_cols+1)
    if 'bwc.numsols' in params.parameters and int(params.parameters['bwc.numsols']) > bwc_mn:
        bwc_mn = int(params.parameters['bwc.numsols'])

    CPUBINDING_CONF_FILE = params.files['CPUBINDING_CONF_FILE']

    # given the bwc sequence length, set the checkpointing interval
    # close to its square root -- actually slightly larger because
    # I like it better with fewer checkpints.
    L = MC.ncols() / bwc_mn + MC.ncols() / bwc_mn + 32
    interval = min(2**int(log(L,2)/2 + 3), L)

    bwc_executable = "linalg/bwc/bwc.pl"

    bwc_common_args = [
        "--matrix", os.path.realpath(matrix_file),
        "--prime", params.parameters['e'],
        "--nullspace", "LEFT",
        "--rhs", os.path.realpath(vector_file),
        "--wdir", WORKDIR,
        f"interval={interval}",
        f"mn={bwc_mn}",
        "balancing_options=reorder=columns",
        f"cpubinding={CPUBINDING_CONF_FILE}",
        "verbose_flags=^all-cmdline,^bwc-timing-grids,^all-bwc-dispatch,^bwc-cache-major-info,perl-cmdline,perl-sections,^perl-checks,^bwc-iteration-timings,^bwc-loading-mksol-files",
    ]

    bwc_add_args = []
    if params.mpi:
        timeprint("Run bwc.pl from solve_filtering_plus_bwc with MPI enabled.")
        if not params.bwc_slurm:
            bwc_add_args += [
                "--mpi", params.parameters['bwc.mpi']
            ]

        # Escape spaces in arguments correctly, depending on whether the command is
        # written to a file (slurm), or executed with subprocess (non-slurm).
        if params.bwc_slurm:
            bwc_add_args.append(f"mpi_extra_args='{params.mpi_extra_args}'")
        else:
            bwc_add_args.append(f"mpi_extra_args={params.mpi_extra_args}")

    start_bwc = time.time()

    if params.bwc_slurm:
        bwc_tasks = ["prep", "krylov", "lingen", "mksol", "gather", "cleanup"]

        if params.mpi:
            if params.parameters["bwc.mpi_slurm"] == "auto":
                a, b = map(int, params.parameters["bwc.mpi"].split("x"))
                bwc_mpi_slurm = find_factors_close_to_square_root((a * b) // bwc_mn)
                timeprint(f"Setting bwc.mpi_slurm to automatically determined value of {bwc_mpi_slurm}")
            else:
                bwc_mpi_slurm = params.parameters["bwc.mpi_slurm"]
                timeprint(f"Setting bwc.mpi_slurm to user-defined value of {bwc_mpi_slurm}")

        single_threaded_tasks = ["prep", "lingen", "cleanup"]

        setup_commands_common = [
            params.files['PYTHON'],
            "wait_for_file.py",
            "-m", "30000",
            "-f"
        ]

        # XXX: It seems that parallelizing mksol across 6 nodes
        # is already working. Maybe we already include the SM columns
        # in the matrix?
        # if not params.chars_in_mat:
        #     single_threaded_tasks.append("mksol")

        for bwc_task in bwc_tasks:
            start_task = time.time()

            processes = []
            njobs = bwc_mn if bwc_task not in single_threaded_tasks else 1

            for jobnum in range(njobs):
                task_args = []
                setup_commands = setup_commands_common[::] if bwc_task != "prep" else []

                if bwc_task == "krylov":
                    task_args += [
                        "--ys", f"{jobnum}..{jobnum+1}",
                        "skip_online_checks=1"
                    ]
                    setup_commands += [
                        f"{WORKDIR}aqrels.out.indexed.matrix.*/*",
                        f"{WORKDIR}X",
                        f"{WORKDIR}V{jobnum}-{jobnum+1}.0"
                    ]
                elif bwc_task == "lingen":
                    for jobnum in range(bwc_mn):
                        setup_commands.append(f"{WORKDIR}A{jobnum}-{jobnum+1}.*")
                elif bwc_task == "mksol":
                    task_args += [
                        f"solutions={jobnum}-{jobnum+1}",
                        "skip_online_checks=1"
                    ]
                    setup_commands += [
                        f"{WORKDIR}A0-{bwc_mn}.*",
                        f"{WORKDIR}F.sols{jobnum}-{jobnum+1}.*"
                    ]
                elif bwc_task == "gather":
                    task_args.append(f"solutions={jobnum}-{jobnum+1}")
                    setup_commands.append(f"{WORKDIR}S.sols{jobnum}-{jobnum+1}.*")
                elif bwc_task == "cleanup":
                    for jobnum in range(bwc_mn):
                        setup_commands.append(f"{WORKDIR}K.sols{jobnum}-{jobnum+1}.*")

                if params.mpi:
                    # Running prep with different mpi= argument results in file corruption that
                    # surfaces as the following bug:
                    #
                    #   terminate called after throwing an instance of 'std::runtime_error'
                    #       what():  code BUG() : condition rc == 15 failed at
                    #       *** Error: caught signal "Aborted"
                    #
                    # Likely because krylov will rebuild the needed matrix cache file on every slurm node
                    # and that appears to write to the same file concurrently.
                    if bwc_task != "prep" and njobs == 1:
                        task_args += [
                            "--mpi", params.parameters['bwc.mpi']
                        ]
                    else:
                        task_args += [
                            "--mpi", bwc_mpi_slurm
                        ]

                command_list = [
                    params.files['SAGE'],
                    "cado_nfs_binaries.py",
                    "--cado-build-dir", params.dirs['CADO_BUILD_DIR'],
                    "--thr", params.parameters["bwc.thr"],
                    bwc_executable,
                    bwc_task,
                    *bwc_common_args,
                    *task_args,
                    *bwc_add_args
                ]

                # all ints must be str to pass on command line.
                command_list = list(map(str, command_list))

                processes.append(
                    slurmit(
                        params,
                        " ".join(command_list),
                        params.prefix[:-1] + f"-bwc-{bwc_task}",
                        jobnum,
                        setup_commands=[' '.join(setup_commands) + '\n'],
                        partition=params.bwc_slurm_job_partition,
                        nested_mpi=params.mpi
                    )
                )

            _, cputime_task = slurm_wait(processes, throw_error_on_failed_job=True, params=params)
            overall_cputime.add(cputime_task)

            end_task = time.time()
            params.timing[f"bwc {bwc_task}"] = end_task - start_task
            params.timing[f"bwc {bwc_task}_cputime"] = cputime_task
    else:
        with open(os.path.join(WORKDIR,"bwc.out"),"w") as outfile, \
            open(os.path.join(WORKDIR,"bwc.err"),"w") as errfile:

            print(f"Writing bwc stdout to {outfile.name}")
            print(f"Writing bwc stderr to {errfile.name}")

            CadoNFS(bwc_executable,
                ":complete",
                *bwc_common_args,
                "--solutions", f"0-{bwc_mn}",
                "--thr", params.parameters["bwc.thr.controller"],
                *bwc_add_args,
                capture=outfile,
                stderr=errfile
            )

    end_bwc = time.time()
    params.timing["linalg/bwc/bwc.pl (complete)"] = end_bwc-start_bwc

    timeprint("bwc completed! Starting some sage matrix multiplications.")

    Ze = Integers(e)

    start_makesolmatrix = time.time()
    sol_matrix_rows = []
    for j in range(bwc_mn):
        with open(f"{WORKDIR}/K.sols{j}-{j+1}.0.txt",'r') as f:
            sol_matrix_rows.append([int(line) for line in f])

    sol_matrix = matrix(Ze, sol_matrix_rows, sparse=True)
    end_makesolmatrix = time.time()

    # The optimization only works if all occurring integers can be represented as float64 (i.e., without losing precision).
    # This is the case for all integers up to 2^53.
    do_parallelized_mat_mult = (max(M.ncols(), S_block.ncols()) * e * e).bit_length() <= 53

    if do_parallelized_mat_mult:
        try:
            from sparse_dot_mkl import dot_product_mkl
            import numpy as np
            from scipy.sparse import csr_array, bmat
        except:
            timeprint("To speed up the following matrix multiplications, install the following python packages: sparse_dot_mkl, numpy, scipy.")
            do_parallelized_mat_mult = False

        if do_parallelized_mat_mult:
            timeprint("Use scipy and Intel MKL for parallelized sparse matrix multiplication (interpreting entries in Z(e) as float64 and reducing mod e later; the code ensures this precision is not an issue, i.e., all integers are < 2^53)")

            def convert_sparse_matrix_scipy_to_sage(R, m):
                """
                Convert scipy.sparse array to a sparse sage matrix *without*
                constructing the full dense matrix during conversion.
                Forces the scipy matrix entries to be integers, then converts
                them to ring R (for us, Ze = Integers(e)).
                """
                m_coo = m.astype(int).tocoo()
                row_col_val_dict = {(r, c): v for r, c, v in zip(m_coo.row, m_coo.col, m_coo.data)}
                return matrix(R, *(m.shape), row_col_val_dict, sparse=True)

            def convert_sparse_matrix_sage_to_scipy(m, dtype):
                """
                Convert sparse sage matrix to scipy.sparse.csr_array *without*
                constructing the full dense matrix during conversion.
                Sparse sage vectors of dimension D are converted to a csr_array
                of form 1xD (because that's how the blocks work out with
                scipy.sparse.bmat later on).
                """
                m_dict = m.dict()
                l = len(m_dict)

                rows = [None] * l
                cols = [None] * l
                data = [None] * l
                i = 0

                if hasattr(m, "is_vector"):
                    shape = (1, m.degree())
                    r = 0
                else:
                    shape = m.dimensions()

                for k, d in m_dict.items():
                    if hasattr(m, "is_vector"):
                        c = k
                    else:
                        r, c = k
                    rows[i] = r
                    cols[i] = c
                    data[i] = d
                    i += 1

                return csr_array((data, (rows, cols)), shape=shape, dtype=dtype)

            ### Time assert
            start_assert = time.time()

            timeprint("Convert sparse sage matrices to scipy (sol_matrix, block_matrix = [[M], [SC_filtered_nochars]])")
            timeprint("type of M", type(M))
            sol_matrix_sc = convert_sparse_matrix_sage_to_scipy(sol_matrix, dtype=np.float64)
            SC_filtered_nochars_sc = convert_sparse_matrix_sage_to_scipy(SC_filtered_nochars, dtype=np.float64)
            if params.scipy_matrix:
                # M is already stored as a scipy matrix
                M_sc = M.scipy_M()
            else:
                M_sc = convert_sparse_matrix_sage_to_scipy(M, dtype=np.float64)
            block_matrix_filtered_nochars_sc = bmat([
                [M_sc],
                [SC_filtered_nochars_sc]
            ])

            timeprint("Multiply sol_matrix * block_matrix = sol_times_filtered_nochars")
            sol_times_filtered_nochars_sc = dot_product_mkl(sol_matrix_sc, block_matrix_filtered_nochars_sc)
            sol_times_filtered_nochars = convert_sparse_matrix_scipy_to_sage(Ze, sol_times_filtered_nochars_sc)
            assert sol_times_filtered_nochars.is_zero()

            end_assert = time.time()

            ### Time matonright
            start_matonright = time.time()

            timeprint("Convert sparse sage matrices to scipy (block_matrix_2 = [[S_block], [SC_filtered_chars]])")
            S_block_sc = convert_sparse_matrix_sage_to_scipy(S_block, dtype=np.float64)
            SC_filtered_chars_sc = convert_sparse_matrix_sage_to_scipy(SC_filtered_chars, dtype=np.float64)
            block_matrix_filtered_chars_sc = bmat([
                [S_block_sc],
                [SC_filtered_chars_sc]
            ])

            timeprint("Multiply sol_matrix * block_matrix_2 = small_matrix_on_the_right")
            small_matrix_on_the_right_sc = dot_product_mkl(sol_matrix_sc, block_matrix_filtered_chars_sc)
            small_matrix_on_the_right = convert_sparse_matrix_scipy_to_sage(Ze, small_matrix_on_the_right_sc)

            end_matonright = time.time()
    if not do_parallelized_mat_mult:
        timeprint("Use sage matrix(Integers(e), entries, sparse=True) for multiplication [this is single threaded and can take hours]")

        # sol_matrix * vertical_join([MC, SC_vector_filtered]) decomposes as
        # follows:

        start_assert = time.time()
        assert (sol_matrix * block_matrix(2,1,[M,matrix(SC_filtered_nochars)], sparse=True)).is_zero()
        end_assert = time.time()

        start_matonright = time.time()
        small_matrix_on_the_right = sol_matrix * block_matrix(2,1,[S_block,matrix([SC_filtered_chars])], sparse=True)
        end_matonright = time.time()

    params.timing["making sol_matrix"] = end_makesolmatrix - start_makesolmatrix
    params.timing["assert sol_matrix"] = end_assert - start_assert
    params.timing["small_matrix_on_the_right"] = end_matonright-start_matonright

    timeprint(f"found sol_matrix, {sol_matrix.nrows()}x{sol_matrix.ncols()}")

    # At this point, SC_filtered_chars is the target characters.  We need
    # to compute the characters for each of the rows in the sol_matrix.
    # Then, some combination of sol_matrix should give us
    # SC_filtered_chars, and we use this as our overall solution
    # (something in the nullspace of M|t).  Recall also that S_block is
    # the character portion of the matrix.  C_block we totally ignore,
    # since that is from the pre-filtered target vector.

    # First get the left nullspace of small_matrix_on_the_right
    second_solve = small_matrix_on_the_right
    start_secondsolve = time.time()
    ns2 = second_solve.left_kernel().basis_matrix()
    end_secondsolve = time.time()
    params.timing["second_solve.left_kernel"] = end_secondsolve-start_secondsolve
    timeprint("Started the second solve!")

    timeprint("ns2 nrows = " + str(ns2.nrows()))
    timeprint("ns2 ncols = " + str(ns2.ncols()))

    if ns2.nrows() == 0:
        raise RuntimeError("Cannot solve the small linear system after bwc.")

    start_lastsolve = time.time()
    # find a row combination that reaches -1 on the last coordinate
    assert not (ns2 * sol_matrix[:,-1:]).is_zero(), "The last column of ns2*sol_matrix is zero, which (I think) means that no combination of the solutions bwc found to the non-character part of the matrix will be able to cancel out the target characters"
    combine = (ns2 * sol_matrix[:,-1:]).solve_left(vector([-1]))
    end_lastsolve = time.time()
    params.timing["final solve_left"] = end_lastsolve-start_lastsolve
    timeprint("Completed the final solve_left!")

    start_solf = time.time()
    sol_f = combine * ns2 * sol_matrix[:,:-1]

    # This is the same assert that we have in solve_filtering_plus_sage.
    # It must hold!
    timeprint("type of sol_f:", str(type(sol_f)))
    timeprint("type of MC:", str(type(MC)))
    timeprint("type of SC_vector_filtered:", str(type(SC_vector_filtered)))


    if params.scipy_matrix:

        from sparse_dot_mkl import dot_product_mkl
        import numpy as np
        from scipy.sparse import csr_array, bmat

        def convert_sparse_vector_scipy_to_sage_dense(R, v):
            v_coo = v.astype(int).tocoo()
            row_col_val_dict = {(r, c): v for r, c, v in zip(v_coo.row, v_coo.col, v_coo.data)}
            assert v.shape[0] == 1
            v_len = v.shape[1]
            v_s = vector(R, v_len)  # not sparse

            for k in row_col_val_dict.keys():
                (r,c) = k
                v = row_col_val_dict[k]
                assert r == 0
                v_s[c] = R(v)

            return v_s

        def convert_sage_dense_vector_to_scipy(v, dtype):
            # vector is not sparse!
            l = len(v)
            rows = [None] * l
            cols = [None] * l
            data = [None] * l
            i = 0

            for c in range(len(v)):
                r = 0
                d = int(v[c])
                rows[i] = r
                cols[i] = c
                data[i] = d
                i += 1

            return csr_array((data, (rows, cols)), shape=(1,len(v)), dtype=dtype)

        def convert_sparse_matrix_sage_to_scipy(m, dtype):
            """
            Convert sparse sage matrix to scipy.sparse.csr_array *without*
            constructing the full dense matrix during conversion.
            Sparse sage vectors of dimension D are converted to a csr_array
            of form 1xD (because that's how the blocks work out with
            scipy.sparse.bmat later on).
            """
            m_dict = m.dict()
            l = len(m_dict)

            rows = [None] * l
            cols = [None] * l
            data = [None] * l
            i = 0

            if hasattr(m, "is_vector"):
                shape = (1, m.degree())
                r = 0
            else:
                shape = m.dimensions()

            for k, d in m_dict.items():
                if hasattr(m, "is_vector"):
                    c = k
                else:
                    r, c = k
                rows[i] = r
                cols[i] = c
                data[i] = d
                i += 1

            return csr_array((data, (rows, cols)), shape=shape, dtype=dtype)

        sol_f_scipy = convert_sage_dense_vector_to_scipy(sol_f, dtype=np.float64)
        M_scipy = M.scipy_M()
        S_block_scipy = convert_sparse_matrix_sage_to_scipy(S_block, dtype=np.float64)
        assert M_scipy.shape[0] == S_block_scipy.shape[0]
        MC_scipy = bmat([[
            M_scipy,
            S_block_scipy
        ]])
        leftside_scipy = dot_product_mkl(sol_f_scipy, MC_scipy)
        leftside_sage = convert_sparse_vector_scipy_to_sage_dense(Ze, leftside_scipy)

        timeprint("type of leftside_sage", str(type(leftside_sage)))
        timeprint("type of SC_vector_filtered", str(type(SC_vector_filtered)))
        assert leftside_sage == SC_vector_filtered

    else:
        # Usual path
        leftside = sol_f * MC
        timeprint("type of leftside:", str(type(leftside)))
        assert sol_f * MC == SC_vector_filtered

    end_solf = time.time()
    params.timing["sol_f computations"] = end_solf-start_solf

    start_filtering = time.time()

    row_to_aquery = MM.row_to_aquery
    ab_per_row = MM.ab_per_row(params)

    aquery_to_row = {q:r for r,q in enumerate(row_to_aquery)}

    sol = vector(Integers(e), len(row_to_aquery))
    for row_in_filtered_M, (coeff, abs) in enumerate(zip(sol_f, ab_per_row)):
        for ab in abs.keys():
            c = abs[ab]
            a = ab[0]
            b = ab[1]
            row_in_og_M = aquery_to_row[(a,b)]
            sol[row_in_og_M] += (c * coeff)

    for (coeff, killer_rel) in saved_kr_rels:
        a = int(killer_rel.split(":")[0].split(",")[0], 16)
        b = int(killer_rel.split(":")[0].split(",")[1], 16)
        r = aquery_to_row[(a,b)]
        sol[r] += coeff

    end_filtering = time.time()
    params.timing["recover sol from sol_f"] = end_filtering - start_filtering
    return sol

@timing
def solve_filtering_plus_bwc_chars(params, M, MM, MC, SC_vector, S_block, C_block, use_existing_M_binary=False):

    # Notes: S_block is the character portion of the matrix. C_block is the
    # character portion of the target vector, BEFORE filtering. The BWC-output
    # matrix (the ...matrix.bin file) does NOT include character columns.

    e = params.parameters['e']

    killer_relations = fast_json_load(params.files['KILLER_RELS_DICT'])
    order_of_cancellation = fast_json_load(params.files['KILLER_RELS_ORDER'])

    SC_vector_filtered, saved_kr_rels = apply_filtering_to_target(
        params, M, MM, MC, -SC_vector, killer_relations, order_of_cancellation
    )

    timeprint("Updated target! Created SC_vector_filtered.")
    timeprint("The rhs vector is of length " + str(len(SC_vector_filtered)))

    matrix_wchar_file = params.dirs['TEMP_OUTPUT_DIR'] + "matrix_wchar.bin"

    if not use_existing_M_binary:
        # Writing this file is super slow, so don't repeat if it exists
        with open(matrix_wchar_file, "wb") as f:
            timeprint("Writing matrix file...")
            for i in range(MC.nrows()):
                nz = MC.row(i).nonzero_positions()
                f.write(int.to_bytes(int(len(nz)), length=4, byteorder='little'))
                for j in nz:
                    f.write(int.to_bytes(j, length=4, byteorder='little'))
                    f.write(int.to_bytes(int(MC[i,j]), length=4, byteorder='little', signed=True))

        timeprint("Wrote matrix with chars file!")

    CadoNFS("linalg/bwc/mf_scan2", "-withcoeffs", "-mfile", matrix_wchar_file)

    timeprint(f"Matrix file MM ({type(MM)}) is {MM.get_matrix_file(params)} of size {MM.nrows}x{MM.ncols}")
    timeprint(f"Matrix MC ({type(MC)}) of size {MC.nrows()}x{MC.ncols()}")
    timeprint(f"Matrix M ({type(M)}) of size {M.nrows()}x{M.ncols()}")
    vector_file = params.dirs['TEMP_OUTPUT_DIR'] + "tgt.ascii"

    assert(MC.ncols()==len(SC_vector_filtered))

    write_ascii_vector(vector_file, SC_vector_filtered, e)

    FILES_DIR = params.dirs['TEMP_OUTPUT_DIR']
    WORKDIR = os.path.join(FILES_DIR,"bwc/")
    e = int(params.parameters['e'])
    Path(WORKDIR).mkdir(exist_ok=True)
    make_and_clean(WORKDIR)

    bwc_add_args = []
    if params.mpi:
        timeprint("Run bwc.pl from solve_filtering_plus_bwc_chars with MPI enabled.")
        bwc_add_args += [
            "--mpi", params.parameters['bwc.mpi']
        ]

        # Escape spaces in arguments correctly, depending on whether the command is
        # written to a file (slurm), or executed with subprocess (non-slurm).
        if params.bwc_slurm:
            bwc_add_args.append(f"mpi_extra_args='{params.mpi_extra_args}'")
        else:
            bwc_add_args.append(f"mpi_extra_args={params.mpi_extra_args}")

    if params.parameters['bwc.numsols']:
        bwc_mn = int(params.parameters['bwc.numsols'])
    else:
        bwc_mn = 1

    with open(os.path.join(WORKDIR,"bwc.out"),"w") as outfile, \
         open(os.path.join(WORKDIR,"bwc.err"),"w") as errfile:

        print(f"Writing bwc stdout to {outfile.name}")
        print(f"Writing bwc stderr to {errfile.name}")

        CadoNFS("linalg/bwc/bwc.pl",
            ":complete",
            "--matrix", os.path.realpath(matrix_wchar_file),
            "--prime", params.parameters['e'],
            "--nullspace", "LEFT",
            "--rhs", os.path.realpath(vector_file),
            "--wdir", WORKDIR,
            "--solutions", f"0-{bwc_mn}",
            f"mn={bwc_mn}",
            "balancing_options=reorder=columns",
            "verbose_flags=^all-cmdline,^bwc-timing-grids,^all-bwc-dispatch,^bwc-cache-major-info,perl-cmdline,perl-sections,^perl-checks,^bwc-iteration-timings,^bwc-loading-mksol-files",
            "--thr", params.parameters["bwc.thr.controller"],
            *bwc_add_args,
            capture=outfile,
            stderr=errfile
        )

    Ze = Integers(e)

    found = False
    for j in range(bwc_mn):
        if found:
            break
        with open(f"{WORKDIR}/K.sols{j}-{j+1}.0.txt",'r') as f:
            lines = [int(line) for line in f]
            sol_f = vector(Ze, lines, sparse=True)
            if sol_f[-1] != 0:
                found = True
            else:
                timeprint("Found a solution... with zero coefficient on the target.")

    if sol_f[-1] == 0:
        raise RuntimeError("Bwc didn't find any solutions with nonzero coefficient on the target.")

    timeprint("Found sol_f! (The filtered solution.)")
    timeprint(len(sol_f))

    assert sol_f[:-1] * MC == -1 * sol_f[-1] * SC_vector_filtered

    sol_ff = sol_f[:-1]
    for i in range(len(sol_ff)):
        sol_ff[i] = sol_ff[i] * inverse_mod(-1 * sol_f[-1], e)

    row_to_aquery = MM.row_to_aquery
    ab_per_row = MM.ab_per_row(params)

    aquery_to_row = {q:r for r,q in enumerate(row_to_aquery)}

    sol = vector(Integers(e), len(row_to_aquery))
    for row_in_filtered_M, (coeff, abs) in enumerate(zip(sol_ff, ab_per_row)):
        for ab in abs.keys():
            c = abs[ab]
            a = ab[0]
            b = ab[1]
            row_in_og_M = aquery_to_row[(a,b)]
            sol[row_in_og_M] += (c * coeff)

    for (coeff, killer_rel) in saved_kr_rels:
        a = int(killer_rel.split(":")[0].split(",")[0], 16)
        b = int(killer_rel.split(":")[0].split(",")[1], 16)
        r = aquery_to_row[(a,b)]
        sol[r] += coeff

    return sol

@timing
def CRT_R_alpha(params, current_num_primes, extra_abks, abm_list):
    e = params.parameters['e']
    N = params.poly.N
    m = params.poly.m
    f = params.poly.f[1]
    ZN = Integers(N)
    K = params.poly.K[1]
    alpha = K.gen()

    R_m = ZN(0)
    R_alpha = 0

    ts = time.time()

    if params.mpi:
        from mpi4py.futures import MPIPoolExecutor
        from mpi4py.futures import wait as mpi_wait

        with MPIPoolExecutor(max_workers=params.parameters['mpi.thr']) as executor:
            proccount = executor.num_workers
            timeprint(f"Running find_good_primes on {proccount} workers.")

            futures = [executor.submit(find_good_primes,f,e,num=current_num_primes/proccount,start=2**41+ZZ.random_element(2**41)) for i in range(proccount)]
            mpi_wait(futures)
            prime_list = [p for future in futures for p in future.result()]
            prime_list = list(set(prime_list))

            timeprint("prime_gen took",time.time()-ts)
            ts = time.time()

            proccount = executor.num_workers
            timeprint(f"Running ethroot_p on {proccount} workers.")

            futures = [executor.submit(ethroot_p,abm_list,f,e,N,p) for p in prime_list]
            mpi_wait(futures)
            local_solns = [future.result() for future in futures]
    else:
        cpucount = multiprocessing.cpu_count()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(find_good_primes,f,e,num=current_num_primes/cpucount,start=2**41+ZZ.random_element(2**41)) for i in range(cpucount)]
            concurrent.futures.wait(futures)
            prime_list = [p for future in futures for p in future.result()]
            prime_list = list(set(prime_list))

        timeprint("prime_gen took",time.time()-ts)

        ts = time.time()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = [executor.submit(ethroot_p,abm_list,f,e,N,p) for p in prime_list]
            concurrent.futures.wait(futures)
            local_solns = [future.result() for future in futures]

    timeprint("local_solns/ethroot_gen took",time.time()-ts)
    ps, coeffs_mod_ps = zip(*local_solns)
    ps = list(ps)
    prodps = prod(ps)
    coeffs_mod_ps = list(coeffs_mod_ps)
    coeffs_by_position = list(itertools.zip_longest(*coeffs_mod_ps, fillvalue=0)) # [[x^0 coeff mod each p], [x^1 coeff mod each p], ...]
    coeffs_over_z = [crt(list(coeffs_i), ps) for coeffs_i in coeffs_by_position]
    coeffs_over_z = [(c - prodps) if c >= prodps//2 else c for c in coeffs_over_z]

    cand_R_alpha = ZZ['x'](coeffs_over_z)(alpha)
    cand_R_alpha /= prod([(a-b*alpha)**k for (a,b,k) in extra_abks])

    return cand_R_alpha

def CRT_R_m(params, *args):
    N = params.poly.N
    m = params.poly.m
    ZN = Integers(N)
    cand_R_alpha = CRT_R_alpha(params, *args)

    return cand_R_alpha.polynomial().change_ring(ZN)(m)

def handle_very_large_special_q(params, q, side, working_pfx, rho=None):
    # This is part of descent initialization, if we are in a parameter regime
    # where we have special-q's over ~100 bits.

    og_polyfile = params.files['POLYFILE']
    og_poly = CadoPolyFile(og_polyfile); og_poly.read()

    m = og_poly.m
    f = og_poly.f[1]
    g = og_poly.f[0]

    q_nickname = str(q)[0:10]
    mod_polyfile = working_pfx + "largeq." + q_nickname + ".poly"
    mod_fbfile = working_pfx + "largeq." + q_nickname + ".fb"
    las_output = working_pfx + "largeq." + q_nickname + ".rels.out"

    if rho is None and side == 0:
        rho = ZZ(g.roots(GF(q))[0][0])
    else:
        raise NotImplementedError("handle_very_large_special_q on side=1")

    newg, newf, coeff = transform_polys_by_q(g, f, q, rho, side)

    rr = Integer(newg.resultant(newf))

    assert(newg.degree() == 1)
    c0 = Integer(newg.list()[0])
    c1 = Integer(newg.list()[1])

    if gcd(c1, rr) > 1:
        print("Hm! The new modulus rr has some nontrivial factors so we will try to get rid of them.")
        major_message("Not actually sure if this is allowed... If something breaks this could be why!")

    extra = gcd(c1, rr)
    rr = Integer(rr / extra)
    if rr < 0:
        rr = rr * -1

    looking = True
    while looking:
        res = ecmfactor(rr, 2**20)
        if res[0]:
            print("found another divisor of rr")
            major_message("Not actually sure if this is allowed... If something breaks this could be why!")
            rr = rr / res[1]
        else:
            looking = False

    new_shared_root = -1 * c0 * inverse_mod(c1, rr) % rr

    print("newg", str(newg))
    print("newf", str(newf))
    print("new_shared_root", str(new_shared_root))
    print("rr", str(rr))

    assert(newg(new_shared_root) % rr == 0)
    assert(newf(new_shared_root) % rr == 0)

    if True:
        write_polyfile(
            params,'hardcoded',rr, params.parameters['e'],
            newf.degree(), mod_polyfile, hardcoded_m=new_shared_root, hardcoded_polys=[newg,newf]
        )

        polyinfo = CadoPolyFile(mod_polyfile)
        polyinfo.read()
        minskew = float(params.parameters.get('desc.lq.skewmin', 0.2))

        if float(polyinfo.skewness) < minskew:
            # Throw it out.
            # Practically, these end up discarding every q.
            timeprint("Throwing out because of skewness " + str(polyinfo.skewness))
            sys.exit(1)

        CadoNFS("sieve/makefb",
                "-poly", 'POLY',
                "-out", 'FB',
                "-lim", str(params.parameters.get('desc.lq.lim', 2000000000)),
                "-maxbits", str(params.parameters.get('desc.lq.maxbits', 16)),
                "-side", "1",
                "-t", str(params.parameters.get('nthreads', 8)),
                inputs={ 'POLY': mod_polyfile, },
                outputs={
                    'FB': mod_fbfile,
                    }
                )

        lim0 = params.parameters.get('desc.lq.lim', 2000000000)
        lim1 = params.parameters.get('desc.lq.lim', 2000000000)
        bkthresh = int(round(min(lim0, lim1)/2))

        CadoNFS("sieve/las",
            "-poly", 'POLY',
            "-fb1", 'FB1',
            "-out", 'OUTFILE',
            "-lim0", str(lim0),
            "-lim1", str(lim1),
            "-lpb0", str(params.parameters.get('desc.lq.lpb0', 80)),
            "-lpb1", str(params.parameters.get('desc.lq.lpb1', 90)),
            "-mfb0", str(params.parameters.get('desc.lq.mfb0', 80)),
            "-mfb1", str(params.parameters.get('desc.lq.mfb1', 300)),
            "-ncurves0", str(params.parameters.get('desc.lq.ncurves0', 20)),
            "-ncurves1", str(params.parameters.get('desc.lq.ncurves1', 150)),
            "-I", str(params.parameters.get('desc.lq.I_sieving', 16)),
            "-sqside", "0",
            "-nq", str(params.parameters.get('desc.lq.nq', 1000)),
            "-q0", str(params.parameters.get('desc.lq.q0', 30)),
            "-t", str(params.parameters.get('desc.lq.nthreads', 8)),
            "-exit-early", "2",
            "--never-discard",
            "-bkthresh1", str(bkthresh),
            "-bkmult", "1.216",
            "--memory-margin", params.parameters.get('las.memory_margin', 20),
            outputs={'OUTFILE': las_output},
            inputs={
                'FB1': mod_fbfile,
                'POLY': mod_polyfile,
            }
            )

    winner = ''
    with open(las_output, "r") as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            else:
                winner = line
                break

    print("winning relation", str(winner))
    ijf = winner.strip().split(':')
    i = ijf[0].split(',')[0]
    i = int(i, 10)
    j = ijf[0].split(',')[1]
    j = int(j, 10)
    fac0 = [ int(x, 16) for x in ijf[1].split(',') ]
    fac1 = [ int(x, 16) for x in ijf[2].split(',') ]

    print("fac0", str(fac0))
    print("fac1", str(fac1))

    a0 = coeff[0][0]
    b0 = coeff[0][1]
    a1 = coeff[1][0]
    b1 = coeff[1][1]
    a = a0*i+a1*j
    b = b0*i+b1*j
    if b < 0:
        a = -a
        b = -b

    return a, b, fac0, fac1, winner, newg, newf, rr, new_shared_root

def call_descent_large_init(target,params,seedval):
    seedval = int(seedval)
    completed_init = subprocess.run([
        "time", "-p",
        params.files['SAGE'],
        "descent_large_init.py",
        "--params", params.files['PARAMS'],
        "--seed", str(seedval)
    ], stderr=subprocess.PIPE, text=True)
    overall_cputime.add(extract_time(completed_init.stderr))

    print(f"yay finished with seedval {seedval}")
    print("returncode",completed_init.returncode)

    # TODO: check for errors, have returncodes?

    return None

@timing
def call_descent_large_q_slurm(params, initdatafile):
    if initdatafile == 'ALL':
        # A special argument to say that we want to launch jobs for all qs,
        # in all initdata files matching the expected regex.
        filelist = glob.glob(params.dirs['DESC'] + 'ecminit.*.initdata')
    else:
        filelist = [ initdatafile ]

    processes = []
    jobnum = 0

    for initdata in filelist:

        try:
            with open(initdata, "r") as fp:
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
            print(f"Uh oh! There was a problem parsing the init_data in {initdata}.")
            raise ex
            sys.exit(1)
            os._exit(1)

        large_q = []
        for ff in u_fac + v_fac:
            ffi = Integer(ff)
            if ffi.nbits() > params.parameters.get('LARGEQ', 90):
                large_q.append(ffi)

        # TODO: Have some guardrail for the number of slurm jobs?
        # But it's usually < 10, for each initdata

        for lq in large_q:
            jobnum += 1
            command_list = [
                params.files['SAGE'],
                "descent_large_q_helper.py",
                "--params", params.files['PARAMS'],
                "--existing-init-data", initdata,
                "--lq", str(lq)
            ]
            processes.append(slurmit(params, " ".join(command_list), params.prefix[:-1]+"-largeq", jobnum))

    finished_processes, cputime_slurm = slurm_wait(processes)
    overall_cputime.add(cputime_slurm)
    print(f"Finished running {jobnum} large qs!")

    # TODO:
    # It's pretty tedious to check on the status of this, and even to know when to stop.
    # What we want is, for at least one initdata, all the associated special-q slurm jobs should
    # have completed successfully.
    # But a lot of the special-q jobs fail because the modulus is very smooth. We might want a
    # couple successful runs also because there can still be problems in linear algebra or with
    # unlinked ideals.
    # The jobs that succeed tend to be very fast, and some jobs that fail tend to hang for a
    # long time.
    # Alternatively, this can be run on one initdata file at a time and then it will at least be
    # obvious when one succeeds.

    good_descent_inits = []

    for initdata in filelist:
        with open(initdata, "r") as fp:
            init_dict = json.load(fp)

        u_fac = init_dict['u_fac']
        v_fac = init_dict['v_fac']
        largeq_rels_file = init_dict['largeq_rels_file']

        if not os.path.exists(largeq_rels_file):
            # no large-q relations found
            continue

        large_q = []
        for ff in u_fac + v_fac:
            ffi = Integer(ff)
            if ffi.nbits() > params.parameters.get('LARGEQ', 90):
                large_q.append(ffi)

        found_rels = ""
        with open(largeq_rels_file, "r") as f:
            for line in f.readlines():
                if "Taking" in line:
                    found_rels += line

        good = True

        for lq in large_q:
            if ("q=" + str(lq)) not in found_rels:
                good = False

        if good:

            todofile = init_dict['todofilename']
            todo_list_0 = []
            todo_list_1 = []

            with open(todofile) as td:
                for line in td.readlines():
                    line = line.strip()
                    lineinfo = line.split()
                    if lineinfo[0] == '0':
                        qsize = Integer(lineinfo[1]).nbits()
                        todo_list_0.append(qsize)
                    elif lineinfo[0] == '1':
                        qsize = Integer(lineinfo[1]).nbits()
                        todo_list_1.append(qsize)

            with open(params.files['GOOD_DESCENT_INIT'], "a") as f:
                f.write(str(initdata) + "\n")
                f.write(f"Total todo qs: {len(todo_list_0) + len(todo_list_1)}\n")
                f.write(f"With {sum(todo_list_0)} side-0 bits and {sum(todo_list_1)} side-1 bits.\n")
                todo_list_0.sort()
                todo_list_1.sort()
                f.write("0: " + str(todo_list_0) + "\n")
                f.write("1: " + str(todo_list_1) + "\n")
                f.write("\n")

    return True

def do_cleanup_large_q(params):
    filelist = glob.glob(params.dirs['DESC'] + 'ecminit.*.initdata')
    good_descent_inits = []

    for initdata in filelist:
        with open(initdata, "r") as fp:
            init_dict = json.load(fp)

        u_fac = init_dict['u_fac']
        v_fac = init_dict['v_fac']
        largeq_rels_file = init_dict['largeq_rels_file']

        if not os.path.exists(largeq_rels_file):
            # no large-q relations found
            continue

        large_q = []
        for ff in u_fac + v_fac:
            ffi = Integer(ff)
            if ffi.nbits() > params.parameters.get('LARGEQ', 90):
                large_q.append(ffi)

        found_rels = ""
        with open(largeq_rels_file, "r") as f:
            for line in f.readlines():
                if "Taking" in line:
                    found_rels += line

        good = True

        for lq in large_q:
            if ("q=" + str(lq)) not in found_rels:
                good = False

        if good:
            with open(params.files['GOOD_DESCENT_INIT'], "a") as f:
                f.write(str(initdata) + "\n")


@timing
def handle_bottom_special_q(
    params, q, side, working_pfx, rho=None, hintinfo=None, overwrite_LPB0=None, overwrite_LPB1=None):
    # Here the q is very close to LPB{0,1}.

    og_polyfile = params.files['POLYFILE']
    og_poly = CadoPolyFile(og_polyfile); og_poly.read()

    f = og_poly.f[1]
    g = og_poly.f[0]

    if side == 0:
        q_nickname = 'side' + str(side) + '-' + str(q)
    elif side == 1:
        q_nickname = 'side' + str(side) + '-' + str(q) + '-' + str(rho)

    mod_polyfile = working_pfx + "botq." + q_nickname + ".poly"
    mod_fbfile = working_pfx + "botq." + q_nickname + ".fb"
    las_output = working_pfx + "botq." + q_nickname + ".rels.out"

    if rho is None and side == 0:
        rho = ZZ(g.roots(GF(q))[0][0])

    newg, newf, coeff = transform_polys_by_q(g, f, q, rho, side)

    print("coeff from myxgcd")
    for w in coeff:
        for ww in w:
            ww = Integer(ww)
            print(str(ww) + " which is " + str(ww.nbits()) + " bits")

    rr = Integer(newg.resultant(newf))

    assert(newg.degree() == 1)
    c0 = Integer(newg.list()[0])
    c1 = Integer(newg.list()[1])

    if gcd(c1, rr) > 1:
        print("Hm! The new modulus rr has some nontrivial factors so we will try to get rid of them.")
        major_message("Not actually sure if this is allowed... If something breaks this could be why!")

    extra = gcd(c1, rr)
    rr = Integer(rr / extra)
    if rr < 0:
        rr = rr * -1

    looking = True
    while looking:
        res = ecmfactor(rr, 2**20)
        if res[0]:
            print("found another divisor of rr")
            major_message("Not actually sure if this is allowed... If something breaks this could be why!")
            rr = rr / res[1]
        else:
            looking = False

    new_shared_root = -1 * c0 * inverse_mod(c1, rr) % rr

    print("newg", str(newg))
    print("newf", str(newf))
    print("new_shared_root", str(new_shared_root))
    print("rr", str(rr))

    assert(newg(new_shared_root) % rr == 0)
    assert(newf(new_shared_root) % rr == 0)

    if not os.path.exists(mod_fbfile):
        write_polyfile(
            params,'hardcoded',rr, params.parameters['e'],
            newf.degree(), mod_polyfile, hardcoded_m=new_shared_root, hardcoded_polys=[newg,newf]
        )

        lim1 = params.parameters.get('desc.botq.lim1', 2**20)
        lim1 = min(2**31, lim1, 2**int(params.parameters['LPB1']))

        CadoNFS("sieve/makefb",
                "-poly", 'POLY',
                "-out", 'FB',
                "-lim", str(lim1),
                "-maxbits", str(params.parameters.get('desc.botq.maxbits', 16)),
                "-side", "1",
                "-t", str(params.parameters.get('nthreads', 8)),
                inputs={ 'POLY': mod_polyfile, },
                outputs={
                    'FB': mod_fbfile,
                    }
                )

    if True:

        lim1 = params.parameters.get('desc.botq.lim1', 2**20)
        lim1 = min(2**31, lim1, 2**int(params.parameters['LPB1']))

        lim0 = params.parameters.get('desc.botq.lim0', 2**20)
        lim0 = min(2**31, lim0, 2**int(params.parameters['LPB0']))

        bkthresh_t = int(round(min(lim0, lim1) / 2))
        bkthresh = params.parameters.get('desc.botq.bkthresh1', bkthresh_t)
        bkmult = params.parameters.get('desc.botq.bkmult', '1.365')

        A_used = params.parameters.get('desc.botq.A_sieving', 32)
        if int(A_used) > 32 and 'las.bigA.hwloc_job_binding_policy' in params.parameters:
            t_used = params.parameters['las.bigA.hwloc_job_binding_policy']
        else:
            t_used = params.las_job_binding_policy

        t_used = params.parameters['desc.comp.thr']

        if hintinfo is not None:
            this_lpb0, this_mfb0, this_lpb1, this_mfb1 = hintinfo
        else:
            this_lpb0 = params.parameters['LPB0']
            this_lpb1 = params.parameters['LPB1']
            this_mfb0 = params.parameters.get(f'desc.botq.mfb0.{side}', 80)
            this_mfb1 = params.parameters.get(f'desc.botq.mfb1.{side}', 300)

        if overwrite_LPB0 is not None:
            this_lpb0 = int(overwrite_LPB0)
        if overwrite_LPB1 is not None:
            this_lpb1 = int(overwrite_LPB1)

        CadoNFS("sieve/las",
            "-poly", 'POLY',
            "-fb1", 'FB1',
            "-out", 'OUTFILE',
            "-lim0", str(lim0),
            "-lim1", str(lim1),
            "-lpb0", str(this_lpb0),
            "-lpb1", str(this_lpb1),
            "-mfb0", str(this_mfb0),
            "-mfb1", str(this_mfb1),
            "-ncurves0", str(params.parameters.get(f'desc.botq.ncurves0.{side}', 20)),
            "-ncurves1", str(params.parameters.get(f'desc.botq.ncurves1.{side}', 150)),
            "-A", str(A_used),
            "-sqside", str(params.parameters.get(f'desc.botq.sqside.{side}', 0)),
            # should be easier on side 0 (?), but could be either
            "-nq", str(params.parameters.get('desc.botq.nq', 1000)),
            "-q0", str(params.parameters.get(f'desc.botq.q0.{side}', 30)),
            "-t", t_used,
            "-B", str(16),
            "--adjust-strategy", "2",
            "-exit-early", "2",
            "-bkmult", str(bkmult),
            "-bkthresh1", str(bkthresh),
            "-v",
            "--never-discard",
            "--memory-margin", str(params.parameters.get('desc.botq.memory_margin', 100)),
            outputs={'OUTFILE': las_output},
            inputs={
                'FB1': mod_fbfile,
                'POLY': mod_polyfile,
            }
            )

    winner = ''
    with open(las_output, "r") as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            else:
                winner = line
                break

    print("winning relation", str(winner))
    ijf = winner.strip().split(':')
    i = ijf[0].split(',')[0]
    i = int(i, 10)
    j = ijf[0].split(',')[1]
    j = int(j, 10)
    fac0 = [ int(x, 16) for x in ijf[1].split(',') ]
    fac1 = [ int(x, 16) for x in ijf[2].split(',') ]

    print("fac0", str(fac0))
    print("fac1", str(fac1))

    a0 = coeff[0][0]
    b0 = coeff[0][1]
    a1 = coeff[1][0]
    b1 = coeff[1][1]
    a = a0*i+a1*j
    b = b0*i+b1*j
    if b < 0:
        a = -a
        b = -b

    return a, b, fac0, fac1, winner, newg, newf, rr, new_shared_root


@timing
def handle_bottom_special_q_composites(
    params, q, side, working_pfx, rho=None, strategy='C', overwrite_LPB0=None, overwrite_LPB1=None, overwrite_mfb=None):
    # Here the q is very close to LPB{0,1}.

    assert strategy in ['C','C1','C2']

    if side == 0:
        q_nickname = 'side' + str(side) + '-' + str(q)
    elif side == 1:
        q_nickname = 'side' + str(side) + '-' + str(q) + '-' + str(rho)

    las_output = working_pfx + "compq." + q_nickname + ".rels.out"
    todofile = working_pfx + "compq." + q_nickname + ".todo"

    og_poly = CadoPolyFile(params.files['POLYFILE']); og_poly.read()
    f = og_poly.f[1]
    g = og_poly.f[0]

    starting_q = q

    if strategy in ['C','C1']:
        # Just one prime multiplier
        nprimes = params.parameters.get(f'desc.compq.{side}.nprimes', 150)
        if strategy == 'C':
            primeset = primes_first_n(nprimes)

        elif strategy == 'C1':
            # Allow composite multipliers (yes, the naming is confusing)
            # side 0 will remain the same, but side 1 will have some new logic

            primeset = []
            cand_primes = primes_first_n(50)[10:]
            multipliers = primes_first_n(50)[10:]

            for cc in cand_primes:
                for mm in multipliers:
                    primeset.append( cc*mm )

        if side == 0:
            with open(todofile, "w") as f:
                # Always include 1*q
                root1 = g.roots(Integers(q), multiplicities=False)
                if len(root1) == 1 and root1[0] > 0:
                    f.write(f"0 {q} {root1[0]}\n")

                for qq in primeset:

                    skip_this = False
                    if strategy == 'C1':
                        facs = list(factor(qq))
                        for fac in facs:
                            if fac[1] > 1:
                                skip_this = True

                    if skip_this:
                        continue

                    Zqqq = Integers(qq*q)
                    roots = g.roots(Zqqq, multiplicities=False)
                    if len(roots) == 1 and roots[0] > 0:
                        f.write(f"0 {qq*q} {roots[0]}\n")

        if side == 1 and strategy == 'C':
            rho = Integer(rho)
            with open(todofile, "w") as ff:
                # Always include 1*q
                ff.write(f"1 {q} {rho}\n")

                for qq in primeset:
                    roots = f.roots(GF(qq))
                    if len(roots) > 0:
                        rho_qq = roots[0][0]
                        try:
                            rho_qqq = crt(Integer(rho), Integer(rho_qq), Integer(q), Integer(qq))
                            assert rho_qqq % q == rho
                            assert rho_qqq % qq == rho_qq
                            ff.write(f"1 {qq*q} {rho_qqq}\n")
                        except ValueError:
                            continue    # sometimes crt can't work

        if side == 1 and strategy == 'C1':
            rho = Integer(rho)
            with open(todofile, "w") as ff:
                # Always include 1*q
                ff.write(f"1 {q} {rho}\n")

                for possibly_composite in primeset:
                    try:
                        facs = list(factor(possibly_composite))
                        rhos = [ Integer(rho) ]
                        mods = [ Integer(q) ]
                        skip_this = False
                        for fac in facs:
                            if fac[1] > 1:
                                skip_this = True
                                raise ValueError
                            roots_fac = f.roots(GF(fac[0]))
                            if len(roots_fac) > 0:
                                rho_fac = roots_fac[0][0]
                                rhos.append(Integer(rho_fac))
                                mods.append(Integer(fac[0]))
                            else:
                                skip_this = True
                                raise ValueError

                        if not skip_this:
                            assert Integer(possibly_composite*q) == prod(mods)
                            rho_total = CRT_list(rhos, mods)
                            for ii in range(len(rhos)):
                                assert rho_total % mods[ii] == rhos[ii]

                            ff.write(f"1 {possibly_composite*q} {rho_total}\n")
                    except ValueError:
                        continue

    if strategy == 'C2':
        nprimes = params.parameters.get(f'desc.compq.{side}.nprimes', 150)
        primesize = int(params.parameters.get(f'desc.compq.{side}.randsize', 15))
        wsize = floor(primesize/2)
        zsize = primesize - wsize

        primeset = []
        cand_primes = primes_first_n(300)[240:]
        multipliers = primes_first_n(30)

        for mm in multipliers:
            for cc in cand_primes:
                if mm != cc:
                    primeset.append((mm,cc))

        if False:
            for _ in range(nprimes):
                #primeset.append( (random_prime(2**wsize), random_prime(2**zsize)) )
                w = random_prime(2**wsize)
                z = random_prime(2**zsize)
                if (w,z) not in primeset and (z,w) not in primeset and w != z:
                    primeset.append((w,z))

        if side == 0:
            with open(todofile, "w") as f:
                for (w,z) in primeset:
                    if w == z:
                        continue
                    Zqwz = Integers(q*w*z)
                    roots = g.roots(Zqwz, multiplicities=False)
                    if len(roots) == 1 and roots[0] > 0:
                        f.write(f"0 {q*w*z} {roots[0]}\n")

        if side == 1:
            rho = Integer(rho)
            with open(todofile, "w") as ff:
                for (w,z) in primeset:
                    if w == z:
                        continue
                    roots_w = f.roots(GF(w))
                    roots_z = f.roots(GF(z))
                    if len(roots_w) > 0 and len(roots_z) > 0:
                        rho_w = roots_w[0][0]
                        rho_z = roots_z[0][0]
                        try:
                            rho_qwz = CRT_list(
                                        [Integer(rho),Integer(rho_w),Integer(rho_z)],
                                        [Integer(q),Integer(w),Integer(z)])
                            assert rho_qwz % q == rho
                            assert rho_qwz % w == rho_w
                            assert rho_qwz % z == rho_z
                            ff.write(f"1 {q*w*z} {rho_qwz}\n")
                        except ValueError:
                            continue    # sometimes crt can't work

    lim0 = params.parameters.get(f'desc.compq.{side}.lim0', 2**20)
    lim1 = params.parameters.get(f'desc.compq.{side}.lim1', 2**20)
    mfb0 = params.parameters.get(f'desc.compq.{side}.mfb0', 120)
    mfb1 = params.parameters.get(f'desc.compq.{side}.mfb1', 120)
    ncurves0 = params.parameters.get(f'desc.compq.{side}.ncurves0', 100)
    ncurves1 = params.parameters.get(f'desc.compq.{side}.ncurves1', 100)
    A_s = params.parameters.get(f'desc.compq.{side}.A', 33)

    lim0 = min(lim0, 2**int(params.parameters['LPB0']))
    lim1 = min(lim1, 2**int(params.parameters['LPB1']))
    bkthresh1 = int(round(min(lim0, lim1) / 2))

    bkthresh1 = params.parameters.get(f'desc.compq.{side}.bkthresh1', bkthresh1)

    if overwrite_LPB0 is not None:
        used_lpb0 = int(overwrite_LPB0)
    else:
        used_lpb0 = params.parameters['LPB0']

    if overwrite_LPB1 is not None:
        used_lpb1 = int(overwrite_LPB1)
    else:
        used_lpb1 = params.parameters['LPB1']

    if overwrite_mfb is not None:
        used_mfb0 = int(overwrite_mfb[0])
        used_mfb1 = int(overwrite_mfb[1])
    else:
        used_mfb0 = mfb0
        used_mfb1 = mfb1

    memm = str(params.parameters.get(f'desc.compq.{side}.memory_margin', 100))
    bkmult = params.parameters.get(f'desc.compq.{side}.bkmult', 1.3)

    #if not os.path.exists(las_output):
    # Note: Often requires SUPPORT_LARGE_Q flag
    CadoNFS("sieve/las",
        "-poly", 'POLY',
        "-fb1", 'FB1',
        "-out", 'OUTFILE',
        "-todo", 'TODOFILE',
        "-lim0", str(lim0),
        "-lim1", str(lim1),
        "-lpb0", str(used_lpb0),
        "-lpb1", str(used_lpb1),
        "-mfb0", str(used_mfb0),
        "-mfb1", str(used_mfb1),
        "-ncurves0", str(ncurves0),
        "-ncurves1", str(ncurves1),
        "-A", str(A_s),
        "-t", params.parameters.get('desc.comp.thr', 8),
        "-B", str(16),
        "--adjust-strategy", "2",
        #"-exit-early", "2",        # to handle missing ideals
        "--memory-margin", memm,
        "--allow_compsq",
        "--allow-largesq",
        "-v",
        "-sqside", str(side),
        "-bkmult", str(bkmult),
        "-bkthresh1", str(bkthresh1),
        outputs={'OUTFILE': las_output},
        inputs={
            'FB1': params.files['CAPPED_FBGZ'],
            'POLY': params.files['POLYFILE'],
            'TODOFILE': todofile,
        }
        )

    # TEMPORARY EDIT (along with exit-early change)
    ext_unlinked_ideals = set()
    ext_unlinked_files = [
        # empty for now
    ]
    for filename in ext_unlinked_files:
        with open(filename, 'r') as f:
            for line in f:
                lineinfo = line.strip().split()
                assert lineinfo[0] == '0'       # 1 sided sieving
                p = Integer(lineinfo[1])
                r = Integer(lineinfo[2])
                ext_unlinked_ideals.add( (1,p,r) )

    #with open(params.files['EXT_UNLINKED_TODOS'], "r") as extfile:
    #    for line in extfile.readlines():
    #        # it's a todo file
    #        lineinfo = line.strip().split(" ")
    #        assert lineinfo[0] == '0'
    #        q = Integer(lineinfo[1])
    #        r = Integer(lineinfo[2])
    #        ext_unlinked_ideals.add((1,q,r))

    latest_rel_line = ''
    all_rel_lines = []

    with open(las_output, "r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                latest_rel_line = line
                #break
                all_rel_lines.append(line)

    if latest_rel_line != '':
        for rel_line in all_rel_lines:
            line = rel_line.strip().split(":")
            ab = line[0].split(",")
            a = int(ab[0])
            b = int(ab[1])
            fac0_hex = line[1].split(",")
            fac1_hex = line[2].split(",")

            a = Integer(a)
            b = Integer(b)

            fac0_int = [ Integer(y,16) for y in fac0_hex ]
            fac1_int = [ Integer(y,16) for y in fac1_hex ]

            good_rel = True

            for alg_fac in fac1_int:

                if int(alg_fac) == int(starting_q):
                    continue

                if gcd(b,alg_fac) != 1:
                    r = alg_fac
                else:
                    r = (a * inverse_mod(b,alg_fac)) % alg_fac

                if (1,alg_fac,r) in ext_unlinked_ideals:
                    # We have an unlinked ideal in our factorization
                    good_rel = False
                    timeprint(f"Seeing a relation that uses an unlinked ideal: {1},{alg_fac},{r}")
                    timeprint(f"relation: {rel_line}")
                if alg_fac > 2**used_lpb1:
                    # shouldn't happen, really
                    good_rel = False

            if good_rel:
                return a, b, fac0_hex, fac1_hex
            else:
                continue

    else:
        return -1, -1, [], []
