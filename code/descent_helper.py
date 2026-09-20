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
import os, sys
from helpers import silent_remove, masked_target, construct_S, CadoExplainRenumberFile, truncate_S, sanity_check_descent_outfile, generate_or_load_target
from contextlib import redirect_stdout,redirect_stderr
from candy import print_command_line,warning_message,major_message,error_message
from cado_sage import CadoPolyFile
from cado_nfs_binaries import CadoNFS,CadoNFSBinaries
import time
import descent_ecm_utils

from run import LuckySplitException


# XXX This is not unified with Params in run.py, which is unfortunate.
# This type definition is quite like argparse.Namespace, by the way.
class Params(dict):
    __getattr__ = dict.get

class ConstructSError(Exception):
    def __init__(self, missed, *args, **kwargs):
        super().__init__(*args, **kwargs)
    def __str__(self):
        return "Error in constructing S"

if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='descent_helper.py')
    topparser.add_argument('--params',dest='params')
    topparser.add_argument('--seed',dest='seedval')
    topparser.add_argument('--existing-init-data',dest='existing_init_data',required=False)
    topargs = topparser.parse_args()

    params = Params(json.loads(open(topargs.params,'r').read()))

    e = params.parameters['e']
    N = params.parameters['N']
    ZN = Integers(N)

    target = generate_or_load_target(params)

    # the initial u,v are both about half the modulus bits, and we'll ask
    # that they factor into primes that are smaller than that by a factor
    # u. The rule of thumb is that this means that we'll have to do
    # 1/dickman_rho(u)^2 trials before we find a winning candidate.

    ratio = params.parameters.get('initial_descent_smoothness_ratio', 2.25)
    ratio = float(ratio)
    initial_smoothness_maxbits = params.parameters['MODULUS_BITS'] / 2 / ratio

    #
    # with seed(topargs.seedval):
    #     for spin in itertools.count():
    #         mask, h, u, v = masked_target(target, e)

    #         print(ZZ(u)/ZZ(v))
    #         special_qs = [p for p,k in factor(ZZ(u)/ZZ(v)) if p > params.BOUNDR]

    #         if Integer(max(special_qs, default=0)).ndigits(2) < initial_smoothness_maxbits:
    #             break

    inside_parser = argparse.ArgumentParser(prog='scripts/descent.py')
    inside_parser.add_argument("--target",
                               help="Element whose DL is wanted",
                               type=str,
                               required=True)

    descent.GeneralClass.declare_args(inside_parser)
    descent.DescentUpperClass.declare_args(inside_parser)
    descent.DescentMiddleClass.declare_args(inside_parser)
    descent.DescentLowerClass.declare_args(inside_parser)

    cpu_count = multiprocessing.cpu_count()

    inside_inputs = [
        "--poly", params.files['POLYFILE'],
        "--fb1", params.files['FBFILE'],
        "--lim1", params.BOUNDA,
        "--lim0", params.BOUNDR,
        "--I", params.parameters['I_sieving'],
        "--B", 16,
        "--lpb0", params.parameters['LPB0'],
        "--lpb1", params.parameters['LPB1'],
        "--mfb0", params.parameters['desc.mfb0'],
        "--mfb1", params.parameters['desc.mfb1'],
        "--descent-hint", params.files['HINTFILE'],
        # XXX: '--threads' is passed down to las_descent,
        # it uses CADO_NFS_MAX_THREADS instead (set below).
        "--threads", cpu_count,
        "--cadobindir", params.dirs['CADO_BUILD_DIR'],
        "--datadir", params.dirs['DESC'],
        "--ell", params.parameters['N'],
        "--renumber", params.files['RENUMBERFILE'],
        "--prefix", "descent",
        "--no-logs",
    ]

    # All parameters of descent_upper_class can be passed here.

    par = params.parameters
    init_mfb = par.get('descent_init.mfb', int(par['MODULUS_BITS'] / 2 - 60))
    init_I   = par.get('descent_init.I',  14)
    init_lim = par.get('descent_init.lim', 2**26)
    init_lpb = par.get('descent_init.lpb', int(initial_smoothness_maxbits))
    init_tkewness = par.get('descent_init.tkewness', int(2**min(30, init_lpb-1)))

    # This is a local change to the environment variable and should not be reflected outside this script.
    if 'desc.thr' in par:
        timeprint(f"Setting environment variable CADO_NFS_MAX_THREADS to desc.thr from config file (used to set '-t' binding policy of las_descent).")
        os.environ["CADO_NFS_MAX_THREADS"] = str(par['desc.thr'])
    elif params.has_hwloc:
        timeprint(f"Setting environment variable CADO_NFS_MAX_THREADS=auto (used to set '-t auto' binding policy of las_descent)")
        os.environ["CADO_NFS_MAX_THREADS"] = "auto"
    else:
        timeprint(f"Setting environment variable CADO_NFS_MAX_THREADS={cpu_count} (used to determine threads of las_descent)")
        os.environ["CADO_NFS_MAX_THREADS"] = str(cpu_count)

    #
    # XXX it might make sense to also pass --init-lim. Currently we use
    # the default value of 2^26, which might be too large for small
    # examples, and perhaps too small for large ones...
    more_inside_inputs = [
        "--target", target,
        "--init-mfb", init_mfb,
        "--init-I", init_I,
        "--init-lpb", init_lpb,
        "--init-lim", init_lim,
        "--init-tkewness", init_tkewness
    ]

    # TODO: what is descent_upper_class's `--slaves`? Is it useful?

    args = inside_parser.parse_args([str(c)
                                     for c in inside_inputs + more_inside_inputs])
    general = descent.GeneralClass(args)

    timeprint("Args passed to CADO descent scripts:", args)

    init = descent.DescentUpperClass(general, args)

    assert target == general.target()

    init_starttime = time.time()

    if topargs.existing_init_data is not None and 'init' in topargs.existing_init_data:
        # Should be a legitimate file containing:
        # todofile (string filename)
        # u (Integer)
        # v (Integer)
        # u_fac (list of Integers)
        # v_fac (list of Integers)
        # mask (Integer)

        try:
            with open(topargs.existing_init_data, "r") as fp:
                init_dict = json.load(fp)

            todofilename = init_dict['todofilename']
            u = Integer(init_dict['u'])
            v = Integer(init_dict['v'])
            u_fac = init_dict['u_fac']
            v_fac = init_dict['v_fac']
            mask = Integer(init_dict['mask'])
            firstrelsfile = None    # Note this is always None anyway
            largeq_rels_file = init_dict['largeq_rels_file']

            init_data = (todofilename,
                        [u, v, u_fac, v_fac],
                        None,
                        mask)
        except Exception as ex:
            print(f"Uh oh! There was a problem parsing the init_data in {topargs.existing_init_data}.")
            raise ex
            sys.exit(1)
            os._exit(1)
    else:
        init_data = init.do_descent_for_real(int(target),
                                            topargs.seedval,
                                            randomize_multiplicatively=int(e))
        largeq_rels_file = None

    init_endtime = time.time()
    descent_init_time = round(init_endtime-init_starttime)

    # I don't know why this check isn't working.
    if not init_data:
        raise RuntimeError("init.do_descent_for_real returned None, quitting")

    todofile, uv_fac, firstrelsfile, mask = init_data

    if uv_fac is None:
        timeprint(init_data)
        raise RuntimeError(f"init.do_descent_for_real(target={target}, seed={topargs.seedval}) failed for some reason, quitting")

    todofile, (u, v, u_fac, v_fac), firstrelsfile, mask = init_data

    mask = ZN(mask)
    u = ZZ(u)
    v = ZZ(v)
    u_fac = Factorization([(ZZ(p),1) for p in u_fac])
    v_fac = Factorization([(ZZ(p),1) for p in v_fac])

    h = (pow(mask, e, N) * target) % N

    print("target", str(target))
    print("mask", str(mask))
    print("h", str(h))
    print("N", str(N))
    print("e",str(e))

    assert h == ZN(u/v)

    fac = Factorization(u_fac) / Factorization(v_fac)

    special_qs = [p for p,k in fac if p > params.BOUNDR]

    u_prefix = str(u)[:40]
    v_prefix = str(v)[:40]
    DESCENT_PREFIX = f"desc.{u_prefix}.{v_prefix}"

    for ext in ["todo", "err", "out", "badideals", "badidealinfo",
                "descent.tgt.middle.rels", "json"]:
        path = params.dirs['DESC']+DESCENT_PREFIX+'.'+ext
        if os.path.exists(path):
            warning_message(f"Removing {path}")
            os.unlink(path)

    path = params.files['TGT_INFO']+"."+DESCENT_PREFIX
    if os.path.islink(path):
        warning_message(f"Removing {path}")
        os.unlink(path)

    output_name = params.dirs['DESC'] + DESCENT_PREFIX + ".descent.tgt.middle.rels"
    descent_info = dict(tgt=target, mask=mask, h=h, u=u,
                        v=v,lucky=False,DRELS_FILE=output_name,init_time=descent_init_time,middle_time=0)

    if len(special_qs) == 0:
        major_message("Got lucky with the initial split! No descent necessary.")
        del descent_info['DRELS_FILE']
        descent_info['lucky'] = True
        json_custom.dump(descent_info,
                         open(params.dirs['DESC']+DESCENT_PREFIX + ".json","w"),
                         indent=True)
        os.symlink(os.path.realpath(params.dirs['DESC']+DESCENT_PREFIX +
                                    ".json"),
                   params.files['TGT_INFO']+"."+DESCENT_PREFIX)
        timeprint("Success")
        timeprint(descent_info)
        u_fac = factor(Integer(u))
        v_fac = factor(Integer(v))

        json_custom.dump(descent_info,
                         open(params.dirs['DESC']+DESCENT_PREFIX + ".json","w"),
                         indent=True)

        # In reality we can get rid of the exception now, I think.
        # raise LuckySplitException(mask, u_fac, v_fac, u, v, h)
        sys.exit(0)

    json_custom.dump(descent_info,
                     open(params.dirs['DESC']+DESCENT_PREFIX + ".json","w"),
                     indent=True)

    timeprint(f"Initializing descent with seed={topargs.seedval}"
          f" tgt=(u,v): {target}={u}/{v}")
    # f" (target found after {spin} factorization attemps)"

    with open(params.dirs['DESC']+DESCENT_PREFIX+'.todo', "w") as file:
        for q in reversed(sorted(special_qs)):
            # 0 for the rational side
            print(f"0 {q}", file=file)

    sq_desc = " ".join([f"{q.ndigits(2)}@0" for q in special_qs])
    len_sq = len(special_qs)

    if topargs.existing_init_data is not None:
        with open(todofile, 'r') as f:
            sq_desc = [x.strip() for x in f.readlines()]
            len_sq = len(sq_desc)

    timeprint(f"Descent initialized {len_sq} special-q's: {sq_desc}")

    # prepare these files so that cado doesn't recompute them
    os.symlink(os.path.realpath(params.files['POLYFILE']) + ".badideals",
               params.dirs['DESC'] + DESCENT_PREFIX + ".badideals")
    os.symlink(os.path.realpath(params.files['POLYFILE']) + ".badidealinfo",
               params.dirs['DESC'] + DESCENT_PREFIX + ".badidealinfo")

    parser = argparse.ArgumentParser(description="Descent sieving")
    parser.add_argument("--target",
                        help="Element whose DL is wanted",
                        type=str,
                        required=True)
    descent.GeneralClass.declare_args(parser)
    descent.DescentMiddleClass.declare_args(parser)


    timeprint("output_name:", output_name)
    silent_remove(output_name)
    silent_remove(output_name+".cond")

    inputs = [
        "--poly", params.files['POLYFILE'],
        "--fb1", params.files['FBFILE'],
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
        "--descent-max-increase-A", params.parameters.get('descent_max_increase_A', 10),
        "--descent-max-increase-lpb", params.parameters.get('descent_max_increase_lpb', 4),
        "--memory-margin", params.parameters.get('desc.memory_margin', 20),
        "--cadobindir", params.dirs['CADO_BUILD_DIR'],
        "--prefix", DESCENT_PREFIX,
        "--datadir", params.dirs['DESC'],
        "--ell", params.parameters['N'],
        "--renumber", params.files['RENUMBERFILE'],
        "--no-logs",
        "--target", "tgt"
    ]

    if topargs.existing_init_data is not None:
        real_todofile = todofile
    else:
        real_todofile = params.dirs['DESC']+DESCENT_PREFIX+".todo"

    with open(real_todofile, "r") as tdf:
        for line in tdf.readlines():
            hopefully_prime = Integer(line.split()[1])
            if not is_prime(hopefully_prime):
                err = f"The todofile {real_todofile} includes a composite q {hopefully_prime}. Sad!"
                raise RuntimeError(err)

    descent_middle_time = 0
    descent_middle_start = time.time()
    outfile = params.dirs['DESC']+DESCENT_PREFIX
    with open(outfile+".out","w") as fout, open(outfile+".err","w") as ferr:
        fout.reconfigure(line_buffering=True)
        ferr.reconfigure(line_buffering=True)
        with redirect_stdout(fout):
            with redirect_stderr(ferr):
                args = parser.parse_args([str(c) for c in inputs])
                general = descent.GeneralClass(args)
                middle = descent.DescentMiddleClass(general, args)
                try:
                    if topargs.existing_init_data is not None:
                        relsfile = middle.do_descent(todofile)
                    else:
                        relsfile = middle.do_descent(params.dirs['DESC']+DESCENT_PREFIX+".todo")
                    descent_middle_time = round(time.time() - descent_middle_start)
                except descent.FailedDescent as ex:
                    warning_message(str(ex))
                    timeprint("Failed")
                    sys.exit(1)

    if largeq_rels_file is not None:
        # During initialization, we collected relations for the very large qs,
        # which need to be accounted for in constructing S, and so on.
        with open(relsfile, "a") as outfile:
            with open(largeq_rels_file, "r") as infile:
                for line in infile.readlines():
                    outfile.write("\n" + line)

    try:
        if topargs.existing_init_data is not None:
            used_todofile = todofile
        else:
            used_todofile = params.dirs['DESC']+DESCENT_PREFIX+".todo"

        # Doing horrible things to fake the params object from a json exportable object
        POLYFILE = params.files['POLYFILE']
        poly = CadoPolyFile(POLYFILE); poly.read()
        params.poly = poly

        sanity_check_descent_outfile(poly.f[0], relsfile, todofile=used_todofile)

        RENUMFILE = params.files['EXPLAIN_RENUMBER_FILE']
        renum = CadoExplainRenumberFile(poly, RENUMFILE) ; renum.read()
        params.R = renum
        CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])

        S_list, S_alg_vector, S_rat_vector = construct_S(params, relsfile, relsfile+".cond.indexed.0", u, v)
        truncate_S(params, S_list, S_alg_vector)

        descent_info['middle_time'] = descent_middle_time
        os.symlink(os.path.realpath(params.dirs['DESC']+DESCENT_PREFIX +
                                    ".json"),
                   params.files['TGT_INFO']+"."+DESCENT_PREFIX)
        timeprint("Success")
        timeprint(descent_info)
        json_custom.dump(descent_info,
                         open(params.dirs['DESC']+DESCENT_PREFIX + ".json","w"),
                         indent=True)
    except RuntimeError as ex:
        if re.match(r"Taken line missing", str(ex)):
            warning_message(str(ex))
            timeprint("Failed")
            exit(1)
        elif re.match(r"constructing S", str(ex)):
            warning_message(str(ex))
            timeprint("Failed")
            exit(1)
        else:
            timeprint("Failed")
            raise ex
    except Exception as ex:
        error_message("Got exception", ex)
        timeprint("Failed")
        raise ex
