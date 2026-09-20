from sage.all import *

import argparse
import json
import time
import re
import os

from cado_nfs_binaries import CadoNFS, CadoNFSBinaries
from cado_sage import CadoPolyFile
from run import parse_config

def estimate_algebraic_query_sieving_rels_per_q(params, polyfile, fbfile, nsamples=1024):
    print("===== Algebraic query sieving unique-rels-per-q estimate =====")
    print("----- Relevant parameters -----")
    print("Number of random samples:", nsamples)

    poly = CadoPolyFile(polyfile); poly.read()

    LPB1_queries = params['LPB1_queries']
    BOUNDA_queries = 2**LPB1_queries

    # The below four lines (setting special-q bounds q0 and q1) come from call_algebraic_query_sieving in helpers.py
    c0 = QQ(params["algebraic_query_sieving.q0_ratio"])
    c1 = 1 #4
    q0 = floor(c0*BOUNDA_queries)
    q1 = c1*BOUNDA_queries

    # The following comes from do_algebraic_query_sieving in helpers.py
    #AQRELS_FILE = params.files['AQRELS_FILE']+"."+jobnum
    #FBFILE = params.files['FBFILE']
    #POLYFILE = params.files['POLYFILE']
    #LPB0 = params.parameters['LPB0']

    print("algebraic_query_sieving.B:", params['algebraic_query_sieving.B'])
    print("algebraic_query_sieving.A:", params['algebraic_query_sieving.A'])
    #print("A_sieving (used as default if algebraic_query_sieving.A is unset):", params['A_sieving'])
    print("BOUNDA_queries:", BOUNDA_queries)
    print("LPB1_queries:", LPB1_queries)
    print("sieve.mfb1:", params['sieve.mfb1'])
    print("sieve.powlim", params['sieve.powlim'])
    print("algebraic_query_sieving.q0_ratio (magic constant c0)", params["algebraic_query_sieving.q0_ratio"])
    print("q0 (determined from BOUNDA_queries and q0_ratio):", q0)
    print("q1 (determined from BOUNDA_queries):", q1)

    print()
    print("----- Running las ----")
    print()
    start = time.time()
    output = CadoNFS("sieve/las",
            "--memory-margin", "80",
            "-sqside", 0,
            "-B", params['algebraic_query_sieving.B'],
            "-A", params['algebraic_query_sieving.A'],
            "--adjust-strategy", params['algebraic_query_sieving.adjust_strategy'],
            "-q0", q0,
            "-q1", q1,
            "-skew", poly.skewness,
            "-lpb0", LPB1_queries,
            "-mfb0", params['sieve.mfb1'],
            "-powlim", params['sieve.powlim'],
            "-poly", 'POLYFILE',
            "-bkmult", "1s:1.1",
            "-fb0", 'FB1',
            #"-out", 'OUT', # not used by this part of the estimator
            "-lim0", BOUNDA_queries,
            "-t", params["las.hwloc_job_binding_policy"],
            # the below lines are just for this part of the estimator
            "-random-sample", nsamples,
            "-dup",
            "-dup-qmin", q0,
            #outputs={'OUT': outfile}, # not used by this part of the estimator
            inputs={
                'FB1': fbfile,
                'POLYFILE': polyfile + ".only-side1",
            },
            capture=True # makes it return stdout
            )
    finish = time.time()
    print("...finished running las")

    print()
    print("----- Results -----")
    print()

    # The last line of output is of the form
    # "# Total xxx reports ..."
    # where xxx is the total number of relations after dedup
    output = output[output.rindex(b'#'):]
    assert output.startswith(b"# Total "), f"las output wasn't in expected format; last line was: {output}"
    num_relations = int(output.split(b" ")[2])
    relations_per_q = num_relations / nsamples
    total_special_qs = float(log_integral(q1) - log_integral(q0))
    estimated_total_relations = relations_per_q * total_special_qs

    print(f"Number of relations after dedup: {num_relations}")
    print(f"Relations per special-q: {relations_per_q}")
    print(f"Estimated total number of relations over full q-range: {estimated_total_relations}")
    print()
    min_aqrel_count = float(0.9 * 2**LPB1_queries / (LPB1_queries * log(2.)))
    print(f"\"Enough\" relations is:", min_aqrel_count)
    if estimated_total_relations < min_aqrel_count:
        print("----> This A_sieving is too small to give enough relations")
    elif estimated_total_relations > 2 * min_aqrel_count:
        print("----> This A_sieving gives more than enough relations; you should probably make it smaller")
    else:
        print("----> This A_sieving seems reasonable, but try decreasing it a little more to see if that works")
    print()

    return estimated_total_relations

def estimate_algebraic_query_sieving_time(params, polyfile, fbfile, tempdir, nsamples=1024):
    print("===== Algebraic query sieving time estimate =====")
    print("----- Relevant parameters -----")
    print("Number of random samples:", nsamples)

    poly = CadoPolyFile(polyfile); poly.read()

    LPB1_queries = params['LPB1_queries']
    BOUNDA_queries = 2**LPB1_queries

    # The below four lines (setting special-q bounds q0 and q1) come from call_algebraic_query_sieving in helpers.py
    c0 = QQ(params["algebraic_query_sieving.q0_ratio"])
    c1 = 1 #4
    q0 = floor(c0*BOUNDA_queries)
    q1 = c1*BOUNDA_queries

    # The following comes from do_algebraic_query_sieving in helpers.py
    #BOUNDA_queries = params.BOUNDA_queries
    #AQRELS_FILE = params.files['AQRELS_FILE']+"."+jobnum
    #FBFILE = params.files['FBFILE']
    #POLYFILE = params.files['POLYFILE']
    #LPB0 = params.parameters['LPB0']
    #LPB1_queries = params.parameters['LPB1_queries']

    # The las output isn't super useful to us, we just care about timing here
    outfile = f"{tempdir}las-estimate-aqrels-{time.time()}.out"

    print("algebraic_query_sieving.B:", params['algebraic_query_sieving.B'])
    print("algebraic_query_sieving.A:", params['algebraic_query_sieving.A'])
    #print("A_sieving (used as default if algebraic_query_sieving.A is unset):", params['A_sieving'])
    print("BOUNDA_queries:", BOUNDA_queries)
    print("LPB1_queries:", LPB1_queries)
    print("sieve.mfb1:", params['sieve.mfb1'])
    print("sieve.powlim", params['sieve.powlim'])
    print("algebraic_query_sieving.q0_ratio (magic constant c0)", params["algebraic_query_sieving.q0_ratio"])
    print("q0 (determined from BOUNDA_queries and q0_ratio):", q0)
    print("q1 (determined from BOUNDA_queries):", q1)

    print()
    print("----- Running las ----")
    print()
    start = time.time()
    CadoNFS("sieve/las",
            "--memory-margin", "80",
            "-sqside", 0,
            "-B", params['algebraic_query_sieving.B'],
            "-A", params['algebraic_query_sieving.A'],
            "--adjust-strategy", params['algebraic_query_sieving.adjust_strategy'],
            "-q0", q0,
            "-q1", q1,
            "-skew", poly.skewness,
            "-lpb0", LPB1_queries,
            "-mfb0", params['sieve.mfb1'],
            "-powlim", params['sieve.powlim'],
            "-poly", 'POLYFILE',
            "-bkmult", "1s:1.1",
            "-fb0", 'FB1',
            "-out", 'OUT',
            "-lim0", BOUNDA_queries,
            "-t", params["las.hwloc_job_binding_policy"],
            "-random-sample", nsamples,
            outputs={'OUT': outfile},
            inputs={
                'FB1': fbfile,
                'POLYFILE': polyfile + ".only-side1",
            }
            )
    finish = time.time()
    print("...finished running las")

    print()
    print("----- Results -----")
    print()

    num_special_qs = float(log_integral(q1) - log_integral(q0))
    core_seconds_per_q = float((finish - start)*88/nsamples) # 88 is the number of cores on our machines, TODO get this programattically
    print(f"Sieved {nsamples} special-q's in wall-clock time {finish - start:.2f}sec")
    print(f"Assuming 88 cores, this gives {core_seconds_per_q} core-seconds per special q")
    print(f"Estimated time for sieving entire range: {num_special_qs * core_seconds_per_q / 3600:.2f} core-hours")
    return

def estimate_fb_extension_sieving_rels(params, polyfile, fbfile, nsamples=1024):
    print("===== fb extension sieving #relations estimate =====")
    print("----- Relevant parameters -----")
    print("Number of random samples:", nsamples)

    poly = CadoPolyFile(polyfile); poly.read()

    LPB1 = params['LPB1']
    q1 = BOUNDA = 2**LPB1
    LPB1_queries = params['LPB1_queries']
    q0 = BOUNDA_queries = 2**LPB1_queries

    B = params.get('extension_sieving.B', 16)
    A = params.get('extension_sieving.A', params['A_sieving'])

    mfb0 = params.get('extension_sieving.mfb',
                    params['sieve.mfb1'])
    lim0 = params.get('extension_sieving.lim',
                                           BOUNDA_queries)

    print("B:", B)
    print("A:", A)
    print("lbp0:", LPB1_queries)
    print("mfb0:", mfb0)
    print("lim0:", lim0)
    print("sieve.powlim", params['sieve.powlim'])
    print("q0:", q0)
    print("q1:", q1)

    print("\n----- Running las ----\n")
    output = CadoNFS("sieve/las",
        "-sqside", 0,
        "--memory-margin", "80",
        "-B", B,
        "-A", A,
        "-skew", poly.skewness,
        "-lpb0", LPB1_queries,
        "-mfb0", mfb0,
        "-powlim", params['sieve.powlim'],
        "-poly", 'POLY',
        "-fb0", 'FB1',
        "-lim0", lim0,
        "--allow-largesq",
        "-q0", q0,
        "-q1", q1,
        "--never-discard",
        "-t", params.get("las.hwloc_job_binding_policy", "auto"),
        "-sync",
        # the below lines are just for this part of the estimator
        "-random-sample", nsamples,
        "-dup",
        "-dup-qmin", q0,
        "-exit-early", 1,
        inputs={
            'FB1': fbfile,
            'POLY': polyfile + ".only-side1",
        },
        capture=True # makes it return stdout
    )

    print("\n----- Results -----\n")

    output = output[output.rindex(b'#'):]
    assert output.startswith(b"# Total "), f"las output wasn't in expected format; last line was: {output}"
    num_relations = int(output.split(b" ")[2])
    relations_per_q = num_relations / nsamples
    total_special_qs = float(log_integral(q1) - log_integral(q0))
    estimated_total_relations = relations_per_q * total_special_qs

    print(f"Number of relations after dedup: {num_relations}")
    print(f"Relations per special-q: {relations_per_q}")
    print(f"Estimated total number of relations over full q-range: {estimated_total_relations}")
    return estimated_total_relations

if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='estimate2.py')
    topparser.add_argument('--config', dest='config', required=True)
    topparser.add_argument('--locations', dest='locations', required=True)
    topparser.add_argument('--coarse', action="store_true", help="Coarse estimate: use only 100 randomly sampled special q's instead of 1024")
    topparser.add_argument('-A', '--A', dest="A_sieving", type=int,
            help="Override the value of A_sieving (aka algebraic_query_sieving.A)")
    topparser.add_argument('--powlim', dest="powlim", type=int, help="Override the value of sieve.powlim")
    topparser.add_argument('--mfb1', dest="mfb1", type=int, help="Override the value of sieve.mfb1")
    topparser.add_argument('--fbfile', dest="fbfile", type=str, help="Specify a fbfile")
    topparser.add_argument('--polyfile', dest="polyfile", type=str, help="Specify a polyfile")
    topparser.add_argument('--extrels', action='store_true', help="compute expected number of extrels instead of aqrels.")
    topargs = topparser.parse_args()
    nsamples = 100 if topargs.coarse else 1024

    with open(topargs.locations, "r") as locations:
        for l in locations.readlines():
            if re.search(r"^#", l):
                continue
            if m := re.match(r"^(\w+)=(\S+)\s*$", l.strip()):
                var,value = m.groups()
                if os.environ.get(var) is not None:
                    print(f"Using {var} from environment")
                else:
                    print(f"Using {var} from config file {topargs.locations}")
                    os.environ[var] = value

    parameters = parse_config(topargs.config)
    nbits = parameters["MODULUS_BITS"]
    CADO_BUILD_DIR=os.environ['CADO_BUILD_DIR']
    TEMP_OUTPUT_DIR=os.environ['TEMP_OUTPUT_DIR'] + f"n{nbits}/"

    if topargs.fbfile:
        fbfile = topargs.fbfile
    else:
        fbfile = TEMP_OUTPUT_DIR + "fb.gz"

    if topargs.polyfile:
        polyfile = topargs.polyfile
    else:
        polyfile = TEMP_OUTPUT_DIR + "f.poly"

    CadoNFSBinaries().set_build_dir(CADO_BUILD_DIR)

    if topargs.A_sieving:
        print("Overriding algebraic_query_sieving.A to", topargs.A_sieving)
        parameters['algebraic_query_sieving.A'] = topargs.A_sieving
    if topargs.powlim:
        print("Overriding sieve.powlim to", topargs.powlim)
        parameters['sieve.powlim'] = topargs.powlim
    if topargs.mfb1:
        print("Overriding sieve.mfb1 to", topargs.mfb1)
        parameters['sieve.mfb1'] = topargs.mfb1

    if 'algebraic_query_sieving.A' not in parameters.keys():
        parameters['algebraic_query_sieving.A'] = parameters['A_sieving']

    if 'algebraic_query_sieving.B' not in parameters.keys():
        parameters['algebraic_query_sieving.B'] = 16

    if topargs.extrels:
        estimate_fb_extension_sieving_rels(parameters, polyfile, fbfile, nsamples)
    else:
        estimate_algebraic_query_sieving_rels_per_q(parameters, polyfile, fbfile, nsamples)
        estimate_algebraic_query_sieving_time(parameters, polyfile, fbfile, TEMP_OUTPUT_DIR, nsamples)
