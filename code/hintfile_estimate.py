import re
import os
import time
import subprocess
import argparse
from random import randint
from run import parse_config
from candy import print_command_line, major_message
from sage.all import *

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='hintfile_estimate.py', description='help write good hintfiles')
    parser.add_argument('--side', dest='side', default="1", required=False)
    parser.add_argument('--bitsize', dest='bitsize', required=True)
    parser.add_argument('--config', dest='config', required=True)
    parser.add_argument('--hintfile', dest='hintfile', required=True)
    parser.add_argument('--locations', dest='locations', required=True)
    parser.add_argument('--nbits', dest='nbits', required=True)
    parser.add_argument('--nq', dest='nq', required=True)
    parser.add_argument('--lpb0', dest='lpb0_custom', required=False, default=0)
    parser.add_argument('--lpb1', dest='lpb1_custom', required=False, default=0)
    parser.add_argument('--polyfile', dest='polyfile', required=False)
    parser.add_argument('--fbfile', dest='fbfile', required=False)
    args = parser.parse_args()

    with open(args.locations, "r") as locations:
        for l in locations.readlines():
            if re.search(r"^#", l):
                continue
            if m := re.match(r"^(\w+)=(\S+)\s*$", l.strip()):
                var,value = m.groups()
                if os.environ.get(var) is not None:
                    print(f"Using {var} from environment")
                else:
                    print(f"Using {var} from config file {args.locations}")
                    os.environ[var] = value

    CADO_BUILD_DIR=os.environ['CADO_BUILD_DIR']
    TEMP_OUTPUT_DIR=os.environ['TEMP_OUTPUT_DIR'] + f"n{args.nbits}/"

    las_descent_cmd = CADO_BUILD_DIR + "sieve/las_descent"

    if args.polyfile:
        poly_file = args.polyfile
    else:
        poly_file = TEMP_OUTPUT_DIR + "f.poly"

    if args.fbfile:
        fb_file = args.fbfile
    else:
        fb_file = TEMP_OUTPUT_DIR + "capped.fb.gz"

    current_hintfile = args.hintfile
    seed = randint(0, 65535)

    parameters = parse_config(args.config)

    lim0 = min(2**31, 2**parameters['LPB0'], parameters['desc.lim'])
    lim1 = min(2**31, 2**parameters['LPB1'], parameters['desc.lim'])

    if int(args.lpb0_custom) > 1:
        lpb0 = int(args.lpb0_custom)
    else:
        lpb0 = parameters.get('LAS_DESCENT_UNTIL_LPB0', 'LPB0')

    if int(args.lpb1_custom) > 1:
        lpb1 = int(args.lpb1_custom)
    else:
        lpb1 = parameters.get('LAS_DESCENT_UNTIL_LPB1', 'LPB1')

    descent_middle_cmd = [
        las_descent_cmd,
        "--recursive-descent",
        "--never-discard",
        '--allow-largesq',
        '--adjust-strategy', '2',
        '--fb1', fb_file,
        '-poly', poly_file,
        '--descent-max-increase-A', str(parameters.get('descent_max_increase_A', 2)),
        '--descent-max-increase-lpb', str(parameters.get('descent_max_increase_lpb', 6)),
        '--descent-hint-table', current_hintfile,
        '--I', str(parameters['I_sieving']),
        '--lim0', str(lim0),
        '--lim1', str(lim1),
        '--lpb0', str(lpb0),
        '--mfb0', str(parameters['desc.mfb0']),
        '--lpb1', str(lpb1),
        '--mfb1', str(parameters['desc.mfb1']),
        '-t', str(parameters['desc.thr']),
        '--B', '16',
        "-q0", str(2**(int(args.bitsize)-1)),
        "-q1", str(2**int(args.bitsize)),
        "-random-sample", str(args.nq),
        "-sqside", str(args.side),
        "-seed", str(seed),
        "-v",
        "-skew", str(parameters['poly.skew']),
        '--bkthresh1', str(parameters.get('desc.bkthresh1', min(2**31, 2**lpb0, 2**lpb1))),
        "-bkmult", parameters.get('desc.bkmult', "1s:1.1")
    ]

    print_command_line(*descent_middle_cmd)
    start = time.time()
    print("launching subprocess...")
    process = subprocess.Popen(descent_middle_cmd)
    process.wait()
    end = time.time()
    print("Time:", str(round(end-start, 5)))
