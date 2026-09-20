#/usr/bin/env sage

from sage.all import *
import argparse
import json
from math import ceil
from functools import reduce
from operator import ior
from os import cpu_count
from time import time
from datetime import datetime
from misc_tools import fast_json_dump


def fast_union(s1, s2):
    s1.update(s2)
    return s1

def timestamp(ts=None):
    if not ts:
        ts = time()
    return datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

def timeprint(*args):
    print(f"{timestamp()}:", *args)

def query_worker(d, N, lines):
    rdict = dict()
    for line in lines:
        a = Integer(line)
        ad = a.powermod(d, N)
        rdict[str(a)] = str(ad)
    return rdict

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='oracle.py',description='eth root computation')
    parser.add_argument('-d', '--d', dest='d', required=True)
    parser.add_argument('-N', '--N', dest='N', required=True)
    parser.add_argument('--in', dest='infile', required=True)
    parser.add_argument('--out', dest='outfile', required=True)
    parser.add_argument('--mpi', action='store_true')
    parser.add_argument('--thr', type=int, default=cpu_count())
    args = parser.parse_args()

    d = Integer(args.d)
    N = Integer(args.N)

    with open(args.infile) as infile:
        lines = infile.readlines()

    nlines = len(lines)
    proccount = args.thr
    batch_size = ceil(nlines / proccount)

    # XXX: deactivated because, for some reason, this doesn't work consistently.
    if False:
    # if args.mpi:
        from mpi4py.futures import MPIPoolExecutor

        with MPIPoolExecutor(max_workers=proccount) as executor:
            timeprint(f"Running sage oracle on {proccount} workers with MPI (WARNING: doesn't seem to be fully stable and throws PMIX_ERROR error; but seems to work correctly despite the PMIX_ERROR error, if it terminates).")
            rdict_futures = [executor.submit(query_worker, d, N, lines[i * batch_size : (i+1) * batch_size]) for i in range(proccount)]
    else:
        from concurrent.futures import ProcessPoolExecutor
        from concurrent.futures import wait as concurrent_wait

        with ProcessPoolExecutor(max_workers=proccount) as executor:
            timeprint(f"Running sage oracle on {proccount} threads on single machine.")
            rdict_futures = [executor.submit(query_worker, d, N, lines[i * batch_size : (i+1) * batch_size]) for i in range(proccount)]
            concurrent_wait(rdict_futures)

    timeprint("Finished oracle queries")

    rdict = reduce(fast_union, [future.result() for future in rdict_futures], {})
    fast_json_dump(rdict, args.outfile)
