#!/usr/bin/env sage

from sage.all import *
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
from contextlib import redirect_stdout,redirect_stderr
from candy import print_command_line,warning_message,major_message,error_message
from cado_sage import CadoPolyFile
from cado_nfs_binaries import CadoNFS,CadoNFSBinaries
import time
import functools


if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='sm_cado_sage_helper.py')
    topparser.add_argument('--infile',dest='infile',required=True)
    topparser.add_argument('--outfile',dest='outfile',required=True)
    topparser.add_argument('--e',dest='e',required=True)
    topparser.add_argument('--poly',dest='poly',required=True)
    topargs = topparser.parse_args()

    e = Integer(topargs.e)
    poly = CadoPolyFile(topargs.poly); poly.read()
    f = poly.f[1]
    K = poly.K[1]
    alpha = K.gen()
    Kw = poly.nt[1]
    sm_maps = Kw.schirokauer_maps(e)

    outfile = open(topargs.outfile, "w")

    with open(topargs.infile, "r") as infile:
        for line in infile.readlines():
            line = line.strip()
            ab = line.split(",")
            a = Integer(ab[0])
            b = Integer(ab[1])
            sm_vec = vector(
                        Integers(e),
                        sum([s(a-b*alpha).list() for s in sm_maps], []))

            sm_str = ",".join( [str(sm_vec[i]) for i in range(len(sm_vec))] )
            outfile.write(sm_str + "\n")

    outfile.close()
