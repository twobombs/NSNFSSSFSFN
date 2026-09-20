#!/usr/bin/env sage
from sage.all import *
import argparse
import json
from cado_nfs_binaries import CadoNFS, CadoNFSBinaries

class Params(dict):
    __getattr__ = dict.get

if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='polyselect_helper.py')
    topparser.add_argument('--params',dest='params')
    topparser.add_argument('--admin',dest='admin')
    topparser.add_argument('--admax',dest='admax')
    topparser.add_argument('--outfile',dest='outfile')
    topargs = topparser.parse_args()

    # Doing horrible things to fake the params object from a json exportable object
    params = Params(json.loads(open(topargs.params,'r').read()))
    CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])

    output = CadoNFS("polyselect/polyselect",
                     "-P", params.parameters['poly.P'],
                     "-N", params.parameters['N'],
                     "-degree", params.parameters['POLY_DEG'],
                     "-t", params.polyselect_nthreads_or_auto,
                     "-admin", topargs.admin,
                     "-admax", topargs.admax,
                     "-incr", params.parameters['poly.incr'],
                     "-nq", params.parameters['poly.nq'],
                     "-keep", params.parameters['poly.keep'],
                     capture=True)

    firstline = True
    lines = iter(output.decode('utf-8').split("\n"))
    with open(topargs.outfile, "w") as f:
        for line in lines:
            if not line.startswith("#") or ("lognorm" in line and "raw" not in line and "optimized" not in line):
                if not firstline:
                    f.write("\n")
                else:
                    firstline = False
                f.write(line)
