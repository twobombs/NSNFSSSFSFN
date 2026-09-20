#!/usr/bin/env sage
from sage.all import *
import argparse
import json
from cado_nfs_binaries import CadoNFS, CadoNFSBinaries

class Params(dict):
    __getattr__ = dict.get

if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='poly_ropt_helper.py')
    topparser.add_argument('--params',dest='params')
    topparser.add_argument('--infile',dest='infile')
    topparser.add_argument('--outfile',dest='outfile')
    topargs = topparser.parse_args()

    # Doing horrible things to fake the params object from a json exportable object
    params = Params(json.loads(open(topargs.params,'r').read()))
    CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])

    output = CadoNFS("polyselect/polyselect_ropt",
                     "-t",params.polyselect_nthreads_or_auto,
                     "-inputpolys",topargs.infile,
                     "-area", params.parameters['poly.area'],
                     "-Bf",params.parameters['poly.Bf'],
                     "-Bg",params.parameters['poly.Bg'],
                     "-ropteffort",params.parameters.get('poly.ropteffort',5),
                     capture=True)

    with open(topargs.outfile,"w") as outfile:
        outfile.write(output.decode('utf-8'))
