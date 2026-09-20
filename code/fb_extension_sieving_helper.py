#!/usr/bin/env sage

from sage.all import *
import argparse
import json
from cado_nfs_binaries import CadoNFS,CadoNFSBinaries
from cado_sage import CadoPolyFile
from relations import strip_rational_part_of_relations
from helpers import do_fb_extension_sieving

class Params(dict):
    __getattr__ = dict.get

if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='fb_extension_sieving_helper.py')
    topparser.add_argument('--params',dest='params')
    topparser.add_argument('--jobnum',dest='jobnum')
    topparser.add_argument('--q0',dest='q0')
    topparser.add_argument('--q1',dest='q1')
    topargs = topparser.parse_args()

    # Doing horrible things to fake the params object from a json exportable object
    params = Params(json.loads(open(topargs.params,'r').read()))
    POLYFILE = params.files['POLYFILE']
    poly = CadoPolyFile(POLYFILE); poly.read()  
    params.poly = poly
    CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])
    
    do_fb_extension_sieving(params,jobnum=topargs.jobnum,q0=topargs.q0,q1=topargs.q1)
