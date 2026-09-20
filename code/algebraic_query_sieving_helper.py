#!/usr/bin/env sage

from sage.all import *
import argparse
import json
from cado_nfs_binaries import CadoNFS,CadoNFSBinaries
from cado_sage import CadoPolyFile
from relations import strip_rational_part_of_relations
from helpers import do_algebraic_query_sieving

class Params(dict):
    __getattr__ = dict.get

# def do_algebraic_query_sieving(params,jobnum,q0,q1):
#     BOUNDA_queries = params.BOUNDA_queries
#     AQRELS_FILE_jobnum = params.files['AQRELS_FILE']+"."+jobnum
#     FBFILE = params.files['FBFILE']
#     POLYFILE = params.files['POLYFILE']
#     poly = CadoPolyFile(POLYFILE); poly.read()
#     LPB0 = params.parameters['LPB0']
#     LPB1_queries = params.parameters['LPB1_queries']
#     A_sieving = params.parameters['A_sieving']

#     c0 = 1/4 # 1/2 # Magic constants
#     c1 = 1 #4

#     CadoNFSBinaries().set_build_dir(params.dirs["CADO_BUILD_DIR"])

#     CadoNFS("sieve/las",
#             "-sqside 1",
#             "-A", A_sieving,
#             "-q0", q0,
#             "-q1", q1,
#             "-skew", poly.skewness,
#             "-lpb0 64", # "-lpb0", str(LPB0),
#             "-lpb1", LPB1_queries,
#             "-mfb0 64", # "-mfb0 15",
#             "-mfb1", params.parameters['sieve.mfb1'],
#             "-powlim", params.parameters['sieve.powlim'],
#             "-poly", 'POLY',
#             "-fb1", 'FB1',
#             "-out", 'AQRELS',
#             "-lim1", params.BOUNDA_queries,
#             "-lim0", 2**params.parameters['alg.loglim0'], #2**31, #params.BOUNDR,
#             "-t", params.las_job_binding_policy,
#             outputs={'AQRELS': AQRELS_FILE_jobnum},
#             inputs={
#                 'FB1': FBFILE,
#                 'POLY': POLYFILE,
#             }
#             )

#     strip_rational_part_of_relations(AQRELS_FILE_jobnum)

#     return


if __name__=='__main__':
    topparser = argparse.ArgumentParser(prog='algebraic_query_sieving_helper.py')
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

    do_algebraic_query_sieving(params,jobnum=topargs.jobnum,q0=topargs.q0,q1=topargs.q1)
