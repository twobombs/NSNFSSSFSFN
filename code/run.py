from sage.all import *
from crt_ethroot import *
from helpers import *
import montgomery_ethroot
import bwc_helpers
import sys
import os
import re
import time
import subprocess
import json
import json_custom
import functools
import argparse
import configparser
import atexit
from pathlib import Path
from collections import OrderedDict, namedtuple
from misc_tools import find_factors_close_to_square_root, fast_persistent_save, fast_persistent_load
import padic_eth_root
import hybrid_root_step
import descent_large_las
import tocfile

try:
    from pebble import ProcessPool
except ImportError as e:
    from concurrent.futures import ProcessPoolExecutor as ProcessPool

from candy import major_message, error_message
from candy import OK, NOK, HURRAH

from cado.scripts import descent

from cado_nfs_binaries import CadoNFS, CadoNFSBinaries

from cado.tests.sagemath import cado_sage
from cado_sage import CadoIdealsDebugFile

x = polygen(QQ, 'x')

class LuckySplitException(Exception): pass

def kill_subprocesses():
    for process in multiprocessing.active_children():
        if process.poll() is None:
            timeprint("killing", process.pid)
            os.kill(process.pid,signal.SIGTERM)
            os.killpg(os.getpgid(process.pid),signal.SIGTERM)
            process.terminate()
            process.wait()

atexit.register(kill_subprocesses)

class Logger(object):
    def __init__(self, params, isstdout=True):
        if isstdout:
            self.terminal = sys.stdout
            self.log = open(params.files['LOGFILE'], "a")
        else:
            self.terminal = sys.stderr
            self.log = open(params.files['ERRFILE'], "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

@timing
def gen_N(params):
    """
    given a bitsize params.parameters['MODULUS_BITS'] and a public
    exponent e params.parameters['e'], compute a public key modulus N
    together with the matching secret exponent d.
    """
    parameters = params.parameters
    MODULUS_BITS = parameters['MODULUS_BITS']
    e = parameters['e']
    prime_bits = int(MODULUS_BITS/2)

    done = False
    while(not done):
        p = random_prime(2**prime_bits-1, False, 2**(prime_bits-1))
        q = random_prime(2**prime_bits-1, False, 2**(prime_bits-1))
        phi = (p-1)*(q-1)
        if gcd(e, phi) == 1:
            done = True

    N = p*q
    d = Integer(inverse_mod(e,phi))
    print("N =", N)
    print("d =", d)
    return N,d

@timing
def run_polysel(params, already_selected=False):
    parameters = params.parameters
    POLY_DEG = parameters['POLY_DEG']

    if 'N' not in parameters:
        timeprint("N not found in config, generating random N")
        parameters['N'], parameters['d'] = gen_N(params)
    else:
        timeprint("Found parameter N=",parameters['N'],"in config")

    N = parameters['N']
    e = Integer(parameters['e'])

    if not already_selected:
        if params.cadopoly:
            write_polyfile(params, 'cado', N, e, POLY_DEG, params.files['POLYFILE'])
        else:
            write_polyfile(params, 'custom', N, e, POLY_DEG, params.files['POLYFILE'])

    major_message("----Polynomial Selection----")
    print(params.poly)
    print("Shared root mod N:", params.poly.m)

    CadoNFS("utils/numbertheory_tool",
            "-poly", 'POLY',
            "-badideals", 'BADIDEALS',
            "-badidealinfo", 'BADIDEALINFO',
            "-ell", e,
            inputs={
                'POLY': params.files['POLYFILE'],
                },
            outputs={
                'BADIDEALS': params.files['POLYFILE'] + ".badideals",
                'BADIDEALINFO': params.files['POLYFILE'] + ".badidealinfo",
                }
            )
    params.save_to_file()

@timing
def run_precomp(params, nopolysel=False):
    parameters = params.parameters
    BOUNDA = params.BOUNDA
    BOUNDA_queries = params.BOUNDA_queries
    BOUNDR = params.BOUNDR
    I = params.I_sieving
    POLY_DEG = parameters['POLY_DEG']
    LPB0 = parameters['LPB0']
    LPB1 = parameters['LPB1']
    LPB1_queries = parameters['LPB1_queries']

    if nopolysel:
        timeprint("Not running polysel, using existing f.poly")

        CadoNFS("utils/numbertheory_tool",
                "-poly", 'POLY',
                "-badideals", 'BADIDEALS',
                "-badidealinfo", 'BADIDEALINFO',
                "-ell", str(parameters['e']),
                inputs={
                    'POLY': params.files['POLYFILE'],
                    },
                outputs={
                    'BADIDEALS': params.files['POLYFILE'] + ".badideals",
                    'BADIDEALINFO': params.files['POLYFILE'] + ".badidealinfo",
                    }
                )
    else:
        run_polysel(params)

    major_message("----Algebraic Sieving----")

    if LPB1 > 32:
        # Make a second, smaller fb.gz
        # (I believe we do not need any of the associated freerel, renum files)

        capped_lim = min(2**32, BOUNDA_queries)

        major_message(f"LPB1={LPB1} is large so we are making capped.fb.gz")
        CadoNFS("sieve/makefb",
                "-poly", 'POLY',
                "-out", 'FB',
                "-lim", str(capped_lim),
                "-t", params.nthreads,
                inputs={ 'POLY': params.files['POLYFILE'], },
                outputs={
                    'FB': params.files['CAPPED_FBGZ'],
                    }
                )
        major_message("Finished making the capped factor base!")

    CadoNFS("sieve/makefb",
            "-poly", 'POLY',
            "-out", 'FB',
            "-lim", BOUNDA,
            "-t", params.nthreads,
            inputs={ 'POLY': params.files['POLYFILE'], },
            outputs={
                'FB': params.files['FBFILE'],
                }
            )
    major_message("Finished making the factor base!")

    CadoNFS("sieve/freerel",
            "-poly", 'POLY',
            "-renumber", 'RENUMBER',
            "-lpb0", LPB0,
            "-lpb1", LPB1,
            "-out", 'FREEREL',
            "-dl",
            inputs={ 'POLY': params.files['POLYFILE'], },
            outputs={
                'RENUMBER': params.files['RENUMBERFILE'],
                'FREEREL': params.files['FREERELFILE'],
                }
            )

    CadoNFS("misc/debug_renumber",
            "-poly", 'POLY',
            "-renumber", 'RENUMBER',
            "-dl",
            inputs={
                'POLY': params.files['POLYFILE'],
                'RENUMBER': params.files['RENUMBERFILE'],
                },
            capture=open(params.files['DEBUG_RENUMBER_FILE'], 'w')
            )

    # This one is actually more machine readable than the debug_renumber
    # format, and also contains more information.
    CadoNFS("misc/explain_indexed_relation",
            "-poly", 'POLY',
            "-renumber", 'RENUMBER',
            "-raw", "-all",
            "-skip-ideal-checks",
            "-dl",
            inputs={
                'POLY': params.files['POLYFILE'],
                'RENUMBER': params.files['RENUMBERFILE'],
                },
            capture=open(params.files['EXPLAIN_RENUMBER_FILE'], 'w')
            )
    major_message("Finished renumbering!")

    #do_algebraic_query_sieving(params)
    call_algebraic_query_sieving(params)
    major_message("Finished query sieving!")

    #call_fb_extension_sieving(params)
    #major_message("Finished extension sieving!")

    MM = filter_relation_file(params)
    major_message("matrix:", MM)
    fast_persistent_save(MM, params.dirs['TEMP_OUTPUT_DIR']+"MM.sobj")

    if params.cado_nfs_filter:
        build_killer_rels_dict(params, MM)

    call_fb_extension_sieving(params)
    major_message("Finished extension sieving!")

    if 'n1024' not in params.files['RENUMBERFILE']:
        timeprint("parse_fb_extension_relations...")
        per_q = parse_fb_extension_relations(params.files["EXTRELS_FILE"], (1, params.BOUNDA_queries, params.BOUNDA))
        timeprint("convert_to_indexed_relations for extrels...")
        indexed_relations_file = big_convert_to_indexed_relation(per_q.values(),
                                                   params,
                                                   params.files["EXTRELS_INDEXED"])


@timing
def run_queries(params, steps=["rqueries", "aqueries", "extqueries"]):
    if steps == ["queries"]:
        steps = ["rqueries", "aqueries", "extqueries"]

    if "rqueries" in steps:
        major_message("----Rational Queries----")
        timeprint("Using rational bound " + str(params.BOUNDR) + "...")
        num_rqueries = do_rational_queries(params)
        timeprint(f"Completed {num_rqueries} rational queries!")

    if "aqueries" in steps:
        major_message("----Algebraic Queries----")
        timeprint(f"Using algebraic bound {params.BOUNDA_queries}...")
        num_aqueries = do_algebraic_queries(params, "algebraic")
        timeprint(f"Completed {num_aqueries} algebraic queries!")

    if "extqueries" in steps:
        timeprint(f"Starting extension queries.")
        num_ext_queries = do_algebraic_queries(params, "extension")
        timeprint(f"Completed {num_ext_queries} extension queries!")

@timing
def run_descent(target, params, seedval=None):
    # Return a target info struct, or None on failure

    # XXX This is dead code! See call_descents instead

    #files = params.files
    parameters = params.parameters

    e = parameters['e']
    N = params.poly.N

    ZN = Integers(N)

    with seed(seedval):
        mask, h, u, v = masked_target(ZN(target), e)

    timeprint(f"Initializing descent with seed={seedval}, tgt=(u,v): {target}={u},{v}")

    prefix = f"desc.{u}.{v}"

    u_fac = factor(Integer(u))
    v_fac = factor(Integer(v))
    special_qs = [p for p,k in factor(ZZ(u)/ZZ(v)) if p > params.BOUNDR]

    target_info = dict(tgt=target, mask=mask, h=h, u=u, v=v)

    target_info['lucky'] = len(special_qs) == 0

    if target_info['lucky']:
        major_message("Got lucky with the initial split! No descent necessary.")
        return target_info

    with open(params.dirs['DESC']+prefix+'.todo', "w") as file:
        for q in special_qs:
            # 0 for the rational side
            print(f"0 {q}", file=file)

    timeprint("Descent initialized " + str(len(special_qs)) + " special-q's on the rational side.")

    try:
        target_info['DRELS_FILE'] = do_descent_sieving(params, prefix)

        return target_info
    except RuntimeError as ex:
        if re.match(r"Failed descents", str(ex)):
            warning_message(str(ex))
            return None
        elif re.match(r"Taken line missing", str(ex)):
            warning_message(str(ex))
            return None
        else:
            raise ex
    except Exception as ex:
        error_message("Got exception", ex)
        raise ex

@timing
def run_validate_descent(params, use_target_info_file=None):
    assert use_target_info_file is not None

    target_info = json.load(open(use_target_info_file))
    u = ZZ(target_info['u'])
    v = ZZ(target_info['v'])

    uv_fac = get_uv_fac(u, v, target_info)

    DRELS_FILE = target_info['DRELS_FILE']
    DRELS_INDEXED = DRELS_FILE + ".cond.indexed"

    timeprint("Making sure construct_S works...")
    S_list, S_alg_vector, S_rat_vector = construct_S(
        params, DRELS_FILE, DRELS_INDEXED, u, v, already_indexed=False, partial_R=True, uv_fac=uv_fac
    )
    timeprint("Done with construct_S")

    timeprint("Making sure run_S_sanity_checks works...")
    run_S_sanity_checks(params, S_list, S_alg_vector, S_rat_vector, u, v, partial_R=True)
    timeprint("Done with run_S_sanity_checks")

    timeprint(f"S_alg_vector nonzero positions: {len(S_alg_vector.nonzero_positions())}")
    timeprint(f"S_rat_vector nonzero positions: {len(S_rat_vector.nonzero_positions())}")
    timeprint(f"len of S_list: {len(S_list)}")

    timeprint("Making sure truncate_S works...")
    T_list, ST_list, ST_alg_vector = truncate_S(params, S_list, S_alg_vector, partial_R=True)

    timeprint("Running the truncation assertion...")
    assert all([ST_list[k] == S_list.get(k,0) + T_list.get(k,0) for k in list(S_list.keys()) + list(T_list.keys())])

    timeprint("Running run_ST_sanity_checks...")
    run_ST_sanity_checks(params, S_list, T_list, ST_alg_vector, S_rat_vector, u, v, partial_R=True)

    timeprint("Everything looks good!")
    return


@timing
def run_linalg(params, use_target_info_file=None):
    # We've had these files for a while already, we could have done this
    # filtering earlier, especially _before the descent_
    # Note: MM depends only on aqrels, NOT on the seed or descent
    mm_path = params.dirs['TEMP_OUTPUT_DIR'] + "MM.sobj"

    try:
        import mr4mp
    except ModuleNotFoundError:
        timeprint("WARNING: Install Python package 'mr4mp' for faster mapreduce operation")

    if ('n1024' in params.files['RENUMBERFILE']) or (params.overwrite_MC):
        # MM.sobj too big!
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
    else:
        timeprint(f"Loading {mm_path} from disk")
        MM = fast_persistent_load(mm_path)
        timeprint(f"Finished loading")

    if use_target_info_file is None:
        target_info = params.target_info #json.load(open(params.files['TGT_INFO']))
        seed_ten = None
        seed_pfx = ""
    else:
        target_info = json.load(open(use_target_info_file))
        seed_ten = str(target_info['seed'])[:10]
        seed_pfx = seed_ten + "-"

    u = ZZ(target_info['u'])
    v = ZZ(target_info['v'])

    uv_fac = get_uv_fac(u, v, target_info)

    #rdict = json.load(cat_or_zcat(params.files['RQUERIES_FILE']))
    #if len(rdict) != params.R.number_of_rational_queries():
    #    raise RuntimeError(f"Weird. {len(rdict)} queries in the json file, {params.R.number_of_rational_queries()} are #expected in the renumber table ({len(params.R._ideals)} ideals)")

    DRELS_FILE = target_info['DRELS_FILE']
    DRELS_INDEXED = DRELS_FILE + ".cond.indexed"

    # depends only on matrix
    MC_saved_file = params.dirs['TEMP_OUTPUT_DIR']+"MC.sobj"
    S_block_saved_file = params.dirs['TEMP_OUTPUT_DIR']+"S_block.sobj"
    # depends on seed and target
    S_list_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"S_list.sobj"
    S_alg_vector_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"S_alg_vector.sobj"
    S_rat_vector_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"S_rat_vector.sobj"
    T_list_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"T_list.sobj"
    ST_list_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"ST_list.sobj"
    ST_alg_vector_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"ST_alg_vector.sobj"
    SC_vector_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"SC_vector.sobj"
    C_block_saved_file = params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"C_block.sobj"

    check_files_exist = [
        MC_saved_file, S_block_saved_file,
        S_list_saved_file, S_alg_vector_saved_file,
        S_rat_vector_saved_file, T_list_saved_file,
        ST_list_saved_file, ST_alg_vector_saved_file,
        SC_vector_saved_file, C_block_saved_file
    ]

    skip_make_linalg_system = True  # it's very slow, so if we can, skip as much as possible
    if params.overwrite_MC or params.overwrite_SC:
        skip_make_linalg_system = False

    for check_file in check_files_exist:
        if not os.path.exists(check_file):
            skip_make_linalg_system = False

    if skip_make_linalg_system:
        major_message("Skipping make_linalg_system, since the files all exist.")
        MC = fast_persistent_load(MC_saved_file)
        S_block = fast_persistent_load(S_block_saved_file)
        S_list = fast_persistent_load(S_list_saved_file)
        S_alg_vector = fast_persistent_load(S_alg_vector_saved_file)
        S_rat_vector = fast_persistent_load(S_rat_vector_saved_file)
        T_list = fast_persistent_load(T_list_saved_file)
        ST_list = fast_persistent_load(ST_list_saved_file)
        ST_alg_vector = fast_persistent_load(ST_alg_vector_saved_file)
        SC_vector = fast_persistent_load(SC_vector_saved_file)
        C_block = fast_persistent_load(C_block_saved_file)

    else:
        if 'n1024' in params.files['RENUMBERFILE']:
            partial_R=True
        else:
            partial_R=False

        timeprint("Start make_linalg_system")

        S_list, S_alg_vector, S_rat_vector = construct_S(
            params, DRELS_FILE, DRELS_INDEXED, u, v, already_indexed=False, partial_R=partial_R, uv_fac=uv_fac
        )

        if 'n1024' not in params.files['RENUMBERFILE']:
            # still WIP
            run_S_sanity_checks(params, S_list, S_alg_vector, S_rat_vector, u, v, partial_R=partial_R)

        timeprint(f"S_alg_vector nonzero positions: {len(S_alg_vector.nonzero_positions())}")
        timeprint(f"S_rat_vector nonzero positions: {len(S_rat_vector.nonzero_positions())}")
        timeprint(f"len of S_list: {len(S_list)}")

        T_list, ST_list, ST_alg_vector = truncate_S(params, S_list, S_alg_vector, partial_R=partial_R)

        # well it's not _exactly_ the concatenation, because if a key k is in
        # both S_list and T_list, we have ST_list[k] = S_list[k] + T_list[k]
        # assert ST_list == S_list | T_list
        if 'n1024' not in params.files['RENUMBERFILE']:
            assert all([ST_list[k] == S_list.get(k,0) + T_list.get(k,0) for k in list(S_list.keys()) + list(T_list.keys())])

        # items from T list are not smooth, they're just extra queries.
        if 'n1024' not in params.files['RENUMBERFILE']:
            run_ST_sanity_checks(params, S_list, T_list, ST_alg_vector, S_rat_vector, u, v, partial_R=partial_R)

        MC, SC_vector, S_block, C_block = make_linalg_system(params,
                                                             MM,
                                                             ST_list,
                                                             ST_alg_vector)

        if 'n1024' not in params.files['RENUMBERFILE']:
            # too expensive currently
            # for n1024

            # Files don't exist, so we need to write them
            # OR we want to overwrite them.
            if (not os.path.exists(MC_saved_file)) or params.overwrite_MC:
                major_message("Overwriting/writing to MC and S_block.")
                fast_persistent_save(MC, MC_saved_file)
                fast_persistent_save(S_block, S_block_saved_file)

            # If we have gone through with all of make_linalg_system,
            # the whole point was to get a new SC vector, so we write that information.
            major_message("Overwriting/writing to the SC vector files.")
            fast_persistent_save(S_list, S_list_saved_file)
            fast_persistent_save(S_alg_vector, S_alg_vector_saved_file)
            fast_persistent_save(S_rat_vector, S_rat_vector_saved_file)
            fast_persistent_save(T_list, T_list_saved_file)
            fast_persistent_save(ST_list, ST_list_saved_file)
            fast_persistent_save(ST_alg_vector, ST_alg_vector_saved_file)
            fast_persistent_save(SC_vector, SC_vector_saved_file)
            fast_persistent_save(C_block, C_block_saved_file)

    # Matrix M : num_rows x num_cols
    # Vector S_alg_vector : num_cols
    # Want : coefficient for each row of M (each algebraic query)
    # [1 x r] [r x c] = [1 x c]

    # Here, the objects are confusingly named:
    # MC: M || S_block, where M is the main matrix, and S_block is the character block
    # S_block: Character block of the matrix
    # SC_vector: ST_alg_vector || C_block, so it is the target vector and the characters
    # C_block: Character segment of the target vector

    # MC is used for all/any target vectors, so should always be loaded if it exists.
    # S_block as well.
    # SC_vector and C_block do depend on the specified seed and target vector.

    # pull these back from our abstraction layer in order to interface
    # with the rest of the code.
    M = MM.matrix()

    # M is MM.matrix(), and MC is M=MM.matrix() + the character block

    row_to_aquery = MM.row_to_aquery
    indexed_relations_file = MM.indexed_relations_filename(params)

    major_message("----Linear Algebra----")
    timeprint(f"The matrix M is {MC.nrows()}x{MC.ncols()}",
              f"({MC.ncols()-M.ncols()} columns are characters)")
    timeprint(f"The SC vector is of length {len(SC_vector)}")

    @timing
    def solve_sage(params=params):
        sol = MC.solve_left(-SC_vector)
        assert(M.nrows()==len(sol))
        return sol

    if params.sage_linalg and not params.cado_nfs_filter:
        timeprint("Solving with Sage solve_left")
        sol = solve_sage()

    elif params.sage_linalg and params.cado_nfs_filter:
        timeprint("Solving with -C (full filtering) and -S (sage solve_left)")
        # sol should be a solution for target=-SC_vector
        sol = solve_filtering_plus_sage(params, M, MM, MC, -SC_vector)
    else:

        bwc_failed = None
        try:
            if type(MM) is LinearAlgebraMatrix_IndexedFileOnly:
                assert(not params.cado_nfs_filter)
                sol = bwc_helpers.solve_system_allatonce(params, MC, SC_vector)
                assert(M.nrows()==len(sol))
            else:
                assert(params.cado_nfs_filter)
                timeprint("Solving with -C (full filtering) and bwc")
                if params.chars_in_mat:
                    if params.use_intermediate:
                        sol = solve_filtering_plus_bwc_chars(
                            params, M, MM, MC, SC_vector, S_block, C_block, use_existing_M_binary=True
                        )
                    else:
                        sol = solve_filtering_plus_bwc_chars(params, M, MM, MC, SC_vector, S_block, C_block)
                else:
                    sol = solve_filtering_plus_bwc(params, M, MM, MC, SC_vector, S_block, C_block, seed_ten)
                #sol = bwc_helpers.solve_system_bwc_from_filtered(params,
                #                                                 MM, ST_alg_vector,
                #                                                 S_block, C_block,
                #                                                 ST_list)
                #sys.exit(0)
                #pass

        except FileNotFoundError as ex:
            error_message("bwc failed", NOK)
            bwc_failed = ex

        if bwc_failed is not None:
            # raise bwc_failed
            if params.cado_nfs_filter:
                error_message("Fallback: -C (full filtering) and -S (sage solve_left)")
                sol = solve_filtering_plus_sage(params, M, MM, MC, -SC_vector)
            else:
                error_message("Fallback: Sage solve_left")
                sol = solve_sage()

    # Make the aqrels table-of-contents if it doesn't already exist, so that we have fast random access to the indexed relations
    timeprint("Making indexed relations table-of-contents file...")
    tocfile.maketoc(indexed_relations_file, indexed_relations_file + ".toc")
    timeprint("Done making indexed relations table-of-contents file")

    run_sol_sanity_checks(params, sol, ST_alg_vector, indexed_relations_file)
    timeprint("Found a solution! We have sol*M = ST mod e.")
    try:
        # I'm 99% sure that sol is a sage vector. It certainly is when I test without slurm. On the off chance the slurm version returns a list or a tuple instead of a sage vector, I'm adding a try/except just to be 100% sure it won't crash from us calling sol.sparse_vector()
        toreturn = LinalgOutput(sol=sol.sparse_vector(), ST_list=ST_list, ST_alg_vector=ST_alg_vector, T_list=T_list, row_to_aquery=row_to_aquery, S_rat_vector=S_rat_vector, indexed_relations_file=indexed_relations_file)
    except:
        toreturn = LinalgOutput(sol=sol, ST_list=ST_list, ST_alg_vector=ST_alg_vector, T_list=T_list, row_to_aquery=row_to_aquery, S_rat_vector=S_rat_vector, indexed_relations_file=indexed_relations_file)
    if seed_ten is None:
        fast_persistent_save(toreturn, params.dirs['TEMP_OUTPUT_DIR']+"linalgoutput.sobj")
    else:
        # It's possible fast_persistent_save won't work for n1024,
        # so at least save the things that are harder to recover otherwise.
        if 'n1024' in params.files['RENUMBERFILE']:
            try:
                sol_sparse = sol.sparse_vector()
                fast_persistent_save(sol_sparse, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-sol_sparse.sobj")
                fast_persistent_save(ST_list, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-ST_list.sobj")
                fast_persistent_save(T_list, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-T_list.sobj")
                fast_persistent_save(ST_alg_vector, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-ST_alg_vector.sobj")
                fast_persistent_save(S_rat_vector, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-S_rat_vector.sobj")
                # we can get row_to_aquery and indexed_relations_file
                # pretty easily otherwise.

                fast_persistent_save(toreturn, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-linalgoutput.sobj")
            except:
                pass

        else:
            fast_persistent_save(toreturn, params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-linalgoutput.sobj")
    return toreturn

@timing
def serialize_ttplus_ttminus(params):
        # Serialize TTPlus and TTMinus to a file
        timeprint("Loading linalgoutput.sobj")
        if len(params.use_descent_init_file) > 2:
            # want to use a specific descent output
            existing_target = params.use_descent_init_file

            with open(params.use_descent_init_file, "r") as fp:
                init_dict = json.load(fp)

            seed_ten = str(init_dict['seed'])[:10]
            linalg_output = fast_persistent_load(params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-linalgoutput.sobj")
        else:
            existing_target = None
            with open(params.files['TGT_INFO'], "r") as fp:
                params.target_info = json.load(fp)
            linalg_output = fast_persistent_load(params.dirs['TEMP_OUTPUT_DIR']+"linalgoutput.sobj")

        timeprint("Calling sol.lift_centered()")
        sol = linalg_output.sol.lift_centered()

        timeprint("Making TT")
        TT = [(linalg_output.row_to_aquery[i],s) for i,s in enumerate(sol)]
        TT += list(linalg_output.ST_list.items())

        timeprint("Making TTplus and TTminus")
        TTplus  = [((a,b),k)  for (a,b),k in TT if k>0]
        TTminus = [((a,b),-k) for (a,b),k in TT if k<0]

        timeprint("Serializing TTplus and TTminus")

        out_prefix = params.dirs['TEMP_OUTPUT_DIR']
        fast_persistent_save(TTplus, out_prefix + "TTplus.sobj")
        fast_persistent_save(TTminus, out_prefix + "TTminus.sobj")

        # also serialize in a way that we can read through it instead of loading the whole thing into memory
        timeprint("Also serializing TTplus and TTminus in text format")
        with open(out_prefix + "TTplus.txt", "w") as f:
            for (a,b),k in TTplus:
                print(a,b,k,file=f)
        with open(out_prefix + "TTminus.txt", "w") as f:
            for (a,b),k in TTminus:
                print(a,b,k,file=f)
        timeprint("Done serializing TTplus and TTminus")

@timing
def run_eth_root(params, linalg_output, use_target_info_file=None):
    sol = linalg_output.sol
    ST_list = linalg_output.ST_list
    T_list = linalg_output.T_list
    row_to_aquery = linalg_output.row_to_aquery
    S_rat_vector = linalg_output.S_rat_vector

    e = params.parameters['e']
    N = params.poly.N
    m = params.poly.m
    f = params.poly.f[1]
    ZN = Integers(N)
    K = params.poly.K[1]
    alpha = K.gen()

    R_m = ZN(0)
    R_alpha = 0

    if use_target_info_file is None:
        seed_pfx = ""
    else:
        target_info = json.load(open(use_target_info_file))
        seed_ten = str(target_info['seed'])[:10]
        seed_pfx = seed_ten + "-"

    major_message("----Starting R^d----")
    if params.montgomery_root or params.montgomery_parallel:

        if (not params.overwrite_SC) and os.path.exists(params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"R_alpha.sobj"):
            major_message("Reading R(alpha) from existing file...")
            R_alpha = fast_persistent_load(params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"R_alpha.sobj")
        else:
            R_alpha = montgomery_ethroot.eth_root_montgomery(params, linalg_output)
            fast_persistent_save(R_alpha, params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"R_alpha.sobj")

        R_m = R_alpha.polynomial().change_ring(ZN)(m)
        # Recompute U_m with lift_centered coeffs
        S_m = prod(ZN(a - b*m)**mult for ((a,b),mult) in ST_list.items())
        U_m = prod(
            ZN(a - b * m)**mult
            for ((a,b),mult) in (
                (row_to_aquery[i], sol[i].lift_centered())
                for i in sol.nonzero_positions()
            )
        )
        if ZN(power_mod(Integer(R_m), e, N)) == -ZN(S_m * U_m):
            R_m = -R_m
            R_alpha = -R_alpha
        assert ZN(power_mod(Integer(R_m), e, N)) == ZN(S_m * U_m), "R_m != S_m U_m mod N. Could be bug in eth_root_montgomery; try rerunning root step without -M and see whether that also fails"
        timeprint("Huzzah! R_m = S_m U_m mod N!")
        # This is the computation we don't want to do because it's
        # expensive, but nevertheless, it's good to have the code that
        # does it!
        if False:
            S_alpha = prod((a - b*alpha)**mult for ((a,b),mult) in ST_list.items())
            U_alpha = prod(
                (a - b * alpha)**mult
                for ((a,b),mult) in (
                    (row_to_aquery[i], sol[i].lift_centered())
                    for i in sol.nonzero_positions()
                )
            )
            assert S_alpha * U_alpha / R_alpha^e == 1

    elif params.padic_gamma_fac:

        if 'n1024' in params.files['RENUMBERFILE']:
            # special preprocessing case
            extqueries_dict = json.load(cat_or_zcat('EXTPRIMES.dict'))
            # preprocessed file, not quite json
            rqueries_dict = dict()
            with open('RATPRIMES.dict', 'r') as subset_file:
                for line in subset_file:
                    line = line.strip()
                    key = line.split(":")[0]    # str
                    val = line.split(":")[1].strip()
                    assert is_prime( int(key) )
                    rqueries_dict[key] = val
        else:
            extqueries_dict = json.load(cat_or_zcat(params.files['EXT_QUERIES_FILE']))
            rqueries_dict = json.load(cat_or_zcat(params.files['RQUERIES_FILE']))

        # This is a special case.
        # We have the gamma factorization stored, but it's too big to compute prod()
        # so we need to recover delta a little differently.
        gamma_fac = fast_persistent_load(params.dirs['TEMP_OUTPUT_DIR']+"gamma-fac.sobj")
        if params.hybrid_root:
            if os.path.exists(params.dirs['TEMP_OUTPUT_DIR']+"delta.sobj"):
                delta = fast_persistent_load(params.dirs['TEMP_OUTPUT_DIR']+"delta.sobj")
            else:
                if not os.path.exists(params.dirs['TEMP_OUTPUT_DIR']+"TTplus.sobj") or not os.path.exists(params.dirs['TEMP_OUTPUT_DIR']+"TTminus.sobj"):
                    serialize_ttplus_ttminus(params);
                delta = hybrid_root_step.run_hybrid_root_step(11000, params.dirs['TEMP_OUTPUT_DIR']) # XXX replace the 11000 with the number of bits to reconstruct
        else:
            delta = padic_eth_root.padic_eth_root(params,
                                                  linalg_output,
                                                  lift_centered=True,
                                                  pre_multiply_root=None,
                                                  pre_multiply_gamma_fac=gamma_fac)
        if not os.path.exists(params.dirs['TEMP_OUTPUT_DIR']+"delta.sobj"):
            delta = delta(alpha)
            fast_persistent_save(delta, params.dirs['TEMP_OUTPUT_DIR']+"delta.sobj")
        # R(alpha) = gamma*delta but we only have gamma_fac
        # delta should be small enough to deal with normally
        R_m = delta.polynomial().change_ring(ZN)(m)
        for fac in gamma_fac:
            nf_elt = fac[0]
            exp = fac[1]
            poly_m = nf_elt.polynomial().change_ring(ZN)(m)
            R_m = R_m * ZN(poly_m ** exp)
        # Recompute U_m with lift_centered coeffs
        S_m = prod(ZN(a - b*m)**mult for ((a,b),mult) in ST_list.items())
        U_m = prod(
            ZN(a - b * m)**mult
            for ((a,b),mult) in (
                (row_to_aquery[i], sol[i].lift_centered())
                for i in sol.nonzero_positions()
            )
        )
        if ZN(power_mod(Integer(R_m), e, N)) == -ZN(S_m * U_m):
            R_m = -R_m
        assert ZN(power_mod(Integer(R_m), e, N)) == ZN(S_m * U_m), "R_m != S_m U_m mod N. Could be bug in eth_root_montgomery; try rerunning root step without -M and see whether that also fails"
        timeprint("Huzzah! R_m = S_m U_m mod N!")
        # We don't get R_alpha in this branch, only R_m

    elif params.padic_root:
        r = padic_eth_root.padic_eth_root(params, linalg_output)
        R_m = r(ZN(m))
    else:
        S_m, U_m, ethroot_gen, extra_abks, abm_list = get_RSU_m(ST_list, sol, row_to_aquery, f, e, N, m)

        START_ETHROOT_PRIMES = params.parameters['START_ETHROOT_PRIMES']
        MAX_ETHROOT_PRIMES = params.parameters['MAX_ETHROOT_PRIMES']
        current_num_primes = START_ETHROOT_PRIMES

        while(current_num_primes <= MAX_ETHROOT_PRIMES):
            cand_R_m = Integer(CRT_R_m(params,current_num_primes,extra_abks,abm_list))

            # R**e = S*U
            if ZN(power_mod(cand_R_m, e, N)) == ZN(S_m * U_m):
                R_m = cand_R_m
                break
            else:
                current_num_primes += 1000
                timeprint("Retrying the eth-root step with " + str(current_num_primes) + " primes...")

        if current_num_primes > MAX_ETHROOT_PRIMES:
            raise RuntimeError("Retried too many times; giving up.")

    timeprint("Starting reconstruction of U")
    ts = time.time()
    # U(m)^d is computed from the algebraic queries

    aqueries_dict = json.load(cat_or_zcat(params.files['AQUERIES_FILE']))

    aqueries_num_accesses = 0
    extqueries_num_accesses = 0
    rqueries_num_accesses = 0

    Um_d = ZN(1)
    #for i in range(len(sol)):
    for i in sol.nonzero_positions():
        a,b = row_to_aquery[i]
        if params.montgomery_root or params.montgomery_parallel or params.padic_root or params.padic_gamma_fac:
            exponent = sol[i].lift_centered()
        else:
            exponent = Integer(sol[i])
        query_inp = str(a-b*m)
        abm_d = Integer(1)
        if query_inp in aqueries_dict.keys():
            abm_d = Integer(aqueries_dict[query_inp])
            aqueries_num_accesses += 1
        elif query_inp in extqueries_dict.keys():
            abm_d = Integer(extqueries_dict[query_inp])
            extqueries_num_accesses += 1
        else:
            raise RuntimeError("Uh oh!"
                               " While computing U(m)^d"
                               " we encountered relation"
                               f" number #{i} (a,b)={(a,b)}"
                               " which is not in the query logs.")
        Um_d = ZN(Um_d * power_mod(abm_d, exponent, N))
    endtime = time.time()
    timeprint("Finished U; took",endtime-ts)
    params.timing['U reconstruction'] = endtime - ts

    # If we have done things correctly, we have ST_list = S_list +
    # T_list, to match the notation in the paper. The thing that is
    # smooth on the rational side is u/v*S(m).  We handle T(m) by using
    # the fact that T(alpha) is smooth on the algebraic side, when
    # involving both the query and extension factor bases.
    # So, we compute (u/v*S(m))^d using the rational queries.
    # And we compute T(m)^d using the algebraic+extension queries.

    if not params.padic_gamma_fac:
        extqueries_dict = json.load(cat_or_zcat(params.files['EXT_QUERIES_FILE']))
        rqueries_dict = json.load(cat_or_zcat(params.files['RQUERIES_FILE']))

    # find the e-th root of the leading coefficient of the rational
    # polynomial

    g1_d = ZN(1)
    for p,k in params.poly.f[0].leading_coefficient().factor(limit=params.BOUNDR):
        if str(p) not in rqueries_dict:
            raise RuntimeError("Uh oh!"
                               " While computing lc(g)^d"
                               f" we encountered {p}"
                               " which is not in the query logs.")
        g1_d *= ZN(rqueries_dict[str(p)]) ** k
        rqueries_num_accesses += 1

    # I have:
    # -ZN(u/v * S_mq * prod([params.poly.f[0][1]**k for (a,b),k in S_list.items()])) == prod([ZN(Primes().unrank(i))**v for i,v in sorted(S_rat_vector.dict().items())])

    # -ZN(u/v * S_mq) == prod([ZN(Primes().unrank(i))**v for i,v in sorted(S_rat_vector.dict().items())]) / ZN(g1)**sum([k for (a,b),k in S_list.items()])

    starttime = time.time()
    uvSm_d = ZN(1)
    P = Primes()
    for i,exponent in S_rat_vector.dict().items():
        p_i = P.unrank(i)
        if str(p_i) not in rqueries_dict.keys():
            raise RuntimeError("Uh oh!"
                               " While computing u/v*S(m)^d"
                               f" we encountered the rational prime {p_i}"
                               " which is not in the query logs.")
        p_i_d = Integer(rqueries_dict[str(p_i)])
        uvSm_d = ZN(uvSm_d * power_mod(p_i_d, exponent, N))
        rqueries_num_accesses += 1

    # This is important when we have a non monic rational polynomial
    sum_ST = sum([k for (a,b),k in ST_list.items()])
    sum_T  = sum([k for (a,b),k in T_list.items()])
    sum_S  = sum_ST - sum_T
    uvSm_d /= g1_d ** sum_S

    Tm_d = ZN(1)
    for (a,b),exponent in T_list.items():
        query_inp = str(a-b*m)
        abm_d = Integer(1)
        if query_inp in aqueries_dict.keys():
            abm_d = Integer(aqueries_dict[query_inp])
            aqueries_num_accesses += 1
        elif query_inp in extqueries_dict.keys():
            abm_d = Integer(extqueries_dict[query_inp])
            extqueries_num_accesses += 1
        else:
            raise RuntimeError("Uh oh!"
                               " While computing T(m)^d"
                               " we encountered relation"
                               f" number #{i} (a,b)={(a,b)}"
                               " which is not in the query logs.")
        Tm_d = ZN(Tm_d * power_mod(abm_d, exponent, N))
    endtime = time.time()
    params.timing["u/v*S(m)^d"] = endtime-starttime
    timeprint("u/v*S(m)^d took", endtime-starttime)

    if use_target_info_file is None:
        target_info = params.target_info #json.load(open(params.files['TGT_INFO']))
    else:
        target_info = json.load(open(use_target_info_file))

    mask = target_info['mask']
    found = ZN(inverse_mod(Integer(mask), N) * inverse_mod(Integer(R_m), N) * uvSm_d * Um_d * Tm_d)

    fast_persistent_save(found, params.dirs['TEMP_OUTPUT_DIR']+seed_pfx+"found.sobj")
    timeprint("found", str(found))

    timeprint(f"aqueries_num_accesses: {aqueries_num_accesses}")
    timeprint(f"extqueries_num_accesses: {extqueries_num_accesses}")
    timeprint(f"rqueries_num_accesses: {rqueries_num_accesses}")

    return found

@timing
def run_indiv(params, existing_init_data=None):

    # TODO: Add an option to skip descent init,
    # and pull that info from a file. (And don't overwrite the target.)
    # Pending the rest of the new descent init code.

    Path(params.dirs['DESC']).mkdir(exist_ok=True)

    N = params.poly.N
    ZN = Integers(N)

    target = generate_or_load_target(params)

    # The individual computation

    major_message("----Descents----")

    ratio = params.parameters.get('initial_descent_smoothness_ratio', 2.25)
    ratio = float(ratio)
    initial_smoothness_maxbits = params.parameters['MODULUS_BITS'] / 2 / ratio
    initial_smoothness_maxbits = int(initial_smoothness_maxbits)

    for b in range(min(params.parameters['LPB0'],
                       params.parameters['LPB1']),
                   initial_smoothness_maxbits+1):
        raw_norms = []
        for side in range(2):
            I = params.I_sieving
            f = params.poly.f[side]
            s = params.poly.skewness
            d = f.degree()
            f_offset = max([log(abs(f[i])*s**(i-d/2),2) for i in range(d+1)])
            raw_norms.append((I+b/2)*d + f_offset)
        for side in range(2):
            n = [x.round() for x in raw_norms]
            n[side] -= b
            print(f"With I={I}, a {b}@{side} special-q"
                  f" will have norms of {n[0]} and {n[1]} bits")

    write_hintfile(params.parameters['desc.hintfile'],params.files['HINTFILE'], params.BOUNDR, params.BOUNDA, params.I_sieving)

    if params.descent_slurm:
        target_info = call_descents_slurm(target, params, existing_init_data)
    else:
        target_info = call_descents(target, params)

    if target_info is None:
        raise RuntimeError("Uh oh! Retried descent too many times.")

    if target_info['lucky']:
        timeprint("Lucky split, computing solution directly")
        mask = ZN(target_info['mask'])
        u = ZZ(target_info['u'])
        v = ZZ(target_info['v'])
        # tgt^d = u^d / (v^d * mask)
        # u^d and v^d come from queries
        rqueries_dict = json.load(cat_or_zcat(params.files['RQUERIES_FILE']))
        rquery = lambda p : ZN(rqueries_dict[str(p)])
        found = prod([rquery(p)**mult for (p, mult) in factor(u/v)]) / mask
        return found

    return True

@timing
def run_descent_ecm_init(params):
    # Run and save the init_data to a file, to be consumed by run_indiv.

    Path(params.dirs['DESC']).mkdir(exist_ok=True)

    N = params.poly.N
    ZN = Integers(N)

    params.save_to_file()
    target = generate_or_load_target(params)

    if (not params.slurm) and (not params.descent_slurm):
        seedval = int(time.time()) * 10**6
        call_descent_large_init(target,params,seedval)
        # this will produce only one output file

    else:
        # We would hope to only need one job.
        # However unfortunate errors sometimes come up which cannot be
        # solved by rerandomizing, like if the "new" polynomials come with
        # so smooth of a modulus that inverses can't be taken.
        # In that case we hope that a different seed has better luck.
        numjobs = params.parameters.get('desc_ecm_init.numjobs', 8)
        processes = []

        for jobnum in range(numjobs):
            seed_j = int(time.time() + jobnum) * 10**6
            command_list = [
                params.files['SAGE'],
                "descent_large_init.py",
                "--params", params.files['PARAMS'],
                "--seed", str(seed_j)
            ]
            processes.append(slurmit(params, " ".join(command_list), params.prefix[:-1]+"-ecminit", jobnum))

        finished_processes, cputime_slurm = slurm_wait(processes)
        overall_cputime.add(cputime_slurm)
        print(f"Finished running {numjobs} large descent inits!")

    return True


def check_solution(params,found):
    found_e = power_mod(found, params.parameters['e'], params.poly.N)
    timeprint("Found^e:", found_e, -found_e)
    target = generate_or_load_target(params)
    timeprint("Target:", target)
    ok = target in [found_e, -found_e]
    timeprint(f"all good! {HURRAH}" if ok else f"NOK NOK NOK {NOK}")
    assert ok


def parse_config(configfilename):
    """
    parse the .config file, probably from the config/ subdirectory. It's
    interpreted in the good old syntax of ms-dos ini files. Most
    arguments are retained as strings, except integer arguments, which
    are converted to ints.
    """
    config = configparser.ConfigParser()
    # preserve case
    config.optionxform=str
    config.read(configfilename)

    parameters = {}
    for key in config['DEFAULT']:
        if config['DEFAULT'][key].isnumeric():
            parameters[key] = config['DEFAULT'].getint(key)
        else:
            parameters[key] = config['DEFAULT'][key]
    return(parameters)

class Params(object):
    def __init__(self, args):
        locationsfilename = args.locations

        locations_dict = {}
        with open(locationsfilename, "r") as locations:
            for l in locations.readlines():
                if re.search(r"^#", l):
                    continue
                if m := re.match(r"^(\w+)=(\S+)\s*$", l.strip()):
                    var,value = m.groups()

                    if var in ['SLURM_JOB_PARTITION', 'BWC_SLURM_JOB_PARTITION']:
                        locations_dict[var] = value
                        continue

                    if os.environ.get(var) is not None:
                        timeprint(f"Using {var} from environment")
                    else:
                        timeprint(f"Using {var} from config file {locationsfilename}")
                        os.environ[var] = value

        CADO_BUILD_DIR          = os.environ['CADO_BUILD_DIR']
        TEMP_OUTPUT_DIR         = os.environ['TEMP_OUTPUT_DIR']
        SAGE                    = os.environ.get('SAGE', 'sage')
        PYTHON                  = os.environ.get('PYTHON', 'python3')
        SLURM_EXCLUDE           = os.environ.get('SLURM_EXCLUDE','')
        SLURM_JOB_PARTITION     = locations_dict.get('SLURM_JOB_PARTITION', os.environ.get('SLURM_JOB_PARTITION', 'hiprio'))
        BWC_SLURM_JOB_PARTITION = locations_dict.get('BWC_SLURM_JOB_PARTITION', os.environ.get('BWC_SLURM_JOB_PARTITION', 'hiprio-88-cores'))

        CadoNFSBinaries().set_build_dir(CADO_BUILD_DIR)

        self.cado_nfs_filter = True #args.cado_nfs_filter
        self.sage_linalg = args.sage_linalg
        self.crt_root = args.crt_root
        self.padic_root = args.padic_root
        self.padic_gamma_fac = args.padic_gamma_fac
        self.hybrid_root = args.hybrid_root
        self.montgomery_root = args.montgomery_root or ((not self.crt_root) and (not args.padic_root))
        if self.padic_gamma_fac:
            self.montgomery_root = False
        self.debug = args.debug
        self.sm_alg = args.sm_alg
        self.slurm = args.slurm
        self.slurm_exclude = SLURM_EXCLUDE
        self.slurm_job_partition = SLURM_JOB_PARTITION
        self.bwc_slurm_job_partition = BWC_SLURM_JOB_PARTITION
        self.cadopoly = True #args.cadopoly
        self.oracle = args.oracle
        self.chars_in_mat = args.chars_in_mat
        self.descent_slurm = args.descent_slurm
        #self.descent_ecm_init = args.descent_ecm_init
        self.bwc_slurm = args.bwc_slurm
        self.save_intermediate = args.save_intermediate
        self.use_intermediate = args.use_intermediate
        self.overwrite_MC = args.overwrite_MC
        self.overwrite_SC = args.overwrite_SC
        self.extra_prefix = args.extra_prefix
        self.montgomery_parallel = args.montgomery_parallel
        self.montgomery_new_renumber = args.montgomery_new_renumber
        self.use_descent_init_file = args.use_descent_init_file
        self.todofile_glob = args.todofile_glob
        self.descent_strategy = args.descent_strategy
        self.descent_counter = args.descent_counter
        self.given_A = args.given_A
        self.given_mfb1 = args.given_mfb1
        self.scipy_matrix = args.scipy_matrix

        self.given_prec = args.given_prec
        self.given_mnb = args.given_mnb
        self.given_lub = args.given_lub
        self.given_fin = args.given_fin
        self.given_t = args.given_t
        self.given_pub = args.given_pub

        self.dirs = dict()
        self.files = dict()
        self.timing = dict()

        self.dirs['CADO_BUILD_DIR'] = CADO_BUILD_DIR

        ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
        self.dirs['ROOT_DIR'] = ROOT_DIR
        PARAMS_DIR = os.path.join(ROOT_DIR, "parameters/")
        self.dirs['PARAMS_DIR'] = PARAMS_DIR

        self.parameters = parse_config(args.config)

        self.mpi = args.mpi
        self.mpi_extra_args = args.mpi_extra_args
        if self.mpi:
            found_mpi_exec = False
            sp = subprocess.run(["ompi_info"], stdout=subprocess.PIPE)
            if sp.returncode == 0:
                m = re.search(rb'\s+Prefix: ([^\s]+).*', sp.stdout)
                if m:
                    mpi_prefix = m.groups()[0].decode()
                    mpi_mpirun_bin = f"{mpi_prefix}/bin/mpirun"

                    if os.path.isfile(mpi_mpirun_bin):
                        self.mpi_prefix     = mpi_prefix
                        self.mpi_mpirun_bin = mpi_mpirun_bin
                        found_mpi_exec      = True

            if not found_mpi_exec:
                timeprint("WARNING: Cannot find MPI executable, disabling MPI!")
                self.mpi = False

            if "mpi.thr" not in self.parameters or self.parameters["mpi.thr"] == "auto":
                if "SLURM_NPROCS" in os.environ:
                    # Rely on Slurm to get the number of processors available (across
                    # whatever partition was allocated to Slurm)
                    self.parameters["mpi.thr"] = int(os.environ["SLURM_NPROCS"]) - 1
                else:
                    # Default to the number of local cores
                    try:
                        self.parameters["mpi.thr"] = len(os.sched_getaffinity(0))
                    except AttributeError:
                        # Getting the number of processors is weird on macos
                        self.parameters["mpi.thr"] = 1
            timeprint(f"Using {self.parameters['mpi.thr']} MPI processes.")

            if "bwc.mpi" not in self.parameters or self.parameters["bwc.mpi"] == "auto":
                if "SLURM_NNODES" in os.environ:
                    self.parameters["bwc.mpi"] = find_factors_close_to_square_root(
                        int(os.environ["SLURM_NNODES"])
                    )
                else:
                    timeprint("WARNING: could not find number of MPI nodes, either set 'bwc.mpi' to <a>x<b> for a*b nodes or use 'salloc' and make sure SLURM_NNODES is set.")
                    self.mpi = False

            if "bwc.mpi_slurm" not in self.parameters and self.mpi and self.bwc_slurm:
                timeprint("Parameter 'bwc.mpi_slurm' not set, trying to auto detect.")
                self.parameters["bwc.mpi_slurm"] = "auto"

        if "bwc.thr" not in self.parameters:
            self.parameters["bwc.thr"] = "auto"

        if not "bwc.thr.controller" in self.parameters or self.parameters["bwc.thr.controller"] == "auto":
            if "SLURM_NTASKS_PER_NODE" in os.environ:
                self.parameters["bwc.thr.controller"] = find_factors_close_to_square_root(
                    int(os.environ["SLURM_NTASKS_PER_NODE"]) // 2
                )
            else:
                # Assume all machines have the same number of CPUs as the current one.
                self.parameters["bwc.thr.controller"] = find_factors_close_to_square_root(
                    os.cpu_count() // 2
                )

        cadopoly_in_params = self.parameters.get('cadopoly')
        if cadopoly_in_params is not None:
            self.cadopoly = bool(cadopoly_in_params)

        timeprint("parameters:",self.parameters)

        if m:=re.search(r"(\w+)\.config$", args.config):
            self.prefix = m.group(1) + "/"
        elif (nbits := self.parameters['MODULUS_BITS']) is not None:
            self.prefix = f"n{nbits}/"
        else:
            self.prefix = ""

        temp = TEMP_OUTPUT_DIR + self.prefix
        self.dirs['TEMP_OUTPUT_DIR'] = temp

        timeprint(f"Ensuring {temp} exists")
        Path(temp).mkdir(exist_ok=True)

        self.dirs['DUP'] = temp+"0/"
        self.dirs['BWC'] = temp+"bwc/"
        self.dirs['BWCBINDIR'] = CADO_BUILD_DIR+"linalg/bwc"
        self.dirs['DESC'] = temp+"desc/"
        self.dirs['SLURMJOBS'] = temp+"slurmjobs/"
        Path(self.dirs['SLURMJOBS']).mkdir(exist_ok=True)
        self.files['SAGE'] = SAGE
        self.files['PYTHON'] = PYTHON
        self.files['LOGFILE'] = temp + "stdoutlog.log"
        self.files['ERRFILE'] = temp + "stderrlog.log"

        self.files['POLYSELECT'] = temp + "polyselect"
        self.files['POLYSELECTOPT'] = temp + "polyselect_ropt"
        self.files['POLYFILE'] = temp + "f.poly"
        self.files['FBFILE'] = temp + "fb.gz"
        self.files['RENUMBERFILE'] = temp + "renumber.gz"
        self.files['HINTFILE'] = temp + "ignore.hint"
        self.files['DEBUG_RENUMBER_FILE'] = temp + "renumber.info"
        self.files['EXPLAIN_RENUMBER_FILE'] = temp + "renumber.explained"
        self.files['RQUERIES_TODO'] = temp + "rqueries.todo"
        self.files['RQUERIES_FILE'] = temp + "rqueries.json"
        self.files['AQUERIES_TODO'] = temp + "aqueries.todo"
        self.files['AQUERIES_FILE'] = temp + "aqueries.json"
        self.files['AQRELS_FILE'] = temp + "aqrels.out"
        self.files['AQRELS_INDEXED'] = temp + "aqrels.out.indexed"
        self.files['TGT_INFO'] = temp + "tgt.json"
        self.files['EXTRELS_FILE'] = temp + "extrels.out"
        self.files['EXTRELS_INDEXED'] = temp + "extrels.out.indexed"
        self.files['EXT_QUERIES_FILE'] = temp + "extqueries.json"
        self.files['EXT_QUERIES_TODO'] = temp + "extqueries.todo"
        self.files['FREERELFILE'] = temp + "freerel"
        self.files['KILLER_RELS_DICT'] = temp + "killer_rels.json"
        self.files['KILLER_RELS_ORDER'] = temp + "cancellation_order.json"
        self.files['PARAMS'] = temp + "params.json"
        # note: unlinked_ideals are within the algebraic factor base (bounda_queries)
        self.files['UNLINKED_IDEALS'] = temp + "unlinked_ideals.out"
        self.files['UNLINKED_IDEALS_TODOS'] = temp + "unlinked_ideals.todo"
        self.files['EXT_UNLINKED_IDEALS'] = temp + "ext_unlinked_ideals.out"
        self.files['EXT_UNLINKED_TODOS'] = temp + "ext_unlinked_ideals.todo"
        # note:
        # Filtering also leaves some ideals unlinked, even if we have a relation for them.
        # That remains a mystery. However we store them separately from the above ideals
        # which are unlinked due to incomplete sieving.
        self.files['FILT_UNLINKED_IDEALS'] = temp + "filt_unlinked_ideals.out"
        self.files['FILT_UNLINKED_TODOS'] = temp + "filt_unlinked_ideals.todo"
        self.files['POLYSEL_SELECT_PFX'] = temp + "select/polyselect.out"
        self.files['THE_TARGET'] = temp + "thetarget"
        self.files['GOOD_DESCENT_INIT'] = temp + "good_descent_inits"

        self.files['CPUBINDING_CONF_FILE'] = PARAMS_DIR + "cpubinding.conf"

        Path(temp + "select").mkdir(exist_ok=True)
        Path(temp + "algrels").mkdir(exist_ok=True)

        if not os.path.isfile(CADO_BUILD_DIR + "misc/debug_renumber"):
            raise RuntimeError("Missing requirement: Need to make debug_renumber in the cado-nfs directory.\n")


        self.BOUNDA = 2**self.parameters['LPB1'] # LPB1
        self.BOUNDA_queries = 2**self.parameters['LPB1_queries']
        self.BOUNDR = 2**self.parameters['LPB0'] # LPB0
        self.I_sieving = self.parameters['I_sieving']

        if int(self.parameters['LPB1']) > 32:
            # We need a second, capped factor base file to give to las programs.
            # To be used only for sieving primes, not restricting what can show up in factorizations.
            # NOT used for renumbering.
            self.files['CAPPED_FBGZ'] = temp + "capped.fb.gz"
        else:
            self.files['CAPPED_FBGZ'] = self.files['FBFILE']

        # LPB1 algebraic smoothness bound for linear algebra
        # LPB1_queries algebraic smoothness bound for descent and queries
        # LPB0 rational smoothness bound
        # A_sieving precomputation las
        # I_sieving las_descent

        self.has_hwloc = False
        if os.uname().sysname == 'Linux':
            sp = subprocess.run(["readelf", "-d",
                                 CADO_BUILD_DIR + "sieve/las"],
                                stdout=subprocess.PIPE)
            if sp.returncode == 0 and re.search(rb"libhwloc.so", sp.stdout):
                self.has_hwloc = True
        elif os.uname().sysname == 'Darwin':
            sp = subprocess.run(["dyld_info", "-dependents",
                                 CADO_BUILD_DIR + "sieve/las"],
                                stdout=subprocess.PIPE)
            if sp.returncode == 0 and re.search(rb"libhwloc.so", sp.stdout):
                self.has_hwloc = True

        self.nthreads = multiprocessing.cpu_count()

        if self.has_hwloc:
            self.polyselect_nthreads_or_auto = "auto"
            self.las_job_binding_policy = self.parameters["las.hwloc_job_binding_policy"]
        else:
            self.polyselect_nthreads_or_auto = self.nthreads
            self.las_job_binding_policy = self.nthreads

        # Superseded by self.las_job_binding_policy now.
        del self.parameters["las.hwloc_job_binding_policy"]

        major_message(f"multiprocessing: -t {self.nthreads},",
                      f"for polyselect: -t {self.polyselect_nthreads_or_auto}",
                      f"and for las: -t {self.las_job_binding_policy}")


    @functools.cached_property
    def poly(self):
        filename = self.files['POLYFILE']
        major_message(f"Reading poly from {filename}")
        poly = CadoPolyFile(filename)
        poly.read()
        for f, K in zip(poly.f, poly.K):
            B = 10**7
            D = f.discriminant()
            primes = [p for p,e in D.factor(limit=B) if e>1]
            timeprint(f"Computing maximal order of {K} with prime list {primes}")
            OK = K.maximal_order(v=primes, assume_maximal=True)
        return poly

    @functools.cached_property
    def R(self):
        """
        of course it only makes sense to use (call) this property once
        the EXPLAIN_RENUMBER_FILE file is created, which gets done in
        run_precomp
        """
        filename = self.files['EXPLAIN_RENUMBER_FILE']
        major_message(f"Reading renumber table info from {filename}")
        timeprint("Initializing R...")
        time0 = time.time()
        R = CadoExplainRenumberFile(self.poly, filename)
        time1 = time.time()
        timeprint("Beginning to run R.read()...")
        R.read()
        time2 = time.time()
        timeprint("Woo! Finished reading R.")

        self.timing['init_renumberfile'] = time1 - time0
        self.timing['read_renumberfile'] = time2 - time1

        return R

    def save_to_file(self):
        self_dict = {
            'files':                        self.files,
            'parameters':                   self.parameters,
            'dirs':                         self.dirs,
            'polyselect_nthreads_or_auto':  self.polyselect_nthreads_or_auto,
            'las_job_binding_policy':       self.las_job_binding_policy,
            'BOUNDA':                       self.BOUNDA,
            'BOUNDA_queries':               self.BOUNDA_queries,
            'BOUNDR':                       self.BOUNDR,
            'timing':                       self.timing,
            'has_hwloc':                    self.has_hwloc
        }
        if hasattr(self,'target'):
            self_dict['target'] = self.target

        with open(self.files['PARAMS'], "w") as fp:
            json_custom.dump(self_dict, fp, indent=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='run.py', description='eth root computation')
    parser.add_argument('-l','--locations', dest='locations', default="locations.config", required=False)
    parser.add_argument(dest='config', default="config/n60.config")
    parser.add_argument('step',nargs='?', choices=[
        'info', 'gen',
        'precomp', 'precomp_nopolysel',
        'polysel_select', 'polysel_ropt', 'polysel_candidates', 'polysel_alreadysel', 'polysel_only',
        'algebraic_sieving_only', 'extra_algebraic', 'extra_granular_algebraic',
        'extension_only', 'extra_granular_extension',
        'extra_granular_filter', 'filter_only',
        'descent_large_init', 'descent_large_q', 'cleanup_large_q', 'descent_large_las_start', 'descent_large_las_rebalance', 'descent_large_las_todo',
        'descent_large_las_extract', 'descent_large_las_finish', 'descent_bottom',
        'validate_descent',
        'indiv','indiv_only',
        'queries', 'rqueries', 'aqueries', 'extqueries',
        'linalg', 'linalg_nobwc',
        'root', 'until_root',
        'serialize_ttplus_ttminus',
        'all',
        'update_alg_unlinked', 'update_ext_unlinked'
    ], default='all')
    #parser.add_argument('-C','--cado-nfs-filter', action='store_true')
    parser.add_argument('--sm-alg', dest='sm_alg', choices=['sm_simple','sm_append','cado_sage'], required=False, default='sm_append')
    parser.add_argument('-S', '--sage-linalg', action='store_true')
    parser.add_argument('-M', '--montgomery-root', action='store_true')
    parser.add_argument('-Mll', '--montgomery-parallel', dest='montgomery_parallel', action='store_true')
    parser.add_argument('--crt-root', action='store_true')
    parser.add_argument('-P', '--padic-root', action='store_true')
    parser.add_argument('--padic-gamma-fac', dest='padic_gamma_fac', action='store_true')
    parser.add_argument('--hybrid-root', dest='hybrid_root', action='store_true', help="Do hybrid p-adic+CRT root. Must be used with --padic-gamma-fac")
    parser.add_argument('--slurm', dest='slurm',action='store_true', help="Distribute some computation across machines using Slurm.")
    parser.add_argument('--mpi', dest='mpi',action='store_true', help="Run with MPI support (requires the corresponding setup, see README).")
    parser.add_argument('--mpi-extra-args', default="--mca usnic ucx --map-by core --bind-to core:overload-allowed",
                        help="Specify additional arguments for mpirun/mpiexec commands. Passed on to bwc.")
    parser.add_argument('--debug', action='store_true', help="Enable extra debugging output.")
    parser.add_argument('--cadopoly', action='store_true')
    parser.add_argument('--chars-in-mat', dest='chars_in_mat',action='store_true')
    parser.add_argument('--scipy-matrix', dest='scipy_matrix',action='store_true')
    parser.add_argument('--descent-slurm', dest='descent_slurm',action='store_true')
    parser.add_argument('--use-descent-init-file', dest='use_descent_init_file', default="")
    parser.add_argument('--descent-counter', dest='descent_counter')
    parser.add_argument('--desc-todofile-glob', dest='todofile_glob', default='')
    parser.add_argument('--bwc-slurm', action='store_true')
    parser.add_argument('--oracle', choices=['sage', 'luna_k6', 'remote_luna_S750'], required=False, default='sage')
    parser.add_argument('--save-intermediate',dest='save_intermediate', action='store_true')
    parser.add_argument('--use-intermediate',dest='use_intermediate', action='store_true')
    parser.add_argument('--overwrite-MC',dest='overwrite_MC',action='store_true')
    parser.add_argument('--overwrite-SC',dest='overwrite_SC',action='store_true')
    parser.add_argument('--extra-prefix', dest='extra_prefix', default="")
    parser.add_argument("--given-A", dest='given_A', default=-1)
    parser.add_argument("--given-mfb1", dest='given_mfb1', default=-1)
    parser.add_argument('--fix-paths', action='store_true',
                        help="set this flag if you continue a run and the temporary ouput folder (TEMP_OUTPUT_DIR) was moved.")
    parser.add_argument("--descent-strategy", dest='descent_strategy', default="")
    parser.add_argument('--montgomery-new-renumber',dest='montgomery_new_renumber',action='store_true', help="Use random-access renumber table in montgomery (need to manually create necessary files first, see comment in random_access_renumber.py)")
    parser.add_argument("--given-prec", dest='given_prec', default=-1)
    parser.add_argument("--given-mnb", dest='given_mnb', default=-1)
    parser.add_argument("--given-lub", dest='given_lub', default=-1)
    parser.add_argument("--given-fin", dest='given_fin', default=-1)
    parser.add_argument("--given-t", dest='given_t', default=-1)
    parser.add_argument("--given-pub", dest='given_pub', default=-1)
    args = parser.parse_args()

    params = Params(args)

    FILES_TO_PATCH = ["PARAMS", "TGT_INFO"]

    if args.fix_paths:
        for fname in FILES_TO_PATCH:
            with open(params.files[fname], "r+") as fp:
                content = fp.read()
                patched_content, nsubs = re.subn(rf"\/.*\/{params.prefix}",
                    params.dirs['TEMP_OUTPUT_DIR'],
                    content,
                    flags=re.M)

                if nsubs > 0:
                    fp.seek(0)
                    fp.write(patched_content)




    slurm_nnodes = os.environ.get("SLURM_NNODES", "n/a")
    slurm_ntasks_per_node = os.environ.get("SLURM_NTASKS_PER_NODE", "n/a")

    config_info = OrderedDict([
        ("step", args.step),
        ("cado-nfs-filter", True),
        ("sm-alg", args.sm_alg),
        ("sage-linalg", args.sage_linalg),
        ("montgomery-root", args.montgomery_root),
        ("padic-root", args.padic_root),
        ("crt-root",args.crt_root),
        ("MPI", args.mpi),
        ("mpi_extra_args", args.mpi_extra_args),
        ("slurm", args.slurm),
        ("descent_slurm", args.descent_slurm),
        ("bwc_slurm", args.bwc_slurm),
        ("slurm_nnodes", slurm_nnodes),
        ("slurm_ntasks_per_node", slurm_ntasks_per_node),
        ("debug", args.debug),
    ])

    timing_info = OrderedDict([("run_precomp",0),
                               ("polysel_higheffort_select",1),
                               ("polysel_higheffort_ropt",1),
                               ("polysel_higheffort_candidates",1),
                               ("run_polysel",1),
                               ("write_polyfile",2),
                               ("call_algebraic_query_sieving",1),
                               ("aqrels_slurm_sieving",2),
                               ("aqrels_extra_slurm_sieving",2),
                               ("call_fb_extension_sieving",1),
                               ("extension_slurm_sieving",2),
                               ("call_todo_sieving",2),
                               ("filter_relation_file",1),
                               ("swap_parts_of_relations",2),
                               ("build_killer_rels_dict",1),
                               ("init_renumberfile",1),
                               ("read_renumberfile",1),
                               #("R.number_of_fb_valuations",1),
                               ("update_collection_of_unlinked_ideals",1),
                               ("run_queries",0),
                               ("do_rational_queries",1),
                               ("do_algebraic_queries",1),
                               ("run_indiv",0),
                               ("call_descents",1),
                               ("call_descents_slurm",1),
                               ("descent_init",2),
                               ("descent_middle",2),
                               ("run_descent_ecm_init",2),
                               ("call_descent_large_q_slurm",2),
                               ("start_descent_large_las",2),
                               ("rebalance_descent_large_las",2),
                               ("todofile_descent_large_las",2),
                               ("finish_descent_large_las",2),
                               ("descent_rock_bottom",2),
                               ("extract_outstanding_qs",2),
                               ("handle_bottom_special_q",2),
                               ("handle_bottom_special_q_composites",2),
                               ("run_linalg",0),
                               ("run_linalg_nobwc",0),
                               ("construct_S",1),
                               ("run_S_sanity_checks",1),
                               ("run_ST_sanity_checks",1),
                               ("truncate_S",1),
                               ("algebraic_consistency_check",2),
                               ("parse_fb_extension_relations",2),
                               ("convert_to_indexed_relation",2),
                               ("make_linalg_system",1),
                               #("sm_simple_init",2),
                               #("sm_append_init",2),
                               #("cado_sage_init",2),
                               ("sm_init",2),
                               ("compute_sm_block_for_matrix",2),
                               ("maps_from_more_ab",2),
                               ("C_block",2),
                               #("variants",2),
                               ("block_matrix",2),
                               ("SC_vector",2),
                               ("solve_filtering_plus_sage",1),
                               ("solve_filtering_plus_bwc",1),
                               ("apply_filtering_to_target",2),
                               ("killer relations maps_from_more_ab", 3),
                               ("write_tgt_vector",2),
                               ("linalg/bwc/bwc.pl (complete)", 2),
                               ("bwc prep", 3),
                               ("bwc krylov", 3),
                               ("bwc lingen", 3),
                               ("bwc mksol", 3),
                               ("bwc gather", 3),
                               ("bwc cleanup", 3),
                               ("making sol_matrix",2),
                               ("assert sol_matrix",2),
                               ("small_matrix_on_the_right",2),
                               ("second_solve.left_kernel",2),
                               ("final solve_left",2),
                               ("sol_f computations",2),
                               ("recover sol from sol_f",2),
                               ("solve_filtering_plus_bwc_chars",1),
                               ("solve_sage",1),
                               ("run_sol_sanity_checks",1),
                               ("run_eth_root",0),
                               ("CRT_R_alpha",1),
                               ("eth_root_montgomery",1),
                               ("sol_times_M",2),
                               ("accumulate",2),
                               ("status",2),
                               ("padic_eth_root",1),
                               ("U reconstruction",1),
                               ("u/v*S(m)^d",1)])


    def pretty_print(field):
        indent = timing_info[field]
        if field not in params.timing:
            print("\t"*indent + f"{field}:","Not run")
        else:
            timing = f"{params.timing[field]:9.2f}s"
            field_cputime = f"{field}_cputime"
            if field_cputime in params.timing:
                timing += f" (cputime: {params.timing[field_cputime]:.2f}s)"
            print("\t"*indent + f"{field}:", timing)

    def print_config_info():
        major_message("----Config Summary-----")
        for k, v in config_info.items():
            print(f"{k}: {v}")

    def print_timing_info(start=None, end=None):
        print_config_info()
        fields = list(timing_info.keys())
        if start: start = fields.index(start)
        if end: end = fields.index(end)
        major_message("----Timing Summary-----")
        for field in fields[start:end]:
            pretty_print(field)

    def print_query_info(ps):
        total = 0
        query_fnames = {
            'Rational': 'RQUERIES_TODO',
            'Algebraic': 'AQUERIES_TODO',
            'Extension': 'EXT_QUERIES_TODO'
        }
        print("Queries:")
        for query_type, query_fname in query_fnames.items():
            if query_fname in ps.files and os.path.isfile(ps.files[query_fname]):
                query_file = ps.files[query_fname]
                nqueries = sum(1 for _ in open(query_file))
                total += nqueries
                print(f"\t{query_type}: {nqueries:,}")
            else:
                print(f"\t{query_type}: n/a")
        print(f"\tTotal: {total:,}")

    sys.stdout = Logger(params)
    sys.stderr = Logger(params, isstdout=False)

    params.save_to_file()

    if args.step == "info":
        print(params)
        sys.exit(0)

    if args.step == "update_alg_unlinked":
        update_collection_of_unlinked_ideals(params, 'alg')
        sys.exit(0)

    if args.step == "update_ext_unlinked":
        update_collection_of_unlinked_ideals(params, 'ext')
        sys.exit(0)

    if args.step == "gen":
        gen_N(params)
        sys.exit(0)

    if args.step == "polysel_select":
        polysel_higheffort_select(params)
        print_timing_info(start='polysel_higheffort_select', end='polysel_higheffort_ropt')
        sys.exit(0)

    if args.step == "polysel_ropt":
        polysel_higheffort_ropt(params)
        print_timing_info(start='polysel_higheffort_ropt', end='polysel_higheffort_candidates')
        sys.exit(0)

    if args.step == "polysel_candidates":
        polysel_higheffort_candidates(params, ncandidates=100)
        print_timing_info(start='polysel_higheffort_candidates', end='run_polysel')
        sys.exit(0)

    if args.step == "precomp":
        run_precomp(params)
        print_timing_info(end='run_queries')
        sys.exit(0)

    if args.step == "polysel_only":
        run_polysel(params)
        print_timing_info(end='call_algebraic_query_sieving')
        sys.exit(0)

    if args.step == "polysel_alreadysel":
        run_polysel(params, already_selected=True)
        print_timing_info(end='call_algebraic_query_sieving')
        sys.exit(0)

    if args.step == "precomp_nopolysel":
        run_precomp(params, nopolysel=True)
        print_timing_info(end='run_queries')
        sys.exit(0)

    if args.step == "algebraic_sieving_only":
        call_algebraic_query_sieving(params)
        major_message("Finished query sieving!")
        sys.exit(0)

    if args.step == "extension_only":
        call_fb_extension_sieving(params)
        major_message("Finished extension sieving!")
        sys.exit(0)

    if args.step == "extra_algebraic":
        assert len(params.extra_prefix) > 0, "Need to specify a file prefix for extra sieving!"
        call_algebraic_query_sieving(params, is_extra=True)
        print_timing_info(start='call_algebraic_query_sieving',end='call_fb_extension_sieving')
        sys.exit(0)

    if args.step == "extra_granular_algebraic":
        assert len(params.extra_prefix) > 1
        assert int(params.given_A) > 1
        assert int(params.given_mfb1) > 1
        timeprint("Doing extra granular algebraic sieving.")
        major_message("MAKE SURE THE UNLINKED IDEAS HAVE BEEN UPDATED (via update_alg_unlinked).")
        call_extra_sieving_granular(params,
                                    params.extra_prefix,
                                    int(params.given_A), int(params.given_mfb1),
                                    is_alg_or_ext=True)
        timeprint("Done with extra granular algebraic sieving.")
        major_message("PLEASE CALL update_alg_unlinked.")
        sys.exit(0)

    if args.step == "extra_granular_filter":
        assert len(params.extra_prefix) > 1
        assert int(params.given_A) > 1
        assert int(params.given_mfb1) > 1
        timeprint("Doing extra granular algebraic (filter) sieving.")
        major_message("MAKE SURE THE UNLINKED IDEAS HAVE BEEN UPDATED (probably via filter_only).")
        call_extra_sieving_granular(params,
                                    params.extra_prefix,
                                    int(params.given_A), int(params.given_mfb1),
                                    is_alg_or_ext=True, is_filter=True)
        timeprint("Done with extra granular algebraic (filter) sieving.")
        sys.exit(0)

    if args.step == "extra_granular_extension":
        assert len(params.extra_prefix) > 1
        assert int(params.given_A) > 1
        assert int(params.given_mfb1) > 1
        timeprint("Doing extra granular extension sieving.")
        major_message("MAKE SURE THE UNLINKED IDEAS HAVE BEEN UPDATED (via update_ext_unlinked).")
        call_extra_sieving_granular(params,
                                    params.extra_prefix,
                                    int(params.given_A), int(params.given_mfb1),
                                    is_alg_or_ext=False)
        timeprint("Done with extra granular extension sieving.")
        major_message("PLEASE CALL update_ext_unlinked.")
        sys.exit(0)

    if args.step == "filter_only":
        MM = filter_relation_file(params)
        major_message("matrix:", MM)

        if 'n1024' not in params.files['RENUMBERFILE']:
            fast_persistent_save(MM, params.dirs['TEMP_OUTPUT_DIR']+"MM.sobj")

        if params.cado_nfs_filter:
            build_killer_rels_dict(params, MM)
        print_timing_info(start='filter_relation_file',end='run_queries')
        sys.exit(0)

    if args.step in ["queries", "rqueries", "aqueries", "extqueries"]:
        run_queries(params, [args.step])
        print_timing_info(start='run_queries',end='run_indiv')
        print_query_info(params)
        sys.exit(0)

    if args.step == "descent_large_init":
        result = run_descent_ecm_init(params)
        # TODO: timing info
        sys.exit(0)

    if args.step == "descent_large_q":
        assert 'slurm.numjobs' in params.parameters
        result = call_descent_large_q_slurm(params, params.use_descent_init_file)
        # TODO: timing info
        sys.exit(0)

    if args.step == "cleanup_large_q":
        do_cleanup_large_q(params)
        sys.exit(0)

    if args.step == "descent_large_las_start":
        assert len(params.use_descent_init_file) > 5
        descent_large_las.start_descent_large_las(params, params.use_descent_init_file)
        sys.exit(0)

    if args.step == "descent_large_las_rebalance":
        assert len(params.use_descent_init_file) > 5
        with open(params.use_descent_init_file, "r") as fp:
            init_dict = json.load(fp)

        seed_ten = str(init_dict['seed'])[:10]
        descent_large_las.rebalance_descent_large_las(params, seed_ten, params.descent_counter)
        sys.exit(0)

    if args.step == "descent_large_las_todo":
        assert len(params.todofile_glob) > 5
        descent_large_las.todofile_descent_large_las(params, params.todofile_glob)
        sys.exit(0)

    if args.step == "descent_large_las_finish":
        assert len(params.use_descent_init_file) > 5
        descent_large_las.finish_descent_large_las(params, params.use_descent_init_file)
        sys.exit(0)

    if args.step == "descent_bottom":
        assert len(params.use_descent_init_file) > 5
        assert params.descent_strategy in ['P','C','C1','C2']
        descent_large_las.descent_rock_bottom(params, params.use_descent_init_file, params.descent_strategy)
        sys.exit(0)

    if args.step == "descent_large_las_extract":
        assert len(params.use_descent_init_file) > 5
        # a lot is assumed here because it's a very particular use case
        with open(params.use_descent_init_file, "r") as fp:
            init_dict = json.load(fp)
        seed_ten = str(init_dict['seed'])[:10]
        working_dir = params.dirs['DESC'] + 'seed' + seed_ten + "/"
        assumed_file = working_dir + "desc.total.rels"
        desired_lpb0 = params.parameters['LPB0']
        desired_lpb1 = params.parameters['LPB1']

        # desired_lpb0 = params.parameters['LAS_DESCENT_UNTIL_LPB0']
        # desired_lpb1 = params.parameters['LAS_DESCENT_UNTIL_LPB1']

        desired_lpb0 = 36
        desired_lpb1 = 35

        descent_large_las.extract_outstanding_qs(params, desired_lpb0, desired_lpb1, assumed_file)
        sys.exit(0)

    if args.step == "indiv":
        if len(params.use_descent_init_file) > 2:
            # it's an existing file
            existing_init = params.use_descent_init_file
        else:
            existing_init = None
        result = run_indiv(params, existing_init_data=existing_init)
        if result is True:
            linalg_output = run_linalg(params)
            found = run_eth_root(params, linalg_output)
        elif result is not None:
            found = result
        else:
            timeprint("Descent failed")
            print_timing_info(start='run_indiv')
            sys.exit(1)
        check_solution(params,found)
        print_timing_info(start='run_indiv')
        print_query_info(params)
        sys.exit(0)

    if args.step == "indiv_only":
        if len(params.use_descent_init_file) > 2:
            # it's an existing file
            existing_init = params.use_descent_init_file
        else:
            existing_init = None
        result = run_indiv(params, existing_init_data=existing_init)
        if result is True or result is not None:
            timeprint("Descent succeeded")
        else:
            timeprint("Descent failed")
        print_timing_info(start='run_indiv',end='run_linalg')
        print_query_info(params)
        sys.exit(0)

    if args.step == "validate_descent":
        run_validate_descent(params, use_target_info_file=params.use_descent_init_file)
        sys.exit(0)

    if args.step == "linalg":
        if len(params.use_descent_init_file) > 2:
            # want to use a specific descent output
            existing_target = params.use_descent_init_file
        else:
            existing_target = None
            with open(params.files['TGT_INFO'], "r") as fp:
                params.target_info = json.load(fp)
        try:
            linalg_output = run_linalg(params, use_target_info_file=existing_target)
            found = run_eth_root(params, linalg_output, use_target_info_file=existing_target)
            check_solution(params,found)
        except Exception as ex:
            print_timing_info(start='run_linalg')
            raise ex
            sys.exit(1)
            os._exit(1)
        check_solution(params,found)
        print_timing_info(start='run_linalg')
        print_query_info(params)

    if args.step == "linalg_nobwc":
        # We are really only running this for the 1024 computation,
        # so some of this may be geared to work in that setup only.
        assert len(params.use_descent_init_file) > 2
        existing_target = params.use_descent_init_file
        linalg_output = run_linalg_nobwc(params, use_target_info_file=existing_target)
        found = run_eth_root(params, linalg_output, use_target_info_file=existing_target)
        check_solution(params, found)
        print_timing_info(start='run_linalg_nobwc')

    if args.step == "root":
        if len(params.use_descent_init_file) > 2:
            # want to use a specific descent output
            existing_target = params.use_descent_init_file

            with open(params.use_descent_init_file, "r") as fp:
                init_dict = json.load(fp)

            seed_ten = str(init_dict['seed'])[:10]
            linalg_output = fast_persistent_load(params.dirs['TEMP_OUTPUT_DIR']+seed_ten+"-linalgoutput.sobj")
        else:
            existing_target = None
            with open(params.files['TGT_INFO'], "r") as fp:
                params.target_info = json.load(fp)
            linalg_output = fast_persistent_load(params.dirs['TEMP_OUTPUT_DIR']+"linalgoutput.sobj")

        try:
            found = run_eth_root(params,linalg_output, use_target_info_file=existing_target)
        except Exception as ex:
            print_timing_info(start="run_eth_root")
            raise ex
            sys.exit(1)
            os._exit(1)
        check_solution(params,found)
        print_timing_info(start="run_eth_root")
        print_query_info(params)

    if args.step == "serialize_ttplus_ttminus":
        serialize_ttplus_ttminus(params)
        print_timing_info()
        print_query_info(params)

    if args.step == "until_root":
        try:
            run_precomp(params)
            run_queries(params)
            result = run_indiv(params)
            if result is True:
                linalg_output = run_linalg(params)
            elif result is not None:
                found = result
            else:
                timeprint("Descent failed")
                print_timing_info()
                sys.exit(1)
                os._exit(1)
        except Exception as ex:
            print_timing_info()
            raise ex
            sys.exit(1)
            os._exit(1)
        print_timing_info()
        print_query_info(params)
        sys.exit(0)
        os._exit(0)

    if args.step == "all":
        # Else, run all steps
        try:
            run_precomp(params)
            result = run_indiv(params)
            if result is True:
                linalg_output = run_linalg(params)
                run_queries(params)     # need extrels.out.indexed to exist
                found = run_eth_root(params, linalg_output)
            elif result is not None:
                found = result
            else:
                timeprint("Descent failed")
                print_timing_info()
                sys.exit(1)
                os._exit(1)
            check_solution(params,found)
        except Exception as ex:
            print_timing_info()
            raise ex
            sys.exit(1)
            os._exit(1)
        print_timing_info()
        print_query_info(params)
        sys.exit(0)
        os._exit(0)
