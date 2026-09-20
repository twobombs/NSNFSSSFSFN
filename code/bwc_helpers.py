from sage.all import *
from tempfile import TemporaryDirectory
import subprocess # for calling bwc binaries
from pathlib import Path
import shutil
import os
#from helpers import timing
from time import time
from cado_nfs_binaries import CadoNFS
from timing import overall_cputime, extract_time, extract_time_from_file, timing

def write_binary_matrix(filename, M, BWCBINDIR):
    """
    Write a matrix M to a file in the binary format bwc.pl expects. Also write the auxiliary files bwc.pl needs.
    Note: When you want to solve things of the form x*M = -y (rather than M*x = -y) you should write M transpose, not M.
    """
    if not M.is_square():
        print("WARNING: writing a non-square matrix. bwc.pl currently fails on these.")
    with open(filename, "wb") as f:
        for i in range(M.nrows()):
            nz = M.row(i).nonzero_positions()
            f.write(int.to_bytes(int(len(nz)), length=4, byteorder='little'))
            for j in nz:
                f.write(int.to_bytes(j, length=4, byteorder='little'))
                f.write(int.to_bytes(int(M[i,j]), length=4, byteorder='little', signed=True))
    p = subprocess.run([
        "time", "-p",
        BWCBINDIR + "/mf_scan2",
        "-withcoeffs",
        "-mfile", filename
    ], stderr=subprocess.PIPE, text=True)
    overall_cputime.add(extract_time(p.stderr))


def write_ascii_vector(filename, vector, modulus):
    with open(filename, 'w') as f:
        f.write(f"{len(vector)} 1 {modulus}\n") # nrows ncols modulus
        for coeff in vector:
            f.write(f"{coeff % modulus}\n")

def make_and_clean(WORKDIR):
    shutil.rmtree(WORKDIR)
    Path(WORKDIR).mkdir(exist_ok=True)

@timing
def solve_linalg_system(params, matrix_filename, vector_filename, m=4, n=4, left=False):
    """
    Given matrix M and vector y, solve for x such that M*x = -y (mod the given modulus)
    i.e., find linear combination of *columns* of M that give -y. (Note the minus sign!)
    For now, bwc.pl fails when M has more columns than rows, but works for square M.

    Note: m and n are parameters for the Block-Wiedemann algorithm. Not sure the optimal values.
    """
    BWCBINDIR = params.dirs['BWCBINDIR']
    thr = params.parameters["bwc.thr.controller"]
    modulus = params.parameters["e"]
    #with TemporaryDirectory(prefix="bwc_") as WORKDIR: # XXX change to delete=False to preserve the workdir # For some reason files are disappearing from the temp directory before they can be read.
    WORKDIR = params.dirs['BWC']
    Path(WORKDIR).mkdir(exist_ok=True)
    make_and_clean(WORKDIR)
    if left:
        left_str = "nullspace=LEFT"
    else:
        left_str = ""

    bwc_cmd = [
        "time", "-p",
        BWCBINDIR + "/bwc.pl",
        ":complete",
        f"prime={modulus}",
        f"matrix={os.path.realpath(matrix_filename)}",
        f"rhs={os.path.realpath(vector_filename)}",
        f"wdir={WORKDIR}",
        f"m={m}", f"n={n}",
        left_str,
        "balancing_options=reorder=columns",
        "save_submatrices=1",
        "verbose_flags=^all-cmdline,^bwc-timing-grids",
        f"thr={thr}"
    ]

    if params.mpi:
        print("Run bwc.pl from solve_linalg_system with MPI enabled.")

        bwc_cmd += [
            f"mpi={params.parameters['bwc.mpi']}",
            "mpi_extra_args='--mca usnic ucx'"
        ]

    print(f"Running bwc command: {' '.join(bwc_cmd)}")

    with open(os.path.join(WORKDIR,"bwc.out"),"w") as outfile, \
         open(os.path.join(WORKDIR,"bwc.err"),"w") as errfile:
        print(f"Writing bwc stdout to {outfile.name}")
        print(f"Writing bwc stderr to {errfile.name}")

        subprocess.call(bwc_cmd,stdout=outfile,stderr=errfile)
        cputime += extract_time_from_file(errfile)

        # output gets written as ascii to {WORKDIR}/K.sols0-1.0.txt
        with open(f"{WORKDIR}K.sols0-1.0.txt",'r') as f:
            # One coefficient per line, with an extra "1" coefficient at the end
            sol = [int(line) for line in f]
            assert sol[-1] != 0, f"bwc.pl returned a solution with RHS coefficient 0. (Hopefully this will stop happening once we run dup1 and dup2 on M before running bwc)"
            scalar = sol[-1]
            sol = vector(Zmod(modulus), sol[:-1])
            sol /= Zmod(modulus)(scalar)
            return sol

def remove_redundant_rows(M):
    """
    Return a largest subset of rows of M that are linearly independent.
    i.e., if M is rank d, return d linearly independent rows of M.
    Also return the indices of said rows.

    This function is used for working around the issue where bwc fails on non-square matrices.
    """
    out = None
    out_indices = []
    for i in range(M.nrows()):
        row = M.row(i)
        if out is None:
            if row:
                out = matrix([row], sparse=True)
                out_indices.append(i)
        else:
            tmp = out.stack(row)
            if tmp.rank() == tmp.nrows():
                out = tmp
                out_indices.append(i)
        if out and out.rank() == M.rank():
            return out, out_indices
    raise ValueError("M is zero (otherwise this should be unreachable)")

def remove_redundant_columns(M):
    """ Same as remove_redundant_rows, but for columns """
    M2, indices = remove_redundant_rows(M.T)
    return M2.T, indices

@timing
def squarify(M,params):
    """
    Remove zero columns and remove most redundant rows.
    The result isn't necessarily full rank (or even square)
    but it's relatively fast and empirically it seems to be good enough for bwc.pl to succeed.
    """
    # FIRST, remove zero columns and make an initial matrix
    col_indices = []
    row_indices = set()
    print(f"Computing rank of M ({repr(M)})...")
    Mrank = M.rank()
    print(f"M is rank {Mrank}")
    for i in range(M.ncols()):
        col = M.column(i)
        if col != 0:
            col_indices.append(i)
            for j in col.nonzero_positions():
                if j not in row_indices:
                    row_indices.add(j)
                    break
    print(f"Number of nonzero columns: {len(col_indices)}")
    print(f"Number of rows so far: {len(row_indices)}")
    if len(row_indices) < Mrank:
        for i in range(M.nrows()):
            if len(row_indices) >= Mrank: break
            if i not in row_indices: row_indices.add(i)
    row_indices = sorted(list(row_indices))
    print("Making Mnew submatrix")
    Mnew = M.matrix_from_rows_and_columns(row_indices, col_indices)
    print("removed {} rows among the first {}".format(
          len(set(range(row_indices[-1]+1))-set(row_indices)),
          row_indices[-1]+1))
    print("removed cols:",
          set(range(M.ncols()))-set(col_indices))
    # Maybe this is already good enough?
    return Mnew, row_indices, col_indices

def solve_system_allatonce(params, M, v, ensure_full_rank=False, **kwargs):
    """
    Solve linear system x*M = -v mod e, doing so all at once: start from a sage matrix and vector rather than from files.
    Include workaround for non-square matrix.
    This function will eventually go away once non-square matrices are okay and once we separate matrix file creation into the precomp phase
    """
    e = params.parameters['e']
    FILES_DIR = params.dirs['TEMP_OUTPUT_DIR']
    BWCBINDIR = params.dirs['CADO_BUILD_DIR']+"linalg/bwc"
    params.dirs['BWCBINDIR'] = BWCBINDIR

    if ensure_full_rank:
        ### First: work around the issue with non-square M

        print("Remove redudant rows starting")
        ts = time()
        # We're solving x*M = -v
        M_fewerrows, row_indices = remove_redundant_rows(M)
        print("Redundant rows",time()-ts)

        recover_x = lambda x_smaller : vector(Zmod(e), M.nrows(), {ri: xi for (ri,xi) in zip(row_indices, x_smaller)})
        # Now we're solving x_smaller * M_fewerrows = -v

        # check that the above code is correct:
        ##xsmaller_test = vector(Zmod(e), [randint(0,e-1) for _ in range(M_fewerrows.nrows())])
        ##assert xsmaller_test * M_fewerrows == recover_x(xsmaller_test) * M
        print("Remove redundant columns starting")
        ts = time()
        M_square, col_indices = remove_redundant_columns(M_fewerrows)
        v_smaller = vector(Zmod(e), [v[i] for i in col_indices])
        # Now we're solving x_smaller * M_square = -v_smaller
        print("Redundant columns",time()-ts)

        assert M_square.is_square()
        assert not M_square.is_singular()
        ### Hooray, we've worked around the non-square M issue
    else:
        print("Squarify starting")
        ts = time()
        M_square, row_indices, col_indices = squarify(M,params)
        print("Squarify",time()-ts)
        recover_x = lambda x_smaller : vector(Zmod(e), M.nrows(), {ri: xi for (ri,xi) in zip(row_indices, x_smaller)})
        v_smaller = vector(Zmod(e), [v[i] for i in col_indices])

    MATRIX_FILE = FILES_DIR + "/M.bin"
    RHS_FILE = FILES_DIR + "/rhs.txt"

    ts = time()
    write_binary_matrix(MATRIX_FILE, M_square, BWCBINDIR)
    write_ascii_vector(RHS_FILE, v_smaller, e)
    print("Writing matrix files",time()-ts)

    print("Calling solve_linalg_system")
    xsmaller = solve_linalg_system(params, MATRIX_FILE, RHS_FILE, **kwargs)
    if not vector(Zmod(e), xsmaller * M_square) == -v_smaller:
        print("*** Linalg error ***")
        print(f"xsmaller * M_square: {list(xsmaller * M_square)[:20]}...")
        print(); print(); print()
        print(f"-vsmaller: {list(-v_smaller)[:20]}...")
        print()
        print(f"Difference has hamming weight {vector(Zmod(e), xsmaller * M_square + v_smaller).hamming_weight()}")
        print(); print(); print()
        raise AssertionError("Linalg error 1. bwc gave a wrong solution.")

    if ensure_full_rank and not vector(Zmod(e), xsmaller * M_fewerrows) == -vector(Zmod(e), v):
        print("*** Linalg error ***")
        print(f"xsmaller * M_fewerrows: {list(xsmaller * M_fewerrows)[:20]}...")
        print()
        print(f"-v: {list(-vector(Zmod(e), v))[:20]}...")
        print()
        print(f"Difference has hamming weight {(vector(Zmod(e), xsmaller * M_fewerrows) + vector(Zmod(e), v)).hamming_weight()}")
        print()
        print("Target vector in row span of M_fewerrows:", vector(Zmod(e), v) in M_fewerrows.row_space())
        print("Target vector in row span of M:", vector(Zmod(e), v) in M.row_space())
        print(); print(); print()
        raise AssertionError("Linalg error 2. Target vector might not be in span of M.")

    x = recover_x(xsmaller)
    assert len(x) == M.nrows()
    if not vector(Zmod(e), x * M) == -vector(Zmod(e), v):
        print("*** Linalg error ***")
        print(f"x * M: {list(x * M)[:20]}...")
        print(); print(); print()
        print(f"-v: {list(-vector(Zmod(e), v))[:20]}...")
        print(); print(); print()
        print((x*M+v).sparse_vector().dict())
        raise AssertionError("Linalg error 3")

    return x


@timing
def solve_system_bwc_from_filtered(params,
                                   MM, ST_alg_vector,
                                   S_block, C_block,
                                   ST_list):

    # We default to m=n=6, because we want as least as many as the rank
    # of the unit group, which is less than len(C_block).
    m = MM.M.params.m
    n = MM.M.params.n

    FILES_DIR = params.dirs['TEMP_OUTPUT_DIR']
    WORKDIR = os.path.join(FILES_DIR,"bwc/")
    RHS_FILE = os.path.join(FILES_DIR, "rhs.txt")
    e = int(params.parameters['e'])
    Path(WORKDIR).mkdir(exist_ok=True)
    make_and_clean(WORKDIR)

    rhs = vector([ST_alg_vector[xj] for xj in MM.column_expand_map()])

    write_ascii_vector(RHS_FILE, rhs, e)

    # XXX: untested code path for mpi=True
    bwc_add_args = []
    if params.mpi:
        print("Run bwc.pl from unfiltered with MPI enabled.")
        bwc_add_args += [
            f"--mpi {params.parameters['bwc.mpi']}",
            "mpi_extra_args='--mca usnic ucx'"
        ]

    with open(os.path.join(WORKDIR,"bwc.out"),"w") as outfile, \
         open(os.path.join(WORKDIR,"bwc.err"),"w") as errfile:

        print(f"Writing bwc stdout to {outfile.name}")
        print(f"Writing bwc stderr to {errfile.name}")

        CadoNFS("linalg/bwc/bwc.pl",
                ":complete",
                "--matrix", os.path.realpath(MM.get_matrix_file(params)),
                "--prime", params.parameters['e'],
                "--nullspace", "LEFT",
                "--rhs", os.path.realpath(RHS_FILE),
                "--wdir", WORKDIR,
                "--solutions", f"0-{len(C_block)}",
                f"m={m}", f"n={n}",
                "balancing_options=reorder=columns",
                "verbose_flags=^all-cmdline,^bwc-timing-grids",
                "--thr", params.parameters["bwc.thr.controller"],
                *bwc_add_args,
                capture=outfile,
                stderr=errfile
            )

    Ze = MM.base_ring()

    sol_matrix_rows = []
    for j in range(len(C_block)):
        with open(f"{WORKDIR}/K.sols{j}-{j+1}.0.txt",'r') as f:
            sol_matrix_rows.append([int(line) for line in f])

    # we now have a space of solutions to v * vjoin(M,RHS) = 0. Check which
    # of those happen to also be in the left nullspace of vjoin(S_block,
    # C_block).

    sol_matrix = matrix(Ze, sol_matrix_rows)
    SC = block_matrix(2,1,[S_block, matrix([C_block])])
    print(sol_matrix * SC)
    print((sol_matrix * SC).rank())
    ## nh, nv = (int(c) for c in params.parameters["bwc.thr.controller"].split('x'))
    ## # sol_matrix = matrix(Integers(67), [vector([int(x) for x in open(f"data/bwc/K.sols{j}-{j+1}.0.txt")]) for j in range(5)])
    ## from cado_sage import bwc
    ## par = bwc.BwcParameters(m=m,n=n,p=ZZ(e),nullspace='left')
    ## B = bwc.BwcMatrix(par, "data/n60.aqrels.out.indexed.matrix.bin")
    ## B.read()
    ## if not os.path.exists(f'data/n60.aqrels.out.indexed.matrix.{nh}x{nv}'):
    ##     os.symlink(f'bwc/n60.aqrels.out.indexed.matrix.{nh}x{nv}',
    ##                f'data/n60.aqrels.out.indexed.matrix.{nh}x{nv}')
    ## B.fetch_balancing(nh, nv)
    assert sol_matrix * block_matrix(2, 1, [MM.M.M, matrix([rhs])]) == 0
    print(sol_matrix[:,-1:])


    # I'm confused. I _think_ that there _has_ to be a solution, and so
    # this matrix can't be full rank. And yet, it is.

    raise NotImplementedError("TBC")

    # rhs = [Ze(c) for c in open('data/rhs.txt').readlines()[1:]]
    # S = matrix(Integers(67), [vector([int(x) for x in open(f"data/bwc/K.sols{j}-{j+1}.0.txt")]) for j in range(5)])




if __name__ == "__main__":
    print("Testing with a random dense 30-by-30 matrix mod 65537...")
    M = random_matrix(Zmod(65537), 30, 30)
    s = vector(Zmod(65537), [randint(0,65536) for _ in range(30)])
    target = -s * M
    with TemporaryDirectory(prefix="bwctest_", delete=False) as WORKDIR:
        print(f"Using directory {WORKDIR}/ as workdir")
        #print("Writing matrix...")
        #write_binary_matrix(WORKDIR + "/MT.bin", M.T, BWCBINDIR=BWCBINDIR)
        #print("Writing target vector...")
        #write_ascii_vector(WORKDIR + "/vec.txt", target, 65537)
        #print("Solving...")
        #sol = solve_linalg_system(WORKDIR+"/MT.bin", WORKDIR+"/vec.txt", 65537)
        #assert sol * M == -target, "solution incorrect"
        #print("Success!")

        sol = solve_system_allatonce(params, M, target, 65537, WORKDIR)
        assert sol * M == -target, "solution incorrect"
        print("Success!")
