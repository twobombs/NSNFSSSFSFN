from sage.all import *
import helpers
from helpers import major_message
from timing import *
from collections import defaultdict
from sage.misc.persist import SagePickler, SageUnpickler

from functools import partial

from cado_sage import CadoPolyFile
from cado_sage import CadoIndexFile
from cado_sage import CadoNumberTheory

from candy import major_message, warning_message, error_message
import montgomery_reduction_new
from montgomery_reduction_new import CadoMontgomeryReductionProcess, AccumulateLog

from misc_tools import fast_persistent_save, fast_persistent_load

from tocfile import RandomAccessIndexedRelations
from random_access_renumber import RandomAccessRenumberTable

import padic_eth_root
import concurrent.futures
import glob


def sanity_check_eth_root_input(params, linalg_output):
    N = params.poly.N
    target = params.target
    ZN = Integers(N)
    target_info = params.target_info
    mask = ZN(target_info['mask'])
    u = target_info['u']
    v = target_info['v']
    e = params.parameters['e']
    assert params.target  * mask**e == ZN(u/v)

    drels = target_info['DRELS_FILE']

    fac_uv = u.factor() / v.factor()

    epsilon = fac_uv.unit()

    rat_primes_queried = [(p,k) for p,k in fac_uv if p < params.BOUNDR]
    rat_primes_descended = [(p,k) for p,k in fac_uv if p >= params.BOUNDR]


def column_to_sage_ideal_parallelizable_new(col_index, explain_renumber_filename, polyfile):
    global my_cadopoly
    # Reading a CadoPolyFile takes around a tenth of a second, so we only want to make a new CadoPolyFile once per worker process instead of once per function call
    if 'my_cadopoly' not in globals():
        my_cadopoly = CadoPolyFile(polyfile)
        my_cadopoly.read()
    with RandomAccessRenumberTable(explain_renumber_filename, my_cadopoly) as R:
        ideal = R.column_to_sage_ideal(col_index)
    # Calling hash(ideal) here caches the hash of the ideal, which is somewhat slow to compute and will be needed later when using the ideal as a dict key.
    # This way it gets done in parallel instead of singlethreaded
    hash(ideal)
    return SagePickler.dumps(ideal)

def column_to_sage_ideal_parallelizable(col_index, ideal, has_merged_J, poly, side_hint=None):
    """
    Parallelization-friendly version of CadoExplainRenumberFile.column_to_sage_ideal
    from helpers.py
    """
    parser, side, Idata, _col_index = ideal

    if parser == 'J' and has_merged_J:
        # This case is special, really.
        assert side_hint is not None    # what can we do?
        assert len(side) == len(_col_index)
        assert _col_index == (_col_index[0],)*len(side)
        side = side[side_hint]
        _col_index = _col_index[side_hint]

    assert _col_index == col_index

    # XXX this assert is a bit excessive in full generality.
    # Perhaps we'd like to assert side == side_hint, at most.
    assert side == 1

    K = poly.K
    J = poly.nt.J()
    OK = poly.nt.maximal_orders()

    if parser == 'J':
        I = J[side]
    elif parser == 'rat':
        (p,) = Idata
        I = OK[side].fractional_ideal(p)
    elif parser == 'proj':
        (p,) = Idata
        I = OK[side].fractional_ideal(p) + J[side]
    elif parser == 'easy':
        (p, r) = Idata
        I = OK[side].fractional_ideal(p, K[side].gen() - r) * J[side]
    elif parser == 'generic':
        (p, theta) = Idata
        I = OK[side].fractional_ideal(p, theta)
    else:
        raise AssertionError("Unknown parser:", parser)
    return SagePickler.dumps(I)


def compute_logmap_parallelizable(abm, *, alpha, real_roots, complex_roots, prec):
    """
    Parallelizable version of LogMap.
    abm: tuple ((a,b),m), where we want the log of (a - b*alpha)**m
    real_roots, complex_roots: roots of defining polynomial f,
        as a list of elements of RealIntervalField / ComplexIntervalField of the desired precisionn
    """
    (a,b), m = abm
    g = (a - b*alpha).polynomial()
    R = RealIntervalField(prec)
    return vector([R(abs(g(root))).log() for root in real_roots] +
                  [2 * R(abs(g(root))).log() for root in complex_roots]) * m

# Create this object so that I don't have to copy-paste the complete
# boilerplate all the time

class eth_root_montgomery_object(object):
    def __init__(self, params, linalg_output):
        self.params = params
        self.linalg_output = linalg_output
        if not linalg_output.sol.is_sparse():
            timeprint("WARNING: linalg_output.sol was stored as a dense vector; sparse would be better")
        timeprint("calling sol.sparse_vector().lift_centered()")
        self.sol           = linalg_output.sol.sparse_vector().lift_centered()
        timeprint("done with lift_centered")
        self.ST_list       = linalg_output.ST_list
        self.ST_alg_vector = linalg_output.ST_alg_vector
        self.T_list        = linalg_output.T_list
        self.row_to_aquery = linalg_output.row_to_aquery
        self.S_rat_vector  = linalg_output.S_rat_vector
        # self.idx_rel_file  = linalg_output.indexed_relations_file
        self.idx_rel_file  = params.files['AQRELS_FILE'] + ".indexed"

        self.e = params.parameters['e']
        self.poly = params.poly
        self.N = params.poly.N
        self.m = params.poly.m
        self.ZN = Integers(self.N)
        self.side = 1
        self.f = params.poly.f[self.side]
        self.K = params.poly.K[self.side]
        self.nt = params.poly.nt[self.side]
        self.alpha = self.K.gen()
        self.d = self.K.degree()
        if not params.montgomery_new_renumber:
            self.R = params.R

        self.extra_checks = False

        self.OK = self.K.ring_of_integers()

        # these are computed by functions below
        self.U_list = None
        self.sol_M = None
        self.ideals = None
        self.valuations = None
        self.big_power_valuations = None

        # self.big_power is only computed when self.extra_checks == True,
        # which is just for debugging (and most of what it does is
        # horribly expensive anyway)
        self.big_power = None
        timeprint("eth_root_montgomery_object init done")

    # Here we compute sol * M, not reducing mod e
    # This code is basically just run_sol_sanity_checks from helpers.py, but not reducing mod e
    @timing
    def sol_times_M(self):
        timeprint("Computing sol*M...")
        timeprint(f"Length of sol (and thus relevant length of indexed_relations) is {len(self.sol)}")
        timeprint(f"Number of nonzero elements of sol is {len(self.sol.nonzero_positions())}")
        assert self.sol.is_sparse()
        self.sol_M = vector(ZZ, len(self.ST_alg_vector), sparse=True)
        with RandomAccessIndexedRelations(self.idx_rel_file) as rels:
            if self.params.montgomery_new_renumber:
                with RandomAccessRenumberTable(self.params.files['EXPLAIN_RENUMBER_FILE'], self.poly) as R:
                    for pos in self.sol.nonzero_positions():
                        sol_r = self.sol[pos]
                        irel = rels[pos]
                        for ii in irel.indices:
                            col = R.renumber_to_column(ii)
                            assert col >= 0
                            self.sol_M[col] += sol_r
            else:
                for pos in self.sol.nonzero_positions():
                    sol_r = self.sol[pos]
                    irel = rels[pos]
                    for ii in irel.indices:
                        col = self.params.R.renumber_to_column(ii)
                        assert col >= 0
                        self.sol_M[col] += sol_r
        e = self.e
        timeprint("sol*M assert")
        assert (self.sol_M + self.ST_alg_vector).change_ring(Integers(e)).is_zero()
        return self.sol_M

    def compute_target_ideals_and_valuations(self):
        """
        based on ST_alg_vector and sol_M, compute the target ideals and
        target valuations.
        """
        # We want the eth root of ST U, which corresponds to the vector
        # ST_alg_vector + sol_M. We take its vector of valuations
        V = self.ST_alg_vector.sparse_vector().change_ring(ZZ)
        W = V + self.sol_M
        assert all(W[x] % self.e == 0 for x in W.nonzero_positions()), "ST_alg_vector + sol_M not divisible by e, I must have constructed the target wrong or flipped a sign somewhere"
        # We actually want the valuations of the eth root, not the eth
        # power, so divide the vector by e
        W = (W / self.e).change_ring(ZZ)

        # We now need the ideals that correspond to the entries of the target_valuations vector (i.e., that correspond to the columns of our matrix)
        # Look up the ideals corresponding to each nonzero column
        timeprint("Looking up relevant sage ideals...")
        nz = set(W.nonzero_positions())
        # always include J
        nz.add(0)
        nz = sorted(list(nz))

        # TODO: MPI deactivated because it runs into the following error when ran across machines:
        # PMIX ERROR: PMIX_ERROR in file ../../../../3rd-party/prrte/src/prted/pmix/pmix_server_dyn.c at line 1112
        if False: #self.params.mpi:
            from mpi4py.futures import MPIPoolExecutor
            max_workers = self.params.parameters["mpi.thr"]
            ProcessPool = MPIPoolExecutor
        else:
            max_workers = self.params.nthreads
            ProcessPool = concurrent.futures.ProcessPoolExecutor

        max_workers = 11

        if self.params.montgomery_new_renumber:
            timeprint(f"Starting parallel column-to-sage-ideal ({max_workers} max workers) WITH random access renumber table")
            explain_renumber_filename = self.params.files['EXPLAIN_RENUMBER_FILE']
            polyfile = self.params.files['POLYFILE']
            with ProcessPool(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(
                        column_to_sage_ideal_parallelizable_new,
                        col,
                        explain_renumber_filename,
                        polyfile
                    )
                    for col in nz
                ]
            self.ideals = [SageUnpickler.loads(future.result()) for future in futures]
            timeprint(f"Finished parallel column-to-sage-ideal")
        else:
            timeprint(f"Starting parallel column-to-sage-ideal ({max_workers} max workers) WITHOUT random access renumber table")
            R = self.params.R
            timeprint("Finished loading R")
            with ProcessPool(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(
                        column_to_sage_ideal_parallelizable,
                        col,
                        R._ideals[R.column_to_renumber(col)],
                        R.has_merged_J,
                        R.poly,
                        self.side
                    )
                    for col in nz
                ]
                timeprint("Launched the futures")
            timeprint("Starting to load self.ideals")
            self.ideals = [SageUnpickler.loads(future.result()) for future in futures]
            timeprint("Finished loading self.ideals")
            timeprint(f"Finished parallel column-to-sage-ideal")

        self.valuations = vector([W[i] for i in nz])
        self.big_power_valuations = self.e * self.valuations

        timeprint(f"len(ideals) = len(valuations) = {len(nz)} (ideals that never appear are not counted)."
              f" (len(target vector) was {len(W)})")


    def column_to_ideal(self, col):
        return self.params.R.column_to_sage_ideal(col,
                                                  side_hint=self.side)

    def compute_U_list(self):
        """
        return the list of (a,b) pairs, together with their exponent,
        that participate in the solution of the linear system
        """
        self.U_list = [(self.row_to_aquery[i], self.sol[i]) for i in
                       self.sol.nonzero_positions()]
        return self.U_list

    def compute_big_power(self):
        timeprint("Doing a bunch of really incredibly slow debugging checks")
        # print(ST_list)
        # print(len(ST_alg_vector))
        # print(len(ST_alg_vector.nonzero_positions()))
        ST_ideal_from_vector = prod(
                self.column_to_ideal(col)**self.ST_alg_vector[col]
                for col in self.ST_alg_vector.nonzero_positions())

        ideals_dict = {self.column_to_ideal(col):col
                       for col in range(self.params.R.number_of_algebraic_columns())}
        # print(ST_ideal_from_vector)
        big_ST = prod((a - b*self.alpha)**m
                      for ((a,b),m) in self.ST_list.items())
        ST_ideal_from_list = self.OK.fractional_ideal(big_ST)
        # print(ST_ideal_from_list)
        if ST_ideal_from_list == ST_ideal_from_vector:
            timeprint("ST_ideal_from_list and ST_ideal_from_vector agree!")
        else:
            I_ST = ST_ideal_from_list / ST_ideal_from_vector
            timeprint([(ideals_dict[aa],k) for aa,k in I_ST.factor()])

        # print(U_list)
        # print(len(sol_M))
        # print(sol_M)
        # print(len(sol_M.nonzero_positions()))
        # print(sol)
        U_ideal_from_vector = prod(
                self.column_to_ideal(col)**self.sol_M[col]
                for col in self.sol_M.nonzero_positions())
        # print(U_ideal_from_vector)
        big_U = prod((a - b*self.alpha)**m for ((a,b),m) in self.U_list)
        U_ideal_from_list = self.OK.fractional_ideal(big_U)
        # print(U_ideal_from_list)
        if U_ideal_from_list == U_ideal_from_vector:
            timeprint("U_ideal_from_list and U_ideal_from_vector agree!")
        else:
            I_U = U_ideal_from_list / U_ideal_from_vector
            timeprint([(ideals_dict[aa],k) for aa,k in I_U.factor()])

        self.big_power = big_ST * big_U

    def check_power(self):
        q = ZZ(1)
        while True:
            try:
                while not q.is_prime(proof=False) or not self.f.roots(GF(q)):
                    q += self.e
                r = self.f.roots(GF(q), multiplicities=False)[0]
                z = prod([(a - b*r)**m for ((a,b),m) in self.ST_list.items()])
                z *= prod([(a - b*r)**m for ((a,b),m) in self.U_list])
                # trigger a zero division error if z == 0
                z = 1 / z
                timeprint(f"Doing a quick {self.e}-th power check"
                      f" mod q={q}=1+{self.e}*{(q-1)//self.e}")
                print("q", str(q))
                print("r", str(r))
                print("self.e", str(self.e))
                print("z", str(z))
                ret = z **((q-1)//self.e)
                print("ret", str(ret))
                return z **((q-1)//self.e)
            except ZeroDivisionError:
                q += self.e

    def check_invariant(self, MM):
        timeprint("Doing consistency check on MM (LONG!)")
        for i,I in enumerate(self.ideals):
            if i == 0 and not self.f.is_monic():
                continue
            else:
                # vb = big_power.valuation(I)
                vb = self.big_power_valuations[i]
                vm = MM.current.get(I,0)
                vg = MM.OK.fractional_ideal(MM.accumulated.prod()).valuation(I)
                if vb == (vm+vg)*self.e:
                    continue
                else:
                    timeprint(f"Weird valuation at ideal {i} [[{I}]]: {vb/e},{vm}")
                    assert False
        timeprint("Doing consistency check on MM: OK")

@timing
def eth_root_montgomery(params, linalg_output):
    """
    Here are all the relevant algebraic numbers and vectors and whatnot:

    ST: an algebraic number that comes from descent. Smooth over our factor base. We get it in two forms:
    ST_alg_vector: vector giving the valuation of ST at each ideal (as an integer vector, not mod e). Output by truncate_S. This is the target vector for linalg (except for character columns and reduction mod e).
    ST_list: list of (a,b,m) pairs where prod (a - b x)^m = ST.

    sol: vector over Zmod(e) that is the output of linalg. sol * M = -ST mod e.
    sol * M  corresponds to  U(alpha) = prod_i (x_i - y_i alpha)^m_i.  These (x - y alpha) are all in our query base (and correspond to rows of M). The (x,y,m) here are NOT ST_list!

    The value we want to take the eth root of is ST(alpha) U(alpha), which is guaranteed to be an eth power in the number field.

    To compute the eth root, we call CadoMontgomeryReductionProcess(poly, side, ideals, valuations, log_embeddings)
    Here ideals and valuations are lists or vectors such that zip(ideals, valuations) gives the ideal factorization of **the eth root** of STU.
    i.e., at each ideal:
        take the valuation of ST at that ideal (i.e, take ST_alg_vector)
        add the valuation of U at that ideal (i.e., add sol * M where we **don't** reduce mod e and don't include character columns)
        the result is divisible by e; **divide the result by e**
        then give this to CadoMontgomeryReductionProcess

    log_embeddings: floating-point vector that is the log of the complex embeddings of the eth root we're taking.
        We compute it as follows: we have some product expression for our target
            STU(alpha) = prod_i (c_i - d_i alpha)^m_i
        that we get by combining ST_list and the (x,y,m) from sol*M.
        Then we just take the sum of the log-embedding map (using CadoNumberTheory) of each term.

    Once we have ideals, valuations, and log_embeddings, the rest is just copied from montgomery.sage.

    Eventually we get an algebraic number R such that R(alpha)**e = STU(alpha)
    """

    MTY = eth_root_montgomery_object(params, linalg_output)
    sol_M = MTY.sol_times_M()
    MTY.compute_target_ideals_and_valuations()

    # Finally, we will need to compute the log-embeddings. We might need
    # to do so several times if we need to restart with increased
    # precision.  To prepare for this computation, we need a product
    # expression for STU.  We already have a product expression for ST
    # (namely, ST_list).  For U (which corresponds to sol*M), we look up
    # the (a,b) pair for each row of M and use the exponent in sol

    U_list = MTY.compute_U_list()

    # We now have all the parameters we need, the rest just follows
    # cado/sqrt/montgomery.sage

    if MTY.extra_checks:
        MTY.compute_big_power()

    assert MTY.check_power() == 1

    f = MTY.f
    nt = MTY.nt
    alpha = MTY.alpha
    e = MTY.e

    # FIXME J_valuation_inconsistency
    # to be DE-activated someday!
    if not f.is_monic():
        J = MTY.ideals[0]
        assert MTY.ideals[0].norm() == f.leading_coefficient().abs()
        MTY.valuations[0] *= -1

    major_message("Starting montgomery reduction process")

    prec = params.parameters['MONTGOMERY_ROOT_PRECISION']
    mnb = int(params.parameters.get('montgomery.max_norm_bits', 1000))
    loop_until_nbits = int(params.parameters.get('montgomery.loop_until_nbits', 10))
    final_nbits = int(params.parameters.get('montgomery.final_nbits', 8))
    parallel_until_nbits = int(params.parameters.get('montgomery.parallel_until_nbits', 200000))

    # If cmd arguments were given we overwrite the above
    if int(params.given_prec) > 0:
        prec = int(params.given_prec)
        major_message(f"Using prec={prec}")
    if int(params.given_mnb) > 0:
        mnb = int(params.given_mnb)
        major_message(f"Using mnb={mnb}")
    if int(params.given_lub) > 0:
        loop_until_nbits = int(params.given_lub)
        major_message(f"Using loop_until_nbits={loop_until_nbits}")
    if int(params.given_fin) > 0:
        final_nbits = int(params.given_fin)
        major_message(f"Using final_nbits={final_nbits}")
    if int(params.given_pub) > 0:
        parallel_until_nbits = int(params.given_pub)
        major_message(f"Using parallel_until_nbits={parallel_until_nbits}")

    while True:
        try:
            timeprint("Computing log-embeddings...")
            timeprint("embedding precision:", prec)
            timeprint("Constructing LogMap...")
            LogMap = nt.LogMap(prec, interval_based=True)
            timeprint("Computing log-embeddings of ST and U...")
            # Now use the log map to compute log-embeddings of STU
            try:
                import mr4mp
                from operator import add
                all_roots = [x[0] for x in f.roots(ComplexIntervalField(prec))]
                real_roots = [x for x,y in all_roots if 0 in y]
                complex_roots = [x for x in all_roots if x.imag() > 0] # actually just one root from each complex-conjugate pair
                log_embeddings = mr4mp.pool().mapreduce(
                    partial(compute_logmap_parallelizable,
                        alpha=alpha,
                        real_roots=real_roots,
                        complex_roots=complex_roots,
                        prec=prec,
                    ),
                    add,
                    list(MTY.ST_list.items())
                )
                log_embeddings += mr4mp.pool().mapreduce(
                    partial(compute_logmap_parallelizable,
                        alpha=alpha,
                        real_roots=real_roots,
                        complex_roots=complex_roots,
                        prec=prec,
                    ),
                    add,
                    MTY.U_list
                )
            except ModuleNotFoundError:
                timeprint("WARNING: Install Python package 'mr4mp' for faster mapreduce operation")
                log_embeddings = sum([LogMap(a - b * alpha)*m
                                      for ((a,b),m) in MTY.ST_list.items()])
                log_embeddings += sum([LogMap(a - b * alpha)*m
                                       for ((a,b),m) in MTY.U_list])
            # Again, we actually want the log-embeddings for the eth-root, not
            # the eth-power, so divide by e
            log_embeddings /= e
            timeprint("Done computing log-embeddings of ST and U.")

            DO_LOGGING = False
            if DO_LOGGING:
                accumulate_logfile = params.dirs['TEMP_OUTPUT_DIR'] + "/mont2/accumulate"
            else:
                accumulate_logfile = None

            MM = CadoMontgomeryReductionProcess(MTY.poly, MTY.side,
                                                MTY.ideals, MTY.valuations,
                                                LogMap, log_embeddings, accumulate_logfile, params=params)
            MM.status()

            MTY.check_invariant(MM)

            nideals_history = defaultdict(int)
            b = infinity

            RECOVER_FROM_LOGS = False
            if RECOVER_FROM_LOGS:
                accumulate_logfile = params.dirs['TEMP_OUTPUT_DIR'] + "/mont/accumulate"
                logfiles = glob.glob(f"{accumulate_logfile}.*.sobj")
                major_message("Recovering progress from existing logfiles, of which we found {len(logfiles}}")

                for lf in logfiles:
                    AL = fast_persistent_load(lf)
                    MM.accumulate(AL.g, AL.num_or_den, hint=AL.hint, dolog=False)
                    MM.status()

                major_message("Done recovering from logs. Continuing on with the usual program.")

            if params.montgomery_parallel:
                # First do a parallel loop, if applicable

                cpucount = int(params.parameters.get('montgomery.cpucount', 2))

                if int(params.given_t) > 0:
                    major_message(f"Using cpucount={params.given_t}")
                    cpucount = int(params.given_t)

                while MM.nplus[1] + MM.nminus[1] > parallel_until_nbits:
                    b = min(round(0.8 * max(MM.nplus[1], MM.nminus[1])), mnb, b)
                    MM.one_internally_parallel_reduction_step(b, cpucount)
                    MM.status()
                    nideals_history[MM.nplus[0]+MM.nminus[0]] += 1
                    if nideals_history[MM.nplus[0]+MM.nminus[0]] > 5:
                        major_message("We seem to be stuck in a loop")
                        break


            # The default (sequential) loop
            while MM.nplus[1] + MM.nminus[1] > loop_until_nbits:
                b = min(round(0.8 * max(MM.nplus[1], MM.nminus[1])), mnb, b)
                MM.one_reduction_step(b, skip_asserts=True, skip_embeddings_reduction=True)
                MM.status()
                nideals_history[MM.nplus[0]+MM.nminus[0]] += 1
                if nideals_history[MM.nplus[0]+MM.nminus[0]] > 5:
                    major_message("We seem to be stuck in a loop")
                    timeprint("Adding some randomness to try to get out...")

                    # reset the counter
                    nideals_history[MM.nplus[0]+MM.nminus[0]] = 0
                    b = min(round(0.8 * max(MM.nplus[1], MM.nminus[1])), mnb, b)
                    MM.one_reduction_step(b, skip_asserts=True,
                                            skip_embeddings_reduction=True, escape_loop=True)
                    MM.status()

            if MM.done():
                break

            # one last round.
            MM.one_reduction_step(final_nbits, skip_asserts=True, skip_embeddings_reduction=True)
            MM.status()

            if MM.done():
                break

            while MM.nminus[0] != 0:
                MM.one_reduction_step(0,
                                      all_ideals_at_once=True,
                                      skip_asserts=True,
                                      skip_embeddings_reduction=True)
                MM.status()
            break
        except montgomery_reduction_new.InsufficientPrecision:
            prec *= 2
            if prec > 2000:
                error_message(f"Bailing out, precision={prec} seems too large")
                raise
            warning_message(f"Starting over with precision={prec}")

    fast_persistent_save(MM.accumulated, params.dirs['TEMP_OUTPUT_DIR']+"gamma-fac.sobj")
    gamma = MM.accumulated.prod()
    fast_persistent_save(gamma, params.dirs['TEMP_OUTPUT_DIR']+"gamma.sobj")

    if not MM.current and max(MM.bound_on_coefficients_of_remaining_part()) == 1:
        timeprint("We're almost certainly done!")
        timeprint("Trying a p-adic root lift just to clear improbable leftovers")

    # Do a round of p-adic root computation with what we have computed so
    # far. This should be really quick. If everything went well, we're
    # actually just computing a root of 1 here, no big deal (so delta=1).
    # But it's not inimaginable that there's still a tiny bit that we
    # need to compute, in which case this is the way to go.
    delta = padic_eth_root.padic_eth_root(params,
                                          linalg_output,
                                          lift_centered=True,
                                          pre_multiply_root=(1/gamma).polynomial())
    delta = delta(alpha)
    timeprint(f"delta = {delta}")
    fast_persistent_save(delta, params.dirs['TEMP_OUTPUT_DIR']+"delta.sobj")

    # then (gamma * delta) is an e-th root of big_power
    if MTY.extra_checks:
        timeprint("Final check: ", MTY.big_power / (gamma*delta)**MTY.e, " (both +1 and -1 are fine)")

    ########## above lines are from cado/sqrt/montgomery.sage ##########
    return gamma*delta
