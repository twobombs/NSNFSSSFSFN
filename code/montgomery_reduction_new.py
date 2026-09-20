"""
This file is copied from cado_sage/montgomery_reduction_process.py,
but fixing the bug that happens when our defining polynomial f has any complex roots.
"""

from sage.modules.free_module_element import vector
from sage.functions.generalized import sign
from sage.functions.other import ceil
from sage.misc.functional import round
from sage.structure.factorization import Factorization
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix, column_matrix
from sage.functions.log import exp, log
from sage.matrix.special import vandermonde
from sage.rings.complex_double import CDF
from sage.rings.complex_mpfr import ComplexField
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RealField
from sage.rings.real_mpfi import RealIntervalField
from sage.misc.misc_c import prod
from sage.misc.functional import sqrt
from sage.functions.other import floor
from sage.functions.other import real, imag
from sage.rings.integer_ring import ZZ
from sage.misc.prandom import randrange
from sage.interfaces.ecm import ecm

from misc_tools import fast_persistent_save, fast_persistent_load
from collections import namedtuple

from helpers import timeprint
from helpers import major_message

from unsortedfactorization import UnsortedFactorization

import math
import copy

from collections import defaultdict
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor as ProcessPool
import multiprocessing
from time import time

from timing import timing

from sage.misc.persist import SagePickler, SageUnpickler

from os import getpid

AccumulateLog = namedtuple('AccumulateLog', ['g', 'num_or_den', 'hint'])

def log_of_big_rational(x):
    xn = x.numerator().ndigits(2)
    xd = x.denominator().ndigits(2)
    return RDF(x.norm().abs() * 2**(xd-xn)).log() + (xn - xd)*RDF(2).log()

def lognorm_cached(I, R):
        """memoized computation of R(I.norm()).log(2)"""
        if not hasattr(I, "_lognorm"):
            I._lognorm = R(I.norm()).log(2)
        return I._lognorm

def lazy_factor(I, trialdiv_limit=10000, pickle=False):
    """
    given an ideal I, lazily factor I in such a way that:
     - all prime ideals below some arbitrary "easy" bound are found by
       trial division
     - ecm is tried for some time, in such a way that factors below 50
       digits or so should all be found fairly easily.
    """

    timeprint(f"lazy_factor start (pid {getpid()})") #XXX
    if pickle:
        I = SageUnpickler.loads(I)
        timeprint(f"lazy_factor done unpickling (pid {getpid()})") #XXX

    if not I.is_integral():
        timeprint(f"lazy_factor splitting N/D (pid {getpid()})") #XXX
        return lazy_factor(I.numerator(), trialdiv_limit) / lazy_factor(I.denominator(), trialdiv_limit)

    timeprint(f"lazy_factor about to call norm (pid {getpid()})") #XXX
    N = ZZ(I.norm())

    prime_factors = []

    factors = []

    timeprint(f"Factoring ideal of {N.abs().ndigits(2)}-bit norm")

    for u in [N]:
        while True:
            p = u.trial_division(trialdiv_limit)
            if p == u:
                break
            u = u // p
            prime_factors.append(p)
        if u > 1:
            factors.append(u)

    prime_bits = sum([p.abs().ndigits(2) for p in prime_factors])
    cofactor_bits = [c.abs().ndigits(2) for c in factors]
    timeprint("After trial division, found"
          f" prime factors of total size {prime_bits} bits,"
          f" and cofactor of {cofactor_bits} bits.")

    timeprint(f"starting ECM (pid {getpid()})")
    for i in range(100):
        if not factors:
            break
        new_factors = []
        new_prime_factors = []
        for x in factors:
            for u in ecm.one_curve(x, factor_digits=15):
                if u == 1:
                    continue
                u = u.perfect_power()[0]
                if u.is_prime(proof=False):
                    new_prime_factors.append(u)
                else:
                    new_factors.append(u)
        # timeprint(f"{i}: {factors} -> {new_prime_factors}, {new_factors}")
        prime_factors += new_prime_factors
        factors = new_factors
    timeprint(f"done with ECM (pid {getpid()})")

    # Now factor by taking the gcd with all the integers we found in the
    # lazy factorization.

    F = UnsortedFactorization([])
    D = defaultdict(int)
    for p in prime_factors + factors:
        D[p] += 1
    for p in set(prime_factors):
        e = N.valuation(p)
        F *= (I + p**e).factor()
    for p in set(factors):
        e = D[p]
        assert e == 1
        F *= (I + p)
    assert F.prod() == I

    if pickle:
        return SagePickler.dumps(F)

    return F


class InsufficientPrecision(Exception):
    pass


class CadoMontgomeryReductionProcess(object):
    def __init__(self, poly, side,
                 ideals, valuations, L, log_embeddings, logfile_path=None, params=None):
        """
        with respect to the given side, and the corresponding matrices
        that give the valuations and log embeddings for that field,
        return an algebraic number that matches these as best as we can.

        L is the logmap with which log_embeddings were computed

        params is a Params object, unneeded except to be able to print aggregate timing information at the end
        """
        timeprint("start initializing self.current")
        # This step is slow because for the hash of an ideal I, sage computes the hnf of the basis for I and hashes that.
        # For n768 it takes about an hour.
        # I'm not sure if there's a better ideal invariant we could use instead of the hnf;
        # we'll be looking up self.current[I'] where I' comes from an ideal factorization, and computing the hnf is the natural way of testing if two ideals are equal.
        # Update: We now parallelize the HNF computation: the ideals are coming from a parallelized column_to_sage_ideal, we can have that call I.hash() before returning (pickling) the ideal, which caches the hash.
        self.current = { I : valuations[i] for i, I in enumerate(ideals)
                        if valuations[i] }
        timeprint("done initializing self.current")
        self.accumulated = UnsortedFactorization([])
        self.log_embeddings = log_embeddings
        self.nplus = (0, 0)
        self.nminus = (0, 0)
        self.K = poly.K[side]
        self.nt = poly.nt[side]
        self.OK = self.K.maximal_order()
        self.L = L
        self.poly = poly
        self.step_number = 0

        self.params = params

        self.accumulate_ctr = 0
        self.logfile_path = logfile_path

        # Ri and Ci _might_ be interval fields. But they don't have to.
        self.Ri = log_embeddings.parent().base_ring()
        self.Ci = self.Ri.complex_field()

        timeprint(log_embeddings.parent())
        timeprint(L.codomain())

        assert log_embeddings.parent() is L.codomain()
        precision = self.Ri.precision()

        self.C = ComplexField(precision)
        self.R = RealField(precision)
        self.precision = precision
        timeprint(f"Using precision {precision}:")
        timeprint(f"  self.C is {self.C}")
        timeprint(f"  self.R is {self.R}")
        timeprint(f"  self.Ci is {self.Ci}")
        timeprint(f"  self.Ri is {self.Ri}")

        f = self.K.defining_polynomial()
        self.V = vandermonde(f.roots(self.Ci, multiplicities=False)).transpose()

        # Computing the discriminant of OK bypasses a sage bug, see https://github.com/sagemath/sage/issues/40770
        timeprint("Compute discbound")
        self.discbound = sqrt(f.discriminant() / self.OK._K.discriminant(self.OK.basis()))

        timeprint("Initialized CadoMontgomeryReductionProcess")

    def divergence_of_embeddings(self, c):
        """
        Given a vector of log embeddings for some element g, write it as
        log(abs(norm(g)))*(1,1,...,1)+delta, with the sum of the
        coordinates of delta equal to zero, and then return delta.

        A priori, we'd like to have delta as small as we can. Or so I
        think. I'm just not too sure about whether the direction
        (1,1,...,1) can indeed be regarded as a privileged one. Sure,
        it's normal to the hyperplane sum(embeddings)==0, but what does
        that give us? Maybe there's sense in choosing another direction?

        Two reasonable contenders:

         - (log(theta_i) for i in range(d)), because an element with
           all-ones coefficients in polynomial representation will
           typically have log embeddings along that direction.

         - (log(theta_i)/skewness for i in range(d)), because if the
           definition polynomial has some skewness, we're not that much
           interested in all-ones coefficients in polynomial
           representation, but rather with skewed coefficients
        """
        K = self.K
        alpha = K.gen()
        V = self.L.codomain()
        r1, r2 = self.K.signature()
        R = V.base_ring()

        all_ones = V([1]*r1+[2]*r2)

        # one of the three options above...

        privileged = all_ones
        # privileged = self.L(alpha)
        # privileged = self.L(alpha)-R(self.poly.skewness).log()*all_ones

        # contribution on the privileged direction
        k = c.dot_product(all_ones) / privileged.dot_product(all_ones)
        return c - k * privileged

    def done(self):
        return not self.current

    def complex_to_real_matrix(self, M):
        return column_matrix(self.R, sum([[col.apply_map(real), col.apply_map(imag)] for col in M.columns()], start=[]))
        # TODO: For some reason, the below leads to worse results and causes bugs where embeddings overflow to infinity. The above seems to work fine, so maybe we'll just do LLL on a 10-dim matrix instead of 5-dim.
        if M == M.conjugate():
            # It's already real!
            return M.change_ring(self.R)
        timeprint("complex skew")
        r1, r2 = self.K.signature()
        assert M.ncols() == r1 + 2*r2, f"Matrix has {M.ncols()} columns, expected r1 + 2*r2 = {r1} + 2*{r2} = {r1+2*r2}"
        assert M[:, :r1] == M[:, :r1].conjugate(), "First r1 columns aren't all real"
        assert M[:, r1::2] == M[:, r1+1::2].conjugate(), "Complex columns aren't in complex-conjugate pairs"
        out = column_matrix(self.R,
            M[:,:r1].columns() +
            sum([[col.apply_map(real), col.apply_map(imag)] for col in M[:, r1::2].columns()], start=[])
        )
        assert out.ncols() == M.ncols(), "bug in constructing matrix: out.ncols != M.ncols"
        return out

    def skewed_LLL(self, M, skew):
        """
        returns an LLL basis where the coefficients of each vector are
        skewed according to the skewness of the number field polynomial

        TODO: what if the skew matrix is a complex matrix? I guess that
        we would have to use the conjugate_transpose, right? I'm slightly
        puzzled that this implies losing the nice structure of a Hankel
        matrix with the Newton sums when we skew by a Vandermonde matrix.

        TODO: we should be accepting an inner product matrix of the
        ambient space, not just a vector of weights. In order to make it
        work with fplll, this requires the computation of the Cholesky
        decomposition.
        """
        # We want to compute (M * skew).LLL() * skew**-1
        # But LLL can only be called on rational or integer matrices, and skew isn't rational
        # Worse yet, skew may be complex!
        # Plan:
        # Start with complex skew matrix
        # Scale up by 1/min(abs(diagonal)) if min(abs(diagonal)) < 1
        # Expand our n-by-n complex matrix into a (wide) n-by-2n real matrix D_real
        # Run LLL on M * D_real, and keep the transformation matrix U
        # Return U * M
        D = skew
        dm = min(abs(x) for x in D.diagonal())
        assert dm.parent() == self.Ri, f"dm is supposed to have parent {self.Ri} but instead has parent {dm.parent()}"
        if dm == 0: raise InsufficientPrecision()
        if dm < 1:
            D = 2 * D / dm
        try:
            Dr = self.complex_to_real_matrix(D).apply_map(round).change_ring(ZZ)
        except (ValueError, OverflowError):
            raise InsufficientPrecision()
        # sage used to have a bug where you couldn't call LLL with transformation=True on a rational matrix, but this bug has been fixed: https://github.com/sagemath/sage/pull/38841
        _, U = (M * Dr).LLL(transformation=True)
        return U * M

    def describe_picked_ideal(self, hint, num_or_den):
        """
        hint is an ideal factorization object
        """
        cumulative_norm_bits = 0
        norm_desc = defaultdict(int)
        for I, e in hint:
            nn = lognorm_cached(I, self.R)
            cumulative_norm_bits += e * nn
            norm_desc[(nn.round(),e)] += 1
        norm_desc = [(f"{nn}" + (f"*{e}" if e != 1 else ""), norm_desc[nn,e])
                     for nn,e in sorted(norm_desc.keys())]
        norm_desc = " + ".join([x if v == 1 else f"({x})*{v}"
                                for x,v in norm_desc])
        timeprint(f"Picked ideal (sign {num_or_den})"
              + f" has {cumulative_norm_bits:.2f}-bit norm ({norm_desc})")

    def pick_an_ideal_product_to_kill(self, num_or_den=1, norm_max_bits=64):
        """
        return an ideal (in factored form) that we will attempt to kill,
        with a bound on the number of bits of its norm
        """
        cumulative_norm_bits = 0
        # i = 0
        # V = vector([0]*len(ideals))
        pool = []
        for I, emax in reversed(self.current.items()):
            if cumulative_norm_bits >= norm_max_bits:
                break
            if sign(emax) != num_or_den:
                continue
            nn = lognorm_cached(I, self.R)
            cap = ceil((norm_max_bits - cumulative_norm_bits) / nn)
            e = min(emax * num_or_den, cap)
            pool.append((I, e))
            # V[i] += e
            cumulative_norm_bits += e * nn
        return UnsortedFactorization(pool)

    def pick_a_few_ideal_products_to_kill(self, count, num_or_den=1, norm_max_bits=64):
        """
        return a few ideals (in factored form) that we will attempt to kill,
        with a bound on the number of bits of each one's norm
        """
        cumulative_norm_bits = 0
        pool = []
        norm_desc = []
        out = []
        for I, emax in reversed(self.current.items()):
            if cumulative_norm_bits >= norm_max_bits:
                out.append(UnsortedFactorization(pool))
                pool = []
                cumulative_norm_bits = 0
                if len(out) >= count:
                    break
                continue
            if sign(emax) != num_or_den:
                continue
            nn = self.R(log(I.norm(), 2))
            cap = ceil((norm_max_bits - cumulative_norm_bits) / nn)
            e = min(emax * num_or_den, cap)
            pool.append((I, e))
            cumulative_norm_bits += e * nn
        if pool:
            out.append(UnsortedFactorization(pool))
        return out

    @classmethod
    def F(self, x, digits=2):
        """
        takes an element of R on interval-based R, and return a formatted
        string with this number of decimal places. So in effect it's the
        same as {x:.2f}, except that we need to handle interval fields
        specially because
        sage.rings.real_mpfi.RealIntervalFieldElement.__format__ doesn't
        support format strings (yet)
        """
        if not hasattr(x, 'center'):
            return f"{x:.{digits}f}"
        else:
            return f"{x}"


    def abort_if_infinity(self, v=None):
        """
        if the base field is an interval field, check if we have one of
        the embeddings that is an infinite interval, and bail out if it
        is the case
        """
        if v is None:
            v = self.log_embeddings
        for c in self.log_embeddings:
            if hasattr(c, 'center'):
                c = c.center()
            if c.is_infinity():
                timeprint(f"log-embeddings are {self.log_embeddings}")
                raise InsufficientPrecision()

    def adjust_log_embeddings(self, adjust):
        nl = self.log_embeddings - adjust
        self.abort_if_infinity(nl)
        self.log_embeddings = nl

    def pvec(self, v):
        F = self.F

        return "({})".format(", ".join([F(x) for x in v.change_ring(RDF)]))

    def printable_log_embeddings(self, ll):
        F = self.F

        lld = self.divergence_of_embeddings(ll)
        return f"{self.pvec(ll)}; divergence: {self.pvec(lld)}"

    @timing
    def status(self):
        """
        Some stats. Note that this incidentally updates nplus and nminus
        """

        timeprint("start status")

        self.abort_if_infinity()

        self.nplus = (0, 0)
        self.nminus = (0, 0)

        timeprint(f"begin loop through self.current ({len(self.current)} items)...")
        for I, v in self.current.items():
            if v > 0:
                self.nplus = (self.nplus[0]+1, self.nplus[1] + v * lognorm_cached(I, self.R))
            else:
                self.nminus = (self.nminus[0]+1, self.nminus[1] - v * lognorm_cached(I, self.R))
        timeprint(f"end loop")

        F = self.F

        major_message(f"num bits {F(self.nplus[1])}, {F(self.nminus[1])}")

        timeprint(f"num bits {F(self.nplus[1])} ({round(self.nplus[0])} ideals)",
              f"denom bits {F(self.nminus[1])} ({round(self.nminus[0])} ideals)")
        r1, r2 = self.K.signature()

        #avg_lognorm = sum([log_of_big_rational(I.norm())*e for I,e in
        #                  self.current.items()])/self.K.degree()
        # XXX Above line is slow (1-2 seconds on n666, out of 5-6 seconds per reduction step),
        # so I commented it out and replaced it with the below one.
        # Why was log_of_big_rational used here but nowhere else?
        # Is there maybe a loss of precision for big rationals with RDF?
        avg_lognorm = sum(lognorm_cached(I, RDF)*e for I,e in self.current.items())/self.K.degree()
        # Regardless, avg_lognorm is only used for the following print statement, so it doesn't matter too much
        timeprint("expected average of log-embeddings (typical):",
              self.pvec(avg_lognorm * vector([1]*r1+[2]*r2)))

        timeprint("Our log embeddings:")
        timeprint("  ", self.printable_log_embeddings(self.log_embeddings))

        if self.nminus[1] < 10:
            denom = prod([I.norm()**-v for I,v in self.current.items() if v < 0])
            timeprint("denominator:", denom)
            timeprint("numerator bounds: ", denom * vector(self.bound_on_coefficients_of_remaining_part()))

    @timing
    def accumulate(self, g, num_or_den, hint=UnsortedFactorization([]), glist=None, hintlist=None, cpucount=-1, dolog=True):
        """
        Take action, and register a new term for the product, updating
        the current state accordingly
        """

        self.accumulate_ctr += 1
        if self.logfile_path is not None and dolog:
            thisfile = f"{self.logfile_path}.{self.accumulate_ctr}.sobj"
            thislog = AccumulateLog(g=g, num_or_den=num_or_den, hint=hint)
            fast_persistent_save(thislog, thisfile)

        self.accumulated *= UnsortedFactorization([(g, num_or_den)])
        F = self.F

        gen_norm = self.R(g.norm().abs()).log(2)

        timeprint(f"Picked gen has {gen_norm:.2f}-bit norm")

        adjust = num_or_den * self.L(g)
        timeprint("Adjusting log embeddings by", self.printable_log_embeddings(adjust))

        self.adjust_log_embeddings(adjust)

        if glist is None:
            timeprint("Dividing fractional_ideal(g) / hint.prod()")
            discovered = (self.OK.fractional_ideal(g) / hint.prod())
            timeprint("Finished dividing fractional_ideal(g) / hint.prod()")
            discovered = lazy_factor(discovered)
            timeprint("Finished lazy_factor")
        else:
            discovered = UnsortedFactorization([])
            ctx = multiprocessing.get_context('fork')
            with ProcessPool(mp_context=ctx, max_workers=cpucount) as executor:
                timeprint(f"Start sending out lazy_factor jobs (max_workers={cpucount})")
                futures = [
                    executor.submit(lazy_factor,
                                    SagePickler.dumps(self.OK.fractional_ideal(glist[i])/hintlist[i].prod()),
                                    10000, True)
                    for i in range(len(glist))
                ]
                timeprint("Done sending out lazy_factor jobs, now waiting for results")
                for future in concurrent.futures.as_completed(futures):
                    discovered = discovered * SageUnpickler.loads(future.result())
                timeprint("All lazy_factor results arrived and combined")

        timeprint("ideal_factorization multiply start")
        ideal_factorization = hint * discovered
        assert isinstance(ideal_factorization, UnsortedFactorization)
        timeprint("ideal_factorization multiply end")

        timeprint("update self.current start")
        for I, e in ideal_factorization:
            # Ibits = round(self.R(I.norm()).log(2))
            if I not in self.current:
                # timeprint(f"Just added a {Ibits:.2f}-bit ideal")
                self.current[I] = 0
            self.current[I] -= num_or_den * e
            if self.current[I] == 0:
                # timeprint(f"Just got rid of a {Ibits:.2f}-bit ideal")
                del self.current[I]
        timeprint("update self.current end")

    def l2norm_on_divergence_hyperplane(self, lld):
        """
        ll is such that the sum of the first r1 coordinates + twice the
        sum of the r2 coordinates is zero. We'll compute the norm of the
        d-dimensional vector that is obtained by duplicating the r2
        complex places
        """
        r1, r2 = self.K.signature()
        sum1 = sum([lld[i]**2 for i in range(r1)])
        sum2 = sum([lld[i]**2 for i in range(r1, r1+r2)])
        return sum1 + sum2/2

    def analyze_generators(self, gens, num_or_den):
        """
        given a list of generators of the target ideal, print the log
        embeddings of each of these, together with the divergence vs the
        multiples of all-ones vectors. Then determine which of these,
        when added or subtracted from the target log embeddings, brings
        us closest to the all-ones line.
        """

        timeprint("Our log embeddings:")
        timeprint("  ", self.printable_log_embeddings(self.log_embeddings))
        ll0d = self.divergence_of_embeddings(self.log_embeddings)

        timeprint("log embeddings of generators:")
        scoreboard = []
        for i,g in enumerate(gens):
            ll = self.L(g)
            lld = self.divergence_of_embeddings(ll)
            reached = self.l2norm_on_divergence_hyperplane(lld-num_or_den*ll0d)

            timeprint("  ", self.printable_log_embeddings(self.L(g)),
                  f"; reaches {RDF(reached):.2f}")
            scoreboard.append((reached, (i, g)))

        best_reached, (i0, g0) = min(scoreboard)
        timeprint(f"Based on divergence minimization, we would choose generator {i0}")

        return i0

    def bound_on_coefficients_of_remaining_part(self):
        """
        This returns absolute bounds on the coefficients of the algebraic
        number that has the log-embeddings currently registered in self.
        If it happens that we also know that this number is an algebraic
        integer (denominator==1), then these can be used to bound
        _integer_ coefficients (well, up to the denominators in the
        maximal order, but this is the idea).

        FIXME (hmmm, it seems very very weird that this works. V would
        need to be positive definite or something like this)
        """
        Vi = self.V**-1
        d = self.K.degree()
        Amax = self.nt.modules_of_embeddings_from_log_embeddings(self.log_embeddings)
        return [floor(sum([abs(Amax[i] * Vi[i,j]) for i in range(d)])) for j in range(d)]

    def one_reduction_step(self,
                           bits=64,
                           skip_asserts=True,
                           single_ideal_at_random=False,
                           all_ideals_at_once=False,
                           skip_embeddings_reduction=False,
                           escape_loop=False
                           ):
        if not self.current:
            # if we have no outstanding ideals, do nothing
            timeprint("cannot reduce further, no ideals left")
            return

        d = self.K.degree()
        s = self.poly.skewness
        f = self.K.defining_polynomial()
        x = f.parent().gen()

        timeprint(f"------------ reduction step {self.step_number} ------------")
        self.step_number += 1

        try:
            if single_ideal_at_random:
                I = list(self.current.keys())[randrange(len(self.current))]
                v = self.current[I]
                num_or_den = 1 if v > 0 else -1
                hint = UnsortedFactorization([(I, v * num_or_den)])
            elif type(all_ideals_at_once) is int:
                # we're cheating on num_or_den, here
                num_or_den = all_ideals_at_once
                hint = UnsortedFactorization([(I, num_or_den*v)
                                      for I,v in self.current.items()
                                      if sign(v) == num_or_den])
            elif all_ideals_at_once:
                # we're cheating on num_or_den, here
                num_or_den = 1 if self.nplus[0] > self.nminus[0] else -1
                hint = UnsortedFactorization([(I, num_or_den*v)
                                      for I,v in self.current.items()])
            else:
                num_or_den = 1 if self.nplus[1] > self.nminus[1] else -1
                hint = self.pick_an_ideal_product_to_kill(num_or_den, bits)
        except OverflowError:
            raise InsufficientPrecision()

        self.describe_picked_ideal(hint, num_or_den)

        timeprint("start computing I = hint.prod()")
        I = hint.prod()
        timeprint("done computing I = hint.prod()")

        # in fact we don't really care if it's integral or not...
        # assert I.is_integral()

        # by doing LLL here, we tolerate the norm of each generator to
        # grow by as much as a constant C_K, which we can compute.
        timeprint("start computing basis for I")
        gens0 = I.basis()
        timeprint("done computing basis for I")
        M0 = matrix([list(c) for c in gens0])

        # There does seem to be a computational advantage in doing the
        # reduction in two steps.
        L2s_norm_of_f = self.R(vector(list(f(s*x) / s**(d / 2))).norm())

        skew0 = diagonal_matrix([(self.Ri(s) ** (i - (d - 1)/2)) for i in range(d)])
        timeprint("starting skewed LLL")
        M1 = self.skewed_LLL(M0, skew0)
        timeprint("done with skewed LLL")

        gens1 = [self.K(list(g)) for g in M1]

        if not skip_asserts:
            # The theory is that the skewed norm of the vector M1[0] is
            # bounded as follows
            bound_on_L2s_norm_of_v = abs(self.R(
                    2**((d - 1) / 4) * abs(M0.determinant())**(1 / d)
                    ))
            L2s_norm_of_v = abs(self.R((M1[0]*skew0).norm()))

            L2s_approximation_ratio = L2s_norm_of_v / bound_on_L2s_norm_of_v

            timeprint("L2s approximation ratio for 1st reduction (expected <= 1):",
                  L2s_approximation_ratio)
            assert L2s_approximation_ratio <= 1

            # we can just ignore the discriminant quotient. It just makes the
            # bound sharper if we happen to know it, that's it.
            alg_norm_quotient = abs(gens1[0].norm() / I.norm())
            bound_on_alg_norm_quotient = self.R(prod([
                2**(d * (d - 1) / 4),
                L2s_norm_of_f**(d - 1),
                self.discbound]))

            # This norm is obtained by Hadamard + Cauchy-Schwarz, and is
            # expected to be very loose.
            alg_norm_approximation_ratio = alg_norm_quotient / bound_on_alg_norm_quotient
            timeprint("alg norm approximation ratio for 1st reduction (expected <= 1):",
                  alg_norm_approximation_ratio)
            assert alg_norm_approximation_ratio <= 1

        new_gens1 = []
        for g in gens1:
            ll = self.L(g)
            if any([x.absolute_diameter().is_infinity() for x in ll]):
                continue
            else:
                new_gens1.append(g)
                #raise InsufficientPrecision()

        if len(new_gens1) == 0:
            timeprint("Oh no. No generators passed the is_infinity check.")
            raise InsufficientPrecision()

        timeprint("Analyzing generators in M1")
        i0 = self.analyze_generators(new_gens1, num_or_den)

        if not skip_embeddings_reduction:
            # We now have alternative, somewhat smaller generators of our
            # ideal. They're finely skewed, which contributes to making the
            # norm a bit smaller than if we had neglected that aspect.

            # Next we want the log embeddings to be small, which we translate
            # into the requirement that the vector that we end up with is as
            # close as we can to the same skewed hypersphere as the one the
            # current embeddings of the target lie on.

            modules = self.nt.modules_of_embeddings_from_log_embeddings(self.log_embeddings)
            V = self.V

            # note that V * V.transpose() is the Hankel matrix whose
            # coefficients are the Newton sums, which we can compute fairly
            # easily via the expansion of (xf'/f) in powers of 1/x

            # alas, we're rather interested by the comparison with the
            # hypersphere, and this kills the nice properties of these
            # expression. The polynomial that gives the target embeddings is
            # unknown anyways. We have to consider a skewed V:

            # if we take:
            # W = V * diagonal_matrix(self.nt.real_representation_of_embeddings(g))**-1
            # then W is so that vector(g.list()) * W is the all-ones vector

            if diagonal_matrix(modules).is_singular():
                raise InsufficientPrecision()

            # but of course, the W that we want is the one that compares to
            # the target number!
            W = V * diagonal_matrix(modules)**-num_or_den
            assert W.base_ring().prec() >= self.precision, f"W.base_ring is {W.base_ring()}, lower precision than {self.precision}"


            # So. At this point, what we have to do is to give a skewed LLL
            # reduction with respect to this matrix W
            M2 = self.skewed_LLL(M1, skew=W)
            gens2 = [self.K(list(g)) for g in M2]

            for g in gens2:
                ll = self.L(g)
                if any([x.absolute_diameter().is_infinity() for x in ll]):
                    raise InsufficientPrecision()

            timeprint("Analyzing generators in M2")
            i0 = self.analyze_generators(gens2, num_or_den)

            timeprint(f"Choosing generator {i0}")
            g = gens2[i0]

            # I don't know how to make sense of this debug print, let's
            # drop it.
            # drift_bits = (vector(g.list()) * W).norm().log(2)
            # F = self.F
            # timeprint(f"Picked gen has {F(drift_bits)}-bit drift")
        else:
            timeprint(f"Choosing generator {i0}")
            g = new_gens1[i0]
            gens2 = []

        if escape_loop:
            # Just pick a random generator
            all_gens = gens1 + gens2
            jj = randrange(0, len(all_gens))
            g = all_gens[jj]

        self.accumulate(g, num_or_den, hint=hint)


    def one_internally_parallel_reduction_step(self, bits=64, cpucount=-1):

        if not self.current:
            # if we have no outstanding ideals, do nothing
            timeprint("cannot reduce further, no ideals left")
            return

        timeprint(f"------------ reduction step {self.step_number} ------------")
        self.step_number += 1

        try:
            num_or_den = 1 if self.nplus[1] > self.nminus[1] else -1
            timeprint("start picking a few ideal products to kill")
            hintlist = self.pick_a_few_ideal_products_to_kill(cpucount, num_or_den, bits)
            timeprint("done picking a few ideal products to kill")
        except OverflowError:
            raise InsufficientPrecision()

        candidate_gens = dict()

        ctx = multiprocessing.get_context('fork')
        timeprint(f"setting up ProcessPool with maxworkers={cpucount}")
        with ProcessPool(mp_context=ctx, max_workers=cpucount) as executor:
            timeprint("start submitting outer_handling_of_hint jobs")
            futures = [
                executor.submit(outer_handling_of_hint, i, SagePickler.dumps(hintlist[i]), num_or_den,
                                SagePickler.dumps(self.K), self.precision, self.poly.skewness)
                for i in range(len(hintlist))
            ]
            timeprint("done submitting outer_handling_of_hint jobs; start waiting for results")
            for future in concurrent.futures.as_completed(futures):
                #idx, gens = future.result()
                idx, gens, I = future.result() # PASSIBACK
                candidate_gens[idx] = SageUnpickler.loads(gens)
                timeprint("start unpickle I")
                I = SageUnpickler.loads(I) # PASSIBACK
                timeprint("end unpickle I")
                hintlist[idx].cache_product(I)
            timeprint("done wating for outer_handling_of_hint results")
        timeprint("done cleaning up ProcessPool")
        for i in range(len(hintlist)):
            for g in candidate_gens[i]:
                ll = self.L(g)
                if any([x.absolute_diameter().is_infinity() for x in ll]):
                    raise InsufficientPrecision()

        timeprint("Analyzing generators over all hints...")

        # Pick one generator from each hint,
        # and combine into a single generator for the whole reduction step,
        # which is passed to accumulate.

        timeprint("start multiplying UnsortedFactorization objects")
        total_hint = prod(hintlist)     # product of UnsortedFactorization objects
        timeprint("done multiplying UnsortedFactorization objects")
        glist = []

        # TODO: actually do something here.
        # Just testing for now.

        total_g = self.K(1)

        for i in range(len(hintlist)):
            #cands = [total_g * x for x in candidate_gens[i]]
            #chosen_idx = self.analyze_generators(cands, num_or_den)
            chosen_idx = 0
            glist.append(candidate_gens[i][chosen_idx])
            total_g = total_g * candidate_gens[i][chosen_idx]
            #total_g = cands[chosen_idx]
            #timeprint(f"Chose generator {chosen_idx} for hint number {i}")

        self.accumulate(total_g, num_or_den, total_hint, glist, hintlist, cpucount)


    def do_many_reduction_steps_in_parallel(self, bits=64, cpucount=2):
        d = self.K.degree()
        s = self.poly.skewness
        f = self.K.defining_polynomial()
        x = f.parent().gen()

        if not self.current:
            # if we have no outstanding ideals, do nothing
            timeprint("cannot reduce further, no ideals left")
            return

        num_or_den = 1 if self.nplus[1] > self.nminus[1] else -1

        timeprint(f"-------- reduction steps {self.step_number} to {self.step_number+cpucount-1} --------")
        self.step_number += cpucount

        ts = time()
        #cpucount = multiprocessing.cpu_count()
        jobs = self.pick_a_few_ideal_products_to_kill(cpucount, num_or_den, bits)
        #timeprint(type(jobs[0]))
        modules = self.nt.modules_of_embeddings_from_log_embeddings(self.log_embeddings)
        ctx = multiprocessing.get_context('fork')

        # Note that find_short_generator does read-only operations on the object
        #cmrp_copy = copy.copy(self)
        #import sage.libs.pari.convert_sage as convert_sage
        with ProcessPool(mp_context=ctx) as executor:
            #futures = [executor.submit(find_short_generator,
            #    hint.prod(), num_or_den, self.K, self.precision, modules, self.V, self.poly.skewness, idx, self) for (idx,hint) in enumerate(jobs)]
            futures = [executor.submit(self.bound_on_coefficients_of_remaining_part, self) for i in range(len(jobs))]


            timeprint("at least submitted futures")
            for future in concurrent.futures.as_completed(futures):
                g, idx = future.result()
                self.accumulate(g, num_or_den, hint=jobs[idx])
                self.status()
        timeprint(" took",time()-ts)


def outer_handling_of_hint(idx, hint, num_or_den, K, precision, skewness):
    timeprint("outer_handling_of_hint")

    hint = SageUnpickler.loads(hint)
    K = SageUnpickler.loads(K)

    d = K.degree()
    Ri = RealIntervalField(precision)

    outer_describe_picked_ideal(hint, num_or_den, precision)
    I = hint.prod()

    gens0 = I.basis()
    M0 = matrix([list(c) for c in gens0])

    skew0 = diagonal_matrix([(Ri(skewness) ** (i - (d - 1)/2)) for i in range(d)])
    M1 = outer_skewed_LLL(M0, skew0, precision)

    gens1 = [K(list(g)) for g in M1]

    #return idx, SagePickler.dumps(gens1)
    return idx, SagePickler.dumps(gens1), SagePickler.dumps(I)
    # Passing I back in addition to gens1 allows cacheing hint.prod() in the main process.
    # A different tradeoff we could make would be to have the main thread do I = OK.ideal(gens1)
    # Less pickling, but a little more compute.
    # For n448 passing I back is slightly faster, but this might be worth revisiting if pickling gets slow.


def outer_describe_picked_ideal(hint, num_or_den, precision):
    """
    hint is an ideal factorization object
    """
    R = RealField(precision)
    cumulative_norm_bits = 0
    norm_desc = defaultdict(int)
    for I, e in hint:
        nn = lognorm_cached(I, R)
        cumulative_norm_bits += e * nn
        norm_desc[(nn.round(),e)] += 1
    norm_desc = [(f"{nn}" + (f"*{e}" if e != 1 else ""), norm_desc[nn,e])
                 for nn,e in sorted(norm_desc.keys())]
    norm_desc = " + ".join([x if v == 1 else f"({x})*{v}"
                            for x,v in norm_desc])
    timeprint(f"Picked ideal (sign {num_or_den})"
          + f" has {cumulative_norm_bits:.2f}-bit norm ({norm_desc})")


def outer_skewed_LLL(M, skew, precision):
    D = skew
    Ri = RealIntervalField(precision)
    dm = min(abs(x) for x in D.diagonal())
    assert dm.parent() == Ri, f"dm is supposed to have parent {Ri} but instead has parent {dm.parent()}"
    if dm == 0: raise InsufficientPrecision()
    if dm < 1:
        D = 2 * D / dm
    try:
        Dr = outer_complex_to_real_matrix(D, precision).apply_map(round).change_ring(ZZ)
    except (ValueError, OverflowError):
        raise InsufficientPrecision()
    # sage used to have a bug where you couldn't call LLL with transformation=True on a rational matrix, but this bug has been fixed: https://github.com/sagemath/sage/pull/38841
    _, U = (M * Dr).LLL(transformation=True)
    return U * M


def outer_complex_to_real_matrix(M, precision):
    R = RealField(precision)
    return column_matrix(R, sum([[col.apply_map(real), col.apply_map(imag)] for col in M.columns()], start=[]))


def find_short_generator(I, num_or_den, K, precision, modules, V, skewness, idx, CMRP):
    ## Non-obvious parameters:
    # modules = self.nt.modules_of_embeddings_from_log_embeddings(self.log_embeddings)
    # V = self.V
    # skewness = self.poly.skewness
    # CMRP = CadoMontgomeryReductionProcess (should do read-only operations)
    timeprint("made it to find_short_generator")

    d = K.degree()
    s = skewness
    f = K.defining_polynomial()
    x = f.parent().gen()

    Ri = RealIntervalField(precision)
    R = RealField(precision)

    # match the one_reduction_step above
    gens0 = I.basis()
    M0 = matrix([list(c) for c in gens0])

    L2s_norm_of_f = self.R(vector(list(f(s*x) / s**(d / 2))).norm())

    skew0 = diagonal_matrix([(Ri(s) ** (i - (d - 1)/2)) for i in range(d)])
    M1 = CMRP.skewed_LLL(M0, skew0, R)

    gens1 = [K(list(g)) for g in M1]

    for g in gens1:
        ll = CMRP.L(g)
        if any([x.absolute_diameter().is_infinity() for x in ll]):
            raise InsufficientPrecision()

    # timeprint("Analyzing generators in M1")
    i0 = CMRP.analyze_generators(gens1, num_or_den)

    if True:
        # We now have alternative, somewhat smaller generators of our
        # ideal. They're finely skewed, which contributes to making the
        # norm a bit smaller than if we had neglected that aspect.

        # Next we want the log embeddings to be small, which we translate
        # into the requirement that the vector that we end up with is as
        # close as we can to the same skewed hypersphere as the one the
        # current embeddings of the target lie on.

        #modules = CMRP.nt.modules_of_embeddings_from_log_embeddings(CMRP.log_embeddings)
        V = CMRP.V

        # note that V * V.transpose() is the Hankel matrix whose
        # coefficients are the Newton sums, which we can compute fairly
        # easily via the expansion of (xf'/f) in powers of 1/x

        # alas, we're rather interested by the comparison with the
        # hypersphere, and this kills the nice properties of these
        # expression. The polynomial that gives the target embeddings is
        # unknown anyways. We have to consider a skewed V:

        # if we take:
        # W = V * diagonal_matrix(self.nt.real_representation_of_embeddings(g))**-1
        # then W is so that vector(g.list()) * W is the all-ones vector

        if diagonal_matrix(modules).is_singular():
            raise InsufficientPrecision()

        # but of course, the W that we want is the one that compares to
        # the target number!
        W = V * diagonal_matrix(modules)**-num_or_den
        assert W.base_ring().prec() >= CMRP.precision, f"W.base_ring is {W.base_ring()}, lower precision than {CMRP.precision}"


        # So. At this point, what we have to do is to give a skewed LLL
        # reduction with respect to this matrix W
        M2 = CMRP.skewed_LLL(M1, W, R)
        gens2 = [K(list(g)) for g in M2]

        for g in gens2:
            ll = CMRP.L(g)
            if any([x.absolute_diameter().is_infinity() for x in ll]):
                raise InsufficientPrecision()

        #timeprint("Analyzing generators in M2")
        i0 = CMRP.analyze_generators(gens2, num_or_den)

        #timeprint(f"Choosing generator {i0}")
        g = gens2[i0]
        return g, i0

        # I don't know how to make sense of this debug print, let's
        # drop it.
        # drift_bits = (vector(g.list()) * W).norm().log(2)
        # F = self.F
        # timeprint(f"Picked gen has {F(drift_bits)}-bit drift")
    else:
        #timeprint(f"Choosing generator {i0}")
        g = gens1[i0]
        return g, i0
