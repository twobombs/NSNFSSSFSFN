# p-adic computation of the e-th root.
#
# We have a product expression for some number field element.
# The product expression is a list of pairs ((a,b),k) and represents
# prod([(a-b*alpha)^k]) in the number field. Some of the k's may be
# negative.
#
# We know that this product is an e-th power in the number field and we
# want to compute its e-th root. All k's in the input are bounded by e.
#
# The input size is bounded by O(A*M*log(e)) where A is the number of
# pairs, M is the max number of bits of the coefficients a,b
#
# The e-th root is of size at most O(A*M).
#
# We want to compute the e-th root by p-adic lifting modulo an inert
# prime. We'll do higher dimensional rational reconstruction.
#
#

from sage.functions.other import factorial
from sage.rings.infinity import infinity
from sage.arith.misc import primes
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.matrix.special import block_matrix
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.arith.misc import inverse_mod, power_mod
from timing import *
from candy import major_message


def find_inert_prime(f, TTplus=[], TTminus=[], avoid=set()):
    for i, p in enumerate(primes(2, infinity)):
        if i > 2 * factorial(f.degree()):
            raise ValueError("Cannot find an inert prime."
                             " Maybe this polynomial has a peculiar Galois group")
        fp = f.change_ring(GF(p))
        if fp.degree() < f.degree():
            continue
        if p in avoid:
            print(f"not using p={p} because"
                  " we listed it among the primes to avoid in the lift")
            continue
        if not fp.is_irreducible():
            continue
        Fpn = GF(p**(f.degree()), 'alpha_p', modulus=f)
        if prod([(a-b*Fpn.gen())**k for (a,b),k in TTplus]) == 0:
            continue
        if prod([(a-b*Fpn.gen())**k for (a,b),k in TTminus]) == 0:
            continue
        return p


def reconstruct_as_short_polynomial_fraction(z):
    """
    Given z an element of a modular quotient ring (with a certain
    characteristic, which is a prime or a prime power, and a certain
    degree), return a univariate integer rational fration that is
    equivalent to z, but with short coefficients
    """

    Zpkn = z.parent()
    n = Zpkn.degree()
    pk = Zpkn.characteristic()
    alpha_pk = Zpkn.gen()
    L = block_matrix(2, 2,
                     [matrix(ZZ,[(z * alpha_pk**i).list()
                                 for i in range(n)]),
                      1, pk, 0
                      ]).LLL()
    # the smaller the better. At some point L[0].norm() should no longer
    # grow.
    # score = float(L[0].change_ring(RDF).norm() / sqrt(pk))
    coeffs = L[0].list()
    ZP = ZZ['x']
    rn = ZP(coeffs[:n])
    rd = ZP(coeffs[n:])
    assert rn(alpha_pk) - z * rd(alpha_pk) == 0
    # We can count the number of digits without having to battle
    # precision issues. This means that we're off in our computation by
    # an ever so slightly bit, but this is totally negligible anyways.
    score = max([1+c.ndigits(2) for c in rn.list() + rd.list()])
    return rn/rd, score


def padic_root_of_algebraic_product(TTplus, TTminus, f, p, e,
                                    pre_multiply_root=None, pre_multiply_gamma_fac=None):
    """
    This is the main part of the root computation

    The main use case is root_fragment==None

    given TTplus and TTminus as lists of pairs ((a,b),k), compute the
    e-th root of the product y/z in the number field
    defined by f, where y is the product of the (a-b*alpha)^k in TTplus,
    and z the product of the (a-b*alpha)^k in TTminus. The calculation is
    done modulo increasing powers of the inert prime p (at step i we work
    modulo p^(2^i)), and we stop when we have found a rational fraction
    that is apparently valid over Z.

    Note that this algorithm will only terminate if the solution indeed
    exists over Z, which is not something we check!

    If pre_multiply_root is given, it's the inverse a part of the root
    that we already know, and we incorporate it in the computation (given
    as a univariate rational fraction over the integers). That is, we
    compute an e-th root of y/z*pre_multiply_root**e
    """
    ZP = ZZ['x']

    d = f.degree()
    # Ri is (Z/p^(2^i))[x]/f(x)
    R0 = GF(p**(f.degree()), 'alpha_p', modulus=f)
    alpha0 = R0.gen()
    y0 = prod([(a-b*alpha0)**k for (a,b),k in TTplus])
    z0 = prod([(a-b*alpha0)**k for (a,b),k in TTminus])
    u0 = 1/z0
    # Our iteration will compute the *inverse* e-th root of y/z
    rf0 = 1
    if pre_multiply_root is not None:
        rf0 = pre_multiply_root.change_ring(R0)(alpha0)

    if pre_multiply_gamma_fac is not None:
        for fac in pre_multiply_gamma_fac:
            nf_elt = fac[0]
            exp = fac[1]
            poly_alpha0 = nf_elt.polynomial().change_ring(R0)(alpha0)
            rf0 *= (1/(poly_alpha0**exp))       # whole product is 1/gamma

    r0 = (z0 / y0 / rf0**e).nth_root(e)

    y,z,u,r,R,pk,i,rf = y0,z0,u0,r0,R0,p,0,rf0

    magic = 0
    while True:
        # invariants:
        assert pk == p**(2**i)
        # R = (Z/p^(2^i))[x]/f(x)
        assert y in R
        assert z in R
        assert u in R
        assert r in R
        assert u * z == 1
        assert r**e * y * rf**e * u == 1

        i += 1
        pk = pk * pk
        R = Integers(pk).extension(f)

        r = ZP(r.list())(R.gen())
        # lift the preinverse of z, too.
        u = ZP(u.list())(R.gen())

        # Newton step
        # We need to recompute the thing we're computing an e-th root
        # for, at larger precision. It smells like a waste of time, but
        # since we're doing it at precisions that grow geometrically, in
        # the end it's still fine: it's linear in the input and output
        # sizes.  The annoying thing is that inverting z is tricky, we
        # need to keep track of various things to make it work.
        y = prod([(a-b*R.gen())**k for (a,b),k in TTplus])
        z = prod([(a-b*R.gen())**k for (a,b),k in TTminus])

        # if u*z = 1+pb, then the higher order inverse of z is u*(1-pb)
        # IOW, a Newton step on u: u becomes u*(1-(u*z-1))
        u = u * (1 - (u * z - 1))

        rf = 1
        if pre_multiply_root is not None:
            rf = pre_multiply_root.change_ring(R)(R.gen())

        if pre_multiply_gamma_fac is not None:
            for fac in pre_multiply_gamma_fac:
                nf_elt = fac[0]
                exp = fac[1]
                poly_R_gen = nf_elt.polynomial().change_ring(R)(R.gen())
                rf *= (1/(poly_R_gen**exp))     # whole product is 1/gamma

        # This is the Newton iteration on the function f(x) = y/z - x^-e
        r += r * (1 - r**e * rf**e * y * u) / e

        fractional, score = reconstruct_as_short_polynomial_fraction(r)

        magic = 2 * d * max(pk.ndigits(2)//2-score, 0)
        timeprint(f"Lifting modulo {p}^(2^{i}): {magic} magic zeros")

        if magic > 100:
            # remember that we have computed the _inverse_ of the e-th
            # root.
            return 1/fractional

@timing
def padic_eth_root(params,
                   linalg_output,
                   lift_centered=True,
                   pre_multiply_root=None,
                   pre_multiply_gamma_fac=None):
    e = params.parameters['e']
    N = params.poly.N
    m = params.poly.m
    f = params.poly.f[1]
    ZN = Integers(N)

    # centered lifts are better because we reach the solution way
    # earlier.
    if lift_centered:
        sol = linalg_output.sol.lift_centered()
    else:
        sol = linalg_output.sol.lift()

    TT = [(linalg_output.row_to_aquery[i],s) for i,s in enumerate(sol)]
    TT += list(linalg_output.ST_list.items())

    TTplus  = [((a,b),k)  for (a,b),k in TT if k>0]
    TTminus = [((a,b),-k) for (a,b),k in TT if k<0]

    SU_m = prod([ZN(a-b*m)**k for (a,b),k in TT])

    p = find_inert_prime(f, TTplus, TTminus, avoid=set([e]))
    print(f"Beginning e-th root computation mod powers of {p}")
    print(f"  TTplus has {len(TTplus)} pairs")
    print(f"  TTminus has {len(TTminus)} pairs")

    if pre_multiply_gamma_fac is not None:
        major_message("Using the saved gamma-fac object...")
        r = padic_root_of_algebraic_product(TTplus, TTminus,
                                            f, p, e,
                                            pre_multiply_gamma_fac=pre_multiply_gamma_fac)
        for fac in pre_multiply_gamma_fac:
            nf_elt = fac[0]
            exp = fac[1]
            poly_m = nf_elt.polynomial().change_ring(ZN)(m)
            #poly_m_inv = inverse_mod(poly_m, N)
            SU_m *= ZN(poly_m)**(-e*exp)

        assert r(ZN(m))**e == SU_m
        return r

    # usual case
    r = padic_root_of_algebraic_product(TTplus, TTminus,
                                        f, p, e,
                                        pre_multiply_root=pre_multiply_root)

    if pre_multiply_root is not None:
        SU_m *= pre_multiply_root.change_ring(ZN)(m)**e

    assert r(ZN(m))**e == SU_m

    return r
