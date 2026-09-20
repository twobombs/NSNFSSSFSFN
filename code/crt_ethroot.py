from sage.all import *
import itertools


"""
eth roots in a number field using the CRT approach.
(Note that this is not what we did in the paper, for that see montgomery_*.py and hybrid_root_*.py.)
"""

"""
Let alpha be the root of f in our number field.
Let P be the root of f mod N. (alpha maps to P).

We want to compute an eth root of S(alpha) = prod_i (a_i + b_i alpha)^(m_i).
Eventually we want an eth root of S(P) mod N.

Assume all the m_i are positive, so that S(alpha) is in Z[alpha], and assume the eth root is also in Z[alpha].

The idea is this: Take everything mod p, where p is inert in our the number field.
We start with Z[alpha], which is an extension of Z, namely Z[x]/(f(x))).
So we map to an extension of F_p, namely F_p[x]/(f(x)).
Since f(x) is irreducible mod p and has degree d, we get F_{p^d}.
Let beta be the root of f(x) in F_{p^d}.

Where does S(alpha) map to? It maps to S(beta). We can compute it explicitly as prod_i (a_i + b_i beta)^(m_i).

Then we take the eth root in F_{p^d}. (If e is relatively prime to p^d-1, it should be unique.)
i.e., find roots in F_{p^d} of (x^e - S(beta))

We now know the coefficients mod p of the eth root of S(alpha).

From here, we have two options:
1) Use Hensel lifting to compute the coefficients mod p^2, p^3, p^4, and so on; or
2) Repeat with other inert primes p, and CRT the results together to compute the coefficients mod p1*p2*p3*...

Here we do the latter.

Either way, eventually we know the coefficients of the eth root of S(alpha) over the integers.
"""

x = polygen(ZZ)

def polyselect(N, d):
    """
    Find a polynomial f of degree d and an integer P such that f(P) = 0 (mod N).
    (In particular, just do base-m.)
    Return (f, P).
    """
    P = Integer(int(floor(Integer(N)**(Integer(1)/d))))
    coeffs = Integer(N).digits(P)
    f = sum(ci * x**i for (i,ci) in enumerate(coeffs))
    assert f(P) % N == 0 and f.degree() == d
    assert f.is_irreducible() # if not, we've already factored N.
    return f, P

def find_good_primes(f, e, num=1000, start=2**41, stop=2**42):
    """Generate primes p satisfying the following:
        1. gcd(p^(deg f) - 1, e) = 1 (so that eth roots are unique in F_p^d)
        2. f is irreducible mod p (i.e., p is inert in the number field defined by f)
    """
    prime_list = []
    for p in primes(start=start, stop=stop):
        if gcd(p-1, e) != 1: continue
        if gcd(p.powermod(f.degree(), e) - 1, e) != 1: print("skip", gcd(pow(p, f.degree(), e) - 1, e)); continue
        if f.change_ring(GF(p)).is_irreducible():
            prime_list.append(p)
            if len(prime_list) > num:
                return prime_list
    return prime_list

def ethroot(abms, f, e, N):
    """
    Find an eth root of S(alpha) in Z[alpha], or at least find its coefficients mod p for a bunch of primes p.
    :param abms: list of (a_i, b_i, m_i), where S(x) = prod_i (a_i + b_i x)^m_i
    :param f: polynomial defining the number field. (i.e., f(alpha) = 0)
    :param e: root we're finding
    :param N: modulus
    Yields (p, [coefficients mod p]) for various primes p
    """
    for p in find_good_primes(f, e):
        assert gcd(p-1, e) == 1
        Fpd = GF(p**f.degree(), modulus=f, name='beta') # note: f must be irreducible mod p
        beta = Fpd.gen()
        assert f(beta) == 0
        Sbeta = prod( (ai + bi * beta) ** mi  for (ai,bi,mi) in abms)
        y = polygen(Fpd, 'y')
        root = (y**e - Sbeta).roots()[0][0]
        # The root, as a polynomial in beta, has coefficients that agree mod p with those of the eth root of S(alpha).
        yield (p, [c.lift() for c in root.list()])

def ethroot_p(abms, f, e, N, p):
    """
    Unyield above function
    """
    assert gcd(p-1, e) == 1
    Fpd = GF(p**f.degree(), modulus=f, name='beta') # note: f must be irreducible mod p
    beta = Fpd.gen()
    assert f(beta) == 0
    Sbeta = prod( (ai + bi * beta) ** mi  for (ai,bi,mi) in abms)
    y = polygen(Fpd, 'y')
    root = (y**e - Sbeta).roots()[0][0]
    # The root, as a polynomial in beta, has coefficients that agree mod p with those of the eth root of S(alpha).
    return (p, [c.lift() for c in root.list()])

def do_test(nbits=1000, e=3):
    N = random_prime(2**(nbits//2)) * random_prime(2**(nbits//2 + 1))
    f, P = polyselect(N, 5)
    print(f)

    Zalpha = ZZ.extension(f, 'alpha')
    alpha = Zalpha.gen(1)

    abms = [(randint(-1000,1000), randint(-1000,1000), randint(1,4)) for _ in range(10)]
    solution = prod((a + b*alpha)**m for (a,b,m) in abms)
    target = prod((a + b*alpha)**(e*m) for (a,b,m) in abms)
    print("Target:", target.list())
    print("Solution:", solution.list())
    # raise to e
    abms_new = [(a,b,e*m) for (a,b,m) in abms]
    #ethroot(abms,f,e,N)
    gen = ethroot(abms_new, f, e, N)
    if True:
        local_solns = [next(gen) for _ in range(100)] # XXX: using 100 primes, but should instead compute a size bound
        ps, coeffs_mod_ps = zip(*local_solns)
        ps = list(ps)
        prodps = prod(ps)
        coeffs_mod_ps = list(coeffs_mod_ps)
        coeffs_by_position = list(itertools.zip_longest(*coeffs_mod_ps, fillvalue=0)) # [[x^0 coeff mod each p], [x^1 coeff mod each p], ...]
        coeffs_over_z = [crt(list(coeffs_i), ps) for coeffs_i in coeffs_by_position]
        coeffs_over_z = [(c - prodps) if c >= prodps//2 else c for c in coeffs_over_z]
        print("Found:", coeffs_over_z)
        if all(foundi == soli for (foundi,soli) in zip(coeffs_over_z, solution.list())):
            print("Success!")
        else:
            print("Failure!")
    else:
        p, coefflist = next(gen)
        print("p:", p)
        print("sol % p:", [c % p for c in solution.list()])
        print("coefflist:", coefflist)


# do_test()
