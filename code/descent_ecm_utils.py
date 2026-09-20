from sage.all import *
from sage.libs.libecm import ecmfactor
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor as ProcessPool
import multiprocessing

def isqrt(n):
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

def myxgcd(a, b, T):
    '''
    I believe this stands for "multiplicative" gcd.
    Returns [[b, x], [a, lastx]]
    such that input_a = b/x = a/lastx mod input_b.
    '''
    assert type(a) is int
    assert type(b) is int
    assert type(T) is int
    # ainit = a
    # binit = b
    bound = isqrt(b * T)
    x = 0
    lastx = 1
    y = 1
    lasty = 0
    while abs(b) > bound:
        q = a // b
        r = a % b
        a = b
        b = r
        newx = lastx - q * x
        lastx = x
        x = newx
        newy = lasty - q * y
        lasty = y
        y = newy
    return [[b, x], [a, lastx]]

def filter_with_ecm(N, B1, ncurves):
    '''
    Return two lists: the "prime" factors of N, and the remaining cofactors.
    '''
    if N.is_pseudoprime():
        return [[N], []]
    current_comp = N
    list_p = []
    list_c = []
    for i in range(ncurves):
        res = ecmfactor(current_comp, B1)
        if res[0]:
            fac = res[1]
            current_comp = current_comp // fac
            if fac.is_pseudoprime():
                list_p.append(fac)
            else:
                list_c.append(fac)
            if current_comp.is_pseudoprime():
                list_p.append(current_comp)
                return (list_p, list_c)
    list_c.append(current_comp)
    return (list_p, list_c)

def filter_both_with_ecm(N1, N2, smooth_bound, cofac_bound, B1, ncurves):
    output = []
    primes = []
    for N in (N1, N2):
        print(f"Checking number {N}")
        outN = filter_with_ecm(N, B1, ncurves)
        for p in outN[0]:
            if p > smooth_bound:
                print(f"  The prime factor {p} is too large ({p.nbits()} bits)")
                return (False, None)
        for c in outN[1]:
            if c > cofac_bound:
                print(f"  The composite factor {c} is too large ({c.nbits()} bits)")
                return (False, None)
        output.append(outN[1])
        primes.append(outN[0]);
    return (True, output, primes)

def try_pair(N1, N2, smooth_bound, cofac_bound, B1, ncurves, line):
    print("***** trying a new pair *****")
    print("ECM filtering:")
    out = filter_both_with_ecm(N1, N2, smooth_bound, cofac_bound, B1, ncurves)
    if not out[0]:
        print("pas glop.")
        return (False, None)
    print("Analyzing remaining composite factors:")
    for side in range(2):
        for c in out[1][side]:
            if c < smooth_bound:
                print(f"  The composite factor {c} is smaller than the smoothness bound")
                out[2][side].append(c)
            else:
                print(f"  Full factoring of {c}")
                ret = c.factor()
                print(str(ret))
                if ret[-1][0] > smooth_bound:
                    print(f"  Prime factor {ret[-1]} is too large")
                    print("pas glop.")
                    return (False, None)
                for fac in ret:
                    out[2][side].append(fac[0])

    print("Youpi!")
    return (True, out[2], line)

def try_file(filename, smooth_bound, cofac_bound, B1, ncurves, nthreads=1):
    '''
    filename should contain a list of (a,b) pairs.
    smooth_bound is the bound for known prime factors.
    cofac_bound is the bound on the composite cofactor (or, unknown factors).
    '''
    list_of_cofacs = []
    survivor_lines = []
    with open(filename, "r") as f:
        # Assuming las "survivor file" output
        for line in f.readlines():
            if '#' in line:
                continue
            abcd = line.split()
            assert(len(abcd)==4)
            cofac_side0 = Integer(abcd[2].strip())
            cofac_side1 = Integer(abcd[3].strip())
            list_of_cofacs.append((cofac_side0, cofac_side1))
            survivor_lines.append(line)

    print(f"Number of candidates is {len(list_of_cofacs)}")
    ctx = multiprocessing.get_context('fork')

    batchsize = 44
    num_processed = 0
    winner = (False, None, "")

    while num_processed < len(list_of_cofacs):
        this_batch = range(num_processed, min(num_processed + batchsize, len(list_of_cofacs)))

        with ProcessPool(mp_context=ctx, max_workers=nthreads) as executor:
            futures = [
                executor.submit(try_pair, list_of_cofacs[i][0], list_of_cofacs[i][1], smooth_bound, cofac_bound, B1, ncurves, survivor_lines[i])
                for i in this_batch
            ]
            try:
                for future in concurrent.futures.as_completed(futures):
                    print("reading a completed future...")
                    try:
                        out = future.result()
                        if out[0] == True:
                            print("Got a winner!")
                            #winner = (True, out)
                            winner = out
                            break
                    except Exception:
                        # One process throwing an exception ruins the above loop
                        continue
            except Exception:
                # This is where the timeout is caught
                continue

            print("broke out of the loop...")
            if winner[0] == True:
                print("trying to shut down...")
                executor.shutdown(wait=False, cancel_futures=True)
                for future in futures:
                    try:
                        future.cancel()
                    except Exception:
                        continue
                return winner

            # Else keep going (?)
            # TODO: could also continue the not_done processes
            num_processed += len(this_batch)

    return False, None, ""

def transform_polys_by_q(g, f, q, rho, side):
    # Given the original polynomial pair, we can create a modified pair that will
    # let us account for one large special-q. This q will not show up in sieving or
    # relations, but we will add it back in to relations under the original polynomials.
    # For our purposes, g is the rational side and f is the algebraic side.

    # q must be prime.
    # rho should be (for side 0) a root g.roots(GF(q)).
    # Our new polynomials newg, newf should have a shared root mod newg.resultant(newf). Note that
    # this root is easy to compute assuming newg is linear (the resultant is often composite).
    # In fact this method only works assuming g, newg are linear.

    assert (side in [0,1])
    coeff = myxgcd(int(rho), int(q), 1)
    a0 = coeff[0][0]
    b0 = coeff[0][1]
    a1 = coeff[1][0]
    b1 = coeff[1][1]

    R = parent(f)
    x = R.gen()
    num = a0*x+a1
    den = b0*x+b1

    ff = f(num/den)*(den**f.degree())
    ffn = ff.numerator()
    gg = g(num/den)*(den**g.degree())
    ggn = gg.numerator()

    if side == 0:
        c = ggn.content()   # gcd of all coefficients
        assert c % q == 0
        newg = R(ggn / q)
        newf = ffn
    else:
        c = ffn.content()   # gcd of all coefficients
        assert c % q == 0
        newg = ggn
        newf = R(ffn / q)

    return (newg, newf, coeff)
