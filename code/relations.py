from sage.all import *
from misc_tools import cat_or_zcat
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.misc.misc_c import prod
from sage.arith.misc import inverse_mod
from sage.rings.finite_rings.integer_mod_ring import Integers
from tempfile import TemporaryDirectory
from pathlib import Path
from cado_nfs_binaries import CadoNFS
from candy import major_message, error_message
import os
import re
from timing import *

class las_relation:
    """
    a "las"-relation is a relation as it is output by the las program,
    which means that we have the following data, colon-separated.
     - comma-separated a,b in decimal, representing a-bX
     - comma-separated primes on side 0, each in hex, and possibly repeated
     - comma-separated primes on side 1 (if there is a side 1).
    """
    def __init__(self, line):
        if line.startswith('#'):
            # we should not reach here: the caller must arrange for the
            # ctor to be called only on non-comment data
            raise KeyError
        ab, *sides = line.strip().split(':')
        a,b = ab.split(',')
        sides = [[int(x, 16) for x in s.split(',')] if s else [] for s in sides]
        self.a = int(a, 10)
        self.b = int(b, 10)
        self.sides = sides

    def __str__(self):
        return ":".join([f"{self.a},{self.b}",
                         *[",".join([f"{x:x}" for x in s]) for s in self.sides]])

    def norm(self, poly, i):
        """
        return the norm on side i
        """
        x = polygen(QQ, 'x')
        phi = self.a - self.b * x
        return poly.f[i].resultant(phi)

    def check(self, poly):
        """
        given a CadoPolyFile argument, check this relation for
        consistency
        """

        # This check is not terribly
        # useful to do: we don't expect las to fail, here.
        for i, s in enumerate(self.sides):
            if not s:
                continue
            try:
                assert(abs(self.norm(poly, i)) == prod(s))
            except AssertionError:
                raise RuntimeError(f"norm check failed on side {i}, norm={self.norm(poly, i)}, rel={self}")

    def pull_large_factors(self, bounds, avoid_sq=None, missing_ideals=None):
        """
        This modifies the relation and removes the primes above the given
        bounds. If avoid_sq is given, the corresponding special-q is also
        removed from the relation.
        The primes above the bounds, on each side, are returned.
        """

        missing_q = set()   # side=1 always.
        if missing_ideals is not None:
            missing_q = missing_ideals

        newsides = []
        rough = []
        for i, s in enumerate(self.sides):
            newsides.append([])
            rough.append([])
            for p in s:
                if (i, p) == avoid_sq:
                    continue
                if p >= bounds[i]:
                    rough[i].append(p)
                    continue

                is_smooth = True

                if i == 1:
                    # check for missing_q
                    if gcd(self.b, p) != 1:
                        r = p
                    else:
                        r = (self.a * inverse_mod(self.b, p)) % p

                    if (1,p,r) in missing_q:
                        is_smooth = False

                if is_smooth:
                    newsides[i].append(p)
                else:
                    rough[i].append(p)

        self.sides = newsides
        return rough

class indexed_relation:
    """
    an "indexed"-relation is a relation as it is output by the dup2,
    purge, or merge programs,
    which means that we have the following data, colon-separated.
     - comma-separated a,b in hexadecimal, representing a-bX
     - comma-separated indices, whose corresponding ideal is to be
       fetched from the renumber table
    """
    def __init__(self, line):
        if line.startswith('#'):
            # we should not reach here: the caller must arrange for the
            # ctor to be called only on non-comment data
            raise KeyError
        ab, indices = line.strip().split(':')
        a,b = ab.split(',')
        self.indices = []
        if indices:
            self.indices = [int(x, 16) for x in indices.split(',')]
        self.a = int(a, 16)
        self.b = int(b, 16)

    def __str__(self):
        return ":".join([f"{self.a:x},{self.b:x}",
                         ",".join([f"{x:x}" for x in self.indices])])



def convert_to_indexed_relation(rels, params, filename=None):
    """
    (ab)Use dup2 to convert these las relations to indexed ones, and return
    the indexed relations.
    """
    if filename is None:
        filename = os.path.join(params.dirs['TEMP_OUTPUT_DIR'], "rels")

    nrels = 0
    with open(filename, 'w') as f:
        for r in rels:
            print(r, file=f)
            nrels += 1

    if rels and nrels:
        CadoNFS("filter/dup2",
                "-poly", 'POLY',
                "-renumber", 'RENUMBER',
                "-nrels", nrels,
                "-dl",
                'IN',
                inputs={
                    'POLY': params.files['POLYFILE'],
                    'RENUMBER': params.files['RENUMBERFILE'],
                    'IN': filename,
                    },
                outputs={'OUT': filename}
                )
    else:
        major_message(f"{filename} contains no relations, leaving it as it is")

    return indexed_relations_from_file(filename)


def big_convert_to_indexed_relation(rels, params, filename=None):
    """
    Some of the a,b may be over 64 bits, which means dup2 will not
    correctly renumber. With the small enough ones, we use dup2.
    With the large a,b we do our best.
    """
    if filename is None:
        filename = os.path.join(params.dirs['TEMP_OUTPUT_DIR'], "rels")

    bigrels = []
    n_smallrels = 0
    file_big = filename + ".big"
    file_small = filename + ".small"
    # combined output will be written to filename

    with open(file_small, 'w') as f:
        for r in rels:
            a = Integer(r.a)
            b = Integer(r.b)
            if a.nbits() < 64 and b.nbits() < 64:
                print(r, file=f)
                n_smallrels += 1
            else:
                bigrels.append(r)

    if n_smallrels > 0:
        CadoNFS("filter/dup2",
                "-poly", 'POLY',
                "-renumber", 'RENUMBER',
                "-nrels", n_smallrels,
                "-dl",
                'IN',
                inputs={
                    'POLY': params.files['POLYFILE'],
                    'RENUMBER': params.files['RENUMBERFILE'],
                    'IN': file_small,
                    },
                outputs={'OUT': file_small}
                )
    else:
        major_message(f"{file_small} contains no relations, leaving it as it is")

    if len(bigrels) > 0:
        alg_abp_to_renum = dict()   # (a, b, alg p) --> renumber index
        alg_abp_to_multiplicity = dict()
        rat_p_to_renum = dict()     # (rat p) --> renumber index
        rat_seen_factors = set()

        for r in bigrels:
            a = Integer(r.a)
            b = Integer(r.b)
            fac0 = r.sides[0]
            fac1 = r.sides[1]
            rat_seen_factors.update(set(fac0))
            for f1 in fac1:
                f1 = Integer(f1)
                if (a,b,f1) not in alg_abp_to_multiplicity.keys():
                    alg_abp_to_multiplicity[(a,b,f1)] = 1
                    alg_abp_to_renum[(a,b,f1)] = None
                else:
                    alg_abp_to_multiplicity[(a,b,f1)] += 1

        # Fill in the dictionary mappings

        nlines = 0
        rat_ordering = list(rat_seen_factors)
        rat_ordering.sort()
        alg_ordering = []

        with open(file_big, 'w') as f:
            # Do all the rational primes in one go (the first line)
            a = 10
            b = 10
            f.write(str(a) + "," + str(b) + ":")
            rfac_hex = [ hex(r)[2:] for r in rat_ordering ]     # don't want leading 0x
            f.write(",".join(rfac_hex))
            f.write(":\n")
            nlines += 1
            uniqueness_mult = 10

            # Next do the algebraic primes. One line for each a-b-p
            for (a,b,p) in alg_abp_to_renum.keys():

                if p.nbits() < 10:
                    # For the "bad" and "exceptional" ideals things seem to be done mod p^k rather than
                    # mod p, where k is usually at most 2, and the size of p is usually < 256 or so.
                    mm = p**3
                elif p.nbits() < 20:
                    mm = p**2
                else:
                    mm = p

                new_a = uniqueness_mult*mm + int(mod(a,mm))      # just want to avoid tiny a
                new_b = uniqueness_mult*mm + int(mod(b,mm))      # just want to avoid tiny b
                multiplicity = alg_abp_to_multiplicity[(a,b,p)]
                f.write(str(new_a) + "," + str(new_b) + ":")
                f.write(":")
                f.write(hex(p)[2:])
                for i in range(multiplicity-1):
                    f.write("," + hex(p)[2:])
                f.write("\n")
                nlines += 1
                uniqueness_mult += 1
                alg_ordering.append((a,b,p))
                print("a", str(a))
                print("b", str(b))
                print("p", str(p))
                print("new_a", str(new_a))
                print("new_b", str(new_b))

        CadoNFS("filter/dup2",
                "-poly", 'POLY',
                "-renumber", 'RENUMBER',
                "-nrels", nlines,
                "-dl",
                'IN',
                inputs={
                    'POLY': params.files['POLYFILE'],
                    'RENUMBER': params.files['RENUMBERFILE'],
                    'IN': file_big,
                    },
                outputs={'OUT': file_big}
                )
    else:
        major_message(f"{file_big} contains no relations, leaving it as it is")

    if os.path.isfile(file_big):
        with open(file_big, "r") as f:
            ratline = f.readline().strip()
            indices = ratline.split(":")[1].split(",")
            assert len(indices) == len(rat_ordering) + 1    # 0 always appended
            indices_int = [int(x,16) for x in indices if int(x,16) != 0]
            indices_int.sort()
            assert len(indices_int) == len(rat_ordering)
            for i in range(len(indices_int)):
                p_i = rat_ordering[i]
                index_p_i = hex(indices_int[i])[2:]
                rat_p_to_renum[p_i] = index_p_i

            # algebraic
            j = 0
            for line in f.readlines():
                ab = line.split(":")[0].split(",")
                a = int(ab[0],16)
                b = int(ab[1],16)
                assert a != 10
                real_a, real_b, real_p = alg_ordering[j]
                print("a", str(a))
                print("b", str(b))
                print("real_a", str(real_a))
                print("real_b", str(real_b))
                print("real_p", str(real_p))

                if real_p.nbits() < 10:
                    mm = real_p**3
                elif real_p.nbits() < 20:
                    mm = real_p**2
                else:
                    mm = real_p

                assert mod(a,mm) == mod(real_a,mm)
                assert mod(b,mm) == mod(real_b,mm)
                indices = line.strip().split(":")[1].split(",")
                assert len(indices) >= 2    # 0 always appended
                indices_int = [int(x,16) for x in indices if int(x,16) != 0]
                assert len(indices_int) >= 1
                # So the thing about bad ideals is that abp maps to multiple ideals
                # over the "bad ideal"
                if len(indices_int) == 1:
                    index = hex(indices_int[0])[2:]
                else:
                    index = ",".join([hex(z)[2:] for z in indices_int])

                # I think the above should account for multiplicity correctly.
                # This comes up with "exceptional ideals" as well, where the number of factors
                # of the given norm matters for mapping back to an ideal. Something about
                # ideals of exceptional inertia...

                alg_abp_to_renum[(real_a,real_b,real_p)] = index
                j += 1

    outfile = open(filename, 'w')

    with open(file_small, 'r') as f:
        for line in f.readlines():
            outfile.write(line)

    for r in bigrels:
        a = Integer(r.a)
        b = Integer(r.b)
        fac0 = r.sides[0]
        fac1 = r.sides[1]
        print("making output")
        print("a", str(a))
        print("b", str(b))

        if a >= 0:
            hexa = hex(a)[2:]
        else:
            hexa = "-" + hex(a)[3:]
        if b >= 0:
            hexb = hex(b)[2:]
        else:
            hexb = "-" + hex(b)[3:]

        outfile.write(str(hexa))
        outfile.write(",")
        outfile.write(str(hexb))
        outfile.write(":")
        renums = []
        for f0 in fac0:
            renums.append(rat_p_to_renum[f0])
        fac1 = list(set(fac1))  # dedupe, since multiplicities should be accounted for in the index
        for f1 in fac1:
            renums.append(alg_abp_to_renum[(a,b,f1)])
        # TODO:
        # I believe one J ideal should be counted in all relations...? Not sure if this matches
        # how we parse generic/proj/easy ideals.
        renums.append("0")
        renum_string = ",".join(renums)
        outfile.write(renum_string)
        outfile.write("\n")

    outfile.close()
    return indexed_relations_from_file(filename)


def las_relations_from_file(filename):
    """
    This function takes a file name, and yields its contents, i.e.,
    relations
    """
    with cat_or_zcat(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            yield las_relation(line)

def indexed_relations_from_file(filename):
    """
    This function takes a file name, and yields its contents, i.e.,
    relations
    """
    with cat_or_zcat(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            yield indexed_relation(line)


def strip_rational_part_of_relations(filename):
    """
    takes a filename containing las relations, and keep only the
    algebraic part of each
    """
    error_message("You should not be calling strip_rational_part_of_relations at all. Use one-side sieving and then swap_parts_of_relations instead")
    with open(filename + ".strip", "w") as outfile:
        for rel in las_relations_from_file(filename):
            rel.sides[0] = []
            print(rel, file=outfile)
    os.rename(filename + ".strip", filename)

@timing
def swap_parts_of_relations(filename, remove_tmp_suffix=False, do_ext_dedup=None, extra_relation_check=None):
    """
    takes a filename containing las relations, and keep only the
    algebraic part of each
    If extra_relation_check is not None, it should be a tuple of:
        (True, renum_info_file)
    Then we assert the two following checks:
        1. All of the factors in the relation are indeed prime
        2. If an ideal is going to be labelled as "projective" it is indeed in the set of allowed
            projective ideals, in the renum_info_file.
    If a relation fails either of the above, we omit it from the output file.
    """
    if extra_relation_check is not None:
        assert extra_relation_check[0]
        renum_info_file = extra_relation_check[1]

        # faster for this specific computation
        if 'n1024' in renum_info_file:
            projective_ideals = [2,3,5,7,13,41,79]

        else:
            # usual case
            projective_ideals = [2]

            with open(renum_info_file, "r") as rf:
                for line in rf:
                    # Looks like:
                    # i=0x10 tab[i]= (0x7,0x7) p=0x7 r=0x7 side 1 proj
                    if line.startswith("#"):
                        continue
                    if 'p=0x8000000b' in line:
                        # TEMPORARY
                        break
                    if 'proj' in line:
                        pstr = line.split()[3]
                        p = int(pstr.split('x')[1], 16)
                        q = Integer(p)
                        projective_ideals.append(q)
                    if 'bad ideal' in line:
                        # an ideal can be both bad AND projective
                        # i=0x5 tab[i]=# bad ideal (number 1/2) above (3,3) on side 1
                        prho = line.split()[7].split(",")
                        p = prho[0][1:]
                        rho = prho[1][:-1]
                        if p == rho:
                            q = Integer(p)
                            projective_ideals.append(q)

        relation_check = True
    else:
        relation_check = False

    if do_ext_dedup is not None:
        assert do_ext_dedup[0]
        BOUNDA_queries = do_ext_dedup[1]
        BOUNDA = do_ext_dedup[2]
    else:
        # check can vacuously fail
        BOUNDA_queries = -1
        BOUNDA = -1

    # If do_ext_dedup=True, we ensure that each "large q" has only one relation.
    # This is specific to extension sieving, where we only want one relation per large q.
    seen_large_qs = set()

    with open(filename + ".strip", "w") as outfile:
        for rel in las_relations_from_file(filename):
            if relation_check:
                a = rel.a
                b = rel.b   # in decimal
                facs = rel.sides[0]     # one-sided relation, in decimal
                good_relation = True
                largeq = None

                for q in facs:
                    if not is_prime(q):
                        # All the factors q should be prime!
                        good_relation = False
                        print(f"Weird! Removed relation ({a},{b}) with composite factor {q}.")
                        break
                    if gcd(b, q) == q and (q not in projective_ideals):
                        # If this happens, q divides b, so this ideal will be treated
                        # as projective, but it is not in the renumber list.
                        good_relation = False
                        print(f"Mysterious! Removed relation ({a},{b}) with faux projective ideal {q}.")
                        break
                    if q > BOUNDA_queries and q < BOUNDA:
                        assert largeq is None    # only one large q
                        if b % q != 0:
                            rho = a*inverse_mod(b,q) % q
                        else:
                            rho = q
                        largeq = (q, rho)

                if do_ext_dedup is not None:
                    assert largeq is not None
                    if largeq in seen_large_qs:
                        good_relation = False
                        print(f"Removed relation ({a},{b}) with duplicate q {largeq}.")
                    else:
                        if good_relation:
                            seen_large_qs.add(largeq)

                if good_relation:
                    rel.sides[1] = rel.sides[0]
                    rel.sides[0] = []
                    print(rel, file=outfile)

            else:
                rel.sides[1] = rel.sides[0]
                rel.sides[0] = []
                print(rel, file=outfile)

    new_filename = filename[:-4] if remove_tmp_suffix else filename
    os.rename(filename + ".strip", new_filename)


def parse_fb_extension_relations(filename, qrange):
    """
    takes a filename containing las relations, and keep only one relation
    per q in the given qrange. Return it as a dictionary.
    The qrange is specified as (side, q0, q1).

    Of course this is used pretty much _only_ in the FB extension
    phase, and it should probably be a subroutine of that code area, in
    fact.
    """

    per_q = dict()

    side, q0, q1 = qrange

    for rel in las_relations_from_file(filename):
        found_a_special_q = False
        for q in rel.sides[side]:
            if q < q0:
                continue
            assert q < q1
            if found_a_special_q:
                raise RuntimeError(f"The following relation has several special-q's in the range {q0}..{q1}. The script probably does not except that.")
            if rel.b % q != 0:
                rho = rel.a*inverse_mod(rel.b,q) % q
                #rho = ZZ(Integers(q)(rel.a)/Integers(q)(rel.b))
            else:
                rho = q
            found_a_special_q = True
            if (side,q,rho) not in per_q:
                per_q[side, q, rho] = rel

    return per_q


def keep_only_one_relation_per_q(filename, qrange, keep_file=False):
    per_q = parse_fb_extension_relations(filename, qrange)

    with open(filename + ".mini", "w") as outfile:
        for qq in sorted(per_q.keys()):
            print(per_q[qq], file=outfile)

    if re.search(r"\.gz$", filename):
        raise NotImplementedError

    if keep_file:
        os.rename(filename, filename + ".old")
    os.rename(filename + ".mini", filename)

    # Return per_q so that missing q's can be found.
    return per_q.keys()
