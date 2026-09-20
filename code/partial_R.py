import sys
import os
from sage.all import *

# This class has the relevant renumber information for our 1024-bit computation,
# specifically for this one descent.
# All the functions here should allow us to call construct_S.

class PartialRenumber1024(object):

    def __init__(self, polynomial):
        self._number_of_rational_queries_36 = 2874398515
        self._number_of_fb_valuations_31 = 105113360
        self._number_of_fb_valuations_35 = 1480227083
        self._renumber_to_rational_prime_index = dict()
        self._rational_prime_to_prime_index = dict()
        self._rational_prime_index_to_prime = dict()
        self._renumber_to_column = dict()
        self._column_to_renumber = dict()
        self._column_to_sage_ideal = dict()
        self._easy_p_r_to_renum = dict()
        self._ideals = dict()   # similar to usual R. dict rather than list. side 1 only. no J.
        self.poly = polynomial

    def number_of_rational_queries(self):
        return self._number_of_rational_queries_36

    def number_of_fb_valuations(self, given_boundA):
        assert given_boundA == 2**31 or given_boundA == 2**35
        if given_boundA == 2**31:
            return self._number_of_fb_valuations_31
        else:
            return self._number_of_fb_valuations_35

    def renumber_to_rational_prime_index(self, index):
        if index not in self._renumber_to_rational_prime_index:
            print(f"WARNING: index={index} not in renumber_to_rational_prime_index.")
            raise KeyError
        return self._renumber_to_rational_prime_index[index]

    def rational_prime_to_prime_index(self, prime):
        if prime not in self._rational_prime_to_prime_index:
            print(f"WARNING: prime={prime} not in rational_prime_to_prime_index.")
            raise KeyError
        return self._rational_prime_to_prime_index[prime]

    def rational_prime_index_to_prime(self, index):
        if index not in self._rational_prime_index_to_prime:
            print(f"WARNING: index={index} not in rational_prime_index_to_prime.")
            raise KeyError
        return self._rational_prime_index_to_prime[index]

    def renumber_to_column(self, index):
        if index not in self._renumber_to_column:
            return -1
        return self._renumber_to_column[index]

    def column_to_renumber(self, col_index):
        if col_index not in self._column_to_renumber:
            print(f"WARNING: col_index={col_index} not in column_to_renumber.")
            raise KeyError
        return self._column_to_renumber[col_index]

    def column_to_sage_ideal(self, col_index, side_hint=None):
        # Needed for run_S_sanity_checks
        if col_index not in self._column_to_sage_ideal:
            print(f"WARNING: col_index={col_index} not in column_to_sage_ideal.")
            raise KeyError
        return self._column_to_sage_ideal[col_index]

    def side_and_index_to_ideal(self, side, side_restricted_index):
        # could be defined for side 0 if needed
        assert side == 1
        if side_restricted_index not in self._ideals:
            print(f"WARNING: col_index={side_restricted_index} not in side_and_index_to_ideal.")
            raise KeyError
        I = self._ideals[side_restricted_index]
        return side, *I

    def easy_p_r_to_renum(self, side, p, r):
        # return the decimal renumber index for this (p,r) ideal
        # this is ONLY called on 'easy' extension ideals
        p = ZZ(p)
        r = ZZ(r)
        assert side == 1
        if (side,p,r) not in self._easy_p_r_to_renum:
            print(f"WARNING: entry=({side},{p},{r}) not in easy_p_r_to_renum.")
            raise KeyError
        return self._easy_p_r_to_renum[(side,p,r)]

    def load(self, filename):
        # lines in filename look like
        # <renumber index> <side-restricted index> <the entire line from renumber.explained>
        with open(filename, 'r') as infile:
            for line in infile:
                line = line.strip()

                renum_index_str, side_restricted_index_str, parser, *data = line.split()
                renum_index = int(renum_index_str)
                side_restricted_index = int(side_restricted_index_str)

                assert parser in ['J', 'rat', 'proj', 'easy', 'generic']
                assert side_restricted_index <= renum_index

                if parser == 'J':
                    assert renum_index == 0
                    assert side_restricted_index == 0
                    self._renumber_to_column[renum_index] = side_restricted_index
                    self._column_to_renumber[side_restricted_index] = renum_index
                    self._column_to_sage_ideal[side_restricted_index] = self.poly.nt.J()[1]

                elif parser == 'rat':
                    side, rat_prime = data
                    side = int(side)
                    assert side == 0
                    rat_prime = int(rat_prime)
                    self._renumber_to_rational_prime_index[renum_index] = side_restricted_index
                    self._rational_prime_to_prime_index[rat_prime] = side_restricted_index
                    self._rational_prime_index_to_prime[side_restricted_index] = rat_prime

                elif parser == 'proj':
                    side, p = data
                    side = int(side)
                    assert side == 1
                    p = ZZ(p)
                    OK = self.poly.nt.maximal_orders()
                    I = OK[side].fractional_ideal(p) + self.poly.nt.J()[side]

                    self._renumber_to_column[renum_index] = side_restricted_index
                    self._column_to_renumber[side_restricted_index] = renum_index
                    self._column_to_sage_ideal[side_restricted_index] = I
                    self._ideals[side_restricted_index] = (p,p)

                elif parser == 'easy':
                    side, p, r = data
                    side = int(side)
                    assert side == 1
                    p = ZZ(p)
                    r = ZZ(r)
                    OK = self.poly.nt.maximal_orders()
                    K = self.poly.K
                    J = self.poly.nt.J()
                    I = OK[side].fractional_ideal(p, K[side].gen() - r) * J[side]

                    self._renumber_to_column[renum_index] = side_restricted_index
                    self._column_to_renumber[side_restricted_index] = renum_index
                    self._column_to_sage_ideal[side_restricted_index] = I
                    self._ideals[side_restricted_index] = (p,r)
                    self._easy_p_r_to_renum[(side,p,r)] = renum_index

                elif parser == 'generic':
                    side, p, denom, *coeffs = data
                    side = int(side)
                    assert side == 1
                    p = ZZ(p)
                    denom = ZZ(denom)
                    OK = self.poly.nt.maximal_orders()
                    K = self.poly.K
                    J = self.poly.nt.J()

                    try:
                        theta = K[side]([ZZ(c) for c in coeffs]) / denom
                    except Exception as e:
                        print(side, p, denom, coeffs)
                        raise e

                    I = OK[side].fractional_ideal(p, theta)
                    self._renumber_to_column[renum_index] = side_restricted_index
                    self._column_to_renumber[side_restricted_index] = renum_index
                    self._column_to_sage_ideal[side_restricted_index] = I
                    self._ideals[side_restricted_index] = (p,theta)

    def do_sanity_checks(self):
        assert self.renumber_to_column(0) == 0
        assert self.column_to_renumber(0) == 0
        assert self.rational_prime_to_prime_index(2) == 0
        assert self.rational_prime_to_prime_index(7) == 3
        assert self.renumber_to_column(14) == -1
        assert self.renumber_to_column(16) == 13
        assert self.renumber_to_rational_prime_index(14) == 2
        assert self.column_to_renumber(8) == 8
        assert self.column_to_renumber(12) == 15
