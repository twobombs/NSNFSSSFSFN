from sage.structure.factorization import Factorization
from sage.structure.sequence import Sequence
from sage.rings.integer import Integer

"""
sage's Factorization objects sort their contents, which is slow and not needed.
It's possible to initialize Factorizations with sort=False, but many methods will
internally create new Factorization objects without sort=False, and thus sort.

So we throw the whole thing out and replace it with a simple hashmap
"""

class UnsortedFactorization(Factorization):
    def __init__(self, x, unit=None):
        """
        Like a Factorization, but simpler and faster for our purposes
        x can either be a dict (which gets used directly, not copied), or a list of pairs (p,e)
        unit is 1 if not specified
        """
        if not isinstance(x,dict):
            x = {p: e for (p,e) in x if e != 0}
        if unit is None:
            unit = Integer(1)
        self.__unit = unit
        self.__x = x

    def __iter__(self):
        return iter(self.__x.items())

    def __len__(self):
        return len(self.__x)

    def unit(self):
        return self.__unit

    def __rmul__(self,left):
        """return left * self, where left is not a factorization"""
        x = self.__x.copy()
        x[left] = x.get(left, 0) + 1
        return UnsortedFactorization(x, unit)

    def __mul__(self, other):
        unit = self.__unit
        
        if isinstance(other, UnsortedFactorization):
            x = {}
            for a in set(self.__x).union(set(other.__x)):
                x[a] = self.__x.get(a,0) + other.__x.get(a,0)
                if x[a] == 0:
                    del x[a]
            unit *= other.__unit
        elif isinstance(other, Factorization):
            otherx = dict(other)
            x = {}
            for a in set(self.__x).union(set(otherx)):
                x[a] = self.__x.get(a,0) + otherx.get(a,0)
                if x[a] == 0:
                    del x[a]
            unit *= other.unit()
        else:
            x = self.__x.copy()
            x[other] = x.get(other, 0) + 1
        return UnsortedFactorization(x, unit)
                
    def __pow__(self, n):
        n = Integer(n)
        if n == 1: return self
        if n == 0: return UnsortedFactorization({})
        return UnsortedFactorization({p: n*e for (p,e) in self.__x.items()}, unit=self.__unit**n)

    def __invert__(self):
        return UnsortedFactorization({p: -e for (p,e) in self.__x.items()}, unit=self.__unit**(-1))

    def __truediv__(self, other):
        unit = self.__unit
        if not isinstance(other, UnsortedFactorization):
            x = self.__x.copy()
            x[other] = x.get(other, 0) - 1
        else:
            x = {}
            for a in set(self.__x).union(set(other.__x)):
                x[a] = self.__x.get(a,0) - other.__x.get(a,0)
                if x[a] == 0:
                    del x[a]
            unit *= other.__unit
        return UnsortedFactorization(x, unit)

    def value(self):
        if not hasattr(self, "_prod"):
            from sage.misc.misc_c import prod
            self._prod = prod([p**e for (p,e) in self.__x.items()], self.__unit)
        return self._prod

    def cache_product(self, knownprod):
        """
        Set the cached value of self.prod() to an already-computed value.

        For parallelized Montgomery reduction, worker processes in a ProcessPool will compute I.prod(), and we want to cache that value in a way the main process can use later.
        The normal caching in I.prod() won't work because the worker has a different copy of I than the main process (they communicate via pickled objects).
        So the worker will include I.prod() in its return value, and the main process will cache it using this method.
        """
        if hasattr(self,"_prod"):
            print("Warning: called cache_product but prod was already computed!")
            assert knownprod == self._prod, "The product we're tring to cache is wrong, or at least disagrees with the value we already cached"
        self._prod = knownprod

    expand = value
    prod = value

    def gcd(self, other):
        if not isinstance(other, UnsortedFactorization):
            raise TypeError("can't take gcd of factorization and non-factorization")
        x = {}
        for a in set(self.__x).intersection(set(other.__x)):
            x[a] = min(self.__x[a], other.__x[a])
        return UnsortedFactorization(x) # unit disappears

    def lcm(self, other):
        if not isinstance(other, UnsortedFactorization):
            raise TypeError("can't take lcm of factorization and non-factorization")
        x = {}
        for a in set(self.__x).union(set(other.__x)):
            x[a] = max(self.__x.get(a,0), other.__x.get(a,0))
        return UnsortedFactorization(x) # unit disappears
