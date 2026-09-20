from sage.all import *

from sparse_dot_mkl import dot_product_mkl
import numpy as np
from scipy.sparse import csr_array, bmat


class M_wrapper_2(object):
    # Wrapper around the expected M sage matrix.
    # We call nrows() and ncols() frequently.
    # Otherwise it's a scipy 64-bit-integer matrix.

    def __init__(self, scipy_matrix):

        self._nrows = scipy_matrix.dimensions()[0]
        self._ncols = scipy_matrix.dimensions()[1]
        self._scipy_M = scipy_matrix

    def nrows(self):
        return self._nrows

    def ncols(self):
        return self._ncols

    def scipy_M(self):
        return self._scipy_M


class MC_wrapper(object):
    # Wrapper around the expected MC sage matrix.
    # We call nrows() and ncols() frequently.
    # Otherwise we can get MC from M and S_block.

    def __init__(self, nrows, ncols):

        self._nrows = nrows
        self._ncols = ncols

    def nrows(self):
        return self._nrows

    def ncols(self):
        return self._ncols
