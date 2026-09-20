from os.path import dirname
import struct

from sage.all import ZZ
from cado_sage import CadoPolyFile

"""
Time/space tradeoff by not loading the whole renumber table into RAM and instead reading entries from disk as needed.
Run this before the montgomery root step: sage random_access_renumber.py input-file
	where input-file is the explain_renumber file (e.g. data/n256/renumber.explained)
	and output is written to e.g. data/n256/renumber.explained.toc, data/n256/renumber_to_column, and data/n256/column_to_renumber
Then run the root step wiht --montgomery-new-renumber
"""

class RandomAccessRenumberTable:
    # Need for montgomery:
    #   renumber_to_column # we can get rid of this if we make ST_alg_vector use renumber indices instead of columns
    #   column_to_sage_ideal # or just renumber index to sage ideal if we do the above
    #   number_of_algebraic_columns
    # Might need elsewhere
    #   column_to_renumber (don't need for montgomery)
    #   lookup _ideal by column (don't need for montgomery)
    #   lookup _ideal by renumber (don't need for montgomery)
    #   various stuff with the rational side, which I'm ignoring for now
    # Load the first 2**31 algebraic into memory for fast access
    # TODO: do side0 index <-> rational prime <-> renumber index
    def __init__(self, explain_renumber_filename, polyfile, tocfilename=None, c2rfilename=None, r2cfilename=None):
        """
        Random-access view of the renumber table. Unlike the regular renumber table, the whole thing isn't loaded into memory.
        Maintains some open file descriptors, so be sure to call close() when finished, or use in a with block. ("with RandomAccessRenumberTable(explain_renumber_filename, polyfile) as R:")
        To create the necessary files from an existing explain_renumber file, use RandomAccessRenumberTable.create(...)
        polyfile can either be a CadoPolyFile or a path to a polyfile (which will then be loaded, taking possibly a tenth of a second or so)

        A single RandomAccessRenumberTable instance can't safely be used in multiple threads, but it's safe for each thread to make its own instance pointing at the same files.
        """
        # The explain renumber file is structured like this:
        # Lines starting with # are comments
        # Ignoring comments, line i corresponds to R._ideals[i], which is "renumber index" i
        # Each ideal has a corresponding side (which is either 0 or 1, except for the J ideal it might sometimes be both)
        # The column is the index among the ideals on side 1
        # So, for montgomery it suffices to be able to map renumber -> column and column -> byte offset
        if not tocfilename: tocfilename = explain_renumber_filename + ".toc"
        if not c2rfilename: c2rfilename = dirname(explain_renumber_filename) + "/column_to_renumber"
        if not r2cfilename: r2cfilename = dirname(explain_renumber_filename)  + "/renumber_to_column"
        self.explainfile = open(explain_renumber_filename, 'rb')
        self.tocfile = open(tocfilename, 'rb')
        self.num_columns = None
        self.c2rfile = open(c2rfilename, 'rb')
        self.r2cfile = open(r2cfilename, 'rb')
        if isinstance(polyfile, str):
            self.poly = CadoPolyFile(polyfile)
            self.poly.read()
        else:
            self.poly = polyfile

    def close(self):
        self.explainfile.close()
        self.tocfile.close()
        self.c2rfile.close()
        self.r2cfile.close()
    def __enter__(self):
        return self
    def __exit__(self, *args, **kwargs):
        self.close()

    @staticmethod
    def create(explain_renumber_filename, tocfilename=None, c2rfilename=None, r2cfilename=None):
        """
        Create the tocfile, c2r, and r2c files from a given explain_renumber file.
        Each of these files will likely be almost the size of the explain_renumber table, but will support random access.

        Format of the tocfile: 
        first uint64: number of algebraic columns
        bytes 8n+8 to 8n+16 (for n >= 0): uint64le byte-offset in explain_renumber file of ideal with renumber index n

        Format of c2rfile:
        bytes 8n to 8n+8 (for n >=0): uint64le renumber index for column number n

        Format of r2cfile:
        bytes 8n to 8n+8 (for n >=0): (signed) int64le column number for renumber index n, or -1 (0xffffffffffffffff) if that renumber index is not on side 1 (and thus has no column number)
        """
        if not tocfilename: tocfilename = explain_renumber_filename + ".toc"
        if not c2rfilename: c2rfilename = dirname(explain_renumber_filename) + "/column_to_renumber"
        if not r2cfilename: r2cfilename = dirname(explain_renumber_filename)  + "/renumber_to_column"

        cur_idx = 0
        cur_side1_idx = 0
        cur_byteoffset = 0
        with open(explain_renumber_filename, 'rb') as explain_file:
            with open(tocfilename, 'wb') as toc_file:
                toc_file.seek(8) # save room for the number of algebraic columns
                with open(c2rfilename, 'wb') as c2r_file:
                    with open(r2cfilename, 'wb') as r2c_file:
                        for line in explain_file:
                            if line[0] == 35: # 35 = ord('#'). This is almost twice as fast as line.startswith(b"#")
                                cur_byteoffset = explain_file.tell()
                                continue

                            parser, *data = line.split()
                            if parser == b"J" and len(data) > 1:
                                # merged J
                                side = tuple(int(s) for s in data)
                            else:
                                side = int(data[0])

                            # renumber index (cur_idx) is at byte offset (cur_byteoffset)
                            toc_file.write(struct.pack("<Q", cur_byteoffset))

                            if side == 1 or (isinstance(side, tuple) and 1 in side):
                                # column (cur_side1_idx) is renumber index (cur_idx)
                                c2r_file.write(struct.pack("<Q", cur_idx))
                                r2c_file.write(struct.pack("<q", cur_side1_idx))
                                cur_side1_idx += 1
                            else:
                                r2c_file.write(b"\xff\xff\xff\xff\xff\xff\xff\xff")
                            cur_idx += 1
                            cur_byteoffset = explain_file.tell()
                # Write the number of algebraic columns to bytes 0-8 of the toc
                toc_file.seek(0)
                toc_file.write(struct.pack("<Q", cur_side1_idx))

    def number_of_algebraic_columns(self):
        if self.num_columns is None:
            self.tocfile.seek(0)
            buf = self.tocfile.read(8)
            self.num_columns = struct.unpack("<Q", buf)[0]
        return self.num_columns

    def renumber_to_offset(self, index):
        self.tocfile.seek(index * 8 + 8)
        offset = self.tocfile.read(8)
        return struct.unpack("<Q", offset)[0]

    def renumber_to_ideal(self, index):
        offset = self.renumber_to_offset(index)
        self.explainfile.seek(offset)
        parser, side, I = parse_ideal(self.explainfile.readline(), self.poly)
        return side, *I

    def renumber_to_sage_ideal(self, index):
        offset = self.renumber_to_offset(index)
        self.explainfile.seek(offset)
        parser, side, I = parse_ideal(self.explainfile.readline(), self.poly)
        return renumber_ideal_to_sage_ideal(parser, side, I, self.poly)

    def column_to_renumber(self, col):
        self.c2rfile.seek(col * 8)
        buf = self.c2rfile.read(8)
        return struct.unpack("<Q", buf)[0]

    def renumber_to_column(self, index):
        self.r2cfile.seek(index * 8)
        buf = self.r2cfile.read(8)
        return struct.unpack("<q", buf)[0]

    def column_to_ideal(self, col):
        return self.renumber_to_ideal(self.column_to_renumber(col))

    def column_to_sage_ideal(self, col):
        return self.renumber_to_sage_ideal(self.column_to_renumber(col))

    def side_and_index_to_ideal(self, side, i):
        # We only ever use this with side=1, which is the same as column_to_ideal.
        if side != 1: raise NotImplementedError("side_and_index_to_ideal only implemented for side 1")
        return self.column_to_ideal(i)


def parse_ideal(line, poly):
    """
    Parse a (non-comment) line of the explain renumber file. Adapted from CadoExplainRenumberFile.read() in helpers.py
    Return (parser, side, I)

    line must be a bytes object, not a string
    """
    parser, *data = line.split()

    if parser == b'J':
        # the "data" field is actually a bit of a lie, b
        has_merged_J = len(data) > 1
        if has_merged_J:
            # We're going to cheat, and not return an ideal
            side = tuple(int(s) for s in data)
            I = tuple([f"J{s}" for s in side])
        else:
            side = int(data[0])
            I = (f"J{side}",)
        #self.index_of_J.append(len(self._ideals))
    elif parser == b'rat':
        side, p = data
        side = int(side)
        p = ZZ(p)
        I = (p,)
    elif parser == b'proj':
        side, p = data
        side = int(side)
        p = ZZ(p)
        I = (p,)
        # we're going to have problems here if we have projective
        # primes.
    elif parser == b'easy':
        side, p, r = data
        side = int(side)
        p = ZZ(p)
        r = ZZ(r)
        I = (p, r)
    elif parser == b'generic':
        side, p, denom, *coeffs = data
        side = int(side)
        p = ZZ(p)
        denom = ZZ(denom)
        try:
            theta = poly.K[side]([ZZ(c) for c in coeffs]) / denom
        except Exception as e:
            print(side, p, denom, coeffs)
            raise e
        I = (p, theta)
    return (parser, side, I)

def renumber_ideal_to_sage_ideal(parser, side, Idata, poly):
        K = poly.K
        J = poly.nt.J()
        OK = poly.nt.maximal_orders()

        if parser == b'J':
            # If we have merged J, the existing column_to_sage_ideal functions use side 1
            if isinstance(side, tuple):
                side = 1
            return J[side]
        elif parser == b'rat':
            (p,) = Idata
            I = OK[side].fractional_ideal(p)
        elif parser == b'proj':
            (p,) = Idata
            I = OK[side].fractional_ideal(p) + J[side]
        elif parser == b'easy':
            (p, r) = Idata
            I = OK[side].fractional_ideal(p, K[side].gen() - r) * J[side]
        elif parser == b'generic':
            (p, theta) = Idata
            I = OK[side].fractional_ideal(p, theta)
        else:
            raise AssertionError("Unknown parser:", parser)
        return I

if __name__ == "__main__":
    from sys import argv
    if len(argv) < 2:
        print("Usage: sage random_access_renumber.py input-file")
        print("    where input-file is the explain_renumber file (e.g. data/n256/renumber.explained)")
        print("    and output is written to e.g. data/n256/renumber.explained.toc, data/n256/renumber_to_column, and data/n256/column_to_renumber")
        exit(1)
    RandomAccessRenumberTable.create(argv[1])
