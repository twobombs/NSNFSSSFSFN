import struct
from relations import indexed_relation
"""
Build a TOC for an indexed relations file, which must be in plain text format instead of gzipped.
The TOC allows random access lookups of indexed relations by giving the byte-offset in the indexed relations file of the relevant line.

In particular the format of the TOC file is one little-endian uint64 per indexed relation, giving the byte offset of that relation in the indexed relations file.

Why not use a different format for the indexed relations file in the first place? Because the cado-nfs tooling expects the indexed relations file to be in the line-oriented format.

ASSUMPTIONS:
 * The indexed relations file is in plain text instead of gzipped
 * The indexed relations file contains only ASCII (we open it in binary instead of text mode so that seeks take O(1) time)
 * The indexed relations file is no more than 2**64 bytes long
 * Comments start with #
"""

class RandomAccessIndexedRelations:
    """
    Indexed relations file supporting array-like access with O(1) lookup time.
    Essentially free to construct (i.e., it doesn't read in the whole file).
    A single instance of RandomAccessIndexedRelations can't safely be used in multiple threads/processes,
    but each thread/process can safely make its own RandomAccessIndexedRelations reading the same file simultaneously.
    Assumes a table of contents has already been created (by running this file or by calling maketoc).

    Example usage:
        with RandomAccessIndexedRelations("aqrels.out.indexed") as rels:
            irel = rels[123]
    """
    def __init__(self, relsfilename, tocfilename=None):
        if tocfilename is None:
            tocfilename = relsfilename + ".toc"
        self.relsfile = open(relsfilename, "rb")
        self.tocfile = open(tocfilename, "rb")
    def __getitem__(self, index):
        self.tocfile.seek(index * 8)
        offset = self.tocfile.read(8)
        if len(offset) != 8:
            raise IndexError("out of range")
        offset = struct.unpack("<Q", offset)[0]
        self.relsfile.seek(offset)
        return indexed_relation(self.relsfile.readline().decode('utf-8'))
    def close(self):
        self.relsfile.close()
        self.tocfile.close()
    def __enter__(self):
        return self
    def __exit__(self, *args, **kwargs):
        self.close()
    def __len__(self):
        raise NotImplementedError("len isn't implemented for RandomAccessIndexedRelations, but it shouldn't be hard to implement if it turns out we need it")
        
def maketoc(infilename, outfilename=None):
    if outfilename is None:
        outfilename = infilename + ".toc"
    with open(infilename, 'rb') as infile:
        with open(outfilename, 'wb') as outfile:
            pos = 0
            for line in infile:
                if line.startswith(b'#'):
                    # Skip comments
                    pos = infile.tell()
                    continue
                outfile.write(struct.pack("<Q", pos))
                pos = infile.tell()

if __name__ == "__main__":
    from sys import argv
    if len(argv) not in [2,3]:
        print(f"Usage: {argv[0]} input-file [output-file]")
        print( "    where input-file is an indexed-relations file (e.g. aqrels.out.indexed)")
        print( "    and output-file defaults to (input-file).toc")
        exit(1)
    infilename = argv[1]
    if len(argv) >= 3:
        outfilename = argv[2]
    else:
        outfilename = argv[1] + ".toc"
    maketoc(infilename, outfilename)
