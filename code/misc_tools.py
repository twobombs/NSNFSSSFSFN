import re
import os
import subprocess
import logging
import json

from math import sqrt, floor, ceil

ncpus = os.cpu_count() // 2

def cat_or_zcat(filename):
    """
    Inspired by `tests/sagemath/cado_sage/tools.py` but faster!

    return an iterable over the file contents, which might
    involve decrypting it.

    Preference order (depending on availability): xopen, zcat.
    """

    try:
        from xopen import xopen
        return xopen(filename, mode="r", threads=ncpus)
    except ModuleNotFoundError:
        logging.warning("Install xopen for faster decompression.")

    if re.search(r"\.gz$", filename):
        # Use pigz binary if available.
        if os.path.isfile("/usr/bin/pigz"):
            return subprocess.Popen([
                    "pigz",
                    "-cd",
                    filename
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE).stdout
        else:
            return subprocess.Popen(["zcat", filename],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE).stdout
    else:
        return open(filename)

def try_xopen_write(filename: str, content):
    if type(content) == str:
        content = content.encode()

    try:
        from xopen import xopen
        with xopen(filename, mode="wb", format="zst", threads=ncpus) as fp:
            fp.write(content)
        return True
    except ModuleNotFoundError:
        logging.warning("Install xopen for faster compression.")

    return False

def try_xopen_load(filename: str):
    try:
        from xopen import xopen
        with xopen(filename, mode="rb", threads=ncpus) as fp:
            return fp.read()
    except ModuleNotFoundError:
        logging.warning("Install xopen for faster decompression.")

    return None

def fast_persistent_save(obj, filename, compress=True):
    """
    Save a sage object `obj` persistently and compressed in file `filename`
    using faster, parallelized compression (if available) than the sage function
    sage.misc.persist.save.

    If compress = True tries first xopen (fastest) then sage.misc.persist.save.

    Else (compress = False), use sage.misc.persist.save with disabled compression.
    """

    import sage.misc.persist
    
    if compress:
        obj_bytes = sage.misc.persist.dumps(obj, compress=False)
        if try_xopen_write(filename, obj_bytes):
            return
        sage.misc.persist.save(obj, filename, protocol=-1)
    else:
        sage.misc.persist.save(obj, filename, protocol=-1, compress=False)

def fast_persistent_load(filename):
    """
    Load the object stored in file `filename` using faster, parallelized
    decompression (if needed) than the sage function sage.misc.persist.load.

    Tries first xopen (fastest), then sage.misc.persist.load.
    """

    import sage.misc.persist
    obj = try_xopen_load(filename)

    if obj == None:
        return sage.misc.persist.load(filename)

    return sage.misc.persist.loads(obj)

def fast_json_dumps(d, *args, decode=True, **kwargs):
    """
    Export dictionary `d` as JSON string, using a faster json package
    than json.dumps (if available).
    """

    try:
        import orjson
        j = orjson.dumps(d, *args, **kwargs)
        if decode:
            j = j.decode()
        return j
    except ModuleNotFoundError:
        logging.warning("Install orjson for faster dumping. Falling back to json.dumps")
        j = json.dumps(d, *args, **kwargs)
        if not decode:
            j = j.encode()
        return j

def fast_json_dump(d, filename, *args, compress="auto", **kwargs):
    """
    Save the dictionary `d` to file `filename`, using a faster json package
    than json.dump (if available) and compression (if the disctionary is large
    or compress = True).
    """

    d_json = fast_json_dumps(d, *args, decode=False, **kwargs)

    if compress == "auto":
        # heuristic boundary.
        compress = len(d) > 5000000

    if compress:
        if try_xopen_write(filename, d_json):
            return
        logging.warning("Falling back to uncompressed storage")

    with open(filename, "w") as fp:
        fp.write(d_json.decode())

def fast_json_loads(d_str, *args, **kwargs):
    """
    Convert JSON string `d_str` to dictionary, using a faster json package
    than json.load (if available).
    """

    try:
        import orjson
        return orjson.loads(d_str, *args, **kwargs)
    except ModuleNotFoundError:
        logging.warning("Install orjson for faster JSON loading. Falling back to json.load.")

    return json.loads(d_str, *args, **kwargs)

def fast_json_load(filename, *args, **kwargs):
    """
    Load a dictionary in json format from file `filename`, using a faster json package
    than json.load (if available) and fast decompression (if necessary).
    """

    d_str = try_xopen_load(filename)

    if d_str == None:
        logging.warning("Falling back to json.load")
        with open(filename, "r") as fp:
            return json.load(fp, *args, **kwargs)

    return fast_json_loads(d_str, *args, **kwargs)

def read_until(in_fp, expected_char, buf):
    """
    Read from in_fp until and including `expected_char`, returning characters read too much in buf.
    """
    buf_size = 1024
    new_data = buf
    while True:
        char_pos = new_data.find(expected_char)
        if char_pos >= 0:
            break
        new_chunk = in_fp.read(buf_size)
        if new_chunk == b"":
            return None, new_data
        new_data += new_chunk
    return new_data[:char_pos+1], new_data[char_pos+1:]

def fast_json_query_load(filename):
    """
    Read queries (key, value) format (int, int), streaming the file manually,
    which is faster for very large files.
    """
    logging.info(f"Starting to parse json from file {filename} using custom parser")
    json_format = re.compile(rb'\s*"(-?\d+)"\s*:\s*"(-?\d+)"[,}]\s*')

    with open(filename, "rb") as fp:
        data, buf = read_until(fp, b'{', b"")
        assert data == b'{'

        while data != None:
            data, buf = read_until(fp, b',', buf)
            to_parse = data if data != None else buf
            m = json_format.search(to_parse)
            if not m:
                logging.error(f"Reached unexpected end of json, ending with: {to_parse}")
                break
            k, v = map(int, m.groups())
            yield (k, v)

        closing_bracket_pos = buf.find(b"}")
        if closing_bracket_pos >= 0:
            buf = buf[closing_bracket_pos+1:]
        if buf.strip() != b"":
            logging.error(f"Unexpected trailing data {buf}")

def find_factors_close_to_square_root(i):
    """
    Given i, find a * b < i such that a, b are close to sqrt(i).
    Returns string "{a}x{b}".
    """
    isqrt = sqrt(i)
    a = floor(isqrt)
    b = ceil(isqrt)
    suboptimal_sol = None

    while a * b != i:
        if a * b > i:
            a -= 1
        elif a * b < i:
            if not suboptimal_sol:
                suboptimal_sol = (a, b)
            b += 1

        # Accept suboptimal solution before it factors diverge too much
        if b - a > isqrt // 2 and suboptimal_sol:
            a, b = suboptimal_sol
            break

    return f"{a}x{b}"
