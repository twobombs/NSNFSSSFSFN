#/usr/bin/env python3

import logging
import re

END_OF_QUEUE        = None
MAX_MODULUS_SIZE    = 1024
# Space for header (2 byte op, 8 byte ID, modulus in hex, and 2 bytes for separators
BUF_SIZE            = 12 + (MAX_MODULUS_SIZE // 4)

OP_CLOSE_CONN   = 0
OP_KEYGEN       = 1
OP_SIGN         = 2
OP_VERIFY       = 3
OP_WRAP         = 4
OP_ANS_KEYGEN   = 101
OP_ANS_SIGN     = 102
OP_ANS_VERIFY   = 103
OP_ANS_WRAP     = 104

def i2h(i):
    h = "{:x}".format(i)
    if len(h) % 2 == 1:
        h = "0" + h
    return h

def parse_query(query):
    """
    Assumed query format is "op:id:hex", where
        'op' is a 1-byte Integer opcode (in hex),
        'id' is a 4-byte Integer (in hex) and
        'hex' is the data that should be signed in hexadecimal.
    """

    try:
        query_str = query.decode()
        m = re.match(r'([0-9a-f]+)\:([0-9a-f]+)\:(.*)', query_str)
        groups = m.groups()

        if len(groups) != 3:
            logging.error(f"Query has {len(groups)} groups, expected 3.")
            return None

        op, ident, data = groups

        if len(op) > 2:
            logging.error(f"Opcode '{op}' larger than 1 byte.")

        if len(ident) > 8:
            logging.error(f"Opcode '{ident}' larger than 4 bytes.")

        return int(op, 16), int(ident, 16), data
    except Exception as e:
        logging.error(f"Failed to parse query '{query}' with the following error:\n{e}")
        return None

def build_query(op: int, ident: int, data):
    """
    Build a query of the format "op:id:hex", where
        'op' is a 1-byte Integer opcode (in hex),
        'id' is a 4-byte Integer (in hex) and
        'hex' is the conversion of integer data to hexadecimal.
    """

    query = "{:02x}".format(op)
    query += ":{:08x}".format(ident)

    if type(data) == int:
        data = i2h(data)
        if len(data) % 2 == 1:
            data = "0" + data
    elif type(data) == bytes:
        data = data.hex()
    elif type(data) != str:
        logging.error(f"'data' argument of 'build_query' must be int, bytes, or str, got {type(data)}")

    return f"{query}:{data}".encode()
