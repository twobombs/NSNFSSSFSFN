#!/usr/bin/env python3

import argparse
import logging
import time
import json

from os import mkdir
from tqdm import tqdm

from oracle_client import OracleClient
from oracle_helpers import *

from misc_tools import fast_json_dump, fast_json_load

class LunaK6HSMOracleClient(OracleClient):
    def __init__(self, host, port, logging) -> None:
        return super().__init__(host, port, logging)

if __name__ == "__main__":
    ## Command-line arguments
    parser = argparse.ArgumentParser(description='Luna K6 HSM signature oracle client')

    parser.add_argument('--reuse_paths', action='store_true',
                        help='Append to existing log files instead of creating a new subdir in <logs_dir>')
    parser.add_argument('--logs_dir', type=str, default='logs',
                        help='Path to directory in which log files will be stored. (The directory must already exist.)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable extra debug output')
    parser.add_argument('--stderr', action='store_true',
                        help='Log also to stderr, in addition to the log file.')
    parser.add_argument('--host', default='XXXX',
                        help='Hostname or IP address of the oracle server')
    parser.add_argument('--port', type=int, default=61338,
                        help='Port that the oracle server is listening on')
    parser.add_argument('--label', type=str,
                        help='Name of the RSA key to use (default = id_oracle_rsa_{bits}_exp_{e})')
    parser.add_argument('--keygen', action="store_true",
                        help='Boolean to instruct server to generate new key pair.')
    parser.add_argument('-e', '--pubexp', type=int, default=65537,
                        help='RSA public exponent')
    parser.add_argument('-N', '--modulus', type=int,
                        help='RSA modulus')
    parser.add_argument('--bits', type=int,
                        help='Number of bits for the RSA modulus')
    parser.add_argument('--in', dest='infile',
                        help='File with raw RSA integers "a" -- one per line -- to sign')
    parser.add_argument('--out', dest='outfile',
                        help='File to output signatures (i.e., a^d mod N) to.')
    parser.add_argument('--cont', action='store_true',
                        help='Continue a previous run, i.e., don\'t query values that are already written to the outfile.')
    parser.add_argument('--export-key', action='store_true',
                        help='Wrap and export a key, must be used together with --wrapping-key-label.')
    parser.add_argument('--wrapping-key-label', type=str,
                        help='label of the key that should be used to wrap the key specified by --label. Use together with --export-key.')
    parser.add_argument('--wrapping-key-path', type=str,
                        help='Path to the RSA secret key of the wrapping key; to decrypt the key exported with --export-key.')
    parser.add_argument('--export-key-out', type=str,
                        help='Path to write the exported and decrypted key to.')

    args = parser.parse_args()

    do_reuse_paths  = args.reuse_paths
    logs_dir        = args.logs_dir
    host            = args.host
    port            = args.port
    label           = args.label
    do_keygen       = args.keygen
    pub_exp         = args.pubexp
    modulus         = args.modulus
    modulus_bits    = args.bits
    infile          = args.infile
    outfile         = args.outfile
    do_continue     = args.cont

    do_export           = args.export_key
    wrapping_key_label  = args.wrapping_key_label
    wrapping_key_path   = args.wrapping_key_path
    export_key_out      = args.export_key_out

    if not do_reuse_paths:
        subdir_name = 'run_' + time.strftime('%Y%m%d-%H%M%S')
        logs_dir += f'/{subdir_name}'
        mkdir(logs_dir)

    log_fpath = f"{logs_dir}/luna_k6_hsm_oracle_client.log"
    log_level = logging.DEBUG if args.verbose else logging.WARNING

    logging.basicConfig(filename=log_fpath, level=log_level)
    print(f"Luna K6 HSM oracle client logging to {log_fpath}")

    if args.stderr:
        logging.getLogger().addHandler(logging.StreamHandler())

    if do_export and wrapping_key_label == None:
        logging.error("Must specify --wrapping-key-label for --export-key.")

    if not label:
        label = f"id_oracle_rsa_{modulus_bits}_exp_{pub_exp}"
        sk_label = f"{label}_secret_key"
        logging.info(f"Using key label '{label}'")
    else:
        logging.info(f"Using existing key '{label}', cannot generate new one.")
        sk_label = label
        do_keygen = False

    client = LunaK6HSMOracleClient(host, port, logging)

    if do_keygen:
        query = client.build_keygen_query(pub_exp, modulus_bits, label)
        N = client.query_single(OP_KEYGEN, query)
        try:
            N = int(N)
            logging.info(f"Generated RSA key pair with modulus N = {N}")
        except Exception as e:
            logging.error(f"Received unexpected response to key generation, expected int(modulus), error:\n{e}")
            raise e

    if infile != None and outfile != None:
        file_mode = "r+" if do_continue else "w"

        with open(outfile, file_mode) as out_fp:
            answers: frozenset = frozenset()

            if do_continue:
                previous_answers_json = out_fp.read()

                if previous_answers_json[-1] == "}":
                    logging.error("Run seems to have already terminated, cannot continue it.")
                    exit(1)

                answers = frozenset(map(int, json.loads(previous_answers_json[:-2] + "\n}").keys()))
            else:
                # Manual json dump so that we can write answers as they come in to avoid
                # losing data in case the application crashes at some time long into the
                # computation.
                out_fp.write("{\n")

            def do_queries(queries, quotients, N, first_line, fp):
                for query, sig in client.query_all(OP_SIGN, queries):
                    a, _ = query.split(":")

                    if not first_line:
                        fp.write(f",\n")
                    else:
                        first_line = False

                    i = int(a, 16)
                    if i in quotients:
                        i += quotients[i] * N
                    fp.write(f"  \"{i}\": \"{sig}\"")

                return first_line

            cnt = 1
            skipped = 0
            first_line = True
            queries: list[str] = []
            quotients: dict = dict()

            with open(infile, "r") as in_fp:
                for line in tqdm(in_fp):
                    # Send query in badges to avoid giant 'queries' data structure.
                    if cnt % 10000 == 0:
                        first_line = do_queries(queries, quotients, modulus, first_line, out_fp)
                        queries = []
                        quotients = dict()

                    i = int(line)
                    if i < 0:
                        if modulus == None:
                            logging.error(f"Signature queries must be postive, but got: {i}.")
                            exit(1)
                        else:
                            quotient, i = divmod(i, modulus)
                            quotients[i] = quotient
                            i = i % modulus
                    if i in answers:
                        logging.info(f"Already queried {i}, skipping")
                        skipped += 1
                        continue

                    queries.append(client.build_sign_query(i, sk_label))
                    cnt += 1

            do_queries(queries, quotients, modulus, first_line, out_fp)
            out_fp.write("\n}")
            logging.info(f"Skipped in total {skipped} queries.")

        # make sure to store queries in compressed form at the end.
        logging.info(f"Storing queries in compressed form.")
        queries = fast_json_load(outfile)
        fast_json_dump(queries, outfile)

    if do_export:
        query = client.build_wrap_query(wrapping_key_label, sk_label)
        wrapped_key_ctxt = client.query_single(OP_WRAP, query)

        if wrapping_key_path != None:
            from Crypto.PublicKey import RSA
            from Crypto.Cipher import PKCS1_v1_5

            with open(wrapping_key_path, "rb") as fp:
                # Assume no PW protection of key
                wrapping_key = RSA.import_key(fp.read())

            cipher = PKCS1_v1_5.new(wrapping_key)
            wrapped_key_ptxt = cipher.decrypt(wrapped_key_ctxt)

            if export_key_out:
                with open(export_key_out, "wb") as fp:
                    fp.write(wrapped_key_ptxt)
            else:
                print(f"Decrypted wrapped key: {wrapped_key_ptxt}")
        else:
            print(f"Received the wrapped key: {wrapped_key.hex()}")
