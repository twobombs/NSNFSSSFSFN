#!/usr/bin/env python3

#
# This script must be run on XXXX (where the physical Luna K6 HSM is installed).
# Run as user sysadmin on XXXX (this requires a lot of custom setup and building pycryptoki from source to work, see /physical_hsms/luna_hsm/README.md):
# $ python3.8 luna_k6_hsm_oracle_server.py --verbose --exec_threads 1

import os
import time
import argparse
import logging
import tempfile
import subprocess

from oracle_server import OracleServer
from luna_k6_hsm_crypto_helpers import *

class LunaK6HSMOracleServer(OracleServer):
    def __init__(self, logs_dir, logging, userpin, *args) -> None:
        self.logs_dir           = logs_dir
        self.logging            = logging
        self.key_info_cache     = dict()
        self.sk_handel_cache    = dict()

        self.logging.info("Starting Luna K6 HSM session...")
        self.session = start_session(userpin=userpin)
        self.logging.info("Luna K6 HSM session started successfully")

        return super().__init__(logging, *args)

    def _get_pk_label(self, label: str):
        return f"{label}_public_key"

    def _get_sk_label(self, label: str):
        return f"{label}_secret_key"

    def _get_key_info(self, label, handel=None):
        if label not in self.key_info_cache:
            if handel == None:
                handels = find_existing_rsa_key(self.session, label)
                if len(handels) != 1:
                    self.logging.error(f"Couldn't find key with label {label}")

                handel = handels[0]

            N, e = get_pub_key_info(self.session, handel)
            modulus_bits = ((N.bit_length() + 7) // 8) * 8

            self.key_info_cache[label] = (N, e, modulus_bits)

        N, e, modulus_bits = self.key_info_cache[label]
        self.logging.info(f"[_get_key_info] Retrieved following info for public key with label {label}: N={N}, e={e}, modulus_bits={modulus_bits}")
        return N, e, modulus_bits

    def _get_rsa_sk_handel(self, rsa_sk_label, rsa_sk_handel=None):
        log_tag = f"[_get_rsa_sk_handel] "
        self.logging.info(f"{log_tag}Retrieve secret key handel for label {rsa_sk_label}")

        if rsa_sk_label not in self.sk_handel_cache:
            if rsa_sk_handel == None:
                rsa_sk_handels = find_existing_rsa_key(self.session, rsa_sk_label)

                if len(rsa_sk_handels) != 1:
                    if len(rsa_sk_handels) == 0:
                        error_msg = f"{log_tag}Secret key with label {rsa_sk_label} does not exist, did you use 'keygen' first?"
                    else:
                        error_msg = f"{log_tag}Invalid HSM state"

                    self.logging.error(error_msg)
                    return None

                rsa_sk_handel = rsa_sk_handels[0]

            self.sk_handel_cache[rsa_sk_label] = rsa_sk_handel

        self.logging.info(f"{log_tag}Retrieved handel {rsa_sk_handel}")
        return self.sk_handel_cache[rsa_sk_label]

    def keygen(self, pub_exp: int, modulus_bits: int, label: str):
        """Generate a new RSA key pair on the Luna K6 HSM"""

        log_tag = f"[{label}, keygen]: "

        if not label:
            label = f"RSA_{modulus_bits}_exp_{pub_exp}"

        rsa_sk_label = self._get_sk_label(label)
        rsa_pk_label = self._get_pk_label(label)

        rsa_sk_handels = find_existing_rsa_key(self.session, rsa_sk_label)

        if len(rsa_sk_handels) == 0:
            self.logging.info(f"{log_tag}Generating RSA key pair with sk/pk labels exp={pub_exp} and bits={modulus_bits} ...")
            rsa_pk_handel, rsa_sk_handel = gen_rsa_keys(self.session, rsa_pk_label, rsa_sk_label, pub_exp, modulus_bits)
        else:
            self.logging.warning(f"{log_tag}RSA key pair for label {label} already exists, not generating new one.")
            rsa_pk_handel = rsa_sk_handels[0]

        N, e, _ = self._get_key_info(rsa_pk_label, handel=rsa_pk_handel)

        self.logging.info(f"{log_tag}The RSA key pair is using the following modulus and exponent:\n\tN = {N}\n\te = {e}")

        return N

    def sign(self, data: bytes, label: str):
        """Sign 'data' using the RSA key specified by 'label'"""

        log_tag = f"[{label}, sign]: "

        N, e, modulus_bits = self._get_key_info(label)
        data_pad = pad(data, modulus_bits)

        rsa_sk_handel = self._get_rsa_sk_handel(label)

        self.logging.info(f"{log_tag}sign the following data with key with label {label}: {data_pad.hex()}")
        sig = sign(self.session, rsa_sk_handel, data_pad)
        self.logging.info(f"{log_tag}Produced signature {sig.hex()}")

        return sig

    def wrap(self, wrapping_key_label: str, wrapped_key_label: str):
        """Export the key with label 'wrapped_key_label' encrypted under key 'wrapping_key_label'"""

        log_tag = f"[{wrapped_key_label}, wrap]: "

        rsa_sk_handel = self._get_rsa_sk_handel(wrapped_key_label)

        wrapping_key_handels = find_existing_rsa_key(self.session, wrapping_key_label)
        if len(wrapping_key_handels) == 0:
            self.logging.error(f"{log_tag}Cannot find wrapping key with label '{wrapping_key_label}'")
        elif len(wrapping_key_handels) > 1:
            self.logging.warning(f"{log_tag}Multiple keys with wrapping key label '{wrapping_key_label}' exist, using first one.")

        wrapping_key_handel = wrapping_key_handels[0]

        self.logging.info(f"{log_tag}wrap key with label {wrapped_key_label} using wrapping key with label {wrapping_key_label}")
        wrapped_key = wrap_key(self.session, rsa_sk_handel, wrapping_key_handel)
        self.logging.info(f"{log_tag}Wrapped key {wrapped_key.hex()}")

        return wrapped_key

if __name__ == "__main__":
    ## Command-line arguments
    parser = argparse.ArgumentParser(description='Luna K6 HSM signature oracle server')

    parser.add_argument('--reuse_paths', action='store_true',
        help='Append to existing log files instead of creating a new subdir in <logs_dir>')
    parser.add_argument('--logs_dir', type=str, default='logs',
        help='Path to directory in which log files will be stored.')
    parser.add_argument('--verbose', action='store_true',
        help='Enable extra debug output')
    parser.add_argument('--host', default='0.0.0.0',
        help='IP address of the interface that the socket should listen on')
    parser.add_argument('--port', type=int, default=61338,
        help='Port that the socket should listen on')
    parser.add_argument('--userpin', type=str,
        help='Secret value generated on the Luna K6 HSM to allow user authentication without physically plugging in the PED.')

    args = parser.parse_args()

    do_reuse_paths      = args.reuse_paths
    logs_dir            = args.logs_dir
    host                = args.host
    port                = args.port
    userpin             = args.userpin

    # The HSM can only handle requests from one thread at a time
    exec_threads = 1

    try:
        os.mkdir(logs_dir)
    except:
        pass

    if not do_reuse_paths:
        subdir_name = 'run_' + time.strftime('%Y%m%d-%H%M%S')
        logs_dir += f'/{subdir_name}'
        os.mkdir(logs_dir)

    log_fpath = f"{logs_dir}/luna_k6_hsm_oracle_server.log"
    log_level = logging.DEBUG if args.verbose else logging.WARNING

    logging.basicConfig(filename=log_fpath, level=log_level)

    print(f"Luna K6 HSM server logging to {log_fpath}")
    oracle = LunaK6HSMOracleServer(logs_dir, logging, userpin, host, port, exec_threads)
    oracle.run()
