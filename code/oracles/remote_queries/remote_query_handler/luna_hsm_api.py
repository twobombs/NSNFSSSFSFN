import argparse
import logging

from logging import Logger

from luna_hsm_operations import *

class LunaHsmApi:
    """
    This class provides an interface to HSM operations (in our case, only signing.)
    """

    def __init__(self, logging: Logger, slot: int, userpin: str) -> None:
        self.logging    = logging
        self.slot       = slot
        self.userpin    = userpin

        self.key_info_cache: dict[str, tuple] = dict()
        self.sk_handel_cache: dict[str, int]  = dict()

        self.logging.info("Starting Luna HSM session...")

        self.session = start_session(self.slot, self.userpin)
        self.logging.info("Luna HSM session started successfully")

    def _get_pk_label(self, label: str):
        return f"{label}_public_key"

    def _get_sk_label(self, label: str):
        return f"{label}_secret_key"

    def _get_key_info(self, label, handel=None):
        log_tag = "[_get_key_info] "
        if label not in self.key_info_cache:
            if handel == None:
                handels = find_existing_rsa_key(self.session, label)
                if len(handels) != 1:
                    self.logging.error(f"{log_tag}Couldn't find key with label {label}, got handels: {handels}.")
                    raise Exception(f"Couldn't find unique result for key with label {label}.")

                handel = handels[0]

            N, e = get_pub_key_info(self.session, handel)
            modulus_bits = ((N.bit_length() + 7) // 8) * 8

            self.key_info_cache[label] = (N, e, modulus_bits)
            self.logging.info(f"{log_tag}Retrieved following info for public key with label {label} (cached for future use): N={N}, e={e}, modulus_bits={modulus_bits}")

        N, e, modulus_bits = self.key_info_cache[label]
        self.logging.debug(f"{log_tag}Retrieved following info for public key with label {label}: N={N}, e={e}, modulus_bits={modulus_bits}")
        return N, e, modulus_bits

    def _get_rsa_sk_handel(self, rsa_sk_label, rsa_sk_handel=None):
        log_tag = "[_get_rsa_sk_handel] "
        self.logging.debug(f"{log_tag}Retrieve secret key handel for label {rsa_sk_label}")

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

        self.logging.debug(f"{log_tag}Retrieved handel {self.sk_handel_cache[rsa_sk_label]}")
        return self.sk_handel_cache[rsa_sk_label]

    def sign(self, data: bytes, label: str):
        """Sign 'data' using the RSA key specified by 'label'"""

        log_tag = f"[{label}, sign]: "

        _, _, modulus_bits = self._get_key_info(label)
        data_pad = pad(data, modulus_bits)

        rsa_sk_handel = self._get_rsa_sk_handel(label)

        self.logging.debug(f"{log_tag}sign the following data with key with label {label}: {data_pad.hex()}")
        sig = sign(self.session, rsa_sk_handel, data_pad)
        self.logging.debug(f"{log_tag}Produced signature {sig.hex()}")

        return sig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Luna HSM API')

    parser.add_argument("label", type=str,
        help='Label of the key that should be used for signing.')
    parser.add_argument("data", type=str,
        help='Data (in hex) to sign.')

    parser.add_argument('--slot', type=int,
        help='Slot in which the HSM is installed.')
    parser.add_argument('--userpin', action='store_true',
        help='Read the secret value from stdin that is generated on the Luna HSM to allow user authentication without physically plugging in the PED.')

    parser.add_argument('--logs_dir', type=str, default='logs',
        help='Path to directory in which log files will be stored. (Must already exist.)')
    parser.add_argument('--verbose', action='store_true',
        help='Enable extra debug output')

    args = parser.parse_args()

    label           = args.label
    data            = bytes.fromhex(args.data)
    logs_dir        = args.logs_dir
    slot            = args.slot
    do_read_userpin = args.userpin

    if do_read_userpin:
        userpin = input(f"Please enter the userpin for your Luna HSM on slot {slot}: ")

    log_fpath = f"{logs_dir}/luna_hsm_api.log"
    log_level = logging.DEBUG if args.verbose else logging.WARNING

    logging.basicConfig(filename=log_fpath, level=log_level)

    print(f"Luna HSM API logging to {log_fpath}")
    oracle = LunaHsmApi(logging.getLogger(), slot, userpin)

    print(f"Signing the following data with key with label {label}: {data.hex()}")
    sig = oracle.sign(data, label)
    print(f"Produced signature: {sig.hex()}")