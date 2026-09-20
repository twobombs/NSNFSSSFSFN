"""
This file contains low-level HSM operations using the pycryptoki
Python package published by Thales. The functions are used by the
LunaHsmApi class to perform the operations on the HSM.
"""


from pycryptoki.default_templates import *
from pycryptoki.defines import *
from pycryptoki.key_generator import *
from pycryptoki.session_management import *
from pycryptoki.object_attr_lookup import *
from pycryptoki.sign_verify import *
from pycryptoki.encryption import *
from pycryptoki.mechanism import *
from pycryptoki.conversions import *

RAW_RSA_MECHANISM = Mechanism(mech_type=CKM_RSA_X_509)

# Must only be called once
c_initialize_ex(CKF_OS_LOCKING_OK)

def pad(data, modulus_bits):
    nbytes = modulus_bits // 8
    return data.rjust(nbytes, b"\x00")

## Wrapper functions for HSM operations
def find_existing_rsa_key(h_session, label):
    template = {CKA_LABEL: label}
    return c_find_objects_ex(h_session, template, 1)

def get_rsa_key(h_session, pk_label, sk_label):
    h_rsa_sks = find_existing_rsa_key(h_session, sk_label)
    h_rsa_pks = find_existing_rsa_key(h_session, pk_label)

    if len(h_rsa_sks) == 1 and len(h_rsa_pks) == 1:
        h_rsa_pk, h_rsa_sk = h_rsa_pks[0], h_rsa_sks[0]
    else:
        err_msg = "Invalid state."
        if len(h_rsa_pks) != 1:
            err_msg += f" Found {len(h_rsa_pks)} public key(s) for label {pk_label}."
        if len(h_rsa_sks) != 1:
            err_msg += f" Found {len(h_rsa_sks)} public key(s) for label {sk_label}."
        raise Exception(err_msg)

    return h_rsa_pk, h_rsa_sk

def start_session(slot, userpin):
    h_session = c_open_session_ex(slot)
    login_ex(h_session, slot, userpin)
    return h_session

def dec(h_session, h_rsa_sk, ctxt):
    return c_decrypt_ex(h_session, h_rsa_sk, ctxt, mechanism=RAW_RSA_MECHANISM)

def sign(h_session, h_rsa_sk, data):
    return dec(h_session, h_rsa_sk, data)

def get_pub_key_info(h_session, h_rsa_pk):
    pub_key_info_template = {
        CKA_MODULUS: None,
        LUNA_ATTR_PUBLIC_EXPONENT: None
    }

    infos = c_get_attribute_value_ex(h_session, h_rsa_pk, pub_key_info_template)

    N = int.from_bytes(bytes.fromhex(infos[LUNA_ATTR_MODULUS].decode()), byteorder="big")
    e = int.from_bytes(bytes.fromhex(infos[LUNA_ATTR_PUBLIC_EXPONENT].decode()), byteorder="big")

    return N, e