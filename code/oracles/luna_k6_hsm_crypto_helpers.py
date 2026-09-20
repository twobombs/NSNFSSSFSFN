from pycryptoki.default_templates import *
from pycryptoki.defines import *
from pycryptoki.key_generator import *
from pycryptoki.session_management import *
from pycryptoki.object_attr_lookup import *
from pycryptoki.sign_verify import *
from pycryptoki.encryption import *
from pycryptoki.mechanism import *
from pycryptoki.conversions import *

RAW_RSA_MECHANISM   = Mechanism(mech_type=CKM_RSA_X_509)

def pad(data, modulus_bits):
    nbytes = modulus_bits // 8
    return data.rjust(nbytes, b"\x00")

## Wrapper functions for HSM operations
def find_existing_rsa_key(h_session, label):
    template = {CKA_LABEL: label}
    return c_find_objects_ex(h_session, template, 1)

def gen_rsa_keys(h_session, pk_label, sk_label, public_exp, modulus_bits):
    """
    Generate a new RSA key pair on the HSM partition
    """

    pub_template = {
        CKA_TOKEN: True,
        CKA_PRIVATE: True,
        CKA_MODIFIABLE: True,
        CKA_ENCRYPT: True,
        CKA_VERIFY: True,
        CKA_WRAP: True,
        CKA_MODULUS_BITS: modulus_bits,
        CKA_PUBLIC_EXPONENT: public_exp,
        CKA_LABEL: pk_label
    }

    priv_template = {
        CKA_TOKEN: True,
        CKA_PRIVATE: True,
        CKA_SENSITIVE: True,
        CKA_MODIFIABLE: True,
        CKA_EXTRACTABLE: True,
        CKA_DECRYPT: True,
        CKA_SIGN: True,
        CKA_UNWRAP: True,
        CKA_LABEL: sk_label
    }

    pub_key, priv_key = c_generate_key_pair_ex(
                            h_session,
                            mechanism=CKM_RSA_X9_31_KEY_PAIR_GEN, # That's what the pycryptoki test cases do too
                            pbkey_template=pub_template,
                            prkey_template=priv_template)

    return pub_key, priv_key

def get_rsa_key(h_session, pk_label, sk_label, public_exp, modulus_bits):
    h_rsa_sks = find_existing_rsa_key(h_session, sk_label)
    h_rsa_pks = find_existing_rsa_key(h_session, pk_label)

    if len(h_rsa_sks) == 0 and len(h_rsa_pks) == 0:
        print("Generate new RSA key pair")
        h_rsa_pk, h_rsa_sk = gen_rsa_keys(h_session, pk_label, sk_label, public_exp, modulus_bits)
    elif len(h_rsa_sks) == 1 and len(h_rsa_pks) == 1:
        print("Use existing RSA key pair")
        h_rsa_pk, h_rsa_sk = h_rsa_pks[0], h_rsa_sks[0]
    else:
        raise Exception("Invalid state, only one key from the specified RSA key pair exists.")

    return h_rsa_pk, h_rsa_sk

def enc(h_session, h_rsa_pk, ptxt):
    return c_encrypt_ex(h_session, h_rsa_pk, ptxt, mechanism=RAW_RSA_MECHANISM)

def dec(h_session, h_rsa_sk, ctxt):
    return c_decrypt_ex(h_session, h_rsa_sk, ctxt, mechanism=RAW_RSA_MECHANISM)

def sign(h_session, h_rsa_sk, data):
    return dec(h_session, h_rsa_sk, data)

def verify(h_session, h_rsa_pk, data, sig):
    return data == enc(h_session, h_rsa_pk, sig)

def get_pub_key_info(h_session, h_rsa_pk):
    pub_key_info_template = {
        CKA_MODULUS: None,
        LUNA_ATTR_PUBLIC_EXPONENT: None
    }

    infos = c_get_attribute_value_ex(h_session, h_rsa_pk, pub_key_info_template)

    N = int.from_bytes(bytes.fromhex(infos[LUNA_ATTR_MODULUS].decode()), byteorder="big")
    e = int.from_bytes(bytes.fromhex(infos[LUNA_ATTR_PUBLIC_EXPONENT].decode()), byteorder="big")

    return N, e

def wrap_key(h_session, h_wrapped_key, h_wrapping_key, mechanism=Mechanism(mech_type=CKM_RSA_PKCS)):
    return c_wrap_key_ex(h_session, h_wrapping_key, h_wrapped_key, mechanism)

def test_enc_dec(h_session, h_rsa_pk, h_rsa_sk):
    data = b"deadbeef" + b"\x00"*120
    assert(len(data) == 128)

    ctxt = enc(h_session, h_rsa_pk, data)
    ptxt = dec(h_session, h_rsa_sk, ctxt)

    print("data", data)
    print("ctxt", ctxt)
    print("ptxt", ptxt)

    assert(ptxt == data)

def test_sign_verify(h_session, h_rsa_pk, h_rsa_sk):
    data = b"deadbeef" + b"\x00"*120
    assert(len(data) == 128)

    sig = sign(h_session, h_rsa_sk, data)
    check = verify(h_session, h_rsa_pk, data, sig)

    print("sig", sig)
    print("data", data)
    print("check", check)

    assert check
