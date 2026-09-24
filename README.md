# Forging 1024-bit RSA Signatures in Nearly SNFS Time

*Alternative title: Nearly SNFS-Speed Signature Forgery Sans Factoring N (NSNFSSSFSFN)*

This repository accompanies the [paper](https://eprint.iacr.org/2026/2131.pdf)
(IACR ePrint 2026/2131) and provides an implementation of a variant of the
number field sieve (NFS). The work demonstrates that an adversary with
*temporary* access to a raw, unpadded RSA signing or decryption oracle can
obtain the *permanent* ability to forge signatures and decrypt ciphertexts.
In effect, the adversary acquires capabilities equivalent to possession of the
private key — forging signatures and decrypting offline — without factoring
the public modulus, and at substantially lower computational cost than
factoring would require.

These results indicate that factoring-based estimates of RSA security may be
overly optimistic and warrant revision. They do not, however, constitute an
immediate operational threat to most deployed RSA.

## Overview

- **The algorithm is not polynomial-time.** It is subexponential — the same
  complexity class as the best known factoring algorithms — but attains the
  faster subexponential regime of the *special* number field sieve (SNFS)
  rather than the *general* number field sieve (GNFS). Factoring a 1024-bit
  modulus is estimated to require 500,000–1,000,000 core-years; executing this
  algorithm on a 1024-bit modulus required approximately **1,380 core-years**.

- **The algorithm is not new.** It was introduced in 2007 by Joux, Naccache,
  and Thomé [[1]](https://eprint.iacr.org/2007/424). This work presents, to our
  knowledge, its first public implementation and large-scale execution. The
  implementation builds substantially on
  [CADO-NFS](https://gitlab.inria.fr/cado-nfs/cado-nfs) [[7]](#references) and
  follows the standard NFS
  pipeline — polynomial selection, relation collection by lattice sieving
  [[3]](https://arxiv.org/abs/2001.10860), linear algebra, and an e-th-root
  step [[6]](https://arxiv.org/abs/2305.17425) generalizing the square-root
  computation of ordinary factoring.

- **The attack requires a raw signing oracle.** Standard RSA signature schemes
  employing PKCS #1 v1.5 or RSA-PSS padding do not expose such an oracle, and
  are therefore not affected. Deployments that may expose a raw oracle include
  blind RSA signatures (e.g. Privacy Pass) and certain HSM APIs.

## General FAQ

**1. Is RSA still in widespread use?**
Yes, particularly for digital signatures (certificates, TLS handshakes,
tokens, OAuth). Key exchange in protocols such as TLS uses ECDH or has
transitioned to ML-KEM; the attack does not apply to these algorithms.

**2. Is there an immediate need to discontinue RSA?**
No. The attack does not present a practical risk for most uses of RSA. Even in
potentially vulnerable scenarios, we estimate the attack cost for 2048-bit RSA
at $2^{90}$. This is below the estimated $2^{112}$ cost of factoring a
2048-bit modulus, yet remains roughly 1000 times the estimated $2^{80}$ cost
of factoring a 1024-bit modulus — and no 1024-bit RSA modulus has been
publicly factored to date. For additional assurance, elliptic-curve schemes
(e.g. ECDSA or Ed25519) do not appear susceptible to this class of attack.

**3. Is there a longer-term case for discontinuing RSA?**
In our view, yes. The attack indicates that RSA key sizes up to 4096 bits do
not meet contemporary security levels. The ongoing migration to post-quantum
cryptography presents an opportunity to retire legacy primitives such as RSA.

**4. What are the implications for 2048-bit blind RSA (e.g. Privacy Pass)?**
The attack applies. We estimate it requires $2^{90}$ work and $2^{43}$ oracle
queries — likely below the intended security margin, though feasible only for
the most capable adversaries. Short-term mitigations include shortening key
epochs and increasing key lengths. In the medium term, protocol designers may
eliminate this vector by incorporating a zero-knowledge proof. In the long
term, we anticipate acceptable post-quantum replacements for constructions
such as blind signatures.

**5. Can this code be compiled to forge 1024-bit RSA signatures directly?**
Only with access to significant computational resources. On our academic-scale
CPU cluster the forgery computation required several months; the total effort
was between the RSA-240 and RSA-250 factoring records
[[2]](https://arxiv.org/abs/2006.06197).

**6. Could AI or GPUs accelerate this implementation?**
Almost certainly. See [`docs/las_gpu_notes.md`](docs/las_gpu_notes.md) and
[`code/gpu/`](code/gpu/) for a first step: a validated OpenCL implementation of
the ECM cofactorization stage, the dominant cost within the siever. Prior work
on ECM and modular arithmetic on GPUs is directly relevant
[[4]](https://arxiv.org/abs/1310.3809)[[5]](https://arxiv.org/abs/2501.03245).

**7. Did you use AI or GPUs?**
No.

**8. I use 2048-bit keys with PKCS or PSS padding for signatures — should I be
concerned?**
No. The attack does not appear feasible in this setting.

**9. How does this differ from the textbook RSA forgery attack?**
This attack is strictly more powerful: temporary signing access is converted
into capabilities equivalent to knowledge of the private key — the ability to
forge arbitrary signatures offline at a later time. The bulk of the cost is a
single precomputation depending only on the public key; subsequent individual
forgeries are comparatively inexpensive.

**10. Why, in September 2026, is an algorithm from 2007 only now being
demonstrated?**
The 1024-bit computation required substantial effort (documented in detail in
the paper). Relatively few researchers work in this area, funding is limited,
and the field increasingly draws students toward post-quantum cryptography
and, more recently, AI. Nevertheless, such computational work is essential to
move these algorithms beyond the theoretical.

**11. What does "SNFS time" mean?**
It refers to the *special* number field sieve — also the running time for
factoring special-form integers such as Mersenne numbers, which is faster than
for generic integers such as well-formed RSA moduli. See the
[Wikipedia article](https://en.wikipedia.org/wiki/Special_number_field_sieve)
for background.

## Special FAQ

**12. I use 1024-bit RSA with modulus
`119761307924183143227323805033445425142683607945593313835302891299066516303843000359184032065280614739228104709307215190162897575231548821264285875476222037118119400262673724895150267064851629266552653543645482302630040586124860255443008330625425740446416414144702538281275665910057429738414800482721215922451`.
What should I do?**
Change your key.

**13. When was the computation performed?**
The 1024-bit computation took several months and concluded on 31 August 2026.

## Repository editions

`NSNFSSSFSFN.bundle.sh` packages the entire repository — code, configurations,
patches, build and run pipelines, CI workflow, data, logs, and the paper — as a
single self-extracting, checksum-verified file:

```bash
bash NSNFSSSFSFN.bundle.sh --list                  # list contents
bash NSNFSSSFSFN.bundle.sh -C NSNFSSSFSFN           # extract
bash NSNFSSSFSFN.bundle.sh -C NSNFSSSFSFN --with-submodule   # also fetch CADO-NFS at the pinned commit
```

Regenerate it after any change with `python3 tools/make_single_file.py`; the
`single-file bundle` GitHub workflow fails if it becomes stale.

`nsnfsssfsfn.py` packages the code alone (`code/` and `patches/`) as a single
Python file, with the embedded sources retained as readable comment lines:

```bash
python3 nsnfsssfsfn.py list                      # or: cat code/run.py
python3 nsnfsssfsfn.py unpack mydir --with-cado  # tree + CADO-NFS at the pinned commit, patched
sage nsnfsssfsfn.py run -l locations.config config/n192.config precomp   # runs code/run.py in ./nsnfsssfsfn
sage nsnfsssfsfn.py exec oracles/sage_oracle.py --help                   # run a tool from memory
```

Regenerate with `python3 tools/make_single_python.py`.

## References

1. A. Joux, D. Naccache, and E. Thomé. *When e-th Roots Become Easier Than
   Factoring.* ASIACRYPT 2007. IACR ePrint 2007/424.
   <https://eprint.iacr.org/2007/424>
2. F. Boudot, P. Gaudry, A. Guillevic, N. Heninger, E. Thomé, and P. Zimmermann.
   *Comparing the Difficulty of Factorization and Discrete Logarithm: A
   240-Digit Experiment.* CRYPTO 2020. arXiv:2006.06197.
   <https://arxiv.org/abs/2006.06197>
3. G. McGuire and O. Robinson. *A New Angle on Lattice Sieving for the Number
   Field Sieve.* 2020. arXiv:2001.10860. <https://arxiv.org/abs/2001.10860>
4. W. Trei. *Efficient Modular Arithmetic for SIMD Devices.* 2013.
   arXiv:1310.3809. <https://arxiv.org/abs/1310.3809>
5. Q. Xiong et al. *gECC: A GPU-based High-throughput Framework for Elliptic
   Curve Cryptography.* 2024. arXiv:2501.03245.
   <https://arxiv.org/abs/2501.03245>
6. O. Bernard, P.-A. Fouque, and A. Lesavourey. *Computing e-th Roots in Number
   Fields.* 2023. arXiv:2305.17425. <https://arxiv.org/abs/2305.17425>
7. The CADO-NFS development team. *CADO-NFS: An Implementation of the Number
   Field Sieve Algorithm.* <https://gitlab.inria.fr/cado-nfs/cado-nfs>
