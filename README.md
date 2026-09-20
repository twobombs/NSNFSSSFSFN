## Forging 1024-bit RSA signatures in nearly SNFS time

*Alternative title: Nearly SNFS-Speed Signature Forgery Sans Factoring N (NSNFSSSFSFN)*

This repository contains a paper (eprint 2026/XXXX) and code implementing a variant of the number field sieve (NFS) algorithm.
It shows that an attacker can use *temporary* access to a raw, unpadded RSA signing/decryption oracle to gain the *permanent* ability to forge signatures / decrypt ciphertexts. In other words, the attacker can steal what is effectively the secret key (in that it can be used to sign/decrypt offline), but without actually factoring the public key, and using much less computation than factoring the public key would have taken.
This demonstrates that factoring-based estimates for RSA security may be too optimistic and should be revised, but likely does not pose an immediate operational threat to most deployed RSA in the real world.

**The algorithm is not polynomial-time.** It's not even close. It's "subexponential-time" which is the same class as the best factoring algorithms. However it manages to be a faster subexponential-time: "special" number field sieve rather than "general" number field sieve.
Factoring a 1024-bit modulus is predicted to take 500,000-1,000,000 core-years. Running this algorithm on a 1024-bit modulus took us 1,380 core-years.

**The algorithm is not new.** It was invented in 2007 by [Joux, Naccache, and Thomé](https://eprint.iacr.org/2007/424). However this is the first public implementation and large-scale run.
Most of the code is not new; it builds on [CADO-NFS](https://gitlab.inria.fr/cado-nfs/cado-nfs).

**The algorithm only works if a raw signing oracle is available.** Most RSA usage in practice (that is, RSA signatures using PKCS\#1v1.5 or RSA-PSS padding) do not expose such an oracle, and thus this attack does not pose a practical risk. Examples of RSA use that do expose such a signing oracle would include blind RSA signatures (e.g. Privacy Pass) or HSM APIs.

### General FAQ

1. **Are people actually still using RSA?** \
Yes, particularly for digital signatures (e.g. certificates, TLS handshakes, tokens, OAuth). Key exchange for protocols like TLS uses ECDH or has transitioned to ML-KEM, and the attack does not apply to these algorithms.
2. **Is there an immediate need to stop using RSA?** \
No. The attack does not pose a practical risk for most uses of RSA. Even in the scenarios that might be vulnerable, we estimate that the attack cost for 2048-bit RSA is $2^{90}$.  This is below the estimated $2^{112}$ cost to factor a 2048-bit RSA key, but is 1000 times more work than the estimated $2^{80}$ cost to factor a 1024-bit RSA key, and nobody has factored a 1024-bit RSA key in public yet. If you're still worried, elliptic curve cryptography (e.g. ECDSA or Ed25519 for signatures) does not appear to be vulnerable to this type of attack.
3. **Is there a non-immediate need to stop using RSA?** \
In our opinion, yes. This attack shows that RSA key sizes up to 4096 bits do not meet modern cryptographic security levels. The ongoing transition to post-quantum cryptography provides an opportunity to move away from legacy cryptography like RSA altogether.
4. **What if I use 2048-bit keys for blind RSA (e.g. Privacy Pass)?** \
The attack does apply. We estimate it would take $2^{90}$ work and $2^{43}$ oracle queries, which is likely below your desired security margin while being infeasible for all but the largest, most sophisticated adversaries to exploit. In the short term, you may consider shortening key epochs and increasing RSA key lengths.  In the medium term, protocol designers could eliminate this attack vector by adding a zero-knowledge proof.  In the long term, we hope that there will be acceptable post-quantum replacements for constructions like blind signatures.
5. **Can I compile this code and immediately forge 1024-bit RSA signatures?** \
Only if you have access to significant computational resources. With our academic-sized CPU cluster it took us months to do the forgery computation.  The amount of computation was somewhere between the RSA-240 and RSA-250 records.
6. **Can an AI/GPUs speed up this implementation?** \
Almost certainly yes.
7. **Did you use AI/GPUs?** \
No.
8. **I use 2048-bit keys with PKCS or PSS padding for signatures. Should I be worried?** \
No, our attack does not seem to be feasible in this case.
9. **Why isn't this just the textbook RSA attack?** \
This attack is more powerful, because the attacker can use the temporary ability to request signatures to steal abilities equivalent to having the secret key, without actually having the secret key. (That is, being able to forge any signature of their choice, offline, in the future.)  Most of the work is a single precomputation that depends only on the public key; after this is finished individual signature forgeries are much more efficient to compute.
10. **Why, in September 2026, am I hearing about an algorithm from 2007?** \
The 1024-bit computation we did took significant effort (see our paper for the hairy details). Relatively few researchers are working in this area, there is little funding, and the job market is more interested in students working in post-quantum cryptography (and now, AI).  However, without this kind of computational work the algorithms would remain theoretical.
11. **What does `SNFS time` mean?** \
"Special" number field sieve. This is also the running time to factor special-form integers like Mersenne numbers, which can be done more quickly than generic integers like well-formed RSA moduli. Read about it on [Wikipedia](https://en.wikipedia.org/wiki/Special_number_field_sieve) if you're curious.

### Special FAQ

12. **I use 1024-bit RSA and my modulus is `119761307924183143227323805033445425142683607945593313835302891299066516303843000359184032065280614739228104709307215190162897575231548821264285875476222037118119400262673724895150267064851629266552653543645482302630040586124860255443008330625425740446416414144702538281275665910057429738414800482721215922451`. What should I do?** \
Change your key.
13. **When did you do this computation?** \
The 1024-bit size took us several months. It finished on August 31, 2026.
