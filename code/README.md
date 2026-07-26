# Exact certificates

All Python calculations use exact SymPy 1.14 arithmetic.  From the repository
root, run

```text
python code/g2_mu6_source_provenance_verify.py
python code/g2_mu6_lowmem_certificate_fast.py
python code/g2_mu6_direct_numerator_verify.py
```

The provenance verifier starts from the normalized elliptic cubic and the
three sections used in the paper.  It computes their cyclic cubic norm,
divides by the fixed line `y-x`, constructs the residual conic determinant
and ordering discriminant, and then verifies both saved Singular inputs
formula by formula and as complete files.  Thus the large saved polynomials
are derived inputs, not unverified fixtures.

These three commands are read-only by default.  The explicit maintenance
options

```text
python code/g2_mu6_source_provenance_verify.py --write-inputs
python code/g2_mu6_lowmem_certificate_fast.py --write-formulas
```

replace generated files only after the exact calculations have completed.
The expected logs are in `code/generated/` and use repository-relative paths.

With Singular installed at `/usr/bin/Singular` in WSL, run

```text
powershell -NoProfile -ExecutionPolicy Bypass -File code/g2_mu6_universal_special123_mod13_audit.ps1
```

The wrapper fails unless the reduction has exactly one nonunit factor of
degree eight and exponent one, the stated point lies on it, and the gradient
is nonzero modulo 13.  Irreducibility over the finite field together with
that smooth rational point certifies geometric integrality of this specified
fiber.  Passage to the characteristic-zero generic family additionally uses
the flatness and openness arguments in the paper.

`g2_mu6_nonsimple_norm_conic_Q.py` regenerates the separate exact sample and
its Singular inputs.  The saved modular log records the finite-field witness
for the nonempty open conditions; it is not used as a numerical proof of a
generic identity.
