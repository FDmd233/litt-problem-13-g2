# Exact certificates

All Python calculations use exact SymPy 1.14 arithmetic.  From the repository
root, run

```text
python code/g2_mu6_source_provenance_verify.py
python code/g2_mu6_lowmem_certificate_fast.py
python code/g2_mu6_direct_numerator_verify.py
python code/g2_mu6_nonsimple_norm_conic_Q.py
```

The provenance verifier starts from the normalized elliptic cubic and the
three sections used in the paper.  It computes their cyclic cubic norm,
divides by the fixed line `y-x`, constructs the residual conic determinant
and ordering discriminant, verifies that the structural factor has
multiplicity one and the residual octic is homogeneous, and checks a Rabin
irreducibility certificate for the finite-field slice.  It then verifies both
saved Singular inputs formula by formula and as complete files.  Thus the
large saved polynomials are derived inputs, not unverified fixtures.

`g2_mu6_nonsimple_norm_conic_Q.py` uses an automatically generated kernel
basis.  It verifies the exact change of basis before identifying its internal
point `(1:0:9)` with the paper's universal-basis point `(10:9:0)`; changing
one coordinate triple without the other is invalid.  Its default Python run
also recomputes and asserts, using exact SymPy arithmetic, the elliptic-curve
discriminant, fixed nonzero minors for the two evaluation ranks, the six conic
coefficients, the vanishing determinant and nonzero rank-two minor, the six
fixed-point line values, both line-section derivative resultants, both branch
resultants, and the pair resultant modulo 13.  These printed residues are the
second finite-field witness listed in the paper's appendix and do not require
Singular to verify.

These four commands are read-only by default.  The explicit maintenance
options

```text
python code/g2_mu6_source_provenance_verify.py --write-inputs
python code/g2_mu6_lowmem_certificate_fast.py --write-formulas
python code/g2_mu6_nonsimple_norm_conic_Q.py --write-inputs
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

With `--write-inputs`, `g2_mu6_nonsimple_norm_conic_Q.py` regenerates the
separate exact sample and its Singular inputs.  Those inputs and the saved
modular log provide a separate Singular rendering of the finite-field
witness; the default Python run performs the underlying assertions directly.
The witness proves nonemptiness of the specified open conditions, not a
generic identity.
