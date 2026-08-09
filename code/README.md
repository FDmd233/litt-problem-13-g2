# Optional exact-arithmetic check

The manuscript now states the polynomial identities and finite-field witnesses used in the pair-incidence argument explicitly. The script in this directory repeats those calculations with exact SymPy arithmetic; it is intended as an independent check, not as a logical input to the proof.

Run:

```text
python code/g2_shifted_norm_certificate.py
```

The script reconstructs the cyclic norm and residual conic from the shifted cover $w^3=y-c$ and verifies:

- `det(M_s) = (u-v) F_8`, with homogeneous residual octic `F_8` and structural factor of multiplicity one;
- the identity `F = A^2 + 3 Delta B^2` for the ordering cover;
- the characteristic-13 irreducibility witness for the residual octic;
- the characteristic-13 admissibility witness for the two residual line factors.

The reference output is stored in [`g2_shifted_norm_certificate.out`](g2_shifted_norm_certificate.out). The passage from the finite-field witnesses to the generic characteristic-zero statements is proved in the manuscript.
