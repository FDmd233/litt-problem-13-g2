# Exact certificate

The proof uses one exact symbolic certificate:

```text
python code/g2_shifted_norm_certificate.py
```

The script starts from the shifted cyclic cover $w^3=y-c$ and the three sections of $H^0(C,3\kappa-D)$ used in the paper. It reconstructs the cyclic norm and residual conic, then verifies:

- `det(M_s) = (u-v) F_8`, with homogeneous residual octic `F_8` and structural factor of multiplicity one;
- the normalized identity `F = A^2 + 3 Delta B^2`, which gives the cyclotomic factor-ordering cover;
- a characteristic-13 geometric-integrality witness for the residual octic;
- a characteristic-13 admissibility witness for the two residual line factors.

All symbolic arithmetic is exact. The reference log is stored in [`g2_shifted_norm_certificate.out`](g2_shifted_norm_certificate.out). The finite-field calculations certify nonempty open conditions; the passage to characteristic zero is the geometric argument in the paper.
