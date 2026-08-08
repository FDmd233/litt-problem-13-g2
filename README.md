# A Realization of G2 for Litt Problem No. 13

This repository contains the revised proof of the $G_2$ case of Litt Problem No. 13. The construction uses the shifted cyclic triple cover

\[
w^3=y-c,\qquad c\ne0,
\]

of an elliptic curve $y^2=x^3+Ax+B$.

> **Revision note (2026-08-09).** The July preprint specialized to $c=0$, where an additional involution leaves a compatibility issue in the surface argument. This revision uses $c\ne0$ and rederives the norm-conic identity and finite-field witnesses. The $G_2$ realization is retained; the fixed-CM/isotrivial Prym claim is withdrawn.

## Paper

- [PDF](output/pdf/Litt_Problem_13_G2.pdf)
- [LaTeX source](paper/main.tex)

## Exact certificate

The calculation used in the pair-incidence argument is reconstructed from the shifted curve in a single exact SymPy script:

- [`code/g2_shifted_norm_certificate.py`](code/g2_shifted_norm_certificate.py)
- [`code/g2_shifted_norm_certificate.out`](code/g2_shifted_norm_certificate.out)

## Publication record

- [Revised release: v1.1.0-revised](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.1.0-revised)
- [First release: v1.0.0-preprint](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.0.0-preprint) — superseded by the revised construction above.
