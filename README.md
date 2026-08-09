# A Realization of G2 for Litt Problem No. 13

This repository contains the corrected manuscript for the $G_2$ case of Litt Problem No. 13. The construction uses the shifted cyclic triple cover

\[
w^3=y-c,\qquad c\ne0,
\]

of an elliptic curve $y^2=x^3+Ax+B$.

> **Correction note (2026-08-09).** In addition to the earlier correction excluding the special case $c=0$, this version corrects the first-cohomology calculation at the simple elliptic singularity. The minimal resolution has $p_g(S)=5$ and $q(S)=4$. The generic-twist argument is rewritten accordingly; the Hodge multiplicities $(2,3,2)$ and the $G_2$ conclusion are unchanged.

## Paper

- [PDF](output/pdf/Litt_Problem_13_G2.pdf)
- [LaTeX source](paper/main.tex)

## Verification

The algebra used in the pair-incidence argument is written explicitly in the manuscript. An exact SymPy script is kept as an optional independent check:

- [`code/g2_shifted_norm_certificate.py`](code/g2_shifted_norm_certificate.py)
- [`code/g2_shifted_norm_certificate.out`](code/g2_shifted_norm_certificate.out)

The script is not a logical input to the proof.

## Publication record

- [Corrected release: v1.2.0-corrected](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.2.0-corrected)
- [Revised release: v1.1.0-revised](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.1.0-revised) — superseded by v1.2.0-corrected.
- [First release: v1.0.0-preprint](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.0.0-preprint) — retained for the record.
