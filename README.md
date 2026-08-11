# A Realization of G2 for Litt Problem No. 13

This repository contains the manuscript for the $G_2$ realization in Litt Problem No. 13. The construction starts from an elliptic curve and a shifted cyclic triple cover

$$
E:\ y^2=x^3+Ax+B,
\qquad
C:\ w^3=y-c,
\quad c\neq 0.
$$

From this cover one obtains a polarized Prym threefold and a theta surface $X$. For general parameters, the intersection cohomology object $\mathrm{IC}_X$ has Euler characteristic $7$, and its full convolution Tannaka group is $G_2$ acting through the seven-dimensional standard representation.

## Manuscript

- [PDF](output/pdf/Litt_Problem_13_G2.pdf)
- [Complete LaTeX source](Litt_Problem_13_G2_complete.tex)
- [Split LaTeX source](paper/main.tex)

## Verification in the manuscript

The pair-incidence splitting is proved intrinsically in Section 4 from the cyclic cubic norm. Appendix A records the remaining exact characteristic-zero and finite-field checks used for the residual octic and the admissible specialization. The proof does not rely on auxiliary verification code.

## Build

With a standard LaTeX installation including `latexmk`, run:

```bash
make
```

The resulting PDF is written to `output/pdf/Litt_Problem_13_G2.pdf`, and its SHA256 checksum to `SHA256SUMS.txt`.

## Release

Current release: [`v1.3.0`](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.3.0)
