# A Realization of G2 for Litt Problem No. 13

This repository contains the manuscript for the $G_2$ realization in Litt Problem No. 13.

The construction starts from an elliptic curve together with a shifted cyclic triple cover,

$$
E:\quad y^2=x^3+Ax+B,
\qquad
C:\quad w^3=y-c,
\qquad c\ne 0.
$$

For every parameter with $c\ne0$, $4A^3+27B^2\ne0$, and $4A^3+27(B-c^2)^2\ne0$, the cover gives a Prym threefold of polarization type $(1,1,3)$ and a normal theta surface with one simple elliptic singularity. Its intersection complex has Euler characteristic $7$ and full convolution Tannaka pair $(G_2,V_7)$. Section 8 proves the extension from general parameters by a simultaneous resolution and nearby cycles. The example $y^2=x^3+2$, $w^3=y-1$ is defined over $\mathbf Q$.

Manuscript source corrected October 21, 2026 (structural-pencil argument). Formal external verification of the entire proof is pending.

## Manuscript

- [PDF](output/pdf/Litt_Problem_13_G2.pdf)
- [Complete LaTeX source](Litt_Problem_13_G2_complete.tex)
- [Legacy split LaTeX source](paper/main.tex) (not the authoritative revision; do not use for the current proof)

## Exact checks

The pair-incidence splitting is proved intrinsically in Section 4 from the cyclic cubic norm. Appendix A records the remaining exact characteristic-zero and finite-field checks for the residual octic and the admissible specialization, including the reductions, remainders, and Bézout identities needed for independent verification.

The [all-parameter specialization verification note](verification/specialization_criterion.md) isolates the semismall nearby-cycles cancellation argument, explains the local smooth simultaneous resolution, and states which other parts of the proof still require independent review.

- [All-parameter factor-ordering identity](verification/verify_generic_factor_ordering.py) verifies the cyclic norm decomposition and distinct rational residual-line intersections **symbolically**, without fixing a characteristic-13 point. The script runs in [GitHub CI](.github/workflows/verify-residual-octic.yml).
- [Geometric two-branch exhaustion audit](verification/geometric_two_branch_exhaustion.md) checks the rationality of the line factors, scheme-level Abel inverse, boundary component exhaustion, and rational global labels, with its earlier geometric prerequisites explicitly recorded.
- [Independent exact symbolic audit](verification/verify_residual_octic.py) reconstructs the norm quotient, the determinant factorization and the characteristic-13 residual octic, and verifies an explicit factorization into two transverse, disjoint, unbranched rational line sections. This is a finite specialization test, not a proof of the generic factor-ordering split; the [CI job](.github/workflows/verify-residual-octic.yml) reruns it.
- [Antisymmetric convolution proof audit](verification/alternating_convolution_audit.md) separately checks the perverse shifts, local-system triviality and divisor-support extraction conditional on the geometric splitting.

## Build

~~~bash
make
~~~

The resulting PDF is written to `output/pdf/Litt_Problem_13_G2.pdf`. `SHA256SUMS.txt` records the checksums of the PDF and the complete TeX source. The complete TeX file is authoritative; the legacy split source may lag behind its proof corrections. After a source revision, the tracked PDF and checksums must be rebuilt with `make` and recommitted before citing them as current.

## Release

Earlier release: [`v1.3.0`](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.3.0). The revised manuscript is on `main`; earlier tags and releases retain their original files.

## Authorship and assistance

The manuscript is signed OpenAI. OpenAI ChatGPT/Codex was used to draft and revise the text, examine the arguments and references, and check the LaTeX compilation. This statement records assistance and does not imply journal review or acceptance.
