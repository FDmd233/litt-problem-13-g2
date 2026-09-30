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

Manuscript revised September 30, 2026.

## Manuscript

- [PDF](output/pdf/Litt_Problem_13_G2.pdf)
- [Complete LaTeX source](Litt_Problem_13_G2_complete.tex)
- [Split LaTeX source](paper/main.tex)

## Exact checks

The pair-incidence splitting is proved intrinsically in Section 4 from the cyclic cubic norm. Appendix A records the remaining exact characteristic-zero and finite-field checks for the residual octic and the admissible specialization, including the reductions, remainders, and Bézout identities needed for independent verification.

## Build

~~~bash
make
~~~

The resulting PDF is written to `output/pdf/Litt_Problem_13_G2.pdf`. `SHA256SUMS.txt` records the checksums of the PDF and the complete TeX source. The complete and split sources contain the same manuscript.

## Release

Earlier release: [`v1.3.0`](https://github.com/FDmd233/litt-problem-13-g2/releases/tag/v1.3.0). The revised manuscript is on `main`; earlier tags and releases retain their original files.

## Authorship and assistance

The manuscript is signed OpenAI. OpenAI ChatGPT/Codex was used to draft and revise the text, examine the arguments and references, and check the LaTeX compilation. This statement records assistance and does not imply journal review or acceptance.
