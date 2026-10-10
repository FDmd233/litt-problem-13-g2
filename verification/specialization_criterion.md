# A semismall nearby-cycles criterion for the all-parameter result

This note isolates the precise specialization argument in Section 8 of
\`Litt_Problem_13_G2_complete.tex\`. It is an audit supplement, not an
independent proof of the generic residual-octic and Tannaka calculations.

## General statement

Let \(T\) be a smooth pointed curve in characteristic zero and
\(\mathcal A/T\) an abelian scheme. Suppose
\(f:\mathcal Z\to\mathcal A\) is a projective morphism with
\(\mathcal Z/T\) **smooth and proper of relative dimension two**,
whose fibrewise images \(X_t\) are integral normal surfaces.
Assume every \(f_t:\mathcal Z_t\to X_t\) is birational and semismall,
an isomorphism off the zero section, and the exceptional fibre has
exactly \(m\) irreducible curve components with nondegenerate
intersection form, with \(m\) independent of \(t\).

The semismall decomposition theorem gives on each geometric fibre
\[
Rf_{t,*}\mathbf Q_{\mathcal Z_t}[2]
  \simeq IC_{X_t}\oplus\delta_{0_t}^{\oplus m}.
\]
By smooth base change for \(\mathcal Z/T\) and proper base change for
\(f\), nearby cycles of the left side identify with the analogous
special-fibre direct image. Nearby cycles preserve the skyscraper
sheaf along the zero section; hence
\[
\psi(IC_{X_{\bar\eta}})\oplus\delta_0^{\oplus m}
  \simeq IC_{X_0}\oplus\delta_0^{\oplus m}.
\]
The special-fibre direct image is semisimple. Krull--Schmidt
cancellation yields
\[
\boxed{\psi(IC_{X_{\bar\eta}})\simeq IC_{X_0}}.
\]
**No simultaneous choice or labeling of exceptional-curve classes
is necessary for this cancellation argument.**

Krämer--Weissauer, *Vanishing theorems for constructible sheaves on
abelian varieties*, Section 14, Lemma 14.1, establishes the
tensor-functor behavior of nearby cycles on abelian schemes in
characteristic zero, including compatibility with convolution and
character twists. Since the specialized object is semisimple,
\[
G(IC_{X_0})\hookrightarrow G(IC_{X_{\bar\eta}})
\]
is a closed embedding. Primary source:
https://ar5iv.labs.arxiv.org/html/1111.4947

## Verification in the cyclic Prym family

Under the paper's assumptions
\[
c\ne0,\quad
4A^3+27B^2\ne0,\quad
4A^3+27(B-c^2)^2\ne0,
\]
the three branch points \(b_i\) of \(w^3=y-c\) are distinct, and
their lifts \(r_i\) define a relative section \(D_0=\sum r_i\)
of the norm fibre \(\mathcal N\).

After étale-local labeling of the branch points, the formal-group
logarithm of the elliptic curve identifies the norm equation near
\(D_0\) with \(t_1+t_2+t_3=0\). At the totally ramified points,
\(t_i=w_i^3u_i(w_i)\) with \(u_i(0)\ne0\). Taking formal cube
roots of the units gives coordinates \(z_i=w_i u_i(w_i)^{1/3}\)
and the **exact relative completed equation**
\[
z_1^3+z_2^3+z_3^3=0.
\]
The blowup of its vertex is smooth, and hence the blowup of the
relative norm fibre along \(D_0\) is smooth over the base at
the singular section. Away from that section, the norm-map
critical-point analysis in Section 3 of the manuscript establishes
relative smoothness: cycles with an unramified point are smooth,
the potentially critical \(2r_i+r_j\) would force
\(b_i=b_k\), and the triple cycle \(3r_i\) has a nonzero linear
\(3e_3\) term in symmetric coordinates.

The resulting projective relative surface
\(\widehat{\mathcal N}\) is therefore smooth. For every fibre,
\(\widehat{\mathcal N}_t\to X_t\) is semismall, with exceptional
elliptic and rational curves of intersection matrix
\[
\begin{pmatrix}-3&1\\1&-1\end{pmatrix},
\qquad\det=2.
\]
Thus \(m=2\), and the general lemma applies to all allowed
parameters.

## Exact boundary of this audit

The argument verifies the *specialization interface*, using the
absolute norm-fibre geometry established earlier in the manuscript.
It does **not** independently certify the generic residual-octic
integrality, the splitting of the factor-ordering cover,
the convolution extraction of \(IC_X\), or the final Hodge-group
comparison. Those must be reviewed separately before calling the
whole manuscript independently verified.
