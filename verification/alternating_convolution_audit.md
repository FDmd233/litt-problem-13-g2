# Independent audit: antisymmetric pair convolution detects the theta IC

**Scope.** This note checks the sheaf-theoretic extraction in Section 7 of
\`Litt_Problem_13_G2_complete.tex\`. It is conditional on the preceding
geometric statement that the pair fibre over a Zariski-dense open of the
smooth theta surface has precisely two geometrically integral
nonexceptional components, globally labelled and exchanged by
transposition. It is **not** an independent proof of that geometric
statement.

Let \(P\) be the Prym threefold, \(X=-X\subset P\) its integral normal
theta surface, \(\rho:S\to X\) the semismall resolution, and
\[
K=IC_X,\qquad \widetilde K=R(j:S\to P)_*\mathbf Q_S[2]
=K\oplus\delta_0.
\]
Let \(\sigma\) swap the two factors of \(S\times S\).

## 1. Exact stalk shifts

The convolution is
\[
\widetilde K*\widetilde K
 =R\widetilde m_*\mathbf Q_{S\times S}[4].
\]
For \(t\in X_{\mathrm{reg}}\) on the pair-splitting locus,
\[
\mathcal H^{-2}(\widetilde K*\widetilde K)_t
 =H^2(\widetilde m^{-1}(t),\mathbf Q).
\]
The four one-dimensional components consist of two genuinely
nonexceptional curves \(A_{1,t},A_{2,t}\) and the structural curves
\(E_{\mathrm{exc}}\times\{s_t\}\),
\(\{s_t\}\times E_{\mathrm{exc}}\), exchanged within each pair.
Top cohomology has one basis vector per irreducible reduced component.
Transposition preserves complex orientation. Therefore the
anti-invariant top cohomology has rank two.

Expanding \(\widetilde K=K\oplus\delta_0\) gives
\[
(\widetilde K*\widetilde K)_-
=(K*K)_-\oplus\bigl((K*\delta_0)\oplus
(\delta_0*K)\bigr)_-,
\]
where the second summand is one copy of \(K\), concentrated in
ordinary degree \(-2\) over \(X_{\mathrm{reg}}\). Consequently
\[
\operatorname{rk}\mathcal H^{-2}((K*K)_-)|_{U}=1.
\]
The geometric shifts are both even (2), so the commutativity
constraint does not introduce a Koszul sign.

## 2. Off-divisor exclusion

For general \(t\in P\setminus X\) the fibre
\(X\cap(t-X)\) is a smooth connected projective curve:
the translated ample divisor is connected on the integral
Cohen--Macaulay surface, and generic smoothness applies outside
the singular points. Its involution \(x\mapsto t-x\) preserves
orientation and acts as the identity on \(H^2\). Hence
\[
\mathcal H^{-2}((K*K)_-)=0
\]
on a nonempty Zariski-open subset of \(P\setminus X\).

## 3. Why the rank-one line cannot be a shifted full-support IC

Use the decomposition theorem over characteristic zero to split
\((K*K)_-\) into shifted simple intersection complexes. At a
general smooth point of the divisor \(X\), a full-support perverse
intersection complex with a local system on its complement is
locally an external product
\[
\mathcal L[2]\boxtimes j_{\!*}\mathcal W[1]
\]
for a punctured normal disk and a normal monodromy operator \(T\).
The middle-extension disk stalk has ordinary cohomology only
in degree \(-1\), with value \(\ker(T-1)\).
Therefore an *unshifted* full-support summand contributes along
\(X\) only in degree \(-3\). For a shifted full-support summand
to contribute in degree \(-2\), the shift must be \([-1]\).
But the same shifted summand contributes a nonzero local system
in degree \(-2\) on a general point of \(P\), contradicting the
off-divisor vanishing. Semisimplicity forbids cancellation between
direct summands.

A proper-support summand containing a generic point of the
two-dimensional \(U\subset X\) must have support precisely \(X\).
To contribute in degree \(-2\) it must be unshifted.
The total rank is one, so its local system is rank one.

## 4. Why the rank-one local system is trivial globally

The two ordered components are globally labelled after spreading
on the base. Their top fundamental classes are therefore constant
and transposition exchanges them. Removing the structural
anti-invariant line leaves a *constant* rank-one local system over
\(U\). Let \(\mathcal V\) be the rank-one local system of the
extracted divisor-supported summand on its smooth dense support
\(X_\mathcal V\subset X_{\mathrm{reg}}\). A nonempty Zariski-open
embedding \(U\cap X_\mathcal V\hookrightarrow X_\mathcal V\)
induces a surjection on topological fundamental groups because
the omitted subset has complex codimension at least one. Thus
the monodromy of \(\mathcal V\) is globally trivial, and the
summand is exactly \(IC_X\), not a twisted version.

After quotienting negligible perverse sheaves, exactness and
\(\chi(K)=7>0\) preserve this summand and yield the genuinely
representation-theoretic consequence
\[
\operatorname{Hom}_{G(K)}\bigl(V,\Lambda^2 V\bigr)\ne0.
\]

## Audit boundary

The argument above checks the signs, degrees, supports,
monodromy, and perverse shift logic **conditional on the
geometry of the globally split pair incidence**. It is not a
substitute for independently re-proving residual-octic
integrality, the rank-two factor ordering cover, or the
boundary exhaustion. The finite/symbolic calculations in
\`verification/verify_residual_octic.py\` address a separate
piece of that geometry.
