# Independent audit of the geometric two-branch exhaustion

This audit revisits Sections 4--7 of `Litt_Problem_13_G2_complete.tex`
and separates the **generic algebraic splitting**, the **scheme-level
incidence conversion**, the **exclusion of boundary components**,
and the **Galois descent of the two global labels**.
The new arbitrary-parameter symbolic check is
[verify_generic_factor_ordering.py](verify_generic_factor_ordering.py).
This is a mathematical proof audit, not a machine-formalized proof.

## 1. The ordering cover is rationally split over the generic residual octic

Let $K_0=\mathbb C(c,\lambda,e_1,e_2,e_3)$;
let $\mathscr R/K_0$ be the geometrically integral
plane octic $F_8(u,v,z)=0$, and $L=K_0(\mathscr R)$.
The conic $Q_s=0$ has symmetric matrix $M_s$
with
$
 \det M_s=(u-v)F_8.
$
For the residual octic, $u-v\ne0$ at its
generic point. If $\operatorname{rank}M_s\le1$
there, its adjugate vanishes and therefore
every first derivative of $\det M_s$ vanishes;
this contradicts generic smoothness of the
geometrically integral characteristic-zero
plane octic. Hence $\operatorname{rank}M_s=2$
generically.

Put $q=-v/z$, $Y_*=q^3$, and
choose a primitive cube root $\omega\in\mathbb C$.
For the section $s=\alpha+\beta w+\gamma w^2$
define $L_j(x)=s(x,\omega^j q)$ after
substituting $Y=Y_*$, $j=0,1,2$.
The exact identity
$
 \operatorname{Nm}(s)
  =(Y-x)Q_s-
     \frac{v^3+z^3Y}{\lambda}F_E(x,Y)
$
holds without imposing any specialization on
$c,\lambda,e_1,e_2,e_3$.
Its restriction to $Y=Y_*$ gives
$
 Q_s(x,Y_*)=-(u-v)L_1(x)L_2(x),
$
since $L_0(x)=(u-v)(x-Y_*)$.
The symbolic script reconstructs $Q_s$ by
exact polynomial division and proves **both**
identities in the rational function field
$
 \mathbb Q(c,\lambda,e_1,e_2,e_3,u,v,z,\omega)/
        (\omega^2+\omega+1).
$
No point over $\mathbb F_{13}$ is substituted.

Write $L_j=a_jx+b_j$. The script also checks
$
 a_1b_2-a_2b_1
 =-\frac{v(2\omega+1)}{z^3}\,R(u,v,z),
$
with
$
 \begin{aligned}
 R={}&e_1uv^2z+2e_1v^3z
       -2e_2uvz^2-e_2v^2z^2\\
   &+e_3uz^3-e_3vz^3+uv^3-v^4.
 \end{aligned}
$
The polynomial $R$ is homogeneous of
degree four in $(u,v,z)$, whereas $F_8$
is irreducible of degree eight, so
$R\not\equiv0$ at the generic point of
$\mathscr R$. The roots $-b_j/a_j$ of
$L_1,L_2$ are therefore distinct there.

The rank-two conic has a unique singular
point $\mathbf s\in\mathbb P^2(L)$:
a nonzero column of
$\operatorname{adj}(M_s)$ gives
its $L$-rational coordinates.
Its two intersections with the
$L$-rational line $Y=Y_*$ are
$
 p_j=[-b_j:a_jY_*:a_j]\in\mathbb P^2(L),
 \qquad j=1,2.
$
They are distinct and cannot equal
$\mathbf s$, since a line through
the singular point meets a union of
two lines with multiplicity two.
The two line components of the conic
are therefore **explicitly**
$
 \ell_j=\overline{\mathbf s\,p_j}
 \quad\text{over }L.
$
This gives two $L$-rational sections of
the intrinsic factor-ordering double cover.
The cover is finite étale on the rank-two
locus, so the rational sections extend
after deleting a closed subset.
This is the algebraic content of
generic factor-ordering splitting.

## 2. Why this produces two genuine Abel incidence branches

Work over the open $\mathscr G$ where:
the three norm lines have transverse,
pairwise disjoint intersection divisors;
all avoid the branch divisor; one simple
zero of the cyclic section lies over each
residual point; and both residual triples
are $h^0=1$ divisors.

The two orderings $(\ell_1,\ell_2)$
give ordered triples $D_1,D_2$ via
the norm identity. This construction
works **scheme-theoretically** because
étale locally above a residual line
the cubic cover splits into three sheets
and the norm product $a_1a_2a_3=0$
has exactly two unit factors; the third
is zero even in a nonreduced base ring.
Transversality makes each zero simple.
The three Cartier divisors
$D,D_1,D_2$ have total degree nine,
equal to the degree of $3\kappa$;
there can be no residual vertical or
nilpotent contribution.

Conversely, an ordered pair
$y_1,y_2$ in the smooth Abel locus
with $x+y_1+y_2=0$ determines its
unique degree-three Abel divisors and
then the unique projective section of
$3\kappa-D$ cutting $D+D_1+D_2$.
Kempf's $h^0=1$ local determinantal
model gives a base-change-compatible
inverse Abel map, so the two constructions
are mutually inverse **as schemes**,
not only on complex points.
Thus on $\mathscr G$, incidence has
two disjoint geometrically integral
branches, exchanged by the involution.

## 3. Why no third one-dimensional branch hides on the boundary

Fix geometric generic $x\in X_{\mathrm{reg}}$ and
write
$
 C_x=X\cap(-x-X),\quad
 C_x^{ne}=C_x\setminus\{0,-x\}.
$
Every $y\in C_x^{ne}$ has a unique
degree-three Abel divisor $D_y$,
and similarly for $-x-y$. Their
sum with the fixed $D_x$ gives
a projective section of
$H^0(C,3\kappa-D_x)$, yielding
a rational map $\Psi:C_x\dashrightarrow\mathbb P^2$.

**Finiteness.** A fixed section $s$
of $3\kappa$ has a fixed effective
degree-nine zero divisor, hence only
finitely many degree-three subdivisors.
Therefore $\Psi$ has finite
fibres wherever defined. Every
one-dimensional component of $C_x^{ne}$
has a one-dimensional image.

**No structural-line component.**
On the pencil $u=v$, $s=h(u+zw)$ and
$
 \operatorname{div}(s)
 =D_x+D_{-x}+D_\kappa,\qquad
 D_\kappa\in|\kappa|.
$
For a general pencil point the two
residual triples project to **distinct,
disjoint lines** on the elliptic curve.
Every degree-three subcycle whose
projection is a line section must use
one complete triple: two points on one
line determine that line and exclude
a third point on the other line.
Hence the only residual Abel classes
are $-x,0$, both removed in $C_x^{ne}$.
Any additional curve component mapped
dominantly to $u=v$ would contradict
this generic statement; a component
collapsed to one parameter would
contradict the finite-fibre property.

**No hidden curve on the residual-octic
boundary.** Every remaining curve
component maps to $F_8=0$.
Geometric integrality of $F_8$
forces its one-dimensional image to
be dense. Consequently it meets the
dense admissible open $\mathscr G$,
where the scheme consists of just
two disjoint geometrically integral
branches. Every such component is
the closure of one of those branches.
There are exactly **two** components;
there is no room for a third curve
supported entirely on the
non-admissible finite subset of the
residual octic.

**Reducedness.** $X$ and its
nontrivial translate are two distinct
Cartier divisors in the smooth Prym
threefold; their intersection is a
pure one-dimensional local complete
intersection, hence Cohen--Macaulay
with no embedded associated points.
It is generically reduced on each of
the two admissible branches.
Therefore the entire geometric
generic pair fibre is reduced.

## 4. Why the two components have rational global labels

Let $\tau$ be the covering automorphism.
The Prym identity $1+\tau+\tau^2=0$
gives a $k(x)$-rational section
$
 x\longmapsto y_+(x)=\tau x
$
of the pair-incidence family, with
partner $\tau^2x$.
If the fibre were singular at this
point for general $x$, the tangent
hyperplanes of $X$ at
$\tau x,\tau^2x$ would coincide.
The Gauss map is $\tau$-equivariant,
so the general Gauss covector would
be projectively fixed by $\tau$.
But on the Prym cotangent space,
$\tau$ has two distinct eigenvalues
with multiplicities one and two;
its projective fixed locus is a point
and a line, a proper subset of
$\mathbb P^2$. This contradicts the
positive conormal Gauss degree six,
which makes the Gauss map dominant.
Thus $y_+(x)$ is a smooth point
of the geometric generic fibre,
on exactly one of the two components.

Galois must preserve the unique
component containing this rational
smooth point; the other is then
also fixed, as there are exactly two.
They both descend to $k(x)$.
Projectivity and generic flatness
spread their closures to flat
subfamilies after shrinking the
integral base. The usual openness
of geometric integrality and
reducedness, followed by removal of
the constructible support of
the residual ideal, yields a dense
open where **all** geometric pair
fibres have exactly these two
labelled nonexceptional components,
scheme-theoretically.

## 5. Audit boundary

The logical geometry in Sections 1--4 has
been checked separately from the
representation-theoretic argument.
The key algebraic factorization and
root-separation formulas in Section 1
now have an **arbitrary-parameter exact
symbolic certificate**, not merely a
special finite-field test.

These deductions continue to rely on
the manuscript's earlier proof of:
the generic residual octic's
geometric integrality (for which a
separate characteristic-13
irreducibility certificate exists);
the unique trigonal pencil
$|\kappa|$; the description of
the theta surface and its Gauss
degree; and the existence of a
nonempty admissible locus.
Those are **explicit dependencies**,
not omitted proof obligations.
The accompanying
[alternating convolution audit](alternating_convolution_audit.md)
then checks the top-cohomology and
perverse-sheaf deduction.

This constitutes a much stronger
scheme-level audit of the geometric
pair-splitting chain, but does not
claim a formal proof assistant
certificate, independent journal
peer review, or verification of
every remaining external
Tannakian/Hodge input.
