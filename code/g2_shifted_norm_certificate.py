"""
Exact certificate for the shifted cyclic triple cover
    E: (Y+c)^2 = (x+c)^2 + lambda*P(x),  C: w^3 = Y.
Normalized chart e1=1.  Verifies:
  (1) det(M_s)=(u-v)F8, deg F8=8, homogeneous;
  (2) the explicit cyclotomic identity F=A^2+3*Delta*B^2;
  (3) a mod-13 geometrically-integral residual-octic witness at
      (c,lambda,e2,e3)=(1,1,0,2);
  (4) a mod-13 admissible rank-two conic witness at (u:v:z)=(3:5:1).
All symbolic arithmetic is exact.
"""
import sympy as sp

x,Y,w,u,v,z,lam,e2,e3,c = sp.symbols('x Y w u v z lam e2 e3 c')

S1 = 1 - 3*e2 + 3*e3
S2 = e2**3 - 3*e2*e3 + 3*e3**2
P = x**3 - S1*x**2 + S2*x - e3**3
G = x**2 + 2*c*x + lam*P                  # Y^2+2cY=G

def redY(expr):
    pp = sp.Poly(sp.expand(expr), Y, domain='EX')
    out = 0
    for (k,), coeff in pp.terms():
        A, B = sp.Integer(1), sp.Integer(0)
        for _ in range(k):
            A, B = sp.expand(B*G), sp.expand(A - 2*c*B)
        out += coeff*(A+B*Y)
    return sp.expand(out)

# Kernel basis in the normalized chart.
s0 = x-Y
s1 = Y-w**2+e2*w-e3
s2 = (x-e3)*w-Y+e2*w**2
s = u*s0+v*s1+z*s2
alpha = u*(x-Y)+v*(Y-e3)-z*Y
beta  = e2*v+z*(x-e3)
gamma = -v+e2*z
Nm = redY(alpha**3 + beta**3*Y + gamma**3*Y**2 - 3*alpha*beta*gamma*Y)

q00,q01,q02,q11,q12,q22 = sp.symbols('q00 q01 q02 q11 q12 q22')
Q = q00+q01*x+q02*Y+q11*x**2+q12*x*Y+q22*Y**2
R = sp.Poly(sp.expand(redY((Y-x)*Q)-Nm), x, Y, domain='EX')
sol = sp.solve([coef for _,coef in R.terms()], [q00,q01,q02,q11,q12,q22], dict=True, simplify=False)[0]
qs = [sol[q] for q in (q00,q01,q02,q11,q12,q22)]
M = sp.Matrix([[qs[0],qs[1]/2,qs[2]/2],
               [qs[1]/2,qs[3],qs[4]/2],
               [qs[2]/2,qs[4]/2,qs[5]]])
detM = sp.factor(M.det())
assert sp.cancel(detM/(u-v)).as_numer_denom()[1] != 0
F8 = sp.cancel(detM/(u-v))
PF8 = sp.Poly(F8,u,v,z)
assert PF8.total_degree() == 8
assert all(sum(mon)==8 for mon,_ in PF8.terms())
assert not sp.Poly(sp.expand(F8.subs(u,v)),v,z,domain='EX').is_zero
print('DETERMINANT_FACTOR=PASS')
print('RESIDUAL_DEGREE=8')
print('RESIDUAL_HOMOGENEOUS=True')
print('STRUCTURAL_FACTOR_MULTIPLICITY=1')

# Normalized v=1 sextic and principal-minor discriminant.
Fraw = sp.cancel(F8.subs(v,1))
Delta = sp.factor((-4*(qs[0]*qs[3]-(qs[1]/2)**2)).subs(v,1))
pF = sp.Poly(Fraw,z)
assert pF.degree() == 6
lc = sp.factor(pF.LC())
F = sp.cancel(Fraw/lc)
coeff = sp.Poly(F,z).all_coeffs()
_,c5,c4,c3,c2,c1,c0 = coeff
Lc = sp.factor(4*lam**2*lc)

L0 = (e2**6*lam**2*(u**2+u+1)
      -3*e2**4*e3*lam**2*(u-1)*(2*u+1)
      +3*e2**3*e3*lam*(1-lam)*(u+1)
      +9*e2**2*e3**2*lam**2*(u-1)**2
      +9*e2*e3**2*lam*(lam-1)*(u-1)
      +3*e3**2*(lam-1)**2)
assert sp.expand(Lc - (L0
      +4*c*(c+e2**3*lam)*(u**2+u+1)
      -6*c*e2*e3*lam*(u-1)*(2*u+1)
      +6*c*e3*(1-lam)*(u+1))) == 0

N20 = e2*lam*(2*e2**3*lam*u+e2**3*lam-6*e2*e3*lam*u+6*e2*e3*lam-3*e3*lam+3*e3)
N10 = -lam*(e2**3*lam*u+2*e2**3*lam+3*e2**2*e3*lam*u**2-6*e2**2*e3*lam*u+3*e2**2*e3*lam-3*e2*e3*lam*u+3*e2*e3*lam-3*e3*lam+3*e3)
N00 = lam*(u-1)*(-e2**3*lam+3*e2*e3*lam*u-3*e2*e3*lam+2*e3*lam-2*e3)
N2 = sp.expand(N20 + 2*c*e2*lam*(2*u+1))
N1 = sp.expand(N10 - 2*c*lam*(u+2))
N0 = sp.expand(N00 - 2*c*lam*(u-1))
b2,b1,b0 = [sp.cancel(N/(2*Lc)) for N in (N2,N1,N0)]
k = sp.cancel(3*Delta)
a2 = sp.cancel(c5/2)
a1 = sp.cancel((c4-a2**2-k*b2**2)/2)
a0 = sp.cancel(c3/2-a2*a1-k*b2*b1)
A = z**3+a2*z**2+a1*z+a0
B = b2*z**2+b1*z+b0
for power, target in [
    (2, c2-(a1**2+2*a2*a0+k*(b1**2+2*b2*b0))),
    (1, c1-(2*a1*a0+2*k*b1*b0)),
    (0, c0-(a0**2+k*b0**2))]:
    num = sp.expand(sp.fraction(sp.together(target))[0])
    assert num == 0, (power, sp.factor(num))
assert sp.expand(sp.fraction(sp.together(F-A**2-k*B**2))[0]) == 0
print('SHIFTED_CYCLOTOMIC_IDENTITY=PASS')
print('F=A^2+3*Delta*B^2')

# ---------- finite-field witnesses ----------
p = 13
subs_par = {c:1,lam:1,e2:0,e3:2}
# The short Weierstrass model is y^2=x^3+2x+5; the branch line y=1
# gives x^3+2x+4.
assert int(sp.discriminant(x**3+2*x+5,x)) % p == 8
assert int(sp.discriminant(x**3+2*x+4,x)) % p == 4
Fsp = sp.Poly(sp.together(F8.subs(subs_par)*4),u,v,z,modulus=p)
H = sp.Poly(Fsp.as_expr().subs(z,1),u,v,modulus=p)
assert H.total_degree()==8 and H.degree(v)==8
assert int(H.coeff_monomial(v**8))%p == 7
h = sp.Poly(H.as_expr().subs(u,4),v,modulus=p)
inv_lc = pow(int(h.LC())%p,-1,p)
h = sp.Poly(h.as_expr()*inv_lc,v,modulus=p)

def powmod_x(exp,modpoly):
    res=sp.Poly(1,v,modulus=p); base=sp.Poly(v,v,modulus=p)
    while exp:
        if exp&1: res=(res*base).rem(modpoly)
        base=(base*base).rem(modpoly); exp//=2
    return res
assert sp.gcd(h,powmod_x(p**4,h)-sp.Poly(v,v,modulus=p)).degree()==0
assert (powmod_x(p**8,h)-sp.Poly(v,v,modulus=p)).rem(h).is_zero
Hu,Hv=sp.diff(H.as_expr(),u),sp.diff(H.as_expr(),v)
def modval(expr,dd): return int(expr.subs(dd))%p
assert modval(H.as_expr(),{u:1,v:6})==0
assert (modval(Hu,{u:1,v:6}),modval(Hv,{u:1,v:6}))==(2,0)
print('MOD13_INTEGRALITY_WITNESS=PASS')
print('ELLIPTIC_DISCRIMINANT=8')
print('BRANCH_DISCRIMINANT=4')
print('RABIN_Q4_GCD=1')
print('RABIN_Q8_REMAINDER=0')
print('SMOOTH_POINT=(1,6), GRADIENT=(2,0)')
print('MONIC_SLICE=',h.as_expr())

# Admissibility witness at [u:v:z]=[3:5:1].
qvals=[modval(q.subs(subs_par),{u:3,v:5,z:1}) for q in qs]
assert int(sp.discriminant(w**3-w**2-2,w)) % p == 1
assert qvals == [5,1,7,7,0,2]
# Exact factorization over F_13:
# Q=(1+9x+Y)(1+12x+3Y), up to scalar 5.
L1=(1,9,1); L2=(1,12,3)
def prodlines(L,M):
    A,B,C=L;D,E,Fc=M
    return (A*D,A*E+B*D,A*Fc+C*D,B*E,B*Fc+C*E,C*Fc)
def normvec(vec):
    for a in vec:
        if a%p:
            ii=pow(a%p,-1,p);return tuple((b*ii)%p for b in vec)
    return tuple(vec)
assert normvec(prodlines(L1,L2)) == normvec(qvals)
# E in shifted coordinates: Y^2+2Y=x^2+2x+P(x).
S1s=7; S2s=12
Rel=sp.Poly(Y**2+2*Y-(x**2+2*x+(x**3-S1s*x**2+S2s*x-8)),x,Y,modulus=p).as_expr()
def line_cubic(L):
    A,Bc,C=L; yi=-(A+Bc*x)*pow(C,-1,p)
    return sp.Poly(sp.expand(Rel.subs(Y,yi)),x,modulus=p)
fixed=(0,-1%p,1) # Y-x=0
fixed_disc=int(sp.discriminant(line_cubic(fixed).as_expr(),x))%p
assert fixed_disc == 4
for L in (L1,L2):
    cub=line_cubic(L)
    assert cub.degree()==3 and int(sp.discriminant(cub.as_expr(),x))%p !=0
    assert L[1]%p !=0  # not a kappa-pencil norm line
    xb=(-L[0]*pow(L[1],-1,p))%p
    assert modval(Rel,{x:xb,Y:0}) != 0  # avoids branch divisor

def intersect(L,M):
    A,Bc,C=L;D,E,Fc=M; det=(Bc*Fc-E*C)%p
    assert det
    xx=(-A*Fc+D*C)*pow(det,-1,p)%p
    yy=(-Bc*D+E*A)*pow(det,-1,p)%p
    return xx,yy
def onE(pt): return modval(Rel,{x:pt[0],Y:pt[1]})==0
assert not onE(intersect(L1,L2))
assert not onE(intersect(fixed,L1))
assert not onE(intersect(fixed,L2))
print('MOD13_ADMISSIBILITY_WITNESS=PASS')
print('SECTION=(3:5:1)')
print('CONIC_COEFFS=',qvals)
print('FIXED_DIVISOR_DISCRIMINANT=1')
print('FIXED_LINE_DISCRIMINANT=',fixed_disc)
print('LINE_INTERSECTIONS=',intersect(L1,L2),intersect(fixed,L1),intersect(fixed,L2))
print('LINES=',L1,L2)
print('LINE_DISCRIMINANTS=',[int(sp.discriminant(line_cubic(L).as_expr(),x))%p for L in (L1,L2)])
