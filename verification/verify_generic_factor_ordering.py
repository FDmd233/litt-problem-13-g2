"""Independent generic cyclic-norm factor-ordering identity.

Check over Q(c, lambda, e1,e2,e3,u,v,z,omega) with omega^2+omega+1=0.
A single exact symbolic identity, not finite-field sampling.
"""
import sympy as S
x,Y,c,lam,e1,e2,e3,u,v,z,omega=S.symbols('x Y c lam e1 e2 e3 u v z omega')
Sigma1=e1**3-3*e1*e2+3*e3
Sigma2=e2**3-3*e1*e2*e3+3*e3**2
PD=x**3-Sigma1*x**2+Sigma2*x-e3**3
FE=(Y+c)**2-(x+c)**2-lam*PD
alpha=u*(x-Y)+v*(Y-e3)-e1*z*Y
beta=e2*v+z*(x-e3)
gamma=-e1*v+e2*z
Nm=S.expand(alpha**3+beta**3*Y+gamma**3*Y**2-3*alpha*beta*gamma*Y)
R=S.expand(lam*Nm+(v**3+z**3*Y)*FE)
quo,rem=S.div(R,Y-x,Y)
assert S.expand(rem)==0
assert S.degree(quo,x)<=2 and S.degree(quo,Y)<=2
Q=S.cancel(quo/lam)
q=-v/z
Ystar=q**3

def simp_omega(expr):
    p=S.Poly(S.expand(expr),omega)
    return S.rem(p,S.Poly(omega**2+omega+1,omega)).as_expr().expand()

def L(w):
    return u*(x-Ystar)+v*(Ystar-e3)-e1*z*Ystar+(e2*v+z*(x-e3))*w+(-e1*v+e2*z)*w**2
L0=S.cancel(L(q))
assert S.factor(L0-(u-v)*(x-Ystar))==0
L1=S.cancel(L(omega*q))
L2=simp_omega(L(omega**2*q).expand().subs(omega**3,1))
A=simp_omega(S.cancel((Q.subs(Y,Ystar)+(u-v)*L1*L2)*lam*z**8).expand())
assert S.expand(A)==0, str(S.factor(A))[:600]
a1=simp_omega(S.expand(L1).coeff(x,1))
b1=simp_omega(S.expand(L1).coeff(x,0))
a2=simp_omega(S.expand(L2).coeff(x,1))
b2=simp_omega(S.expand(L2).coeff(x,0))
Rpoly=e1*u*v**2*z+2*e1*v**3*z-2*e2*u*v*z**2-e2*v**2*z**2+e3*u*z**3-e3*v*z**3+u*v**3-v**4
root_sep=simp_omega(S.expand((a1*b2-a2*b1+v*(2*omega+1)*Rpoly/z**3)*z**6))
assert S.expand(root_sep)==0, str(S.factor(root_sep))[:600]
assert S.Poly(Rpoly,u,v,z).total_degree()==4
print('PASS all-parameter rational-root separation formula of degree 4')
print('PASS generic off-curve norm identity (all c,lam,e1,e2,e3)')
print('PASS generic restrictions Q(x,Y*)=-(u-v)L1(x)L2(x)')
print('PASS no finite-field specialization used')
