"""Independent exact checks for the G2 residual-octic calculation.

Reconstruct the conic by polynomial division of a cyclic norm, rather
than importing the manuscript's six q_ij coefficient expressions.
Only the normalized parameter c=lambda=e1=1, e2=0, e3=2 is checked.
"""
import sympy as s

u,v,z,x,Y,t=s.symbols('u v z x Y t')
E=Y**2+2*Y-x**3+6*x**2-14*x+8
alpha=u*(x-Y)+v*(Y-2)-z*Y
beta=z*(x-2)
gamma=-v
norm=s.expand(alpha**3+beta**3*Y+gamma**3*Y**2-3*alpha*beta*gamma*Y)
quot,rem=s.div(s.Poly(s.expand(norm+(v**3+z**3*Y)*E),Y,x),s.Poly(Y-x,Y,x))
assert rem.is_zero, 'The fixed norm line does not divide exactly'
Q=s.expand(quot.as_expr())
co=lambda i,j:Q.coeff(x,i).coeff(Y,j)
a,b,c,d,e,f=[co(*ij) for ij in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]]
det4=s.expand(4*a*d*f-a*e**2-d*c**2-f*b**2+b*c*e)
F8,rem=s.div(s.Poly(det4,u,v,z),s.Poly(u-v,u,v,z))
assert rem.is_zero and s.Poly(F8.as_expr(),u,v,z).total_degree()==8
assert s.expand(det4.subs({u:1+t,v:1,z:1})-2*t*(24*t**3+56*t**2+33*t+6))==0
H=s.expand(F8.as_expr().subs(z,1))
H_reference=(6*u**3*v**5-u**3*v**4+6*u**3*v**3-u**3*v**2-u**3*v
  -5*u**2*v**6+5*u**2*v**4+3*u**2*v**3+4*u**2
  +5*u*v**7+3*u*v**6-3*u*v**4+3*u*v**3+4*u*v
  -6*v**8-2*v**7+2*v**6+v**5-2*v**4+4*v**2)
assert s.Poly(H-H_reference,u,v,modulus=13).is_zero
h=t**8-3*t**7-2*t**6+3*t**5+4*t**4+4*t**3-3*t**2-5*t-2
assert s.Poly(H.subs({u:4,v:t})-7*h,t,modulus=13).is_zero
assert s.Poly(h,t,modulus=13).is_irreducible
assert int(H.subs({u:1,v:6}))%13==0
assert [int(s.diff(H,w).subs({u:1,v:6}))%13 for w in (u,v)]==[2,0]
M=s.Matrix([[a,b/2,c/2],[b/2,d,e/2],[c/2,e/2,f]])
A=[[int(M[i,j].subs({u:3,v:5,z:1}))%13 for j in range(3)] for i in range(3)]
assert int(det4.subs({u:3,v:5,z:1}))%13==0
assert any((A[i][i]*A[j][j]-A[i][j]*A[j][i])%13 != 0 for i in range(3) for j in range(i+1,3))
print('PASS exact norm division, structural factor, characteristic-zero specialization')
print('PASS characteristic-13 octic, irreducible slice, smooth point, rank-2 conic')
print('Finite specialization only; generic factor-ordering split and Tannaka steps not certified')

# At the F_13 point [u:v:z]=[3:5:1], independently factor the conic.
Q13=s.Poly(Q.subs({u:3,v:5,z:1})-2*(Y+9*x+1)*(Y+4*x+9),x,Y,modulus=13)
assert Q13.is_zero, 'Finite-field rank-two conic failed to split'
lines=[x,-9*x-1,-4*x-9]
section_polys=[s.Poly(E.subs(Y,l),x,modulus=13) for l in lines]
branch=s.Poly(E.subs(Y,0),x,modulus=13)
assert all(s.gcd(g,g.diff()).degree()==0 for g in section_polys)
assert all(s.gcd(g,branch).degree()==0 for g in section_polys)
for i in range(3):
    for j in range(i+1,3):
        xx=s.solve(lines[i]-lines[j],x)[0]
        num,den=s.cancel(E.subs({x:xx,Y:lines[i].subs(x,xx)})).as_numer_denom()
        assert int(den)%13!=0 and int(num)%13!=0
print('PASS char-13 splitting into two rational lines, disjoint transverse sections')
