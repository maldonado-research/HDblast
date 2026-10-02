"""Symbolic local subtraction algebra, no field evolution or numerical experiment."""
import sympy as s
w,r,h,H,d,e,f,z,zz=s.symbols('w r h H d e f z zz')
D=lambda q:s.diff(q,w)*(-H*(w*w-r)/w)+s.diff(q,H)*d+s.diff(q,d)*e+s.diff(q,e)*f+s.diff(q,h)*z+s.diff(q,z)*zz
u=(h-d-2*H**2)/(2*w)-r*(d+3*H**2)/(4*w**3)+s.Rational(5,8)*r**2*H**2/w**5
up=s.expand(H*u+D(u)); upp=s.expand(2*H*up+D(up))
wp=r*H/w; wpp=r*(d+3*H**2)/w-r**2*H**2/w**3
v=-u**2/(2*w)-upp/(4*w**2)+wpp*u/(4*w**3)+3*wp*up/(4*w**3)-3*wp**2*u/(4*w**4)
B=-(w*w+2*r)/3; C=1-B/w**2
sh2=H*wp/w**2+wp**2/(4*w**3)
sh4=H*up/w**2-2*H*wp*u/w**3+wp*up/(2*w**3)-3*wp**2*u/(4*w**4)
rho=w/2+((h+H**2)/w+sh2)/4+(u**2/w-(h+H**2)*u/w**2+sh4)/4
p=(w*w-r)/(6*w)+(C*u+(H**2-h)/w+sh2)/4+(C*v+B*u**2/w**3-(H**2-h)*u/w**2+sh4)/4
q=1/(2*w)-u/(2*w**2)
DDq=(D(D(q))-3*H*D(q)-3*d*q)/2
anomaly=s.expand(DDq-(h+r)*q+rho-3*p)
I={n:s.sqrt(s.pi)*s.gamma(s.Rational(n-3,2))/(4*s.gamma(s.Rational(n,2)))*r**s.Rational(3-n,2) for n in range(5,16,2)}
out=0
for term,coeff in s.Poly(anomaly,h,H,d,e,f,z,zz).terms():
 integ=0
 for piece in s.Add.make_args(s.expand(coeff)):
  fac,pow=piece.as_coeff_exponent(w)
  if -pow not in I: raise RuntimeError((term,coeff,pow))
  integ+=fac*I[-pow]
 integ=s.simplify(integ)
 if integ:
  print(term,'integrand:',s.factor(coeff),'integral*2pi2:',integ)
  out+=integ*s.prod(a**b for a,b in zip([h,H,d,e,f,z,zz],term))
print('Anomaly times 2 pi^2:',s.factor(out))
expected=(58*H**4-14*H**2*d-60*H**2*h-42*H*e+15*H*z-9*d**2-30*d*h-6*f+15*h**2+5*zz)/240
assert s.expand(out-expected)==0
rho4=s.expand((u**2/w-(h+H**2)*u/w**2+sh4)/4)
p4=s.expand((C*v+B*u**2/w**3-(H**2-h)*u/w**2+sh4)/4)
rholog=-rho4.coeff(w,-3)/4
plog=-p4.coeff(w,-3)/4
qlog=-s.expand(-u/(2*w**2)).coeff(w,-3)/4
expected_rholog=h**2/64-h*H**2/32-H*z/32+(-d**2+6*H**2*d+2*H*e)/64
expected_plog=-h**2/64+h*d/48+h*H**2/32+zz/96+H*z/48-3*d**2/64-3*H**2*d/32-H*e/16-f/96
expected_qlog=(h-d-2*H**2)/16
for actual,expect in [(rholog,expected_rholog),(plog,expected_plog),(qlog,expected_qlog)]:
 assert s.expand(actual-expect)==0,(actual,expect)
print('PASS: universal rational-symbol local trace, rho4, p4 and Q2 log-coefficient identities. No numerical field evolution.')
rho0=w/2; p0=(w*w-r)/(6*w); q0=1/(2*w)
rho2=s.expand(((h+H**2)/w+sh2)/4)
p2=s.expand((C*u+(H**2-h)/w+sh2)/4)
q2=s.expand(-u/(2*w**2))
for label,residual in [('Ward0',D(rho0)+3*H*p0),('Ward2',D(rho2)+3*H*p2-z*q0/2),('Ward4',D(rho4)+3*H*p4-z*q2/2)]:
 assert s.expand(residual)==0,(label,s.factor(residual))
print('PASS: all three graded Ward identities for general symbolic FRW jets.')
# At order two the I_1 and I_3 logarithms both contribute.  This
# directly fixes the affine-curvature action coefficients in pressure
# without inferring them by dividing the Ward identity by H.
rho2log=s.expand(r*rho2.coeff(w,-1)/8-rho2.coeff(w,-3)/4)
p2log=s.expand(r*p2.coeff(w,-1)/8-p2.coeff(w,-3)/4)
q0log=s.expand(r*q0.coeff(w,-1)/8-q0.coeff(w,-3)/4)
assert s.expand(rho2log-r*(h-H**2)/32)==0
assert s.expand(p2log-(-r*h/32+r*(2*d+3*H**2)/96))==0
assert s.expand(q0log-r/16)==0
print('PASS: direct order-two rho and pressure, and order-zero Q, PV logarithmic action coefficients.')
print('SymPy version:',s.__version__)
