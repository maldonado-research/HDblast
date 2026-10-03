"""Pure symbolic finite-cutoff Ward audit. No physical source evaluation."""
import sympy as S

v, H = S.symbols('v H', positive=True)
d, dt, dtt = S.symbols('d d_dot d_ddot', real=True)
pi = S.pi
vdot = -H*v*(1-v*v)
Q0 = H**2/(2*pi**2)*(v*v/(2*(1+v))-v**3/48-v**5/16)
C = v**3*(H*H*d-H*dt)/(96*pi**2)
D = lambda f: S.diff(f,v)*vdot+S.diff(f,d)*dt+S.diff(f,dt)*dtt
J5 = v**3/3
J7 = v**3/3-v**5/5
J9 = v**3/3-2*v**5/5+v**7/7
A = ((-6*H*H*d+5*H*dt+dtt)*J5/16
     -(50*H*H*d+10*H*dt)*J7/32+70*H*H*d*J9/64)/(2*pi**2)
defect = S.factor((D(Q0)/2+H*Q0)*d+D(C)+4*H*C+H*A)
print('Q0_K =', Q0)
print('Q0_infinity_minus_Q0_K =', S.factor(H*H/(12*pi**2)-Q0))
print('Q0_dot_K =', S.factor(D(Q0)))
print('finite_comoving_cutoff_cosmic_Ward_defect =', defect)
print('defect_coefficient_d =', S.factor(defect.coeff(d)))
print('defect_at_v1 =', S.simplify(defect.subs(v,1)))

# Independent per-comoving-mode subtraction Ward identity. At a=1 write
# w=sqrt(P^2+r). A physical momentum P redshifts as Pdot=-HP and every
# per-comoving-mode observable has an additional a^-3 normalization.
w, r = S.symbols('w r', positive=True)
sp = dt+2*H*d
spp = dtt+5*H*dt+6*H*H*d
wp = r*H/w
wpp = 3*r*H*H/w-r*r*H*H/w**3
U = -H*H/w-wpp/(4*w*w)+3*wp*wp/(8*w**3)
u = d/(2*w)
up = sp/(2*w)-d*wp/(2*w*w)
upp = spp/(2*w)-sp*wp/w**2-d*wpp/(2*w**2)+d*wp**2/w**3
dv4 = (-U*u/w-upp/(4*w*w)+wpp*u/(4*w**3)
       +3*wp*up/(4*w**3)-3*wp*wp*u/(4*w**4))
j4 = H*up/w**2-2*H*wp*u/w**3+wp*up/(2*w**3)-3*wp*wp*u/(4*w**4)
b = -(w*w+2*r)/3
c = 1-b/w**2
esub = (d/w-H*H*u/w**2+j4)/4
psub = (c*u-d/w+c*dv4+2*b*U*u/w**3+d*U/w**2-H*H*u/w**2+j4)/4
q0sub = 1/(2*w)-U/(2*w*w)
DU = lambda f: -H*(w*w-r)/w*S.diff(f,w)+S.diff(f,d)*dt+S.diff(f,dt)*dtt
subtraction_ward = S.factor((DU(esub)+3*H*psub-q0sub*dt/2).subs(r,2*H*H))
print('per_mode_subtraction_Ward_defect =', subtraction_ward)

# Differentiate the proposed baseline primitive before any endpoint limit.
# sqrt factors cancel analytically after expressing P=v*w.
baseline_derivative = H**2/(2*pi**2)*(
    (v-v*v)/(1-v*v)**2-v*v/(2*(1-v*v))
    -S.Rational(3,8)*v*v+S.Rational(5,16)*v*v*(1-v*v))
baseline_defect = S.factor(S.diff(Q0,v)-baseline_derivative)
print('baseline_antiderivative_defect =', baseline_defect)
baseline_zero = S.simplify(Q0.subs(v,0))
baseline_limit = S.simplify(Q0.subs(v,1)-H*H/(12*pi**2))
print('baseline_zero_endpoint_defect =', baseline_zero)
print('baseline_continuum_endpoint_defect =', baseline_limit)
if any(x != 0 for x in (defect, subtraction_ward, baseline_defect, baseline_zero, baseline_limit)):
    raise RuntimeError('A claimed exact finite-cutoff Ward identity failed')
