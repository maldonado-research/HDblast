# Reconstructed from actual executed tool-call text after the failed run.
# This is not an original on-disk execution snapshot and has not been replayed.
import sympy as S
w,r,H,h,z,zz=S.symbols('w r H h z zz', nonzero=True)
# Unit scale factor only at observation; h=delta x, z=delta x dot, zz=delta x ddot.
sp=z+2*H*h; spp=zz+5*H*z+6*H**2*h
wp=r*H/w; wpp=3*r*H**2/w-r**2*H**2/w**3
U=-H**2/w-wpp/(4*w**2)+3*wp**2/(8*w**3)
u=h/(2*w); up=sp/(2*w)-h*wp/(2*w**2)
upp=spp/(2*w)-sp*wp/w**2-h*wpp/(2*w**2)+h*wp**2/w**3
v=-U*u/w-upp/(4*w**2)+wpp*u/(4*w**3)+3*wp*up/(4*w**3)-3*wp**2*u/(4*w**4)
j=H*up/w**2-2*H*wp*u/w**3+wp*up/(2*w**3)-3*wp**2*u/(4*w**4)
b=-(w**2+2*r)/3; c=1-b/w**2
psub=(c*u-h/w+c*v+2*b*U*u/w**3+h*U/w**2-H**2*u/w**2+j)/4
q=-h/(4*w**3)
qp=-sp/(4*w**3)+H*h/(2*w**3)+3*h*wp/(4*w**4)
qpp=-spp/(4*w**3)+H*sp/w**3-H**2*h/(2*w**3)+3*sp*wp/(2*w**4)-3*H*h*wp/w**4+3*h*wpp/(4*w**4)-3*h*wp**2/w**5
q0=1/(2*w)-U/(2*w**2)
local=S.expand((qpp+H*qp-3*H**2*q)/6-q0*h/6-psub)
print('pressure subtraction mismatch per mode:',S.factor(local))
integ=0
for term in S.Add.make_args(local):
    power=term.as_powers_dict().get(w,0)
    n=-power
    if n<=3: raise RuntimeError(('nonintegrable mismatch', term))
    moment=S.sqrt(S.pi)*S.gamma((n-3)/2)*r**((3-n)/2)/(4*S.gamma(n/2))
    integ+=term/w**power*moment
integ=S.simplify(integ)
expect=(zz+2*H*z-11*H**2*h)/144
print('integral times 2pi²:',integ)
if S.simplify(integ-expect)!=0: raise RuntimeError(('mismatch',S.simplify(integ-expect)))
print('PASS exact direct-pressure local contact coefficient')
