import sympy as s
L,M,a,w,h,h1,h2,h3,h4,h5=s.symbols('L M a w h h1 h2 h3 h4 h5', nonzero=True)
jets=[h,h1,h2,h3,h4,h5]
def D(e):
    return s.expand(s.diff(e,L)*L**2+s.diff(e,M)*2*L*M+s.diff(e,a)*L*a+s.diff(e,w)*L*M/w+sum(s.diff(e,jets[i])*jets[i+1] for i in range(5)))
def clean(e): return s.expand(e)
k2=w*w-M
C=2*L**2
wp=L*M/w; wpp=D(wp)
U=-C/(2*w)-wpp/(4*w*w)+3*wp**2/(8*w**3)
Up=D(U); Upp=D(Up)
V=-U**2/(2*w)-Upp/(4*w*w)+wpp*U/(4*w**3)+3*wp*Up/(4*w**3)-3*wp**2*U/(4*w**4)
b=-k2/3-M;c=1-b/w**2
J2=L*wp/w**2+wp**2/(4*w**3)
J4=L*Up/w**2-2*L*wp*U/w**3+wp*Up/(2*w**3)-3*wp**2*U/(4*w**4)
R=w/2+(L**2/w+J2)/4+(U**2/w-L**2*U/w**2+J4)/4
P=k2/(6*w)+(c*U+L**2/w+J2)/4+(c*V+b*U**2/w**3-L**2*U/w**2+J4)/4
S=1/(2*w)-U/(2*w**2)
z=M*h/w; z1=D(z);z2=D(z1)
u=-(h2+2*L*h1)/(2*w)+C*z/(2*w**2)-z2/(4*w**2)+wpp*z/(2*w**3)+3*wp*z1/(4*w**3)-9*wp**2*z/(8*w**4)
u1=D(u);u2=D(u1)
v=-U*u/w+U**2*z/(2*w**2)-u2/(4*w**2)+Upp*z/(2*w**3)+(z2*U+wpp*u)/(4*w**3)-3*wpp*U*z/(4*w**4)+3*(z1*Up+wp*u1)/(4*w**3)-9*wp*Up*z/(4*w**4)-3*(2*wp*z1*U+wp**2*u)/(4*w**4)+3*wp**2*U*z/w**5
# Independent dual-number directional differentiation of canonical formulas.
eps=s.symbols('eps')
ww,ww1,ww2,uu,uu1,uu2,ll,mm=s.symbols('ww ww1 ww2 uu uu1 uu2 ll mm')
VV=-uu**2/(2*ww)-uu2/(4*ww**2)+ww2*uu/(4*ww**3)+3*ww1*uu1/(4*ww**3)-3*ww1**2*uu/(4*ww**4)
BB=-s.symbols('kk2')/3-mm;CC=1-BB/ww**2
JJ2=ll*ww1/ww**2+ww1**2/(4*ww**3)
JJ4=ll*uu1/ww**2-2*ll*ww1*uu/ww**3+ww1*uu1/(2*ww**3)-3*ww1**2*uu/(4*ww**4)
RR=ww/2+(ll**2/ww+JJ2)/4+(uu**2/ww-ll**2*uu/ww**2+JJ4)/4
PP=s.symbols('kk2')/(6*ww)+(CC*uu+ll**2/ww+JJ2)/4+(CC*VV+BB*uu**2/ww**3-ll**2*uu/ww**2+JJ4)/4
SS=1/(2*ww)-uu/(2*ww**2)
base={ww:w,ww1:wp,ww2:wpp,uu:U,uu1:Up,uu2:Upp,ll:L,mm:M,s.symbols('kk2'):k2}
variations={ww:z,ww1:z1,ww2:z2,uu:u,uu1:u1,uu2:u2,ll:h1,mm:2*M*h}
def metric(f): return s.expand(sum(s.diff(f,t).subs(base)*dv for t,dv in variations.items()))
dR=metric(RR);dP=metric(PP);dS=metric(SS)
# Scalar-source variations at fixed geometry, using sigma as canonical mass source.
sigma=2*M*h-h2-2*L*h1
us=sigma/(2*w);us1=D(us);us2=D(us1)
vs=-U*us/w-us2/(4*w**2)+wpp*us/(4*w**3)+3*wp*us1/(4*w**3)-3*wp**2*us/(4*w**4)
js=L*us1/w**2-2*L*wp*us/w**3+wp*us1/(2*w**3)-3*wp**2*us/(4*w**4)
dRs=(sigma/w-L**2*us/w**2+js)/4
dPs=(c*us-sigma/w+c*vs+2*b*U*us/w**3+sigma*U/w**2-L**2*us/w**2+js)/4
dSs=-sigma/(4*w**3)
contacts={'Q':clean(dSs-dS),'rho':clean(dRs-dR+(h2+4*L*h1)*S/2),'p':clean(dPs-dP-h2*S/2)}
for label,e in contacts.items():
 e=s.expand(e.subs(M,2*L**2)); print(label,s.factor(e))
 integ=0
 for term in s.Add.make_args(e):
  power=term.as_powers_dict().get(w,0);coeff=term/w**power
  if power>=-3: raise ValueError(('not individually convergent',label,term))
  n=-power
  I=s.sqrt(s.pi)*s.gamma((n-3)/2)/(4*s.gamma(n/2))*(2*L**2)**((3-n)/2)
  integ+=coeff*I
 print('integral no 1/(2pi^2):',s.factor(s.simplify(integ)))
print('ward baseline checks')
for label,X,Y in [('0',w/2,k2/(6*w)),('2',(L**2/w+J2)/4,(c*U+L**2/w+J2)/4),('4',(U**2/w-L**2*U/w**2+J4)/4,(c*V+b*U**2/w**3-L**2*U/w**2+J4)/4)]:
 print(label,s.factor((D(X)-L*X+3*L*Y).subs(M,2*L**2)))
print('ward metric checks')
# Physical density a^-4: baseline D R -L R +3L P. Metricvariation gives delta divergence.
print(s.factor((D(dR)-L*dR+3*L*dP-h1*R+3*h1*P).subs(M,2*L**2)))
# Machine-readable exact Laurent inventories, without physical evaluations.
import json
from pathlib import Path
inventory={}
for label,e in contacts.items():
 e=s.expand(e.subs(M,2*L**2))
 coefficients={str(-power):str(s.factor(e.coeff(w,power))) for power in range(-15,-4) if e.coeff(w,power)!=0}
 inventory[label]={'integrand_coefficients_of_omega_minus_n':coefficients}
Path(__file__).with_name('FINITE_K_CONTACT_INVENTORY.json').write_text(json.dumps({'symbols':{'L':'-1/eta','h':'h','h1':'h prime','h2':'h second','h3':'h third','h4':'h fourth','M':'a^2 r = 2L^2'},'measure':'k^2 dk/(2 pi^2)','primitive':'I_n(K,M)=M^((3-n)/2) sum_l (-1)^l binomial((n-5)/2,l) v^(2l+3)/(2l+3), v=K/sqrt(K^2+M)','contacts':inventory},indent=2)+'\n')
