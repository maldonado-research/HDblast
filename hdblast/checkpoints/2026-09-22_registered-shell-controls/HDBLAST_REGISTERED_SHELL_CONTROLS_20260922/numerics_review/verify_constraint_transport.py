#!/usr/bin/env python3
"""Independent exact checks of the continuum constraint transport identities.

Requires SymPy. Does not import or run any HDBLAST evolution solver.
Writes only the adjacent CONSTRAINT_TRANSPORT_CHECKS.json.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

t, z = s.symbols('t z', real=True)
A, B, p = [s.Function(q)(t, z) for q in ('A', 'B', 'phi')]
U = s.Function('U')(p)
At, Az, Bt, Bz, pt, pz = [s.diff(f, x) for f in (A, B, p) for x in (t, z)]
EA = s.diff(A,z,2)-3*At**2+3*Az**2+s.Rational(2,3)*s.exp(2*B)*U
EB = s.diff(B,z,2)+3*At**2-3*Az**2-pt**2/2+pz**2/2-s.exp(2*B)*U/3
EP = s.diff(p,z,2)-3*At*pt+3*Az*pz-s.exp(2*B)*s.diff(U,p)
H = -2*s.exp(2*B)*U+6*At**2+6*At*Bt-12*Az**2+6*Az*Bz-6*s.diff(A,z,2)-pt**2-pz**2
M = -3*s.diff(A,t,z)-3*At*Az+3*At*Bz+3*Az*Bt-pt*pz
evolution = {
    s.diff(A,t,2,z):s.diff(EA,z),
    s.diff(A,t,2):EA,
    s.diff(B,t,2):EB,
    s.diff(p,t,2):EP,
}
def reduced(expr):
    return s.simplify(s.expand(expr.subs(evolution, simultaneous=True).doit()))
checks = []
def check(name, expr, expected_zero=True):
    result = reduced(expr)
    ok = (result == 0) if expected_zero else (result != 0)
    checks.append({'name':name,'pass':bool(ok),'negative_control':not expected_zero})
    if not ok: raise RuntimeError(name+': '+str(result))

check('Hamiltonian transport from direct differentiation',s.diff(H,t)-2*s.diff(M,z)-6*Az*M+3*At*H)
check('Momentum transport from direct differentiation',s.diff(M,t)-s.diff(H,z)/2-s.Rational(3,2)*Az*H+3*At*M)
Cp,Cm=H+2*M,H-2*M
check('Outgoing density follows z=z0-t',s.diff(s.exp(3*A)*Cp,t)-s.diff(s.exp(3*A)*Cp,z))
check('Incoming density follows z=z0+t',s.diff(s.exp(3*A)*Cm,t)+s.diff(s.exp(3*A)*Cm,z))
# Independent negative controls: wrong transport sign, wrong density weight,
# incorrect Einstein-scalar kinetic normalization and wrong source coefficient.
check('reject outgoing characteristic sign swap',s.diff(s.exp(3*A)*Cp,t)+s.diff(s.exp(3*A)*Cp,z),False)
check('reject density e^(2A)',s.diff(s.exp(2*A)*Cp,t)-s.diff(s.exp(2*A)*Cp,z),False)
check('reject Hamiltonian scalar kinetic coefficient',s.diff(H+pt**2/2,t)-2*s.diff(M,z)-6*Az*M+3*At*(H+pt**2/2),False)
check('reject momentum coupling 3 Az H',s.diff(M,t)-s.diff(H,z)/2-3*Az*H+3*At*M,False)

# The junction identity and its velocity derivative, independent of bulk PDE.
sig, sig1, sig2, eb, qB, qp = s.symbols('sigma sigma1 sigma2 eB B_t phi_t')
gA=eb*sig/6; gP=-eb*sig1/2
gtA=eb*(sig*qB+sig1*qp)/6
gtP=-eb*(sig1*qB+sig2*qp)/2
shellM=s.expand(-3*gtA+3*gA*qB-qp*gP)
check('junction implies shell momentum identity',shellM)
check('reject scalar junction sign reversal',-3*gtA+3*gA*qB+qp*gP,False)
check('zero velocity-slope is generally invalid',3*gA*qB-qp*gP,False)

# Identify defects when numerical dissipation is placed in velocity equations.
# a_t=pa,b_t=pb,phi_t=pf and pa_t=EA+QA, pb_t=EB+QB, pf_t=EP+QP.
QA,QB,QP, QAz=s.symbols('Q_A Q_B Q_phi Q_Az')
HtQ=(12*At+6*Bt)*QA+6*At*QB-2*pt*QP
MtQ=-3*QAz+3*(Bz-Az)*QA+3*Az*QB-pz*QP
# These follow directly by partial differentiation in the independent jets.
va,vb,vp,vaz=s.symbols('v_A v_B v_phi v_Az')
Hj=H.xreplace({At:va,Bt:vb,pt:vp})
Mj=M.xreplace({At:va,Bt:vb,pt:vp,s.diff(A,t,z):vaz})
check('momentum-only dissipation source in H transport',
      (s.diff(Hj,va)*QA+s.diff(Hj,vb)*QB+s.diff(Hj,vp)*QP).subs({va:At,vb:Bt,vp:pt})-HtQ)
check('momentum-only dissipation source in M transport',
      (s.diff(Mj,va)*QA+s.diff(Mj,vb)*QB+s.diff(Mj,vp)*QP+s.diff(Mj,vaz)*QAz).subs({va:At,vb:Bt,vp:pt,vaz:s.diff(A,t,z)})-MtQ)

data={'status':'PASS','checks':checks,'positive_controls':sum(not x['negative_control'] for x in checks),
      'negative_controls':sum(x['negative_control'] for x in checks),'sympy_version':s.__version__,
      'definition':'H=2(G_tt-T_tt), M=G_tz-T_tz; evolution has E_zz=H/2, E_ij=0',
      'outgoing_density':'exp(3A)*(H+2M), transported toward decreasing z at unit coordinate speed',
      'incoming_density':'exp(3A)*(H-2M), transported toward increasing z at unit coordinate speed',
      'dissipation_scope':'Identities are for continuum nondissipative PDE; specified momentum-only Q adds source terms in report.',
      'self_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('CONSTRAINT_TRANSPORT_CHECKS.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'status':data['status'],'positive_controls':data['positive_controls'],'negative_controls':data['negative_controls']}))
