"""Independent exact source-free bare-bilinear and Fourier-primitive proof."""
from __future__ import annotations
import json
from pathlib import Path
import sys
import sympy as s
from exact_binary import require


def run():
    k, L, eps = s.symbols('k L eps', real=True, nonzero=True)
    ur, ui, wr, wi = s.symbols('ur ui wr wi', real=True)
    C0, Ck2, Zr, OZi, K2Zr = s.symbols('C0 Ck2 Zr OZi K2Zr', real=True)
    R = ((2*k*k+3*L*L)*ur-k*wi-L*wr)/(2*k*eps)
    P = ((s.Rational(2,3)*k*k-L*L)*ur-k*wi-L*wr)/(2*k*eps)
    vector = {L:L*L,ur:wr,ui:wi,wr:-2*k*wi,wi:2*k*wr}
    derivative_R = sum(s.diff(R, x)*dx for x,dx in vector.items())
    identities = {}
    def identity(name, value):
        require(s.simplify(s.expand(value)) == 0, 'Exact independent identity failed: '+name)
        identities[name] = 'PASS'
    identity('unrestricted_bare_Ward_vector_field', derivative_R-L*(R-3*P))
    cr, ci, dr, di, Er, Ei = s.symbols('cr ci dr di Er Ei',real=True)
    free_ur = cr+dr*Er-di*Ei
    free_wr = -2*k*(di*Er+dr*Ei)
    free_wi = 2*k*(dr*Er-di*Ei)
    zreal, omegazimag = dr*Er-di*Ei,2*k*(di*Er+dr*Ei)
    inventory_R = (2*k*k*cr+3*L*L*(cr+zreal)+L*omegazimag)/(2*k*eps)
    inventory_P = (s.Rational(2,3)*k*k*cr-L*L*(cr+zreal)-s.Rational(4,3)*k*k*zreal+L*omegazimag)/(2*k*eps)
    identity('free_R_grouped_without_Re_c_projection', R.subs({ur:free_ur,wr:free_wr,wi:free_wi})-inventory_R)
    identity('free_P_grouped_without_Re_c_projection', P.subs({ur:free_ur,wr:free_wr,wi:free_wi})-inventory_P)
    moment_R = 2*Ck2+3*L*L*(C0+Zr)+L*OZi
    moment_P = s.Rational(2,3)*Ck2-L*L*(C0+Zr)-s.Rational(4,3)*K2Zr+L*OZi
    moment_F = 6*L**3*(C0+Zr)+4*L*K2Zr-2*L*L*OZi
    moment_G = 3*L*L*(C0+Zr)+L*OZi
    identity('independent_moment_F_inventory',moment_F-L*(moment_R-3*moment_P))
    inventory_vector = {L:L*L,Zr:-OZi,OZi:4*K2Zr}
    dG = sum(s.diff(moment_G,x)*dx for x,dx in inventory_vector.items())
    identity('independent_F_antiderivative',dG-moment_F)
    identity('analytic_free_increment_equals_primitive',moment_R-moment_G-2*Ck2)
    dur, dui, dwr, dwi = s.symbols('dur dui dwr dwi',real=True)
    flow_projection = ((2*k*k+3*L*L)*dur-k*dwi-L*dwr)/(2*k*eps)
    identity('unrestricted_endpoint_defect_density_projection',
             R.subs({ur:ur+dur,ui:ui+dui,wr:wr+dwr,wi:wi+dwi})-R-flow_projection)
    # Real complex products in the canonical direct phase antiderivative.
    La,Lb = s.symbols('La Lb',real=True)
    direct = 3*cr*(Lb**2-La**2)+(dr*(3*(Lb**2*Er-La**2)+2*k*Lb*Ei)
                                      -di*(3*Lb**2*Ei-2*k*(Lb*Er-La)))
    endpoint = 3*Lb**2*(cr+dr*Er-di*Ei)+2*k*Lb*(di*Er+dr*Ei) \
        -(3*La**2*(cr+dr)+2*k*La*di)
    identity('direct_phase_antiderivative_expansion',direct-endpoint)
    mutations = {'drop_Re_c':moment_F-(moment_F-6*L**3*C0),
                 'wrong_Omega_imag_sign':moment_F-(6*L**3*(C0+Zr)+4*L*K2Zr+2*L*L*OZi),
                 'omit_pressure_k2_moment':moment_F-L*(moment_R-3*(moment_P+s.Rational(4,3)*K2Zr))}
    for name, residual in mutations.items():
        require(s.simplify(residual) != 0,'Symbolic mutation undetected: '+name)
    return {'status':'PASS_EXACT_SYMBOLIC_ONLY','physical_values_loaded':False,
            'identities_passed':len(identities),'identities':identities,
            'mutations_rejected':list(mutations),'optimized':not __debug__,
            'sympy_version':s.__version__}


if __name__=='__main__':
    result=run()
    path=Path(sys.argv[1]) if len(sys.argv)>1 else None
    if path:
        require(not path.exists(),'Proof receipt must be fresh')
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
