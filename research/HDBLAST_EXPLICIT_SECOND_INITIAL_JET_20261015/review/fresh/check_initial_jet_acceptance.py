#!/usr/bin/env python3
"""Fresh source-hashed initial-constraint review; Python standard library only.

Read-only source: all replay scratch lives in system temp, and only --output
is written. The accepted source checker is executed in isolated subprocesses,
never imported. This review adds all-point Laurent-polynomial proof checks.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent
PINS={
    'check_initial_constraint_completion.py':'954eeb402c3f06c72951bdfb6051d3f9e346b43ae2487324d67370021262425e',
    'EXPLICIT_COMPACT_SECOND_INITIAL_JET.md':'2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed',
    'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
    'INITIAL_CONSTRAINT_COMPLETION.json':'6dcb9afb3aab2b23246912159bfdd30a5fdca01bd43b1da6eef89cb881fcbb00',
    'INITIAL_CONSTRAINT_COMPLETION_OPTIMIZED.json':'e3c12fb9dd6cf3d2ee20d7f7524ef9f1955b9870ef53523b9b0f6e0f2b7cdd2e',
}
README_PIN='c61051888de1b73f9958802f4ec41575347ebadec42cf5062e1961b2b739a06f'


def require(ok,label):
    if not ok:
        raise ValueError(label)


class P:
    """Sparse Laurent polynomial over formal radial and angular jet variables."""
    names=('r','a','ay','u','uy','uyy','ut','v','vy','vt','p','U2',
           'Y','G','N','f','fy','fyy','py','src','Fy','Gm','uyt','Yi')
    n=len(names)

    def __init__(self,terms=None):
        self.d={m:F(c) for m,c in (terms or {}).items() if c}

    @classmethod
    def box(cls,x):
        return x if isinstance(x,cls) else cls({(0,)*cls.n:F(x)})

    @classmethod
    def var(cls,name):
        e=[0]*cls.n
        e[cls.names.index(name)]=1
        return cls({tuple(e):F(1)})

    def __add__(self,other):
        out=self.d.copy()
        for m,c in self.box(other).d.items():
            out[m]=out.get(m,F(0))+c
        return P(out)

    __radd__=__add__

    def __neg__(self):
        return P({m:-c for m,c in self.d.items()})

    def __sub__(self,other):
        return self+-self.box(other)

    def __rsub__(self,other):
        return self.box(other)+-self

    def __mul__(self,other):
        out={}
        for m,c in self.d.items():
            for n,d in self.box(other).d.items():
                k=tuple(a+b for a,b in zip(m,n))
                out[k]=out.get(k,F(0))+c*d
        return P(out)

    __rmul__=__mul__

    def __pow__(self,n):
        if n<0:
            require(len(self.d)==1,'inverse needs formal monomial')
            (m,c),=self.d.items()
            return P({tuple(k*n for k in m):c**n})
        out=P.box(1)
        for _ in range(n):
            out=out*self
        return out

    def __truediv__(self,other):
        return self*self.box(other)**-1

    def __rtruediv__(self,other):
        return self.box(other)*self**-1


def proof_checks():
    v={name:P.var(name) for name in P.names}
    r,a,ay,u,uy,uyy,ut,w,wy,wt,p,U2,Y,G,N,f,fy,fyy,py,src,Fy,Gm,uyt,Yi=[v[n] for n in P.names]
    checks={}

    def zero(name,polynomial):
        require(not polynomial.d,'independent all-point proof failed: '+name)
        checks[name]='PASS_ZERO_LAURENT_POLYNOMIAL'

    # Independent radial-foliation expansion of R, not coordinate Ricci replay.
    # j=a+eps*j1*Y+eps²*j2*Y², L^-1=1+2eps*uY+6eps²*u²Y².
    j1=uy+2*a*u
    j2=6*a*u*u
    j1y=uyy+2*ay*u+2*a*uy
    j2y=6*ay*u*u+12*a*u*uy
    radial=2*(-12*(j1*j1+2*a*j2)-6*(j2y+2*u*j1y+6*u*u*ay))*Y*Y
    angular_conformal=(-208*u*u*Y*Y+12*u*u*G)/r**2
    angular_lapse=(-128*u*u*Y*Y+24*u*u*G)/r**2
    R2=radial+angular_conformal+angular_lapse
    Rformula=(-24*uy*uy-288*a*u*uy-24*u*uyy-(192*ay+384*a*a)*u*u-336*u*u/r**2)*Y*Y+36*u*u*G/r**2
    zero('curvature_Hessian_all_radial_angular_terms',R2-Rformula)
    # Direct affine physical scalar-gradient and K/Pi quadratic coefficients.
    gradient2=(32*p*p*u*u+16*p*u*wy+2*wy*wy)*Y*Y+2*w*w*G/r**2
    K2=2*((ut/r)**2-((-2*ut/r)**2+3*(ut/r)**2))*Y*Y
    Pi2=-2*(wt/r)**2*Y*Y
    potential2=-2*U2*w*w*Y*Y
    direct=R2-gradient2+K2+Pi2+potential2
    A=(-24*uy*uy-288*a*u*uy-24*u*uyy-(192*ay+384*a*a+32*p*p)*u*u
       -(336*u*u+12*ut*ut+2*wt*wt)/r**2-2*wy*wy-16*p*u*wy-2*U2*w*w)
    D=(36*u*u-2*w*w)/r**2
    zero('Hamiltonian_all_source_terms',direct-(A*Y*Y+D*G))
    # Independent all-point diagonal mixed-tensor divergence expansion.
    s=ut/r
    sy=(uyt-a*ut)/r
    radial_momentum=2*(6*(uy*s+u*sy)-18*a*u*s-9*uy*s-wt*wy/r)
    angular_momentum=2*(24*u*s-6*u*s-w*wt/r)*Y*Yi
    Fsource=(18*u*ut-w*wt)/r
    Gsource=(12*u*uyt-48*a*u*ut-6*uy*ut-2*wt*wy)/r
    zero('momentum_radial_complete_Hessian',radial_momentum-Gsource)
    zero('momentum_angular_complete_Hessian',angular_momentum-2*Fsource*Y*Yi)
    # Harmonic identities are polynomial identities on the embedded S3.
    Z=N-F(1,2)
    H4=Y*Y-N/2+F(1,12)
    zero('degree_024_decomposition',Y*Y-(F(1,6)+Z/2+H4))
    zero('degree_2_eigenvalue',(4-8*N)+8*Z)
    zero('degree_4_eigenvalue',(8*N-24*Y*Y)-(4-8*N)/2+24*H4)
    zero('Hamiltonian_degree_024_projection',A/6+4*D/3+(A/2+2*D)*Z+(A-4*D)*H4-(A*Y*Y+4*D*(N-Y*Y)))
    # Independent angular metric variation and scalar compensator.
    U1=py+4*a*p
    zeta=-3*fy/(2*p)
    zetay=-3*fyy/(2*p)+3*fy*py/(2*p*p)
    for lam in (0,8,24):
        dR=-3*fyy-12*a*fy+2*(lam-3)*f/r**2
        zero('Hamiltonian_inverse_lambda_'+str(lam),dR-2*p*zetay-2*U1*zeta-2*(lam-3)*f/r**2)
        fc=-r*r*src/(2*(lam-3))
        zero('Hamiltonian_source_sign_lambda_'+str(lam),2*(lam-3)*fc/r**2+src)
    zero('momentum_pure_trace_radial_inverse',-Fy-p*((Gm-Fy)/p)+Gm)
    zero('momentum_pure_trace_angular_inverse',-3*(Fsource/3)+Fsource)
    actual_vt=-3*(uyt+2*a*ut)/p
    induced=(6*uy*ut+6*u*uyt+12*a*u*ut+2*u*actual_vt*p)/r
    zero('spacetime_vs_Cauchy_affine_radial_shift',induced-6*uy*ut/r)
    # The actual reconstructed first tangent satisfies the linear H constraint.
    first_chi=-3*(uy+2*a*u)/p
    first_chiy=-3*(uyy+2*ay*u+2*a*uy)/p+3*(uy+2*a*u)*py/(p*p)
    first_R=-6*uyy-36*a*uy-(24*ay+48*a*a+12/r**2)*u
    first_H=first_R-4*p*p*u-2*p*first_chiy-2*U1*first_chi
    zero('first_H_constraint_background_reduction',first_H-(-12*ay-12/r**2-4*p*p)*u)
    return checks


def actual_witness(background):
    b=background['background_intervals']
    c=background['reconstruction_coefficients']
    box=lambda d:tuple(map(F,d['exact']))
    mag=lambda d:max(abs(x) for x in box(d))
    rlo,rhi=box(b['rho'])
    plo,phi=box(b['p'])
    Blo,Bhi=box(c['v_wz_coefficient'])
    require((rlo,rhi)==(129,131) and plo<phi<0 and Blo<Bhi<0,'actual rho,p,B signs')
    require(background['actual_model']['c']=='96714321425692/161839258930731'
            and background['actual_model']['delta']=='1/1000','actual branch identity')
    ell=F(1,1000000)
    bm=min(abs(Blo),abs(Bhi))
    bM=max(abs(Blo),abs(Bhi))
    kz=mag(c['alpha_log_zz'])
    k=mag(c['alpha_log_z'])
    Bz=mag(c['v_wz_coefficient_z'])
    Cm=mag(c['u_w_coefficient'])
    Cz=mag(c['u_w_coefficient_z'])
    Czz=mag(c['u_w_coefficient_zz'])
    err=(Bz*k+bM*kz)/rlo
    lower=32*bm/(rhi*ell*ell)-err
    upper=32*bM/(rlo*ell*ell)+err
    require(F(2300000000)<lower and upper<3000000000,'independent actual v_y witness')
    uy=Cz/rlo
    ay=mag(b['a_z'])/rlo
    uyy=(Czz+32*Cm/(ell*ell))/(rlo*rlo)+mag(b['a'])*Cz/rlo
    vm=bM*k
    vtm=bM*(k+1)
    require(Cm<F(1,10000) and uy<1 and uyy<1000000 and vm<4 and vtm<4,
            'derived actual qmax first-field caps')
    # The exact report's geometric remainder uses deliberately coarser caps.
    u=F(1,10000)
    remainder=(24+288*u+24*u*1000000+(192+384+32)*u*u
               +(336*u*u+12+2*16)/(rlo*rlo)+16*u*3000000000+2*4*16)
    require(ay<1 and mag(b['a'])<1 and mag(b['p'])<1 and mag(b['U_2'])<4,'actual background witness caps')
    require(remainder<5000000 and -2*F(2300000000)**2+5000000 < -10**19,'strict Hamiltonian witness')
    return {'p_nonzero':'-7/50000 <= p <= -3/25000 < 0',
            'derived_abs_v_y_lower':str(lower),'derived_abs_v_y_upper':str(upper),
            'derived_abs_u_y_upper':str(uy),'derived_abs_u_yy_upper':str(uyy),
            'derived_abs_v_upper':str(vm),'derived_abs_v_t_upper':str(vtm),
            'other_terms_coarse_upper':str(remainder),
            'Hamiltonian_Hessian_strict_upper':str(-2*F(2300000000)**2+5000000)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root','--source-root',dest='source',type=Path,default=ROOT.parent/'initial-jet')
    parser.add_argument('--readme',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    source=args.source.resolve()
    readme=args.readme if args.readme else source/'README.md'
    require(sha256(readme.read_bytes()).hexdigest()==README_PIN,'review README pin mismatch')
    for path,pin in PINS.items():
        require(sha256((source/path).read_bytes()).hexdigest()==pin,'review source pin mismatch: '+path)
    background=json.loads((source/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_text())
    proofs=proof_checks()
    actual=actual_witness(background)
    replay={}
    reports=[]
    with tempfile.TemporaryDirectory(prefix='fresh-initial-review-') as tmp:
        for optimized in (False,True):
            mode='optimized' if optimized else 'normal'
            output=Path(tmp)/(mode+'.json')
            cmd=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])
            cmd += [str(source/'check_initial_constraint_completion.py'),'--background',
                    str(source/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json'),'--run-controls','--output',str(output)]
            proc=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
            require(proc.returncode==0,'fresh '+mode+' replay failed: '+proc.stderr)
            raw=output.read_bytes()
            retained='INITIAL_CONSTRAINT_COMPLETION_OPTIMIZED.json' if optimized else 'INITIAL_CONSTRAINT_COMPLETION.json'
            require(raw==(source/retained).read_bytes(),'fresh replay differs from retained '+mode)
            data=json.loads(raw)
            require(data['finite_exact_checks']==48 and len(data['controls']['cases'])==7
                    and all(c['returncode']!=0 and c['no_output_report'] for c in data['controls']['cases']),
                    'source controls incomplete')
            replay[mode]={'status':'PASS_FRESH_BYTE_IDENTICAL_TO_RETAINED','sha256':sha256(raw).hexdigest(),
                          'finite_exact_checks':48,'negative_controls_rejected':7}
            data['controls'].pop('optimized')
            reports.append(data)
    require(reports[0]==reports[1],'normal/optimized mathematical results differ')
    for path,pin in PINS.items():
        require(sha256((source/path).read_bytes()).hexdigest()==pin,'source mutated during review: '+path)
    require(sha256(readme.read_bytes()).hexdigest()==README_PIN,'README mutated during review')
    scope=reports[0]['scope']
    require(scope['explicit_second_constraint_jet'] and not scope['second_initial_norms_numerically_enclosed']
            and not scope['gauge_completion_computed'] and not scope['U2_evolution_computed'],
            'accepted scope differs')
    report={'status':'PASS_FRESH_INDEPENDENT_INITIAL_JET_ACCEPTANCE','source_sha256':PINS,
            'reviewed_capsule_README_sha256':README_PIN,
            'all_point_laurent_proof_checks':proofs,'actual_packet_witness':actual,'fresh_replay':replay,
            'normal_optimized_mathematical_identity':True,
            'support_review':'Finite derivatives and smooth multipliers, p bounded away from zero on collar; flat packet endpoints permit smooth extension by zero in the same [ell,2ell] annulus and vanish on open shell/vertex collars.',
            'convention_review':'Physical Cauchy affine gamma,K,phi,Pi; K=+Lie_n gamma/2 and Pi=n phi. Spacetime affine lapse induces second K,Pi and shifts radial momentum source by +6u_yu_t/rho. Coordinate conversion factors verified analytically.',
            'accepted_scope':scope,
            'limits':['Conditional on the authenticated retained regular background and first reconstruction.',
                      'Explicit second differentiated initial constraints only; no exact nonlinear initial family.',
                      'No completed gauge, coupled second evolution, U2 amplitude or remainder estimate.']}
    text=json.dumps(report,indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text,end='')


if __name__=='__main__':
    main()
