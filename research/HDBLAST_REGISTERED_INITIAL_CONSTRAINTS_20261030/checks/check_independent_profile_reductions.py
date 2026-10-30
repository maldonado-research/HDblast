#!/usr/bin/env python3
"""Independent exact abstract-profile audit; no imports of project producers."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys
import tempfile

NAMES=('ir','a','p','h','Vpp','u','v','vy','E','f','fy','fyy','F','Fy','G')
N=len(NAMES)
class P:
    def __init__(self, terms=0):
        if isinstance(terms,P): self.t=terms.t.copy()
        elif isinstance(terms,dict): self.t={k:F(v) for k,v in terms.items() if v}
        else: self.t={} if not terms else {(0,)*N:F(terms)}
    @staticmethod
    def var(name):
        e=[0]*N;e[NAMES.index(name)]=1
        return P({tuple(e):F(1)})
    def __add__(self,b):
        b=P(b);t=self.t.copy()
        for k,v in b.t.items():t[k]=t.get(k,F(0))+v
        return P(t)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.t.items()})
    def __sub__(self,b):return self+-P(b)
    def __rsub__(self,b):return P(b)+-self
    def __mul__(self,b):
        b=P(b);t={}
        for ka,va in self.t.items():
            for kb,vb in b.t.items():
                k=tuple(a+b for a,b in zip(ka,kb));t[k]=t.get(k,F(0))+va*vb
        return P(t)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            if len(self.t)!=1:raise ValueError('non-monomial inverse')
            k,v=next(iter(self.t.items()));return P({tuple(e*n for e in k):v**n})
        r=P(1)
        for _ in range(n):r=r*self
        return r
    def __truediv__(self,b):return self*P(b)**-1
    def dy(self,rules):
        t=P()
        for k,v in self.t.items():
            for j,e in enumerate(k):
                if e:
                    kk=list(k);kk[j]-=1
                    t=t+P({tuple(kk):v*e})*rules[NAMES[j]]
        return t
    def at_packet_zero(self):
        return P({k:v for k,v in self.t.items() if not any(k[NAMES.index(s)] for s in ('u','v','vy'))})
    def encoded(self):return [[list(k),str(v)] for k,v in sorted(self.t.items())]

PINS={
 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md':'2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed',
 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md':'bc3ae85c3166575ed15b017468a3e9668615e15dfaca34bb2a13afe27ea28f92',
 'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/full23_initial.py':'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1',
 'research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/numeric/evaluate_initial_profiles.py':'d9cb8529e47b9590ec8dc94f2f4484b06dbeb0ddee6c16015feccde3e4440240',
}
MUTATIONS={
 'linear_lapse_laplacian_omitted':'first_physical_hamiltonian',
 'wrong_background_ay':'first_physical_hamiltonian',
 'half_hamiltonian_hessian':'registered_aa_equals_physical_A',
 'omit_angular_gradient_source':'registered_dd_equals_physical_D',
 'wrong_hamiltonian_denominator':'hamiltonian_sector_completion_l0',
 'scalar_compensator_sign':'hamiltonian_differential_right_inverse_l0',
 'spacetime_affine_radial_source':'p_registered_Pi_equals_Gmom_minus_Fmom_y',
 'momentum_compensator_sign':'radial_momentum_generic_completion',
}
def auth(root):
    for rel,want in PINS.items():
        p=root/rel
        if not p.is_file() or p.is_symlink() or hashlib.sha256(p.read_bytes()).hexdigest()!=want:raise ValueError('SOURCE_AUTH:'+rel)

def calculate(mutation):
    V={s:P.var(s) for s in NAMES};ir,a,p,h,Vpp,u,v,vy,E,f,fy,fyy,Fv,Fy,Gv=(V[s] for s in NAMES)
    ay=-ir**2-p**2/3
    if mutation=='wrong_background_ay':ay=-ir**2+p**2/3
    py=p*(h-4*a);uy=-2*a*u-p*v/3
    rules={'ir':-a*ir,'a':ay,'p':py,'h':Vpp-h*h+4*a*h,
                  'u':uy,'v':vy,'E':ir*E,'f':fy,'fy':fyy}
    dy=lambda q:P(q).dy(rules)
    ut=E*u;vt=E*(v-3*ir*u/p);uyt=dy(ut);uyy=dy(uy)
    checks=[]
    def zero(name,r):
        r=P(r)
        if r.t:raise ValueError('NONZERO:'+name+':'+json.dumps(r.encoded(),separators=(',',':')))
        checks.append(name)
    zero('first_radial_reconstruction',p*v+3*(uy+2*a*u))
    zero('first_velocity_reconstruction',p*vt+3*(uyt+2*a*ut))
    H1=-6*uyy-36*a*uy-(24*ay+48*a*a+4*p*p+12*ir*ir)*u-2*p*vy-2*p*h*v
    if mutation=='linear_lapse_laplacian_omitted':H1=H1+32*ir*ir*u
    zero('first_physical_hamiltonian',H1)
    zero('first_physical_radial_momentum',ir*(-3*(uyt+2*a*ut)-p*vt))
    # Physical-Cauchy-affine unsimplified Hessian, inherited analytic eq. (6).
    AH=-24*uy*uy-288*a*u*uy-24*u*uyy-(192*ay+384*a*a+32*p*p)*u*u-ir*ir*(336*u*u+12*ut*ut+2*vt*vt)-2*vy*vy-16*p*u*vy-2*Vpp*v*v
    DH=ir*ir*(36*u*u-2*v*v)
    Ar=-2*vy*vy-8*p*u*vy+8*p*(h+2*a)*u*v+(16*p*p-192*ir*ir)*u*u-(2*Vpp+F(8,3)*p*p)*v*v-ir*ir*(12*ut*ut+2*vt*vt)
    Dr=ir*ir*(36*u*u-2*v*v)
    if mutation=='half_hamiltonian_hessian':AH=AH/2
    if mutation=='omit_angular_gradient_source':DH=P()
    zero('registered_aa_equals_physical_A',Ar-AH)
    zero('registered_dd_equals_physical_D',Dr-DH)
    Fmom=ir*(18*u*ut-v*vt)
    Gmom=ir*(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*vy)
    if mutation=='spacetime_affine_radial_source':Gmom=Gmom+6*ir*uy*ut
    Kreg=E*ir/3*(18*u*u-v*v+3*ir*u*v/p)
    Pireg=E*ir/p*(3*ir*u*vy/p+(3*ir*(h-ir)/p+10*p)*u*v+(30*a-6*ir)*u*u+(2*ir-a)*v*v)
    zero('three_registered_K_equals_physical_Fmom',3*Kreg-Fmom)
    zero('p_registered_Pi_equals_Gmom_minus_Fmom_y',p*Pireg-Gmom+dy(Fmom))
    # Generic f,fy,fyy checks do not use its registered formula to cancel radial derivatives.
    zeta=-3*fy/(2*p)
    if mutation=='scalar_compensator_sign':zeta=-zeta
    for lam in (0,8,24):
        DHlin=-3*fyy-12*a*fy+2*(lam-3)*ir*ir*f-2*p*dy(zeta)-2*p*h*zeta
        zero('hamiltonian_differential_right_inverse_l'+str(lam),DHlin-2*(lam-3)*ir*ir*f)
    sectors=[Ar/6+F(4,3)*Dr,Ar/2+2*Dr,Ar-4*Dr]
    for j,(lam,den) in enumerate(zip((0,8,24),(-6,10,42))):
        used=den+1 if mutation=='wrong_hamiltonian_denominator' and j==0 else den
        fj=-sectors[j]/(ir*ir*used)
        zero('hamiltonian_sector_completion_l'+str(lam),2*(lam-3)*ir*ir*fj+sectors[j])
    q=Fv/3;Pi=(Gv-Fy)/p
    if mutation=='momentum_compensator_sign':Pi=-Pi
    zero('angular_momentum_generic_completion',-3*q+Fv)
    zero('radial_momentum_generic_completion',-Fy-p*Pi+Gv)
    # Full-metric lowering is irrelevant only AFTER C1=0; no QC use is inferred here.
    for label,item in [('AH',AH),('DH',DH),('Fmom',Fmom),('Gmom',Gmom),('Kreg',Kreg),('Pireg',Pireg)]:
        zero('packet_vanishing_ideal_'+label,item.at_packet_zero())
    return checks

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--mutation',choices=MUTATIONS);ap.add_argument('--run-controls',action='store_true');args=ap.parse_args()
    if args.output.exists():raise ValueError('REFUSE_OVERWRITE')
    auth(args.repo_root);checks=calculate(args.mutation);auth(args.repo_root)
    rejected=[]
    if args.run_controls:
        with tempfile.TemporaryDirectory() as d:
            for name in MUTATIONS:
                out=Path(d)/(name+'.json');cmd=[sys.executable,'-I','-S','-B']
                if sys.flags.optimize:cmd+=['-O']
                cmd += [str(Path(__file__).resolve()),'--repo-root',str(args.repo_root),'--mutation',name,'--output',str(out)]
                r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
                if r.returncode==0 or out.exists() or ('NONZERO:'+MUTATIONS[name]+':') not in r.stderr:raise ValueError('CONTROL_NOT_REJECTED:'+name)
                rejected.append(name)
    report={'status':'PASS_INDEPENDENT_ALL_POINT_ABSTRACT_PROFILE_REDUCTIONS','schema':'independent-abstract-profile-reductions-v1','scope':'Exact identities in a free Laurent polynomial ring after the stated actual background and first-profile differential identities; inherited physical Hessian/ADM derivations remain analytic premises. No sampled point or interval selection is used.','source_SHA256':PINS,'checks':checks,'check_count':len(checks),'mutations_rejected':rejected,'all_residual_coefficient_counts':0,'independent_from_producer_imports':True,'full_registered_46_entry_binding_checked_here':False,'global_subsidiary_completion_proved_here':False,'finite_amplitude_or_cosmological_realization':False}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
if __name__=='__main__':main()
