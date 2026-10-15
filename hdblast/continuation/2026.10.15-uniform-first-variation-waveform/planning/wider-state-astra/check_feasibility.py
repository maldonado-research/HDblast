#!/usr/bin/env python3
"""New exact-rational wider-state feasibility audit; imports no frozen producer."""
from dataclasses import dataclass
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
INPUTS={
 'inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json':
 '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
 'inputs/previous/archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md':
 '0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c',
}


def require(condition,reason):
    if not condition:
        raise ValueError(reason)


@dataclass(frozen=True)
class I:
    lo:F
    hi:F
    def __init__(self,lo,hi=None):
        object.__setattr__(self,'lo',F(lo))
        object.__setattr__(self,'hi',F(lo if hi is None else hi))
        require(self.lo<=self.hi,'reversed interval')
    def __add__(self,b):
        b=as_i(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-as_i(b)
    def __rsub__(self,b):return as_i(b)+-self
    def __mul__(self,b):
        b=as_i(b);ends=[x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return I(min(ends),max(ends))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=as_i(b);require(not b.lo<=0<=b.hi,'division by zero-containing interval')
        return self*I(1/b.hi,1/b.lo)
    def __rtruediv__(self,b):return as_i(b)/self
    def __pow__(self,n):
        require(isinstance(n,int) and n>=0,'invalid power')
        answer=I(1)
        for _ in range(n):answer*=self
        return answer
    def pair(self):return [str(self.lo),str(self.hi)]


def as_i(x):return x if isinstance(x,I) else I(x)


def bulk_potential(eta):
    W=F(1,3)+eta**2+eta**3/3
    Wp=2*eta+eta**2;Wpp=2+2*eta
    Up=Wp*(Wpp-F(4,3)*W)
    Upp=Wpp**2+2*Wp-F(4,3)*(Wp**2+W*Wpp)
    return Up,Upp


def majorant(ell,m,M,B,gamma):
    T=3*ell;U=1+T;r=M*T*T/2;s=B*T*T/(2*m)
    den=(1-r)*(1-s)-T*r/m
    require(r<1 and s<1 and den>0,'majorant positivity failed')
    A=U*(1-s+T/m)/den;V=U/(m*den)
    eta=ell*(M*T*A+(1+gamma*T)*V)
    vtail=ell/m**2+ell/m+r*A/m+s*V
    tail=ell*(M*T*A+vtail+gamma*T*V)
    require((1-r)*A-T*V==U and (1-s)*V-r*A/m==U/m,'majorant row mismatch')
    return dict(ell=ell,T=T,m=m,M=M,B=B,gamma=gamma,r=r,s=s,determinant=den,
                A=A,V=V,eta=eta,V_tail=vtail,eta_tail=tail)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,default=ROOT,
                        help='Extracted frozen waveform capsule containing inputs/previous/archive')
    args=parser.parse_args()
    source_root=args.source_root.resolve()
    checks=[];negatives=[]
    def check(name,condition):require(condition,name);checks.append(name)
    def reject(name,fn):
        try:fn()
        except ValueError:negatives.append(name)
        else:raise ValueError('negative control accepted: '+name)
    for name,digest in INPUTS.items():
        check('source_pin:'+name,sha256((source_root/name).read_bytes()).hexdigest()==digest)
    source=json.loads((source_root/'inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
    check('same_background_parameters',source['parameters_exact']=={
        'delta':'1/1000','c':'96714321425692/161839258930731'})
    boxes=source['collar']['first_exit_inclusions']
    outer={name:I(*v['outer_box']['exact']) for name,v in boxes.items()}
    initial={name:I(*v['initial_box']['exact']) for name,v in boxes.items()}
    archived_rhs={name:I(*v['rhs_enclosure']['exact']) for name,v in boxes.items()}
    rho,eta,p,a=(outer[name] for name in ['rho','eta','p','a'])
    Up,Upp=bulk_potential(eta)
    rhs={'rho':-rho**2*a,'eta':-rho*p,'p':-rho*(Up-4*a*p),'a':1/rho+rho*p**2/3}
    for name in rhs:
        check('independent_ODE_interval:'+name,rhs[name]==archived_rhs[name])
        check('initial_strict_inclusion:'+name,outer[name].lo<initial[name].lo<=initial[name].hi<outer[name].hi)

    thresholds={}
    for name,Fbox in rhs.items():
        if Fbox.lo<0:thresholds[name+'_lower']=(initial[name].lo-outer[name].lo)/(-Fbox.lo)
        if Fbox.hi>0:thresholds[name+'_upper']=(outer[name].hi-initial[name].hi)/Fbox.hi
    limit=min(thresholds.values())
    check('limiting_eta_upper_face',limit==thresholds['eta_upper']
          ==F(5669686020285051,123672167032900272500))
    ell=F(1,10**6);Z=8*ell
    def margins(depth):
        return {name:(initial[name].lo+min(F(0),depth*rhs[name].lo)-outer[name].lo,
                      outer[name].hi-initial[name].hi-max(F(0),depth*rhs[name].hi))
                for name in rhs}
    margins_Z=margins(Z)
    check('wider_depth_below_strict_threshold',Z<limit)
    check('all_wider_first_exit_margins_positive',all(v>0 for pair in margins_Z.values() for v in pair))
    check('strict_threshold_is_not_attained_by_this_test',margins(limit)['eta'][1]==0)
    check('concrete_rounded_longer_collar',F(1,25000)<limit)

    h=Up/p
    potential=rho**2*(2*h**2-10*a*h+F(35,4)*a**2-Upp-p**2/6)+I(F(-35,4),0)
    sb=source['boxes']
    shell_eta,shell_rho,shell_p,shell_a=(I(*sb[k]['exact']) for k in ['eta_b','rho_b','p_b','a_b'])
    shell_Up,_=bulk_potential(shell_eta)
    shell_h=shell_Up/shell_p
    shell_B=shell_h-4*shell_a+2+2*shell_eta
    mu=1/(shell_rho*shell_B)
    gamma=shell_rho*(2*shell_h-5*shell_a)/2
    beta=mu*I(F(-25,4),F(5,2))+gamma
    check('same_mu_interval_rederived',mu==I(*sb['mu']['exact']))
    check('coarse_coefficients_hold',mu.lo>F(1,1000) and
          max(abs(potential.lo),abs(potential.hi))<10**6 and
          max(abs(beta.lo),abs(beta.hi))<1000 and max(abs(gamma.lo),abs(gamma.hi))<250)
    check('improved_rounded_coefficients_hold',mu.lo>F(1,500) and
          max(abs(potential.lo),abs(potential.hi))<102000 and
          max(abs(beta.lo),abs(beta.hi))<250 and max(abs(gamma.lo),abs(gamma.hi))<250)

    coarse=majorant(ell,F(1,1000),F(10**6),F(1000),F(250))
    improved=majorant(ell,F(1,500),F(102000),F(250),F(250))
    check('coarse_eta_exact',coarse['eta']==F(2007523995553919,1999981973040500000))
    check('coarse_tail_exact',coarse['eta_tail']==F(9538004522513419,1999981973040500000000))
    check('coarse_eta_rounded',coarse['eta']<F(251,250000))
    check('coarse_tail_rounded',coarse['eta_tail']<F(477,10**8))
    check('coarse_L2_rounded',coarse['eta']**2+coarse['eta_tail']**2<F('0.0010038')**2)
    check('improved_eta_rounded',improved['eta']<F('0.000500684'))
    check('improved_tail_rounded',improved['eta_tail']<F('0.000000932473'))
    check('improved_L2_rounded',improved['eta']**2+improved['eta_tail']**2<F('0.000500685')**2)
    check('strict_linear_causal_margin',Z-3*ell==5*ell>0)
    check('retained_example_nonlinear_spatial_margin',Z-(3*ell+F(5,4)*F(7,2)*ell)==F(5,8)*ell>0)
    check('retained_example_nonlinear_time_margin',F(7,2)*ell-3*ell==ell/2>0)

    reject('reverse_inward_rho_sign',lambda:require(-rhs['rho']==archived_rhs['rho'],'inward sign mismatch'))
    reject('claim_closed_strict_threshold',lambda:require(all(v>0 for pair in margins(limit).values() for v in pair),'strict inclusion fails at threshold'))
    reject('claim_ell_1e5_uses_this_collar',lambda:require(all(v>0 for pair in margins(F(8,10**5)).values() for v in pair),'too-wide collar'))
    reject('forget_scaled_output_factor',lambda:require(coarse['eta']/ell<F(251,250000),'unscaled output exceeds gate'))
    reject('claim_tail_exactly_zero',lambda:require(coarse['eta_tail']==0,'tail is nonzero bound'))

    def encoded(x):
        if isinstance(x,F):return str(x)
        if isinstance(x,I):return x.pair()
        if isinstance(x,dict):return {k:encoded(v) for k,v in x.items()}
        if isinstance(x,(tuple,list)):return [encoded(v) for v in x]
        return x
    receipt={
        'status':'PASS_PLANNING_FEASIBILITY_ONLY',
        'source_pins':INPUTS,'checks':checks,'negative_controls_rejected':negatives,
        'fixed_rectangle_strict_inclusion_depth_supremum':limit,
        'threshold_interpretation':'Every 0<=L<L_star passes the stated strict closed-box inclusion; equality does not. This is not the actual maximal background continuation depth.',
        'directional_depth_thresholds':thresholds,'independently_rederived_rhs':rhs,
        'proposed_state':{'ell':ell,'T':3*ell,'Z_and_L_energy':Z},
        'first_exit_margins_at_proposed_Z':margins_Z,
        'independently_rederived_coefficients':{'W_all_real_t':potential,'mu':mu,'beta_all_real_t':beta,'gamma':gamma},
        'coarse_majorant':coarse,'improved_rounded_majorant':improved,
        'shape_gate_inputs':{'coarse_eta_upper':F(251,250000),'ell':ell,
                            'required_peak_gap':2*F(251,250000)+ell,
                            'conservative_gap_target':F(1,50),
                            'status':'See independent shape certificate for gap proof'},
        'assumptions':[
            'retained analytic branch and shell box at delta=1/1000',
            'exact inward background ODE and smooth continuation in rho>0 rectangle',
            'retained coupled first-variation equation and physical reconstruction',
            'for an exact-family interpretation, localized inverse/support/gauge/evolution hypotheses must hold for the newly specified wider tangent',
        ],
        'not_claimed':['updated published waveform','transfer of narrow-state C interval',
                       'finite nonlinear amplitude radius','observed physical scale',
                       'actual maximal collar depth','execution of archived producers'],
    }
    print(json.dumps(encoded(receipt),indent=2,sort_keys=True))


if __name__=='__main__':main()
