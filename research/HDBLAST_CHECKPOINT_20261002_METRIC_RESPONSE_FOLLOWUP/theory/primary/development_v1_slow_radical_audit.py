#!/usr/bin/env python3
"""Pure formal metric contacts; no physical profile, quadrature, or mode run."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback
import sympy as S

L,W,V=S.symbols('L omega v',positive=True)
H=S.symbols('h0:7',real=True)
K2=W**2-2*L**2


def D(expr):
    return S.expand(L**2*S.diff(expr,L)+2*L**3/W*S.diff(expr,W)
                    +sum(H[i+1]*S.diff(expr,H[i]) for i in range(6)))

def plus(a,b):return tuple(S.expand(a[i]+b[i]) for i in (0,1))
def neg(a):return(-a[0],-a[1])
def sub(a,b):return plus(a,neg(b))
def mul(a,b):return(S.expand(a[0]*b[0]),S.expand(a[0]*b[1]+a[1]*b[0]))
def power(a,n):return(S.expand(a[0]**n),S.expand(n*a[0]**(n-1)*a[1]))
def div(a,b):return mul(a,power(b,-1))
def scale(a,b):return(S.expand(a[0]*b),S.expand(a[1]*b))
def deriv(a):return(D(a[0]),D(a[1]))


def inventory():
    ell=(L,H[1]); omega=(W,2*L**2*H[0]/W)
    T=(2*L**2,2*L*H[1]+H[2]); wp=deriv(omega); wpp=deriv(wp)
    U=plus(plus(scale(div(T,omega),-S.Rational(1,2)),
               scale(div(wpp,power(omega,2)),-S.Rational(1,4))),
           scale(div(power(wp,2),power(omega,3)),S.Rational(3,8)))
    Up=deriv(U); Upp=deriv(Up)
    U4=plus(plus(plus(plus(scale(div(power(U,2),omega),-S.Rational(1,2)),
                         scale(div(Upp,power(omega,2)),-S.Rational(1,4))),
                    scale(div(mul(wpp,U),power(omega,3)),S.Rational(1,4))),
               scale(div(mul(wp,Up),power(omega,3)),S.Rational(3,4))),
           scale(div(mul(power(wp,2),U),power(omega,4)),-S.Rational(3,4)))
    b=(-K2/3-2*L**2,-4*L**2*H[0]); c=sub((1,0),div(b,power(omega,2)))
    J2=plus(div(mul(ell,wp),power(omega,2)),
            scale(div(power(wp,2),power(omega,3)),S.Rational(1,4)))
    J4=plus(plus(plus(div(mul(ell,Up),power(omega,2)),
                     scale(div(mul(mul(ell,wp),U),power(omega,3)),-2)),
                scale(div(mul(wp,Up),power(omega,3)),S.Rational(1,2))),
            scale(div(mul(power(wp,2),U),power(omega,4)),-S.Rational(3,4)))
    r0=scale(omega,S.Rational(1,2)); p0=scale(div((K2,0),omega),S.Rational(1,6))
    r2=scale(plus(div(power(ell,2),omega),J2),S.Rational(1,4))
    p2=scale(plus(plus(mul(c,U),div(power(ell,2),omega)),J2),S.Rational(1,4))
    r4=scale(plus(sub(div(power(U,2),omega),div(mul(power(ell,2),U),power(omega,2))),J4),S.Rational(1,4))
    p4=scale(plus(sub(plus(mul(c,U4),div(mul(b,power(U,2)),power(omega,3))),
                        div(mul(power(ell,2),U),power(omega,2))),J4),S.Rational(1,4))
    Nr=plus(plus(r0,r2),r4); Np=plus(plus(p0,p2),p4)
    Fq=sub(scale(power(omega,-1),S.Rational(1,2)),scale(div(U,power(omega,2)),S.Rational(1,2)))
    g=4*L**2*H[0]-2*L*H[1]-H[2]
    dUm=g/(2*W); dUm1=D(dUm); dUm2=D(dUm1)
    dVm=(-U[0]*dUm/W-dUm2/(4*W**2)+wpp[0]*dUm/(4*W**3)
          +3*wp[0]*dUm1/(4*W**3)-3*wp[0]**2*dUm/(4*W**4))
    dJm=(L*dUm1/W**2-2*L*wp[0]*dUm/W**3
          +wp[0]*dUm1/(2*W**3)-3*wp[0]**2*dUm/(4*W**4))
    dNrm=S.expand((g/W+2*U[0]*dUm/W-g*U[0]/W**2-L**2*dUm/W**2+dJm)/4)
    dNpm=S.expand((c[0]*dUm-g/W+c[0]*dVm+2*b[0]*U[0]*dUm/W**3
                  +g*U[0]/W**2-L**2*dUm/W**2+dJm)/4)
    return {'g':g,'g_prime':D(g),'g_second':D(D(g)), 'U2':U,'U4':U4,'Fq':Fq,
            'Nr':Nr,'Np':Np,'dNrm':dNrm,'dNpm':dNpm,
            'grades':[(r0,p0),(r2,p2),(r4,p4)]}


def moment(n):
    return S.expand((2*L**2)**S.Rational(3-n,2)*sum(
        (-1)**j*S.binomial((n-5)//2,j)*V**(2*j+3)/S.Integer(2*j+3)
        for j in range((n-5)//2+1)))


def integrate_high(expr):
    expr=S.expand(expr)
    return S.expand(sum(expr.coeff(W,-n)*moment(n) for n in range(5,16,2)))


def closed_quantities(data):
    Rlead=L**4*V**2*(2*V+3)/(4*(V+1)**2)
    Plead=-L**4*V**2*(4*V+3)/(12*(V+1)**2)
    R0=S.factor(Rlead-integrate_high(data['Nr'][0]))
    P0=S.factor(Plead-integrate_high(data['Np'][0]))
    Q0=V**2/(2*(1+V))-V**3/48-V**5/16
    Cq=S.factor(integrate_high(data['Fq'][1]+data['g']/(4*W**3)))
    Er=S.factor((L*H[1]+H[2]/4)*L**2*V**2/(1+V)-integrate_high(data['Nr'][1]-data['dNrm']))
    Ep=S.factor(-H[2]*L**2*V**2/(4*(1+V))-integrate_high(data['Np'][1]-data['dNpm']))
    return {'scaled_baseline_rho_times_2pi2':R0,'scaled_baseline_p_times_2pi2':P0,
            'physical_baseline_Q_times_2pi2':Q0,'Cq_times_2pi2':Cq,
            'Er_times_2pi2':Er,'Ep_times_2pi2':Ep}


def audit():
    data=inventory(); closed=closed_quantities(data); checks=[]; mutants=[]
    def equal(name,left,right=0):
        residue=S.factor(S.expand(left-right))
        if residue != 0:raise RuntimeError(name+': '+str(residue))
        checks.append(name)
    def reject(name,residue):
        if S.factor(S.expand(residue)) == 0:raise RuntimeError('Mutation invisible: '+name)
        mutants.append(name)
    for index,(r,p) in enumerate(data['grades']):
        ward=sub(plus(deriv(r),scale(mul((L,H[1]),p),3)),mul((L,H[1]),r))
        equal('fixed_mass_subtraction_grade_'+str(2*index)+'_Ward_baseline',ward[0])
        equal('fixed_mass_subtraction_grade_'+str(2*index)+'_Ward_metric_variation',ward[1])
    # Exact rational k-to-v substitution; derivative reconstructs original integral.
    vv=W*S.sqrt(1-2*L**2/W**2)/W  # overwritten below with separate positive k
    k=S.symbols('k',positive=True); wk=S.sqrt(k*k+2*L**2); vk=k/wk
    for n in range(5,16,2):
        equal('finite_moment_'+str(n)+'_upper_derivative',S.diff(moment(n).subs(V,vk),k),k*k/wk**n)
        equal('finite_moment_'+str(n)+'_zero_endpoint',moment(n).subs(V,0))
    Nr=data['Nr'][0].subs(W,wk);Np=data['Np'][0].subs(W,wk)
    equal('density_baseline_closed_derivative',S.diff(closed['scaled_baseline_rho_times_2pi2'].subs(V,vk),k),k*k*(k/2+3*L**2/(4*k)-Nr))
    equal('pressure_baseline_closed_derivative',S.diff(closed['scaled_baseline_p_times_2pi2'].subs(V,vk),k),k*k*(k/6-L**2/(4*k)-Np))
    equal('variance_baseline_closed_derivative',S.diff(closed['physical_baseline_Q_times_2pi2'].subs(V,vk),k),k*k/L**2*(1/(2*k)-data['Fq'][0].subs(W,wk)))
    for name in ('scaled_baseline_rho_times_2pi2','scaled_baseline_p_times_2pi2','physical_baseline_Q_times_2pi2','Cq_times_2pi2','Er_times_2pi2','Ep_times_2pi2'):
        equal(name+'_zero_endpoint',closed[name].subs(V,0))
    qdiff=S.expand(data['Fq'][1]+data['g']/(4*W**3))
    rdiff=S.expand(data['Nr'][1]-data['dNrm']);pdiff=S.expand(data['Np'][1]-data['dNpm'])
    equal('variance_contact_closed_derivative',S.diff(closed['Cq_times_2pi2'].subs(V,vk),k),k*k*qdiff.subs(W,wk))
    equal('density_contact_closed_derivative',S.diff(closed['Er_times_2pi2'].subs(V,vk),k),k*k*((L*H[1]+H[2]/4)/k-rdiff.subs(W,wk)))
    equal('pressure_contact_closed_derivative',S.diff(closed['Ep_times_2pi2'].subs(V,vk),k),k*k*(-H[2]/(4*k)-pdiff.subs(W,wk)))
    equal('continuum_baseline_density',closed['scaled_baseline_rho_times_2pi2'].subs(V,1),11*L**4/480)
    equal('continuum_baseline_pressure',closed['scaled_baseline_p_times_2pi2'].subs(V,1),-11*L**4/480)
    equal('continuum_metric_variance_contact',closed['Cq_times_2pi2'].subs(V,1),-L**2*H[0]/2+L*H[1]/12+H[2]/24)
    equal('continuum_metric_density_contact',closed['Er_times_2pi2'].subs(V,1),L*(-21*L**2*H[1]+7*L*H[2]+3*H[3])/120)
    equal('continuum_metric_pressure_contact',closed['Ep_times_2pi2'].subs(V,1),-7*L**3*H[1]/360+7*L**2*H[2]/180+L*H[3]/120-H[4]/120)
    Dv=lambda expr:S.expand(L**2*S.diff(expr,L)-L*V*(1-V**2)*S.diff(expr,V)+sum(H[i+1]*S.diff(expr,H[i]) for i in range(6)))
    Rb=closed['scaled_baseline_rho_times_2pi2'];Pb=closed['scaled_baseline_p_times_2pi2']
    equal('renormalized_finite_K_baseline_Ward',Dv(Rb)-L*Rb+3*L*Pb)
    equal('metric_potential_first_derivative',data['g_prime'],8*L**3*H[0]+2*L**2*H[1]-2*L*H[2]-H[3])
    equal('metric_potential_second_derivative',data['g_second'],24*L**4*H[0]+12*L**3*H[1]-H[4])
    reject('drop_baseline_finite_K_Ward_term',3*H[1]*(Rb+Pb))
    reject('omit_full_variance_metric_contact',closed['Cq_times_2pi2'])
    reject('omit_full_density_metric_contact',closed['Er_times_2pi2'])
    reject('omit_full_pressure_metric_contact',closed['Ep_times_2pi2'])
    reject('omit_a_minus_four_measure_variation',4*H[0]*Rb)
    reject('substitute_continuum_density_baseline_at_finite_K',Rb-11*L**4/480)
    reject('substitute_continuum_pressure_baseline_at_finite_K',Pb+11*L**4/480)
    reject('retain_scalar_mass_law_contact_under_fixed_phi',H[0]*closed['physical_baseline_Q_times_2pi2']/4)
    return data,closed,{'passed':True,'exact_identity_count':len(checks),'identities':checks,
                       'symbolic_mutation_count':len(mutants),'mutations_rejected':mutants,
                       'scope':'Exact symbolic identities only. No profile values, mode solutions, physical responses, momentum or time quadratures evaluated.'}


def export(data,closed):
    expressions={name:str(S.factor(expr)) for name,expr in closed.items()}
    expressions.update({'g':str(data['g']),'g_prime':str(data['g_prime']),'g_second':str(data['g_second'])})
    coefficients={}
    for name,expr in {'density_subtraction_geo_minus_mass':data['Nr'][1]-data['dNrm'],
                      'pressure_subtraction_geo_minus_mass':data['Np'][1]-data['dNpm'],
                      'variance_subtraction_geo_minus_mass':data['Fq'][1]+data['g']/(4*W**3),
                      'baseline_scaled_density_subtraction':data['Nr'][0],
                      'baseline_scaled_pressure_subtraction':data['Np'][0]}.items():
        coefficients[name]={str(n):str(S.expand(expr).coeff(W,n)) for n in range(-15,2) if S.expand(expr).coeff(W,n)!=0}
    return {'convention':'omega=sqrt(k^2+2L^2); hN=d^N h/deta^N; all closed expressions omit common1/(2pi^2); Qbaseline is physical, stress baselines are scaled a0^4.',
            'moments':{str(n):str(moment(n)) for n in range(5,16,2)},
            'expressions':expressions,'coefficients':coefficients,
            'mapping':{'qmetric':'qmass - 2*h0*L^2*Q0K - Cq',
                       'Rmetric':'Rmass - 4*h0*R0K + Er',
                       'Pmetric':'Pmass - 4*h0*P0K + Ep',
                       'Jmetric':'qmetric; fixed physical phi gives no scalar mass-law contact'}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--formulas',type=Path)
    args=parser.parse_args()
    if args.output.exists() or args.formulas and args.formulas.exists():raise RuntimeError('Refusing to overwrite development evidence')
    provenance={'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':S.__version__,'python_optimization':sys.flags.optimize}
    try:
        data,closed,result=audit();result['provenance']=provenance
        if args.formulas:
            args.formulas.parent.mkdir(parents=True,exist_ok=True)
            args.formulas.write_text(json.dumps(export(data,closed),indent=2,sort_keys=True)+'\n')
    except BaseException as exc:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps({'passed':False,'exception':repr(exc),'traceback':traceback.format_exc(),'provenance':provenance},indent=2)+'\n')
        raise
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('Pure metric contact algebra PASS:',result['exact_identity_count'],'identities;',result['symbolic_mutation_count'],'mutations.')

if __name__=='__main__':main()
