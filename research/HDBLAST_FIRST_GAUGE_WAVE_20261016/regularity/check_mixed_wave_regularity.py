#!/usr/bin/env python3
"""Independent exact mixed D/N Green identities and a C1 contraction bound.

No project scientific program is imported. The existence and C3 bootstrap
proof, including the limits of finite witnesses, is in MIXED_WAVE_REGULARITY.md.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import subprocess
import sys

PINS={
    'background':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
    'master':'92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df'}
ELL=F(1,10**6)
X=F(3)
CONTROLS={
    'extra_boundary_half':'boundary first-ray row mass is x, not x over two',
    'odd_neumann_source':'Neumann forcing uses even reflection and can have nonzero boundary value',
    'drop_D_second_trace':'Dirichlet second spatial derivative retains minus the shell forcing',
    'drop_N_third_trace':'Neumann third spatial derivative retains minus the shell radial forcing derivative',
    'assume_D_C2_without_corner':'Dirichlet forcing one gives a genuine characteristic second-derivative jump',
    'assume_N_C3_without_corner':'Neumann forcing r gives a genuine characteristic third-derivative jump'}


def need(ok,message):
    if not ok:raise ValueError(message)


def int_ray(m,k,A,B,lo,hi):
    return sum((F(comb(k,j))*A**(k-j)*B**j*(hi**(m+j+1)-lo**(m+j+1))/(m+j+1)
                for j in range(k+1)),F(0))


def pieces(x,r,inside):
    split=x-r if inside else F(0)
    return [(F(0),split,-1),(split,x,1)]


def kernel(mode,m,k,x,r,inside=None):
    inside=x>r if inside is None else inside
    upper=int_ray(m,k+1,r+x,F(-1),F(0),x)
    lower=sum(((sign if mode=='N' else 1)*int_ray(m,k+1,sign*(r-x),F(sign),a,b)
               for a,b,sign in pieces(x,r,inside)),F(0))
    return (upper-lower)/(2*(k+1))


def force(m,k,x,r,ix=0,ir=0):
    if m<ix or k<ir:return F(0)
    return F(factorial(m),factorial(m-ix))*F(factorial(k),factorial(k-ir))*x**(m-ix)*r**(k-ir)


def rays(m,k,x,r,order,odd,inside):
    if order>k:return F(0),F(0)
    coefficient=F(factorial(k),factorial(k-order))
    upper=coefficient*int_ray(m,k-order,r+x,F(-1),F(0),x)
    lower=coefficient*sum(((sign if odd else 1)*int_ray(m,k-order,sign*(r-x),F(sign),a,b)
                           for a,b,sign in pieces(x,r,inside)),F(0))
    return upper,lower


def formula(mode,m,k,x,r,ix,ir,inside=None):
    inside=x>r if inside is None else inside
    if ix+ir==1:
        a,b=rays(m,k,x,r,0,mode=='D',inside)
        return (a+(b if ix else -b))/2
    a,b=rays(m,k,x,r,1,mode=='N',inside)
    trace=force(m,k,x-r,F(0)) if inside and mode=='D' else F(0)
    second={(0,2):(a-b)/2-trace,(1,1):(a+b)/2+trace}
    second[2,0]=force(m,k,x,r)+second[0,2]
    if ix+ir==2:return second[ix,ir]
    a,b=rays(m,k,x,r,2,mode=='D',inside)
    if inside:
        trace=(force(m,k,x-r,F(0),ix=1) if mode=='D' else -force(m,k,x-r,F(0),ir=1))
    else:trace=F(0)
    third={(0,3):(a-b)/2+trace,(1,2):(a+b)/2-trace}
    third[2,1]=force(m,k,x,r,ir=1)+third[0,3]
    third[3,0]=force(m,k,x,r,ix=1)+third[1,2]
    return third[ix,ir]


@lru_cache(None)
def lagrange_weights(degree,order):
    if order>degree:return [F(0)]*(degree+1)
    result=[]
    for j in range(degree+1):
        p=[F(1)];den=F(1)
        for k in range(degree+1):
            if k==j:continue
            new=[F(0)]*(len(p)+1)
            for i,v in enumerate(p):new[i]-=k*v;new[i+1]+=v
            p=new;den*=j-k
        result.append(p[order]*factorial(order)/den)
    return result


def exact_branch_derivative(mode,m,k,x,r,ix,ir,inside):
    degree=m+k+2;step=F(1,100*(degree+1))
    wx=lagrange_weights(degree,ix) if ix else [F(1)]
    wr=lagrange_weights(degree,ir) if ir else [F(1)]
    return sum((a*b*kernel(mode,m,k,x+i*step,r+j*step,inside)
                for i,a in enumerate(wx) for j,b in enumerate(wr)),F(0))/step**(ix+ir)


def reflection_checks():
    derivatives=[(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,1),(1,2),(0,3)]
    witnesses=joins=0
    for mode in ('D','N'):
        for m in range(3):
            for k in range(3):
                for x,r,inside in ((F(2),F(1,2),True),(F(1,2),F(2),False),(F(2),F(0),True)):
                    for ix,ir in derivatives:
                        need(formula(mode,m,k,x,r,ix,ir,inside)==exact_branch_derivative(mode,m,k,x,r,ix,ir,inside),
                             'independent mixed-boundary polynomial derivative')
                        witnesses+=1
                for ix,ir in derivatives:
                    inner=formula(mode,m,k,F(1),F(1),ix,ir,True)
                    outer=formula(mode,m,k,F(1),F(1),ix,ir,False)
                    compatible=True
                    if mode=='D' and ix+ir>=2:compatible=force(m,k,F(0),F(0))==0
                    if mode=='D' and ix+ir>=3:compatible=compatible and force(m,k,F(0),F(0),ix=1)==0
                    if mode=='N' and ix+ir>=3:compatible=force(m,k,F(0),F(0),ir=1)==0
                    if compatible:
                        need(inner==outer,'compatible mixed-boundary characteristic join')
                        joins+=1
    need(formula('D',0,0,F(2),F(0),1,0)==0 and formula('D',0,0,F(2),F(0),0,1)==2,
         'Dirichlet boundary first derivatives')
    need(formula('N',0,0,F(2),F(0),0,1)==0 and formula('N',0,0,F(2),F(0),1,0)==2,
         'Neumann boundary first derivatives')
    return {'physical_monomials_per_mode':9,'independent_derivative_witnesses':witnesses,
            'compatible_characteristic_joins':joins,'modes':['Dirichlet','Neumann']}


def authenticate(path,name):
    raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    need(digest==PINS[name],'pinned actual input: '+name)
    return json.loads(raw)


def actual_C1_certificate(background,master):
    need(background['status']=='PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ','actual background status')
    state=background['domain']
    need((F(state['ell']),F(state['T']),F(state['Z']))==(ELL,X*ELL,8*ELL),'actual whole causal collar')
    b=background['background_intervals'];c=background['reconstruction_coefficients']
    lo=lambda row:F(row['exact'][0])
    hi=lambda row:F(row['exact'][1])
    absolute=lambda row:max(abs(lo(row)),abs(hi(row)))
    need(lo(b['rho'])>=129 and hi(b['rho'])<=131 and absolute(b['a'])<=F(14,125) and
         absolute(b['p'])<=F(7,50000),'retained global rho a p boxes')
    need(absolute(c['u_w_coefficient'])<F(61,10**6) and absolute(c['u_w_coefficient_z'])<F(3,200),
         'retained global C and Cz bounds')
    need(master['scope']['third_master_derivatives_certified'],'actual global master certificate')
    adapter=master['source_bound_adapter']
    need(adapter['ell']==str(ELL) and adapter['domain']=='0<=x<=3, r>=0, x+r<=8','actual master causal triangle')
    jets=list(map(F,adapter['scaled_order_bounds']))
    need(jets[0]<F(1001,1000) and jets[1]<F(16001,1000),'conservative actual master value and first jets')
    rho,a,p,H=F(131),F(14,125),F(7,50000),X*ELL
    ay=F(1,129**2)+p*p/3
    # Raw amplitude order is (xi^y,xi^t,d), not a rescaled vector norm.
    L1rows=[ELL*(3*H+2*rho*rho*a+3*rho*a),
            ELL*(2*a+3*H+5*rho*a),ELL*(5*H+5*rho*a)]
    L0rows=[ELL**2*(8+4*rho*rho*ay+8*rho*rho*a*a+6*rho*rho*a*H+16*rho*rho*a),
            ELL**2*(3+3*H*H+8+16*H),ELL**2*(2*a+2*H+4)]
    L1,L0=max(L1rows),max(L0rows)
    q=X*(L1+X*L0/2)
    need(L1<F('0.003889') and L0<F('3.3e-8') and q<F(3,250),'actual global C1 contraction less than 0.012')
    u=F(61,10**6)*F(1001,1000)
    uy=(F(3,200)*F(1001,1000)+F(61,10**6)*F(16001,1000)/ELL)/129
    source=6*ELL**2*rho*rho*(uy+4*a*u)
    need(source<F(8,10**7),'scaled actual radial source below 8e-7')
    bound=X*source/(1-q)
    need(bound<F(6,2470000),'strict weighted C1 solution bound')
    return {'ell':str(ELL),'X':str(X),'L1':str(L1),'L0':str(L0),'contraction_q':str(q),
            'scaled_source_sup':str(source),'weighted_C1_solution_sup':str(bound),
            'raw_amplitude_value_sup':str(X*bound/2),'raw_scaled_first_derivative_sup':str(bound),
            'coarse_strict_scaled_first_cap':'243/100000000','coarse_strict_value_cap':'911/250000000',
            'display_only':{'q':float(q),'source':float(source),'weighted_C1':float(bound)}}


def inject(name):
    reason=CONTROLS[name]
    if name=='extra_boundary_half':need(formula('D',0,0,F(2),F(0),0,1)/2==2,reason)
    elif name=='odd_neumann_source':need(kernel('D',0,0,F(2),F(0))==kernel('N',0,0,F(2),F(0)),reason)
    elif name=='drop_D_second_trace':need(F(0)==formula('D',0,0,F(2),F(1),0,2),reason)
    elif name=='drop_N_third_trace':need(F(0)==formula('N',0,1,F(2),F(1),0,3),reason)
    elif name=='assume_D_C2_without_corner':need(formula('D',0,0,F(1),F(1),0,2,True)==formula('D',0,0,F(1),F(1),0,2,False),reason)
    elif name=='assume_N_C3_without_corner':need(formula('N',0,1,F(1),F(1),0,3,True)==formula('N',0,1,F(1),F(1),0,3,False),reason)
    raise ValueError('deliberate corruption unexpectedly passed')


def controls(output):
    cases=[]
    for optimized in (False,True):
        for name,reason in CONTROLS.items():
            cmd=[sys.executable,'-I']+(['-O'] if optimized else [])+['-B',str(Path(__file__).resolve()),'--inject',name]
            run=subprocess.run(cmd,text=True,capture_output=True,check=False)
            need(run.returncode==1 and run.stdout=='' and run.stderr=='FAIL: '+reason+'\n','intended reflected control failure: '+name)
            cases.append({'name':name,'optimized':optimized,'status':'EXPECTED_EXPLICIT_FAILURE'})
    result={'status':'PASS','control_runs':len(cases),'cases':cases}
    if output:output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','control_runs':len(cases)}))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--background',type=Path)
    ap.add_argument('--master',type=Path)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--inject',choices=list(CONTROLS))
    ap.add_argument('--run-controls',action='store_true')
    args=ap.parse_args()
    try:
        if args.inject:inject(args.inject)
        if args.run_controls:controls(args.output);return 0
        need(args.background is not None and args.master is not None,'explicit actual background and master inputs required')
        bg=authenticate(args.background,'background');master=authenticate(args.master,'master')
        result={'status':'PASS_INDEPENDENT_MIXED_WAVE_C1_AND_REFLECTION_CERTIFICATE','optimized':not __debug__,
            'source_pins':PINS,'reflection':reflection_checks(),'C1_actual_raw_amplitudes':actual_C1_certificate(bg,master),
            'scope':{'conditional_on_retained_smooth_actual_inputs':True,'C1_solution_exists_by_analytic_contraction':True,
                'C3_solution_exists_by_compatible_C2_coefficient_source_bootstrap':True,
                'C2_C3_numerical_solution_bounds_evaluated':False,'full_second_evolution_solved':False,
                'high_Sobolev_or_nonlinear_IBVP_claim':False}}
        text=json.dumps(result,indent=2)+'\n'
        if args.output:args.output.write_text(text)
        else:print(text,end='')
    except (ValueError,OSError,KeyError) as error:
        print('FAIL: '+str(error),file=sys.stderr);return 1
    return 0


if __name__=='__main__':sys.exit(main())
