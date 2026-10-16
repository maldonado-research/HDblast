#!/usr/bin/env python3
"""Exact C1 contraction, explicit comparator, defect cone and actual point signs."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

def require(ok,msg):
    if not ok:raise ValueError(msg)

def eq(a,b,msg):
    if a!=b:raise ValueError(msg+': '+repr(a)+' != '+repr(b))

def mm(a,b):return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def mv(a,v):return [sum(a[i][j]*v[j] for j in range(len(v))) for i in range(len(a))]

def inverse(a):
    n=len(a);v=[row[:]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for i in range(n):
        pivot=next((j for j in range(i,n) if v[j][i]),None)
        require(pivot is not None,'resolvent matrix singular')
        v[i],v[pivot]=v[pivot],v[i];p=v[i][i];v[i]=[x/p for x in v[i]]
        for j in range(n):
            if j!=i:
                q=v[j][i];v[j]=[v[j][k]-q*v[i][k] for k in range(2*n)]
    inv=[row[n:] for row in v]
    eq(mm(a,inv),[[F(i==j) for j in range(n)] for i in range(n)],'exact resolvent inverse')
    return inv

def cap(x):
    require(x>=0,'nonnegative majorant')
    if not x:return {'exact_upper':'0','strict_simple_upper':'0','display_only':0.0}
    s=F(1)
    while x>=s:s*=10
    while x<s/10:s/=10
    unit=s/1000;c=(x//unit+1)*unit
    require(c>x,'strict upper cap')
    return {'exact_upper':str(x),'strict_simple_upper':str(c),'display_only':float(x)}

def sign(s):return F((s>0)-(s<0))

def green_identity_checks(injection):
    """Free profiles are extended; coefficients and forcing never are."""
    n=0;c=F(7,13)
    Q=lambda s:s**5
    Qp=lambda s:5*s**4
    Qpp=lambda s:20*s**3
    P=lambda s:(s*Q(s) if injection=='wrong_product_parity' else sign(s)*abs(s)**6)
    Pp=lambda s:6*abs(s)**5
    Ppp=lambda s:30*sign(s)*abs(s)**4
    pref=c/(8 if injection=='half_commutator' else 4)
    for x in (F(0),F(1,2),F(1),F(3,2)):
        for r in (F(0),F(1,2),F(1),F(3,2),F(2)):
            sp,sm=r+x,r-x
            value=pref*(P(sp)+P(sm)-r*(Q(sp)+Q(sm)))
            if r==0:eq(value,F(0),'odd product-profile Dirichlet boundary')
            if x==0:
                eq(value,F(0),'comparator zero initial displacement')
                px=pref*(Pp(sp)-Pp(sm)-r*(Qp(sp)-Qp(sm)))
                eq(px,F(0),'comparator zero initial velocity')
            pxx=pref*(Ppp(sp)+Ppp(sm)-r*(Qpp(sp)+Qpp(sm)))
            prr=pref*(Ppp(sp)+Ppp(sm)-2*(Qp(sp)+Qp(sm))-r*(Qpp(sp)+Qpp(sm)))
            eq(pxx-prr,c*(Qp(sp)+Qp(sm))/2,'explicit comparator commutator source')
            n+=1+(r==0)+2*(x==0)
    # Integrate the actual indicator supports for a constant physical source.
    for x,r in ((F(1),F(0)),(F(3),F(1)),(F(1),F(2)),(F(3,2),F(3,2))):
        split=max(F(0),x-r)
        first=(r+x)*split-split*split/2+(x-split)**2 if split else x*x
        image=split*split/2
        dsign=1 if injection=='wrong_D_image' else -1
        nsign=-1 if injection=='wrong_N_image' else 1
        directD=(first+dsign*image)/2;directN=(first+nsign*image)/2
        expectedD=x*x/2 if r>=x else r*x-r*r/2
        eq(directD,expectedD,'Dirichlet kernel row mass')
        eq(directN,x*x/2,'Neumann kernel row mass')
        n+=2
    return n

def special_Q(s):
    a=abs(s)
    if a<=1 or a>=2:return F(0)
    eq(a,F(3,2),'chosen point uses only exact bump peak or zeros')
    return F(1)

def comparator_value(x,r,c):
    sp,sm=r+x,r-x
    return c/4*(sp*special_Q(sp)+sm*special_Q(sm)
                -r*(sign(sp)*special_Q(sp)+sign(sm)*special_Q(sm)))

def solve_caps(coefficient,master):
    require(coefficient['status']=='PASS_WHOLE_COLLAR_GAUGE_COEFFICIENT_AND_SOURCE_CERTIFICATE','coefficient status')
    require(master['status']=='PASS' and master['scope']['third_master_derivatives_certified'] is True,'actual master C3 input')
    ell=F(1,1000000);Xmax=F(3)
    eq(F(coefficient['domain']['ell']),ell,'coefficient registered width')
    eq([F(x) for x in coefficient['domain']['z']],[F(0),8*ell],'coefficient audited whole collar')
    eq([F(x) for x in coefficient['domain']['t']],[F(0),3*ell],'coefficient audited time window')
    eq(F(master['exact']['interior']['ell']),ell,'registered actual master width')
    eq(master['exact']['interior']['scaled_domain']['x_max'],'3','master scaled time domain')
    eq(master['exact']['interior']['scaled_domain']['x_plus_r_max'],8,'master causal triangle')
    A,B,D=([[F(x) for x in row] for row in coefficient['matrix_elementwise_absolute_caps'][key]] for key in ('A','B','D'))
    require(all(x>=0 for mat in (A,B,D) for row in mat for x in row),'nonnegative coefficient caps')
    pslo,pshi=(F(x) for x in coefficient['source_profile_mixed_derivative_intervals']['P_s']['0,0']['exact'])
    p0=max(abs(F(x)) for x in coefficient['source_profile_mixed_derivative_intervals']['P_0']['0,0']['exact'])
    require(0<pslo<=pshi,'positive actual derivative source coefficient')
    ps=max(abs(pslo),abs(pshi));c=6*ell/130;delta=max(abs(pslo-c),abs(pshi-c))
    data=master['exact'];W0=F(data['interior']['U']);W1=F(data['interior']['first'])
    byprime=F(data['boundary']['yprime']);potential=F(data['background']['W'])
    eq(3*ell**2*potential*W0,F(data['interior']['first_error']),'authenticated master derivative Duhamel error')
    rows={};numeric={}
    for X in (F(1,2),F(3,4),F(1),F(2),F(3)):
        M=[[X*(A[i][j]+B[i][j]+X*D[i][j]/2) for j in range(3)] for i in range(3)]
        q=max(sum(row) for row in M);require(q<1,'actual full-window C1 contraction')
        inv=inverse([[F(i==j)-M[i][j] for j in range(3)] for i in range(3)])
        require(all(x>=0 for row in inv for x in row),'positive contraction resolvent')
        boundary=F(0) if X<=1 else byprime
        velocity=ell
        dw=velocity+boundary+X*ell**2*potential*W0
        source_defect=16*delta+ps*dw+p0*W0
        p_value=c*X/2;p_first=c*(F(1,2)+8*X);p_norm=max(c,p_first)
        R=[(source_defect if i==0 else F(0))+(A[i][0]+B[i][0])*p_first+D[i][0]*p_value for i in range(3)]
        error=mv(inv,[X*x for x in R])
        direct=mv(inv,[X*(ps*W1+p0*W0),F(0),F(0)])
        val=[min(X*direct[i]/2,(p_value if i==0 else 0)+X*error[i]/2) for i in range(3)]
        deriv=[min(direct[i],(p_first if i==0 else 0)+error[i]) for i in range(3)]
        norms=[min(direct[i],(p_norm if i==0 else 0)+error[i]) for i in range(3)]
        key=str(X)
        numeric[key]={'error':error,'value':val,'first':deriv,'M':M}
        rows[key]={'X':str(X),'contraction_matrix':[[str(x) for x in row] for row in M],
                   'contraction_row_norm':str(q),'resolvent':[[str(x) for x in row] for row in inv],
                   'master_radial_error_from_displacement_comparator':str(dw),
                   'actual_boundary_derivative_term':str(boundary),'source_defect':str(source_defect),
                   'comparator_value':str(p_value),'comparator_each_first_derivative':str(p_first),
                   'forcing_residual_by_channel':[str(x) for x in R],
                   'weighted_C1_error_by_channel':[cap(x) for x in error],
                   'value_error_by_channel':[cap(X*x/2) for x in error],
                   'actual_amplitude_value_by_channel':[cap(x) for x in val],
                   'actual_each_scaled_first_derivative_by_channel':[cap(x) for x in deriv],
                   'actual_weighted_C1_by_channel':[cap(x) for x in norms],
                   'pure_contraction_weighted_C1_by_channel':[cap(x) for x in direct]}
    points=[]
    for x,r in ((F(1,2),F(1)),(F(1,2),F(2)),(F(3,4),F(3,4)),
                (F(1),F(1,2)),(F(2),F(1,2))):
        p=comparator_value(x,r,c);err=x*numeric[str(x)]['error'][0]/2
        lo,hi=p-err,p+err
        point_sign='positive' if lo>0 else 'negative' if hi<0 else 'not decided'
        require(point_sign!='not decided','actual signed gauge amplitude point certificate')
        convenient={(F(1,2),F(1)):(F(56,10**10),F(6,10**9)),
                    (F(1,2),F(2)):(-F(6,10**9),-F(56,10**10)),
                    (F(3,4),F(3,4)):(F(83,10**10),F(9,10**9))}.get((x,r))
        if convenient:
            require(convenient[0]<lo and hi<convenient[1],'displayed signed point strict enclosure')
        points.append({'scaled_x':str(x),'scaled_r':str(r),'physical_t':str(ell*x),'physical_z':str(ell*r),
                       'comparator_radial_amplitude':str(p),'actual_radial_amplitude_interval':[str(lo),str(hi)],
                       'actual_sign':point_sign,'characteristic_crossing':x==r,
                       'convenient_strict_actual_interval':[str(a) for a in convenient] if convenient else None})
    globalX=rows['3'];g=numeric['3']
    # Physical coordinate derivative conversions; not a background-frame Lie norm.
    rho_min=min(F(x) for x in coefficient['static_z_jets']['rho'][0]['exact'])
    require(rho_min>0,'actual proper radial conversion radius')
    physical={'each_time_derivative':[cap(x/ell) for x in g['first']],
              'each_proper_y_derivative':[cap(x/(ell*rho_min)) for x in g['first']]}
    return {'scaled_channel_order':['f=xi^y/Y','b=xi^t/Y','d=angular_gradient_amplitude'],
            'source_caps':{'P_s_interval':[str(pslo),str(pshi)],'P_0_absolute':str(p0),
                           'freezing_constant':str(c),'P_s_minus_freezing_absolute':str(delta),
                           'actual_master_value':str(W0),'actual_master_each_scaled_first_derivative':str(W1),
                           'actual_boundary_scaled_derivative':str(byprime)},
            'time_window_certificates':rows,'signed_actual_points':points,
            'physical_coordinate_derivative_caps_on_full_window':physical,
            'Picard_error_contract':'From comparator p, n further actual-source Picard steps have error cone M(X)^n E(X); scalar max norm <=q(X)^n max(E(X)).'}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    root=Path(__file__).resolve().parent;ap=argparse.ArgumentParser()
    ap.add_argument('--input-root',type=Path,default=root)
    ap.add_argument('--pins',type=Path,default=root/'INPUT_PINS.json')
    ap.add_argument('--output',type=Path)
    ap.add_argument('--run-controls',action='store_true')
    mutations={'wrong_D_image':'Dirichlet kernel row mass','wrong_N_image':'Neumann kernel row mass',
               'wrong_product_parity':'odd product-profile Dirichlet boundary',
               'half_commutator':'explicit comparator commutator source'}
    ap.add_argument('--inject-error',choices=['none']+list(mutations),default='none')
    a=ap.parse_args();checks=green_identity_checks(a.inject_error)
    pins=json.loads(a.pins.read_text());inputs={};hashes={}
    for key,pin in pins.items():
        path=a.input_root/pin['path'];raw=path.read_bytes()
        require(hashlib.sha256(raw).hexdigest()==pin['sha256'] and len(raw)==pin['bytes'],'gauge solver input authentication: '+key)
        hashes[key]=pin['sha256']
        if path.suffix=='.json':inputs[key]=json.loads(raw)
    require(inputs['coefficients']['input_sha256']['HIGHER_JETS_NORMAL.json']==hashes['master'],'actual master/coefficient dependency link')
    report={'status':'PASS_ACTUAL_GAUGE_WAVE_C1_SOLUTION_AND_EXPLICIT_DEFECT_CERTIFICATE',
            'finite_exact_kernel_and_comparator_checks':checks,'input_sha256':hashes,
            'domain':{'physical':'0<=t<=3ell,z>=0,t+z<=8ell','ell':'1/1000000',
                      'scaled':'0<=x<=3,r>=0,x+r<=8'},
            'boundary_rows':['f(x,0)=0','b_r(x,0)=0','d_r(x,0)=0'],
            'initial_rows':'all amplitudes and their x derivatives zero at x=0',
            'norm':'Raw amplitudes. Per channel N_i=max(2 sup|U_i|/X,sup|U_i_x|,sup|U_i_r|).',
            'scope':{'actual_registered_master_used':True,'causal_C1_fixed_point_constructed':True,
                     'explicit_comparator_evaluated':True,'actual_signed_gauge_amplitude_points':True,
                     'numerical_C2_C3_gauge_bounds':False,'transformed_first_fields_evaluated':False,
                     'transformed_quadratic_sources_evaluated':False,'second_evolution_computed':False,
                     'gauge_amplitude_sign_is_invariant_force':False},
            **solve_caps(inputs['coefficients'],inputs['master'])}
    if a.run_controls:
        cases=[]
        with tempfile.TemporaryDirectory(prefix='actual-gauge-controls-') as td:
            for mutation,reason in mutations.items():
                output=Path(td)/(mutation+'.json');cmd=[sys.executable]
                if sys.flags.optimize:cmd.append('-O')
                cmd += [str(Path(__file__).resolve()),'--input-root',str(a.input_root),'--pins',str(a.pins),
                        '--inject-error',mutation,'--output',str(output)]
                run=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
                require(run.returncode!=0 and ('ValueError: '+reason) in run.stderr and not output.exists(),
                        'gauge solver control intended rejection: '+mutation)
                cases.append({'mutation':mutation,'rejection_prefix':reason,'no_output_report':True})
        report['controls']={'status':'PASS','optimized':bool(sys.flags.optimize),'cases':cases}
    text=json.dumps(report,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
