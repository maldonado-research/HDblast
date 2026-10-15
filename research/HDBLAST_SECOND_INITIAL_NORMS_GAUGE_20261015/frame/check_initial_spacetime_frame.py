#!/usr/bin/env python3
"""Authenticated exact bounds for the completed initial spacetime second 1-jet."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

def require(ok,msg):
    if not ok:raise ValueError(msg)

def eq(a,b,msg):require(a==b,msg+': '+str(a-b))

def arr(n):return [[F(0) for _ in range(n)] for _ in range(n)]

def identity_checks(injection):
    """Coordinate-to-frame conversion with the nonzero boost connection."""
    count=0
    for seed in (271,1927,6469):
        rng=random.Random(seed);rnd=lambda:F(rng.randint(-7,7),rng.randint(1,6))
        rho=F(7,3);a,U,Uy,T,q,f,fy,Pi2,Pi1=[rnd() for _ in range(9)]
        if not a:a=F(2,5)
        if not U:U=F(3,7)
        if not fy:fy=F(5,9)
        if not Pi1:Pi1=F(5,11)
        DU=[rnd() for _ in range(3)];Df=[rnd() for _ in range(3)]
        # Indices 0=time,1=radial,2..4=sphere. Background frame values.
        k=arr(5);k[0][0]=-2*U*U
        for j in (2,3,4):k[j][j]=f
        dt=arr(5)
        dt[0][0]=-8*rho**3*q+44*rho**2*U*T
        dt[0][1]=dt[1][0]=rho**2*(-26*U*Uy-F(3,2)*fy-2*a*U*U-3*a*f)
        if injection=='omit_f_y':dt[0][1]+=F(3,2)*rho**2*fy;dt[1][0]=dt[0][1]
        for j in (2,3,4):
            dt[0][j]=dt[j][0]=rho**2*(22*U*DU[j-2]-Df[j-2]/2)
            dt[j][j]=(2*rho*q+4*U*T)*rho**2
        dt[1][1]=2*rho*q-8*U*T
        scales=[1/rho,F(1),1/rho,1/rho,1/rho]
        omega=arr(5);omega[1][0]=omega[0][1]=a
        partial=[[dt[i][j]*scales[i]*scales[j]/rho for j in range(5)] for i in range(5)]
        clock=(rho if injection=='wrong_normal_clock' else F(1))
        eq(clock*partial[0][0],-8*q+44*U*T/rho,'normal-frame clock conversion')
        connection=-sum(omega[c][0]*k[c][1]+omega[c][1]*k[0][c] for c in range(5))
        if injection=='omit_frame_connection':connection=F(0)
        eq(connection,2*a*U*U,'normal-frame connection correction')
        cov=[[partial[i][j]-sum(omega[c][i]*k[c][j]+omega[c][j]*k[i][c] for c in range(5))
              for j in range(5)] for i in range(5)]
        eq(cov[0][1],-26*U*Uy-F(3,2)*fy-3*a*f,'mixed normal-frame metric formula')
        eq(cov[1][1],2*q-8*U*T/rho,'radial normal-frame metric formula')
        for j in (2,3,4):
            eq(cov[0][j],22*U*DU[j-2]/rho-Df[j-2]/(2*rho),'angular normal-frame metric formula')
            eq(cov[j][j],2*q+4*U*T/rho,'tangent normal-frame metric formula')
            eq(cov[1][j],F(0),'vanishing mixed spatial normal derivative')
        phit=rho*(Pi2+2*U*Pi1)
        scalar=(Pi2+U*Pi1 if injection=='half_scalar_lapse_cross' else Pi2+2*U*Pi1)
        eq(phit/rho,scalar,'scalar normal-frame lapse cross term')
        count+=14
    return count

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def abs_jet(report,name,order=0):
    row=report['z_jets'][name][order]['exact']
    return max(abs(F(x)) for x in row)

def cap(x):
    require(x>=0,'nonnegative frame majorant')
    scale=F(1)
    while x>=scale:scale*=10
    while x<scale/10 and scale>F(1,10**100):scale/=10
    unit=scale/1000;upper=(x//unit+1)*unit
    require(upper>x,'strict frame cap')
    return {'exact_upper':str(x),'strict_simple_upper':str(upper),'display_only':float(x)}

def bound(inputs):
    prod,bg,bump=inputs['producer'],inputs['background'],inputs['bump']
    require(prod['status']=='PASS_ACTUAL_SECOND_INITIAL_VALUE_AND_FIRST_SPATIAL_BOUNDS','physical initial norm certificate status')
    require(bg['status']=='PASS_EXACT_ACTUAL_ANNULUS_STATIC_COEFFICIENT_JETS','actual annulus background status')
    require(bump['status']=='PASS','actual bump status')
    require(all(row['C0_C1_C2_all_zero'] is True for row in inputs['gauge']['positive_checks']),
            'initial harmonic completion identity certificate')
    require(prod['scope']['full_first_spatial_derivatives'] is True,'full physical spatial first derivatives missing')
    ell=F(1,1000000)
    require(prod['domain']['ell']==str(ell) and bg['domain']['ell']==str(ell) and bump['ell']==str(ell),'registered frame width')
    rlo=min(F(x) for x in bg['z_jets']['rho'][0]['exact'])
    require(rlo>0,'positive background frame radius')
    C=abs_jet(bg,'C');Cz=abs_jet(bg,'C',1)
    B=abs_jet(bg,'B');k=abs_jet(bg,'k');a=abs_jet(bg,'a')
    qz=F(bump['dimensionless_Q_derivative_sup'][1])/ell
    eq(qz,F(bump['physical_q_derivative_sup'][1]),'first bump physical scale')
    U=C;T=C;Uy=(Cz+C*qz)/rlo;DU=2*C
    Pi1=B*(qz+k+1)/rlo
    val=lambda key:F(prod['bounds'][key]['exact_upper'])
    gamma=val('gamma2_components');phi=val('phi2');q=val('K2_components');Pi2=val('Pi2')
    gamma_spatial=val('gamma2_full_spatial_covariant_components')
    phi_spatial=val('phi2_full_spatial_derivative')
    fy=val('gamma2_radial_covariant_components')
    sphere_f_frame=val('gamma2_angular_tangent_covariant_components')
    lapse_value=2*U*U
    lapse_radial=4*U*Uy;lapse_angular=4*U*DU/rlo
    spatial=max(gamma_spatial,lapse_radial,lapse_angular)
    normal00=8*q+44*U*T/rlo
    normal0y=26*U*Uy+F(3,2)*fy+3*a*gamma
    normal0a=22*U*DU/rlo+sphere_f_frame/2
    normalyy=2*q+8*U*T/rlo
    normalab=2*q+4*U*T/rlo
    scalar_normal=Pi2+2*U*Pi1
    full_metric=max(gamma,lapse_value)
    full_derivative=max(spatial,normal00,normal0y,normal0a,normalyy,normalab)
    full_scalar_derivative=max(phi_spatial,scalar_normal)
    result={'metric2_components':full_metric,'metric2_full_initial_covariant_derivative_components':full_derivative,
            'scalar2':phi,'scalar2_full_initial_covariant_derivative_components':full_scalar_derivative,
            'metric2_time_time_value':lapse_value,'metric2_spatial_time_time_radial_derivative':lapse_radial,
            'metric2_spatial_time_time_angular_derivative':lapse_angular,
            'metric2_all_spatial_covariant_derivative_components':spatial,
            'metric2_normal_00':normal00,'metric2_normal_0y':normal0y,'metric2_normal_0a':normal0a,
            'metric2_normal_yy':normalyy,'metric2_normal_ab':normalab,'scalar2_normal_derivative':scalar_normal}
    coarse={'metric2_components':F(106,100)*10**24,
            'metric2_full_initial_covariant_derivative_components':F(64,10)*10**29,
            'scalar2':F(513,100)*10**33,
            'scalar2_full_initial_covariant_derivative_components':F(306,100)*10**39}
    for key,value in coarse.items():require(result[key]<value,'convenient full frame theorem cap: '+key)
    return {'bounds':{key:cap(value) for key,value in result.items()},
            'convenient_strict_caps':{key:str(value) for key,value in coarse.items()},
            'first_tangent_majorants':{key:str(value) for key,value in
                    {'U':U,'T':T,'U_y':Uy,'sphere_gradient_U':DU,'Pi1':Pi1,'rho_lower':rlo,'abs_a':a}.items()}}

def main():
    root=Path(__file__).resolve().parent;ap=argparse.ArgumentParser()
    ap.add_argument('--pins',type=Path,default=root/'INPUT_PINS.json')
    ap.add_argument('--input-root',type=Path,default=root)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--run-controls',action='store_true')
    mutations={'wrong_normal_clock':'normal-frame clock conversion',
               'omit_frame_connection':'normal-frame connection correction',
               'omit_f_y':'mixed normal-frame metric formula',
               'half_scalar_lapse_cross':'scalar normal-frame lapse cross term'}
    ap.add_argument('--inject-error',choices=['none']+list(mutations),default='none')
    a=ap.parse_args();count=identity_checks(a.inject_error)
    pins=json.loads(a.pins.read_text());inputs={};hashes={}
    for key,pin in pins.items():
        path=a.input_root/pin['path'];raw=path.read_bytes()
        require(hashlib.sha256(raw).hexdigest()==pin['sha256'] and len(raw)==pin['bytes'],'frame input authentication: '+key)
        hashes[key]=pin['sha256']
        if path.suffix=='.json':inputs[key]=json.loads(raw)
    require(inputs['producer']['input_sha256']['background']==hashes['background']
            and inputs['producer']['input_sha256']['bump']==hashes['bump'],'frame producer dependency links')
    report={'status':'PASS_COMPLETED_INITIAL_SPACETIME_SECOND_FRAME_1JET',
            'finite_coordinate_frame_checks':count,'input_sha256':hashes,
            'domain':{'time':'t=0 only','radial':'ell<=z<=2ell, with smooth zero extension','ell':'1/1000000','angular':'entire unit S3'},
            'norm':'Maximum absolute individual components in fixed background Lorentzian orthonormal frame (rho^-1 partial_t,partial_y,rho^-1 E_i), including the covariant derivative index. Euclidean absolute values of frame coefficients; no positive Lorentzian energy is asserted.',
            'scope':{'initial_relative_harmonic_completion_through_order_two':True,
                     'parameter_second_derivative_multiplies_epsilon_squared_over_two':True,
                     'all_initial_spacetime_first_covariant_directions':True,
                     'C1_of_gauge_time_velocities':False,'time_evolution_computed':False,
                     'high_Sobolev_or_nonlinear_remainder_bound':False},**bound(inputs)}
    if a.run_controls:
        cases=[]
        with tempfile.TemporaryDirectory(prefix='initial-frame-controls-') as td:
            for mutation,reason in mutations.items():
                output=Path(td)/(mutation+'.json');cmd=[sys.executable]
                if sys.flags.optimize:cmd.append('-O')
                cmd += [str(Path(__file__).resolve()),'--pins',str(a.pins),'--input-root',str(a.input_root),
                        '--inject-error',mutation,'--output',str(output)]
                run=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
                require(run.returncode!=0 and ('ValueError: '+reason) in run.stderr and not output.exists(),
                        'frame mutation intended rejection: '+mutation)
                cases.append({'mutation':mutation,'rejection_prefix':reason,'no_output_report':True})
        report['controls']={'status':'PASS','optimized':bool(sys.flags.optimize),'cases':cases}
    text=json.dumps(report,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
