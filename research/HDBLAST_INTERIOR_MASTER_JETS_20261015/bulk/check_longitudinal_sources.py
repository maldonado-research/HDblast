#!/usr/bin/env python3
"""Exact direct-coordinate Hessian checks and actual longitudinal source bounds.

No prior scientific module is imported. The finite generic-jet checks supplement
the analytic identities; the actual source enclosure uses authenticated global
master/background norm certificates, rather than values at those generic jets.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import sys

N=5
I=range(N)
A=range(1,N)

def arr(shape):
    return [arr(shape[1:]) for _ in range(shape[0])] if shape else F(0)

def eq(a,b,label):
    if a!=b:
        raise ValueError(label+': exact residual '+str(a-b))

def require(condition,label):
    if not condition:
        raise ValueError(label)

def mm(a,b):
    return [[sum(a[i][k]*b[k][j] for k in I) for j in I] for i in I]

def direct_path(g,h,dg,dh,ddg,ddh):
    """Independent epsilon-series inverse, connection, Ricci, scalar box."""
    inv0=[[1/g[i][i] if i==j else F(0) for j in I] for i in I]
    ih=mm(mm(inv0,h),inv0)
    inv=[inv0,[[-v for v in row] for row in ih],mm(mm(ih,h),inv0)]
    gd=[dg,dh,arr([N,N,N])]
    gdd=[ddg,ddh,arr([N,N,N,N])]
    dinv=arr([3,N,N,N])
    for k,l,a,b in itertools.product(range(3),I,I,I):
        dinv[k][l][a][b]=-sum(inv[r][a][c]*gd[s][l][c][d]*inv[k-r-s][d][b]
            for r in range(k+1) for s in range(k-r+1) for c in I for d in I)
    gamma,dgamma=arr([3,N,N,N]),arr([3,N,N,N,N])
    for k,a,b,c in itertools.product(range(3),I,I,I):
        gamma[k][a][b][c]=sum(inv[r][a][d]*(gd[k-r][b][c][d]+gd[k-r][c][b][d]-gd[k-r][d][b][c])
            for r in range(k+1) for d in I)/2
        for l in I:
            dgamma[k][l][a][b][c]=sum(
                dinv[r][l][a][d]*(gd[k-r][b][c][d]+gd[k-r][c][b][d]-gd[k-r][d][b][c])
                +inv[r][a][d]*(gdd[k-r][l][b][c][d]+gdd[k-r][l][c][b][d]-gdd[k-r][l][d][b][c])
                for r in range(k+1) for d in I)/2
    r2=arr([N,N])
    for a,b in itertools.product(I,I):
        r2[a][b]=2*(sum(dgamma[2][c][c][a][b]-dgamma[2][b][c][a][c] for c in I)
            +sum(gamma[r][c][c][d]*gamma[2-r][d][a][b]-gamma[r][c][b][d]*gamma[2-r][d][a][c]
                 for r in range(3) for c in I for d in I))
    return inv,gamma,r2

def generic_checks(seed,injection):
    rng=random.Random(seed)
    rnd=lambda:F(rng.randint(-7,7),rng.randint(1,6))
    rho,a,ay=F(7,3),rnd(),rnd()
    g=arr([N,N]);g[0][0]=F(1)
    for j in A:g[j][j]=rho*rho*F(-j if j==1 else j)
    dg,ddg=arr([N,N,N]),arr([N,N,N,N])
    for i in A:
        dg[0][i][i]=2*a*g[i][i]
        ddg[0][0][i][i]=2*(ay+2*a*a)*g[i][i]
    for i in A:
        for j in range(i,N):
            for k in A:
                v=rnd()
                dg[k][i][j]=dg[k][j][i]=v
                for target in ((0,k,i,j),(k,0,i,j),(0,k,j,i),(k,0,j,i)):
                    x,y,z,w=target;ddg[x][y][z][w]=2*a*v
                for l in range(k,N):
                    v2=rnd()
                    for x,y in ((k,l),(l,k)):
                        ddg[x][y][i][j]=ddg[x][y][j][i]=v2
    psi=rnd();ps1=[rnd() for _ in I];ps2=arr([N,N])
    ch=rnd();ch1=[rnd() for _ in I];ch2=arr([N,N])
    for i in I:
        for j in range(i,N):
            ps2[i][j]=ps2[j][i]=rnd()
            ch2[i][j]=ch2[j][i]=rnd()
    h,dh,ddh=arr([N,N]),arr([N,N,N]),arr([N,N,N,N])
    h[0][0]=-4*psi
    for k in I:
        dh[k][0][0]=-4*ps1[k]
        for l in I:ddh[k][l][0][0]=-4*ps2[k][l]
    for i,j in itertools.product(A,A):
        h[i][j]=2*psi*g[i][j]
        for k in I:
            dh[k][i][j]=2*(ps1[k]*g[i][j]+psi*dg[k][i][j])
            for l in I:
                ddh[k][l][i][j]=2*(ps2[k][l]*g[i][j]+ps1[k]*dg[l][i][j]+ps1[l]*dg[k][i][j]+psi*ddg[k][l][i][j])
    inv,gamma,directR=direct_path(g,h,dg,dh,ddg,ddh)
    psH=[[ps2[i][j]-sum(gamma[0][c][i][j]*ps1[c] for c in A) for j in I] for i in I]
    box=sum(inv[0][i][j]*psH[i][j] for i in A for j in A)
    grad2=sum(inv[0][i][j]*ps1[i]*ps1[j] for i in A for j in A)
    radial=8*ps1[0]**2+104*a*psi*ps1[0]+8*psi*ps2[0][0]+48*(ay+4*a*a)*psi**2
    candidate=arr([N,N])
    candidate[0][0]=-8*ps1[0]**2-32*a*psi*ps1[0]+16*psi*ps2[0][0]+16*grad2-8*psi*box
    for i in A:
        candidate[0][i]=candidate[i][0]=12*psi*ps2[0][i]-48*a*psi*ps1[i]
        for j in A:
            candidate[i][j]=12*ps1[i]*ps1[j]+24*psi*psH[i][j]+g[i][j]*(4*(psi*box+grad2)-radial)
    if injection=='half_ricci':candidate=[[v/2 for v in row] for row in candidate]
    if injection=='drop_mixed_radial':
        for i in A:candidate[0][i]+=48*a*psi*ps1[i]
    p,py,U1,U2,U3=rnd(),rnd(),rnd(),rnd(),rnd()
    count=0
    for i in I:
        for j in range(i,N):
            matter=-2*ch1[i]*ch1[j]-F(2,3)*(U2*ch*ch*g[i][j]+2*U1*ch*h[i][j])
            eq(candidate[i][j]+matter,directR[i][j]+matter,'longitudinal Einstein Hessian '+str((i,j)))
            count+=1
    ph1=[p if i==0 else F(0) for i in I]
    ph2=[[py if i==0 and j==0 else F(0) for j in I] for i in I]
    directS=2*sum(inv[2][i][j]*ph2[i][j]+inv[1][i][j]*ch2[i][j]
        -sum(inv[r][i][j]*gamma[2-r][c][i][j]*ph1[c] for r in range(3) for c in I)
        -sum(inv[r][i][j]*gamma[1-r][c][i][j]*ch1[c] for r in range(2) for c in I)
        for i in I for j in I)-U3*ch*ch
    chH=[[ch2[i][j]-sum(gamma[0][c][i][j]*ch1[c] for c in A) for j in I] for i in I]
    chbox=sum(inv[0][i][j]*chH[i][j] for i in A for j in A)
    fullS=(8*psi*(ch2[0][0]+4*a*ch1[0])-4*psi*chbox+12*ps1[0]*ch1[0]
           +48*p*psi*ps1[0]+32*(py+4*a*p)*psi*psi-U3*ch*ch)
    eq(fullS,directS,'full longitudinal scalar Hessian')
    # On the background U1=py+4ap. Enforce the first scalar equation only.
    U1=py+4*a*p
    chyy=-4*a*ch1[0]-chbox+U2*ch-6*p*ps1[0]-4*U1*psi
    full_on=(8*psi*(chyy+4*a*ch1[0])-4*psi*chbox+12*ps1[0]*ch1[0]
             +48*p*psi*ps1[0]+32*U1*psi*psi-U3*ch*ch)
    reduced=-12*psi*chbox+8*U2*psi*ch+12*ps1[0]*ch1[0]-U3*ch*ch
    if injection=='wrong_scalar_wave_factor':reduced+=8*psi*chbox
    eq(full_on,reduced,'on-shell scalar cancellation')
    return count+2

def interval_abs(item):
    return max(abs(F(x)) for x in item['exact'])

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def actual_bounds(background,jet_bounds,injection):
    b=background['background_intervals'];c=background['reconstruction_coefficients']
    rho_min=min(F(x) for x in b['rho']['exact'])
    require(rho_min>0,'positive radius required')
    av=interval_abs(b['a']);ay=interval_abs(b['a_z'])/rho_min
    U1,U2,U3=(interval_abs(b[k]) for k in ('U_1','U_2','U_3'))
    C,Cz,Czz=(interval_abs(c[k]) for k in ('u_w_coefficient','u_w_coefficient_z','u_w_coefficient_zz'))
    B,Bz=(interval_abs(c[k]) for k in ('v_wz_coefficient','v_wz_coefficient_z'))
    k,kz=(interval_abs(c[q]) for q in ('alpha_log_z','alpha_log_zz'))
    ell=F(1,1000000);T=3*ell
    W0,W1,W2,W3=(F(x) for x in jet_bounds)
    if injection=='omit_physical_scale':ell=F(1)
    D1,D2,D3=W1/ell,W2/ell**2,W3/ell**3
    for j,(physical,scaled) in enumerate(((D1,W1),(D2,W2),(D3,W3)),1):
        eq(physical*F(1,1000000)**j,scaled,'physical order-'+str(j)+' derivative scale')
    # c(t)=cosh(t)^(-3/2), |c|<=1, |c'/c|<=3T/2,
    # |c''/c|=|-3/2+15tanh(t)^2/4|<=3/2 on this window.
    kt=3*T/2;ktt=F(3,2)
    u=C*W0
    ut=C*(D1+kt*W0)
    utt=C*(D2+2*kt*D1+ktt*W0)
    uz=C*D1+Cz*W0
    uzz=C*D2+2*Cz*D1+Czz*W0
    uyt=(C*D2+Cz*D1+kt*(C*D1+Cz*W0))/rho_min
    uy=uz/rho_min;uyy=uzz/rho_min**2+av*uz/rho_min
    v=B*(D1+k*W0)
    vt=B*(D2+k*D1+kt*(D1+k*W0))
    vtt=B*(D3+k*D2+2*kt*(D2+k*D1)+ktt*(D1+k*W0))
    vz=Bz*(D1+k*W0)+B*(D2+k*D1+kz*W0)
    vy=vz/rho_min
    gu=max(ut,2*u)/rho_min
    gv=max(vt,2*v)/rho_min
    guy=max(uyt,2*uy)/rho_min
    hu=max(utt,2*(ut+T*u),4*u+T*ut)/rho_min**2
    boxu=(utt+3*T*ut+8*u)/rho_min**2
    normdu=(ut*ut+4*u*u)/rho_min**2
    radial=8*uy*uy+104*av*u*uy+8*u*uyy+48*(ay+4*av*av)*u*u
    syy=(8*uy*uy+32*av*u*uy+16*u*uyy+16*normdu+8*u*boxu
         +2*vy*vy+F(2,3)*U2*v*v+F(16,3)*U1*u*v)
    sya=48*av*u*gu+12*u*guy+2*vy*gv
    sab=(12*gu*gu+24*u*hu+4*(u*boxu+normdu)+radial+2*gv*gv
         +F(2,3)*U2*v*v+F(8,3)*U1*u*v)
    scalar=(12*u*(vtt+3*T*vt+8*v)/rho_min**2+8*U2*u*v+12*uy*vy+U3*v*v)
    result={'u':u,'u_t':ut,'u_tt':utt,'u_y':uy,'u_yy':uyy,'u_yt':uyt,
            'v':v,'v_t':vt,'v_tt':vtt,'v_y':vy,
            'Einstein_nn':syy,'Einstein_na_each':sya,'Einstein_ab_each':sab,'scalar':scalar}
    # Coarse rational displays are accepted only after all exact arithmetic.
    caps={'u':F(1,16000),'u_t':1000,'u_tt':16000000000,'u_y':8,
          'u_yy':1000000,'u_yt':125000000,'v':190000,'v_t':3000000000000,
          'v_tt':80000000000000000000,'v_y':24000000000,
          'Einstein_nn':1200000000000000000000,'Einstein_na_each':1200000000000000000000,
          'Einstein_ab_each':1200000000000000000000,'scalar':10000000000000}
    for key,cap in caps.items():
        require(0<=result[key]<cap,'actual source/reconstruction cap failed: '+key)
    return {key:{'exact_upper':str(value),'strict_simple_upper':str(caps[key]),'display_only':float(value)}
            for key,value in result.items()}

def main():
    ap=argparse.ArgumentParser()
    capsule=Path(__file__).resolve().parents[1]
    ap.add_argument('--background',type=Path,default=capsule/'potential/ACTUAL_POTENTIAL_JET_BOUNDS.json')
    ap.add_argument('--jets',type=Path,default=capsule/'jets/HIGHER_JETS_NORMAL.json')
    ap.add_argument('--pins',type=Path,default=Path(__file__).resolve().parent/'INPUT_PINS.json')
    ap.add_argument('--identities-only',action='store_true')
    ap.add_argument('--run-controls',action='store_true')
    ap.add_argument('--output',type=Path)
    ap.add_argument('--inject-error',choices=['none','half_ricci','drop_mixed_radial','wrong_scalar_wave_factor','omit_physical_scale'],default='none')
    args=ap.parse_args()
    checks=sum(generic_checks(seed,args.inject_error) for seed in (271,1927,6469))
    report={'status':'PASS_EXACT_LONGITUDINAL_SOURCE_IDENTITIES','generic_exact_checks':checks,
            'seeds':[271,1927,6469],
            'scope':{'registered_source_enclosed':False,'all_fifteen_Einstein_components':True,
                     'full_scalar_Hessian':True,'source_representative':'longitudinal, first shell displacement zero',
                     'harmonic_gauge_source_enclosed':False,'constraint_compatibility_proved_by_finite_jets':False,
                     'second_order_solution_computed':False,'gauge_invariant_force_computed':False}}
    if not args.identities_only:
        pins=json.loads(args.pins.read_text())
        for key,path in (('background',args.background),('master_jets',args.jets)):
            require(sha(path)==pins[key]['sha256'] and path.stat().st_size==pins[key]['bytes'],
                    'actual source input authentication failed: '+key)
        background=json.loads(args.background.read_text());jets=json.loads(args.jets.read_text())
        require(background['status']=='PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ','actual potential certificate missing')
        require(jets['status']=='PASS' and jets['scope']['third_master_derivatives_certified'] is True
                and jets['scope']['second_master_derivatives_certified'] is True
                and jets['scope']['Wzz_input_pending'] is False,'actual master C3 certificate missing')
        require(jets['potential_report_sha256']==sha(args.background),'master/background dependency pin mismatch')
        # An explicit adapter contract is intentionally kept small and audited.
        values=jets['source_bound_adapter']['scaled_order_bounds']
        require(jets['source_bound_adapter']['ell']=='1/1000000','master scale mismatch')
        require(jets['source_bound_adapter']['domain']=='0<=x<=3, r>=0, x+r<=8','master domain mismatch')
        require(len(values)==4,'need four global master order norms')
        result=actual_bounds(background,values,args.inject_error)
        report.update(status='PASS_ACTUAL_LONGITUDINAL_BULK_SOURCE_ENVELOPES',bounds=result,
                      input_sha256={'background':sha(args.background),'master_jets':sha(args.jets)},
                      input_paths={'background':pins['background']['path'],'master_jets':pins['master_jets']['path']},
                      domain={'physical':'0<=t<=3ell, z>=0, t+z<=8ell','ell':'1/1000000','angular':'entire unit S3'},
                      master_scaled_order_bounds=values,
                      norm='Absolute individual components in the background warped orthonormal frame (n,e0,e1,e2,e3); scalar absolute value. No sum or positive Einstein energy is implied.')
        report['scope']['registered_source_enclosed']=True
    if args.run_controls:
        controls=[]
        for mutation in ('half_ricci','drop_mixed_radial','wrong_scalar_wave_factor','omit_physical_scale'):
            cmd=[sys.executable]
            if sys.flags.optimize:cmd.append('-O')
            cmd.extend([str(Path(__file__).resolve()),'--background',str(args.background),
                        '--jets',str(args.jets),'--pins',str(args.pins),'--inject-error',mutation])
            result=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
            require(result.returncode!=0 and 'ValueError:' in result.stderr,
                    'corruption was not explicitly rejected: '+mutation)
            controls.append({'mutation':mutation,'returncode':result.returncode,
                             'rejection':result.stderr.split('ValueError:')[-1].strip()})
        report['controls']={'status':'PASS','optimized':bool(sys.flags.optimize),'cases':controls}
    raw=json.dumps(report,indent=2)+'\n'
    if args.output:args.output.write_text(raw)
    print(raw,end='')

if __name__=='__main__':main()
