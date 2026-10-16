#!/usr/bin/env python3
"""Actual master and forced gauge C3 approximation errors; exact rationals."""
import argparse, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
ELL=F(1,10**6); X=F(3); TH=X/2

def need(ok,msg):
 if not ok: raise ValueError(msg)
def load():
 pins=json.loads((ROOT/'INPUT_PINS.json').read_text());out={}
 for key,pin in pins.items():
  raw=(ROOT/pin['path']).read_bytes()
  need(len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256'],'input authentication: '+key)
  out[key]=json.loads(raw)
 return out,pins

def mv(A,v):return [sum((a*b for a,b in zip(row,v)),F(0)) for row in A]
def add(*vs):return [sum(vals,F(0)) for vals in zip(*vs)]
def scale(c,v):return [c*x for x in v]
def inverse(A):
 n=len(A);v=[row[:]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
 for i in range(n):
  need(v[i][i]!=0,'nonzero exact inverse pivot');p=v[i][i];v[i]=[x/p for x in v[i]]
  for j in range(n):
   if i!=j:
    p=v[j][i];v[j]=[a-p*b for a,b in zip(v[j],v[i])]
 out=[row[n:] for row in v]
 need(all(sum(A[i][k]*out[k][j] for k in range(n))==F(i==j) for i in range(n) for j in range(n)),'exact inverse identity')
 need(all(x>=0 for row in out for x in row),'nonnegative resolvent')
 return out

def mag(row):return max(abs(F(x)) for x in row['exact'])
def pack(x):return {'exact_upper':str(x),'display_only':float(x)}
def vector(v):return [pack(x) for x in v]
def run(injection='none'):
 inp,pins=load();master=inp['master'];co=inp['coefficients'];bump=inp['bump']
 need(master['status']=='PASS' and co['status']=='PASS_WHOLE_COLLAR_GAUGE_COEFFICIENT_AND_SOURCE_CERTIFICATE','actual input status')
 need(co['input_sha256']['HIGHER_JETS_NORMAL.json']==pins['master']['sha256'],'correlated master coefficient dependency')
 m=master['exact']['interior'];bd=master['exact']['boundary'];f=F(m['forcing']['f']);Q=list(map(F,bump['dimensionless_Q_derivative_sup']))
 need(F(m['ell'])==ELL and m['scaled_domain']['x_max']=='3' and m['scaled_domain']['x_plus_r_max']==8,'actual scaled domain')
 # Error is relative to the displacement-only odd d'Alembert wave Wbar.
 # These are differentiated integral identities, never derivatives of C0 errors.
 velocity=[X*ELL,ELL,ELL*(Q[1]+ELL),ELL*(Q[2]+2*ELL*Q[1]+ELL**2)]
 lift=list(map(F,[bd['y'],bd['yprime'],m['ysecond'],m['ythird']]))
 kernel=[X**2*f/2,F(m['first_error']),max(map(F,m['second_error'].values())),max(map(F,m['third_error'].values()))]
 if injection=='drop_dynamic_lift':lift=[F(0)]*4
 ew=[velocity[k]+lift[k]+kernel[k] for k in range(4)]
 need(ew[3]>=F(m['ythird']),'dynamic boundary third-jet retained')
 need(ew[2]>=ELL*Q[1]+F(m['ysecond']),'nonzero velocity radial jet retained')
 actual=list(map(F,master['source_bound_adapter']['scaled_order_bounds']))
 # Explicit radial gauge comparison P. Each order bound uses support |s|<=2.
 c=6*ELL/130;pcap=[5*c]+[c*(5*Q[k]+k*Q[k-1]) for k in range(1,4)]
 # Independent endpoint-ray commutator: pure characteristic terms cancel.
 # On each signed branch Qo(s)=sign(s)*s^4; Pe(s)=s^5.
 for sgp,sgm in ((1,1),(1,-1)):
  for xx,rr in ((F(1,2),F(3,4)),(F(3,2),F(1,4))):
   sp,sm=rr+xx,rr-xx
   qp,qm=sgp*4*sp**3,sgm*4*sm**3
   qpp,qmm=sgp*12*sp**2,sgm*12*sm**2
   pp,pm=20*sp**3,20*sm**3
   pxx=c/4*(pp+pm-rr*(qpp+qmm))
   prr=c/4*(pp+pm-2*(qp+qm)-rr*(qpp+qmm))
   rhs=c*(qp+qm)/2
   if injection=='half_comparator':pxx/=2;prr/=2
   need(pxx-prr==rhs,'explicit product-profile commutator')
 mats={}
 for name in ('A','B','D'):
  mats[name]={}
  for i,j in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2)):
   mats[name][i,j]=[[mag(entry[f'{i},{j}']) for entry in row] for row in co['matrix_mixed_derivative_intervals'][name]]
 A,B,D=(mats[k] for k in ('A','B','D'))
 ps=co['source_profile_mixed_derivative_intervals']['P_s'];p0=co['source_profile_mixed_derivative_intervals']['P_0']
 delta={}
 for i,j in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2)):
  row=ps[f'{i},{j}'];delta[i,j]=max(abs(F(v)-(c if i+j==0 else 0)) for v in row['exact'])
 residual={};source_defect={}
 for i,j in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2)):
  sd=F(0);rv=[F(0)]*3
  for a in range(i+1):
   for b in range(j+1):
    mult=math.comb(i,a)*math.comb(j,b);k=i+j-a-b
    sd+=mult*(delta[a,b]*Q[k+1]+mag(ps[f'{a},{b}'])*ew[k+1]+mag(p0[f'{a},{b}'])*actual[k])
    for channel in range(3):
     rv[channel]+=mult*((A[a,b][channel][0]+B[a,b][channel][0])*pcap[k+1]+D[a,b][channel][0]*pcap[k])
  source_defect[i,j]=sd;rv[0]+=sd;residual[i,j]=rv
 M=[[X*(A[0,0][i][j]+B[0,0][i][j]+TH*D[0,0][i][j]) for j in range(3)] for i in range(3)]
 q=max(map(sum,M));need(q<1,'true gauge error contraction');R=inverse([[F(i==j)-M[i][j] for j in range(3)] for i in range(3)])
 L=lambda key:[[A[key][i][j]+B[key][i][j]+TH*D[key][i][j] for j in range(3)] for i in range(3)]
 r0,rx,rr,rxx=(residual[k] for k in ((0,0),(1,0),(0,1),(2,0)))
 n0=mv(R,scale(X,r0));n1=mv(R,add(scale(2,r0),scale(X,add(mv(L((1,0)),n0),rx))))
 h1=add(mv(A[0,0],r0),rx)
 if injection=='drop_initial_offset':h1=rx
 need(h1==add(mv(A[0,0],r0),rx),'stabilized third initial jet A(0) residual(0)+residual_x(0)')
 gw=[max(2*r0[i]/X+2*h1[i],rr[i]+h1[i]) for i in range(3)]
 n2=mv(R,add(gw,scale(X,add(scale(2,mv(L((1,0)),n1)),mv(L((2,0)),n0),rxx))))
 errors={'0,0':scale(TH,n0),'1,0':n0,'0,1':n0,'2,0':[min(n1[i],TH*n2[i]) for i in range(3)],'1,1':n1,'3,0':n2,'2,1':n2}
 errors['0,2']=add(errors['2,0'],mv(A[0,0],n0),mv(B[0,0],n0),mv(D[0,0],errors['0,0']),r0)
 errors['1,2']=add(n2,mv(A[1,0],n0),mv(A[0,0],errors['2,0']),mv(B[1,0],n0),mv(B[0,0],n1),mv(D[1,0],errors['0,0']),mv(D[0,0],n0),rx)
 radial=rr if injection!='drop_radial_residual' else [F(0)]*3
 need(radial==rr,'third radial recovery retains residual_r')
 errors['0,3']=add(n2,mv(A[0,1],n0),mv(A[0,0],n1),mv(B[0,1],n0),mv(B[0,0],errors['0,2']),mv(D[0,1],errors['0,0']),mv(D[0,0],n0),radial)
 return {'status':'PASS_ACTUAL_MASTER_AND_GAUGE_CONTINUOUS_SCALED_C3_APPROXIMATION_ERRORS','domain':{'ell':str(ELL),'scaled':'0<=x<=3,r>=0,x+r<=8','physical':'t=ell*x,z=ell*r','includes_shell_and_reflection_transitions':True},'conditional_premises':'Retained smooth actual background/master/gauge solution and authenticated integral/reflection identities. Input PDE certificates are authenticated premises, not re-proved here.','input_sha256':{k:v['sha256'] for k,v in pins.items()},'comparison_functions':{'master':'Wbar=(Qo(r+x)+Qo(r-x))/2, Qo(s)=sign(s)Q(abs(s))','gauge':'(P,0,0), P=c/4[Pe(r+x)+Pe(r-x)-r(Qo(r+x)+Qo(r-x))], Pe(s)=s*Q(abs(s)), c=6ell/130','initial_velocity_in_master_comparison':'Zero; true ell exp(-ell*r)Q velocity is retained in error, including its first two radial derivatives.'},'master_error_budget_by_scaled_total_order':{'velocity':list(map(str,velocity)),'dynamic_Dirichlet_lift':list(map(str,lift)),'reflected_Duhamel':list(map(str,kernel)),'actual_minus_comparison':[pack(x) for x in ew]},'gauge_comparison_absolute_derivative_caps_by_total_order':list(map(str,pcap)),'source_approximation_defect_mixed_scaled_derivative_caps':{f'{i},{j}':pack(v) for (i,j),v in source_defect.items()},'complete_gauge_residual_mixed_scaled_derivative_caps':{f'{i},{j}':vector(v) for (i,j),v in residual.items()},'gauge_error_resolvent':{'matrix':[[str(v) for v in row] for row in M],'row_norm':str(q),'n0':list(map(str,n0)),'n1':list(map(str,n1)),'n2':list(map(str,n2)),'stabilized_initial_third_jet_majorant':list(map(str,h1))},'actual_minus_gauge_comparison_mixed_scaled_derivative_caps':{k:vector(v) for k,v in errors.items()},'scope':{'actual_approximation_errors_not_only_magnitudes':True,'actual_master_dynamic_boundary_retained':True,'derivatives_obtained_from_differentiated_integral_PDE':True,'numerical_fixed_point_iterates_evaluated':False,'correlated_background_point_approximation':False,'common_gauge_quadratic_sources_evaluated':False,'second_response_evaluated':False,'invariant_force':False}}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--inject-error',default='none',choices=['none','drop_dynamic_lift','half_comparator','drop_initial_offset','drop_radial_residual']);a=p.parse_args()
 need(not a.output.exists(),'refuse overwrite');r=run(a.inject_error);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'status':r['status'],'master_C3_error':[v['display_only'] for v in r['master_error_budget_by_scaled_total_order']['actual_minus_comparison']],'gauge_C3_error_by_channel':[max(row[i]['display_only'] for k,row in r['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'].items() if sum(map(int,k.split(',')))==3) for i in range(3)]},indent=2))
if __name__=='__main__':main()
