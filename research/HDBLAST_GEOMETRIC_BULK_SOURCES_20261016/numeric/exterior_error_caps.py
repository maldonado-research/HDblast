#!/usr/bin/env python3
"""Differentiated error majorants on backward-closed exterior r>=x.

The exact master lift is y(x-r) for x>=r and zero for x<=r.  Its jet
through3 is zero on r>=x (flat compatible corner).  Every predecessor
Duhamel, coefficient, and derivative majorant stays valid on this smaller
domain.  Its backward cones never meet the shell at positive time, so the
same error resolvent estimates apply without a physical boundary term.
"""
import argparse,importlib.util,json,math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('inherited_error_algebra',ROOT/'check_inherited_errors.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
ELL=old.ELL;X=old.X;TH=old.TH;need=old.need

def derive(exterior=True):
 inp,pins=old.load();master=inp['master'];co=inp['coefficients'];bump=inp['bump']
 m=master['exact']['interior'];bd=master['exact']['boundary'];f=F(m['forcing']['f']);Q=list(map(F,bump['dimensionless_Q_derivative_sup']))
 velocity=[X*ELL,ELL,ELL*(Q[1]+ELL),ELL*(Q[2]+2*ELL*Q[1]+ELL**2)]
 lift=[F(0)]*4 if exterior else list(map(F,[bd['y'],bd['yprime'],m['ysecond'],m['ythird']]))
 kernel=[X**2*f/2,F(m['first_error']),max(map(F,m['second_error'].values())),max(map(F,m['third_error'].values()))]
 ew=[velocity[k]+lift[k]+kernel[k] for k in range(4)]
 actual=list(map(F,master['source_bound_adapter']['scaled_order_bounds']))
 c=6*ELL/130;pcap=[5*c]+[c*(5*Q[k]+k*Q[k-1]) for k in range(1,4)]
 mats={name:{key:[[old.mag(entry[f'{key[0]},{key[1]}']) for entry in row] for row in co['matrix_mixed_derivative_intervals'][name]] for key in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))} for name in ('A','B','D')}
 A,B,D=(mats[k] for k in ('A','B','D'))
 ps=co['source_profile_mixed_derivative_intervals']['P_s'];p0=co['source_profile_mixed_derivative_intervals']['P_0']
 delta={key:max(abs(F(v)-(c if sum(key)==0 else 0)) for v in ps[f'{key[0]},{key[1]}']['exact']) for key in A}
 residual={};source_defect={}
 for i,j in A:
  sd=F(0);rv=[F(0)]*3
  for a in range(i+1):
   for b in range(j+1):
    mult=math.comb(i,a)*math.comb(j,b);k=i+j-a-b
    sd+=mult*(delta[a,b]*Q[k+1]+old.mag(ps[f'{a},{b}'])*ew[k+1]+old.mag(p0[f'{a},{b}'])*actual[k])
    for channel in range(3):rv[channel]+=mult*((A[a,b][channel][0]+B[a,b][channel][0])*pcap[k+1]+D[a,b][channel][0]*pcap[k])
  source_defect[i,j]=sd;rv[0]+=sd;residual[i,j]=rv
 M=[[X*(A[0,0][i][j]+B[0,0][i][j]+TH*D[0,0][i][j]) for j in range(3)] for i in range(3)]
 R=old.inverse([[F(i==j)-M[i][j] for j in range(3)] for i in range(3)])
 L=lambda key:[[A[key][i][j]+B[key][i][j]+TH*D[key][i][j] for j in range(3)] for i in range(3)]
 r0,rx,rr,rxx=(residual[k] for k in ((0,0),(1,0),(0,1),(2,0)))
 n0=old.mv(R,old.scale(X,r0))
 n1=old.mv(R,old.add(old.scale(2,r0),old.scale(X,old.add(old.mv(L((1,0)),n0),rx))))
 h1=old.add(old.mv(A[0,0],r0),rx)
 gw=[max(2*r0[i]/X+2*h1[i],rr[i]+h1[i]) for i in range(3)]
 n2=old.mv(R,old.add(gw,old.scale(X,old.add(old.scale(2,old.mv(L((1,0)),n1)),old.mv(L((2,0)),n0),rxx))))
 errors={'0,0':old.scale(TH,n0),'1,0':n0,'0,1':n0,'2,0':[min(n1[i],TH*n2[i]) for i in range(3)],'1,1':n1,'3,0':n2,'2,1':n2}
 errors['0,2']=old.add(errors['2,0'],old.mv(A[0,0],n0),old.mv(B[0,0],n0),old.mv(D[0,0],errors['0,0']),r0)
 errors['1,2']=old.add(n2,old.mv(A[1,0],n0),old.mv(A[0,0],errors['2,0']),old.mv(B[1,0],n0),old.mv(B[0,0],n1),old.mv(D[1,0],errors['0,0']),old.mv(D[0,0],n0),rx)
 errors['0,3']=old.add(n2,old.mv(A[0,1],n0),old.mv(A[0,0],n1),old.mv(B[0,1],n0),old.mv(B[0,0],errors['0,2']),old.mv(D[0,1],errors['0,0']),old.mv(D[0,0],n0),rr)
 return ew,errors,source_defect,residual

def run():
 prior=json.loads((ROOT/'inputs/APPROXIMATION_ERRORS_NORMAL.json').read_text())
 ew0,eu0,sd0,r0=derive(False)
 need(list(map(str,ew0))==[v['exact_upper'] for v in prior['master_error_budget_by_scaled_total_order']['actual_minus_comparison']],'whole-domain master algebra exact replay')
 need({k:list(map(str,v)) for k,v in eu0.items()}=={k:[x['exact_upper'] for x in v] for k,v in prior['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'].items()},'whole-domain gauge algebra exact replay')
 ew,eu,sd,res=derive()
 need(all(a<=b for a,b in zip(ew,ew0)) and all(a<=b for k in eu for a,b in zip(eu[k],eu0[k])),'exterior restriction monotonically improves derivative error caps')
 return {'status':'PASS_BACKWARD_CLOSED_EXTERIOR_ACTUAL_C3_ERRORS','domain':'0<=x<=3,x<=r,x+r<=8',
         'master_actual_minus_comparison_by_total_order':list(map(str,ew)),
         'gauge_actual_minus_comparison_mixed':{k:list(map(str,v)) for k,v in eu.items()},
         'source_defect_mixed':{f'{i},{j}':str(v) for (i,j),v in sd.items()},
         'complete_gauge_residual_mixed':{f'{i},{j}':list(map(str,v)) for (i,j),v in res.items()},
         'proof':'Exact dynamic lift Ly(x,r)=1_{x>=r}y(x-r), with compatible flat corner, has zero jets through3 on r>=x. Every backward characteristic point (s,q) obeys q>=r-x+s>=s, so this restriction is closed under all Duhamel integrals. All other inherited global differentiated estimates remain upper bounds. The retained positive matrix-resolvent derivative algebra and stabilized initial jets apply on this restriction without positive-time shell contact.',
         'baseline_formula_replay':'exact equality for all 4 master and 30 gauge derivative caps',
         'scope':{'actual_solution_errors':True,'no_numerical_Picard_iteration':True,'second_response_evaluated':False}}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();need(not a.output.exists(),'refuse overwrite')
 r=run();a.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
if __name__=='__main__':main()
