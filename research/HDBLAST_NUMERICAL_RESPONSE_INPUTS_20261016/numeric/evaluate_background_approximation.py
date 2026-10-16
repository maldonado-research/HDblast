#!/usr/bin/env python3
"""One shared rational affine background approximator and actual C4 errors."""
import argparse,importlib.util,json,math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('initial_profiles',ROOT/'evaluate_initial_profiles.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
I,J=m.I,m.J;ELL=m.ELL

def derived(r,e,p,a,alpha):
 h=m.poly(m.UP[1],e)/p;k=r*(F(5,2)*a-h);C=alpha/r**2;B=3*alpha/(r**3*p)
 V=r*r*(2*h*h-10*a*h+F(35,4)*a*a-m.poly(m.UP[2],e)-p*p/6)
 return {'rho':r,'eta':e,'p':p,'a':a,'alpha':alpha,'h':h,'k':k,'C':C,'B':B,'static_master_potential':V}

def run(cells):
 m.need(isinstance(cells,int) and cells>0,'nonempty background collar partition');data,pins=m.load();shell=data['shell'];mid={k:sum(map(F,shell['boxes'][k+'_b']['exact']))/2 for k in ('rho','eta','p','a')}
 r,e,p,a=(mid[k] for k in ('rho','eta','p','a'));h=m.poly(m.UP[1],e)/p
 slopes={'rho':-r*r*a,'eta':-r*p,'p':-r*(m.poly(m.UP[1],e)-4*a*p),'a':1/r+r*p*p/3,'alpha':r*(F(5,2)*a-h)};mid['alpha']=F(1)
 errs={};truecaps={};rows=[]
 for j in range(cells):
  rc=I(F(8*j,cells),F(8*(j+1),cells));z=ELL*rc
  rr,ee,pp,aa,hh,CC,BB,kk,U2,E=m.background(rc,data,'none');actual=derived(rr,ee,pp,aa,CC*rr**2)
  # A single common z variable defines all five hat states. Coefficients are
  # evaluated from these shared functions, not independently selected boxes.
  hat={k:J([mid[k]+slopes[k]*z,slopes[k],0,0,0]) for k in mid};approx=derived(*(hat[k] for k in ('rho','eta','p','a','alpha')))
  m.need(approx['rho'].val().lo>0 and approx['p'].val().hi<0,'approximate branch denominators')
  for key in actual:
   if key not in errs:errs[key]=[F(0)]*5;truecaps[key]=[F(0)]*5
   for n in range(5):
    difference=actual[key].c[n]-approx[key].c[n];error=ELL**n*math.factorial(n)*difference.mag();errs[key][n]=max(errs[key][n],error);truecaps[key][n]=max(truecaps[key][n],ELL**n*math.factorial(n)*actual[key].c[n].mag())
  rows.append({'scaled_r':rc.row(),'rho_actual':actual['rho'].val().row(),'p_actual':actual['p'].val().row(),'C_actual':actual['C'].val().row(),'B_actual':actual['B'].val().row()})
 m.need(rows[0]['scaled_r'][0]=='0' and F(rows[-1]['scaled_r'][1])==8 and all(F(rows[k]['scaled_r'][1])==F(rows[k+1]['scaled_r'][0]) for k in range(len(rows)-1)),'gap-free whole background collar cover')
 return {'status':'PASS_CORRELATED_RATIONAL_BACKGROUND_APPROXIMATION_SCALED_C4_ERRORS','domain':{'ell':str(ELL),'scaled_r':['0','8'],'physical_z':['0',str(8*ELL)],'time':'static background'},'definition':{'shared_affine_state':'Xhat(z)=shell_midpoint+z*F(shell_midpoint) for rho,eta,p,a; alphahat(z)=1+z*k(shell_midpoint). Derived coefficients use these same functions.','shell_midpoint':{k:str(v) for k,v in mid.items()},'physical_z_slopes':{k:str(v) for k,v in slopes.items()},'derived':'h=U1(eta)/p,k=rho*(5a/2-h),C=alpha/rho²,B=3alpha/(rho³p),Vstatic=rho²*(2h²-10ah+35a²/4-U2-p²/6)','master_potential':'W(t,z)=Vstatic(z)-35 tanh²(t)/4; temporal term retained exactly in true and hat, so errors cancel.'},'input_sha256':{k:v['sha256'] for k,v in pins.items()},'proof':'At each complete cell, actual shell+z*F(outer) is valid by authenticated first-exit enclosure. Differentiate actual ODE through4 and the stated rational hat functions separately. Their interval difference, multiplied by ell^n*n!, encloses D_r^n(true-hat), uniformly including cell joins. Actual alpha uses its integrated negative k bound and alpha(0)=1. No derivative of point errors is taken.','uniform_actual_minus_approximation_scaled_radial_derivative_caps':{key:[{'exact_upper':str(v),'display_only':float(v)} for v in vals] for key,vals in errs.items()},'actual_scaled_derivative_magnitude_caps':{key:list(map(str,vals)) for key,vals in truecaps.items()},'whole_cells':cells,'cell_enclosures':rows,'scope':{'actual_background_error_certified':True,'common_correlated_functions_defined':True,'coefficient_correlations_do_not_imply_independent_interval_choices':True,'accurate_regular_branch_shell_solve_performed':False,'first_master_or_gauge_refined_using_this_background':False,'second_response_evaluated':False}}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--cells',type=int,default=64);a=p.parse_args();m.need(not a.output.exists(),'refuse overwrite');r=run(a.cells);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'status':r['status'],'value_errors':{k:v[0]['display_only'] for k,v in r['uniform_actual_minus_approximation_scaled_radial_derivative_caps'].items()}},indent=2))
if __name__=='__main__':main()
