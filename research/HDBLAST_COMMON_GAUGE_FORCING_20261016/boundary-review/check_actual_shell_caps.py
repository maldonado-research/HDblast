#!/usr/bin/env python3
"""Independent exact-rational audit of actual common-gauge shell majorants."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


EXPECTED={
 'GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json':'4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53',
 'HIGHER_GAUGE_NORMS_NORMAL.json':'f246bc6840ea36f123c7c53da28a9c04925954187e0f06df25add779ce083c6d',
 'HIGHER_JETS_NORMAL.json':'92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df',
 'COORDINATE_SHELL_SOURCE_ENCLOSURE.json':'03aba506022e8b5fba38280db9b568d9d1c4fb0c5f76d0f33bf096a0ec8c75f1',
}


def require(test,reason):
 if not test: raise ValueError(reason)


def boundinterval(x):
 lo,hi=map(F,x);require(lo<=hi,'ordered interval');return max(abs(lo),abs(hi))


def readcap(x): return F(x['exact_upper'])


def strictcap(x,denominator=1000000):
 return F(x*denominator//1+1,denominator)


def packed(x):
 return {'exact_upper':str(x),'strict_readable_cap':str(strictcap(x))}


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--shell-root',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 require(not args.output.exists(),'refuse overwrite')
 raw=(args.shell_root/'COMMON_HARMONIC_SHELL_NORMAL.json').read_bytes();rep=json.loads(raw)
 inputs={}
 for name,digest in EXPECTED.items():
  data=(args.shell_root/'inputs'/name).read_bytes()
  require(hashlib.sha256(data).hexdigest()==digest,'authenticated retained input '+name)
  require(rep['input_sha256'][name]==digest,'producer input linkage '+name)
  inputs[name]=json.loads(data)
 c=inputs['GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json'];g=inputs['HIGHER_GAUGE_NORMS_NORMAL.json']
 w=inputs['HIGHER_JETS_NORMAL.json'];s=inputs['COORDINATE_SHELL_SOURCE_ENCLOSURE.json']['exact']
 ell=F(1,1000000);tm=3*ell
 require(c['domain']['ell']==str(ell) and g['domain']['ell']==str(ell),'retained ell')
 for name in ('rho','a','p'):
  require(c['first_exit'][name]['shell']['exact']==s[name],'shared actual shell '+name)
 rlo,rhi=map(F,s['rho']);plo,phi=map(F,s['p'])
 require(0<rlo<=rhi and plo<=phi<0,'positive rho and negative p throughout shell box')
 ar=boundinterval(s['a']);pr=boundinterval(s['p']);sp=2*pr;ss=boundinterval(s['sigma2'])
 invr=1/rlo;cm=1/(1-tm*tm/2)
 def gv(i,j,ch):return F(g['scaled_mixed_partial_caps'][str(i)+','+str(j)][ch])/ell**(i+j)
 bb,bt,btt=gv(0,0,1),gv(1,0,1),gv(2,0,1)
 dd,dt,dtt=gv(0,0,2),gv(1,0,2),gv(2,0,2)
 fr,ftr=gv(0,1,0),gv(1,1,0)
 y,yt,Hin=map(lambda key:F(w['exact']['boundary'][key]),('y','yprime','H'))
 Hit=F(w['exact']['interior']['Hprime']);gam=boundinterval(w['exact']['boundary']['gamma'])
 kz=boundinterval(c['static_z_jets']['k'][0]['exact'])
 # alpha(0)=1; vartheta=cosh(t)^(-3/2)<=1. Every input is the actual master.
 B=3/(rlo**3*abs(phi))
 u=min(boundinterval(s['global_envelopes']['u']),y*invr**2)
 ut=(yt/ell+3*tm*y/2)*invr**2
 vv=min(boundinterval(s['global_envelopes']['v']),B*(Hin/ell+(gam+kz)*y))
 vt=B*(Hit/ell**2+(gam+kz)*yt/ell+3*tm*(Hin/ell+(gam+kz)*y)/2)
 h,ht=4*u+2*fr*invr,4*ut+2*ftr*invr
 x,xt=u+bt,ut+btt;z,zt=u+tm*bb,ut+bb+tm*bt
 j,jt=cm*dt+bb,cm*tm*dt+cm*dtt+bt+tm*bb
 tau=8*u+2*bt+6*tm*bb+16*dd
 taut=8*ut+2*btt+6*bb+6*tm*bt+16*dt
 sq=4*x*x+12*z*z+64*z*dd+96*dd*dd+8*j*j
 sqt=8*x*xt+24*z*zt+64*(zt*dd+z*dt)+192*dd*dt+16*j*jt
 amps={'u':u,'u_t':ut,'v':vv,'v_t':vt,'J':h,'J_t':ht,'b':bb,'b_t':bt,
       'b_tt':btt,'d':dd,'d_t':dt,'d_tt':dtt,'X':x,'X_t':xt,'Z':z,'Z_t':zt,
       'V_gauge':j,'V_gauge_t':jt,'tau':tau,'tau_t':taut,
       'abs_trace_A_squared':sq,'abs_trace_A_squared_t':sqt}
 checks=[]
 def compare(name,ours,theirs):
  require(0<=ours<=readcap(theirs),'producer underestimates '+name)
  require(readcap(theirs)<F(theirs['strict_readable_cap']['exact']) or readcap(theirs)==0,
          'strict readable cap '+name)
  checks.append(name)
 for name,val in amps.items():compare('amplitude/'+name,val,rep['bounds']['actual_field_amplitude_caps'][name])
 A=ar*h*h/4+sp*h*vv/6+ss*vv*vv/6
 At=ar*h*ht/2+sp*(ht*vv+h*vt)/6+ss*vv*vt/3
 L,Lt=ar*h+sp*vv/3,ar*ht+sp*vt/3
 components={'H':h,'eta':vv,'m_00':2*x,'m_0a':2*j,'m_aa':2*z+8*dd,
             'm_ab_a_ne_b':8*dd,'tau':tau,'abs_trace_A_squared':sq}
 for name,val in components.items():compare('component/'+name,val,rep['bounds']['actual_field_component_caps'][name])
 normals={'H_y':4*ar*h+2*ar*tau+4*sp*vv/3,'eta_y':pr*h/2+ss*vv/2,
          'nabla_y_m_diagonal':L,'nabla_y_m_offdiagonal':F(0),
          'coordinate_m_y_in_frame_diagonal':2*ar*max(2*x,2*z+8*dd)+L,
          'coordinate_m_y_in_frame_offdiagonal':2*ar*max(2*j,8*dd)}
 for name,val in normals.items():compare('proper-first-normal-field/'+name,val,rep['bounds']['proper_first_normal_field_caps'][name])
 PP=pr*h*h/4+ss*h*vv/2+2*vv*vv
 PPt=pr*h*ht/2+ss*(ht*vv+h*vt)/2+4*vv*vt
 coeff=(A+2*L*x,L*j,A+2*L*z,2*L*dd)
 coefft=(At+2*Lt*x+2*L*xt,Lt*j+L*jt,At+2*Lt*z+2*L*zt,2*Lt*dd+2*L*dt)
 q00,q0a,qaa,qab=coeff[0],coeff[1],coeff[2]+4*coeff[3],4*coeff[3]
 t00,t0a,taa,tab=[invr*v for v in (coefft[0],coefft[1],coefft[2]+4*coefft[3],4*coefft[3])]
 # Metric compatibility cancels the A gamma term from the angular 0a boosts.
 da00=invr*(2*coeff[0]+2*tm*q0a)
 da0a=invr*(8*coeff[1]+2*tm*L*(bt+tm*bb+4*dd))
 daab=invr*(2*coeff[2]+16*coeff[3]+2*tm*q0a)
 daoff=invr*(16*coeff[3]+2*tm*q0a)
 rows={}
 for i in range(4):
  for jj in range(i,4):
   row=str(i)+str(jj)
   if i==jj==0:v,t,d=q00,t00,da00
   elif i==0:v,t,d=q0a,t0a,da0a
   elif i==jj:v,t,d=qaa,taa,daab
   else:v,t,d=qab,tab,daoff
   rows[row]={'value':v,'proper_time_derivative':t,'each_angular_covariant_derivative':d,
              'each_tangential_covariant_derivative':max(t,d)}
 for family,factor,pscale in (('one_copy_outward_residual_hessian',1,1),
                             ('two_copy_outward_residual_hessian',2,2),
                             ('coordinate_physical_forcing',2,1)):
  fam=rep['bounds'][family]
  for row,rs in rows.items():
   for kind,val in rs.items():compare(family+'/'+kind+'/'+row,factor*val,fam['metric'][kind][row])
  for kind,val in {'value':PP,'proper_time_derivative':PPt*invr,
                   'each_angular_derivative':2*PP*invr,'each_tangential_derivative':max(PPt,2*PP)*invr}.items():
   compare(family+'/scalar/'+kind,pscale*val,fam['scalar'][kind])
 # Use the retained producer field caps, rather than smaller independent caps,
 # to make the new gauge-source adapters portable to that authenticated report.
 pp=rep['bounds']['parameters'];cc=rep['bounds']['actual_field_component_caps'];aa=rep['bounds']['actual_field_amplitude_caps']
 ar=F(pp['a_abs']);sp=F(pp['sigma_prime_abs']);ss=F(pp['sigma_second_abs'])
 h=readcap(cc['H']);vv=readcap(cc['eta']);tau=readcap(cc['tau']);sq=readcap(cc['abs_trace_A_squared'])
 ht=readcap(aa['J_t']);vt=readcap(aa['v_t']);taut=readcap(aa['tau_t']);sqt=readcap(aa['abs_trace_A_squared_t'])
 GN=4*ss*vv*vv/3+4*sp*h*vv+10*ar*h*h+8*ar*h*tau+4*ar*sq
 QC=4*ar*h*h+3*ar*h*tau+2*ar*sq+4*sp*h*vv/3+sp*vv*tau/3
 GNt=(8*ss*vv*vt/3+4*sp*(ht*vv+h*vt)+20*ar*h*ht+8*ar*(ht*tau+h*taut)+4*ar*sqt)*invr
 QCt=(8*ar*h*ht+3*ar*(ht*tau+h*taut)+2*ar*sqt+4*sp*(ht*vv+h*vt)/3+sp*(vt*tau+vv*taut)/3)*invr
 require(GN<F(12000000000) and QC<1,'portable normal gauge value caps')
 out={'status':'PASS_INDEPENDENT_ACTUAL_SHELL_MAJORANTS_AND_NORMAL_GAUGE_CAPS',
      'shell_report_sha256':hashlib.sha256(raw).hexdigest(),'retained_input_sha256':EXPECTED,
      'independent_scalar_comparisons':len(checks),'checks':checks,
      'arithmetic':'Exact fractions, no intermediate rounding, no producer code import',
      'new_actual_shell_gauge_caps':{'G_N_value':packed(GN),'Q_C_lower_y_value':packed(QC),
                'G_N_proper_time_derivative':packed(GNt),'Q_C_lower_y_proper_time_derivative':packed(QCt)},
      'portable_caps':{'G_N_value':'12000000000','Q_C_lower_y_value':'1'},
      'derivative_scope':'Proper shell clock derivative of lower normal component in radially parallel background normal frame; no source normal derivative or bulk source C1 bound.',
      'limitations':'Bounds inherit retained continuous master/gauge certificates, which are authenticated rather than all freshly replayed. No U2 solve, inverse norm, nonlinear family or force claim.'}
 args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
 try:main()
 except (ValueError,KeyError) as exc:raise SystemExit(str(exc))
