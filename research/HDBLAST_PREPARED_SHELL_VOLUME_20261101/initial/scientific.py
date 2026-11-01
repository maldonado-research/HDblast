"""Portable science of the unchanged frozen 1024-cell initial campaign."""
from fractions import Fraction as F
import math

def need(ok,msg):
 if not ok:raise ValueError(msg)

def first(m,data,r,tail=False):
 rho,e,p,a,h,C,B,k,U2,E=m.background(r,data,'none')
 q=m.bump(r,data,tail);dy=lambda j:-(1/rho)*j.dz()
 u=C*q;uy=dy(u);uyy=dy(uy);v=B*(q.dz()+k*q);vy=dy(v)
 ut=E*u;vt=E*(v-3*u/(rho*p));s=ut/rho;j=vt/rho
 ap=dy(a)
 QH=(-2*uy**2-24*a*u*uy-2*u*uyy-(16*ap+32*a*a+F(8,3)*p*p)*u*u
     -4*u*u/rho**2-s*s-(j*j+vy*vy+U2*v*v)/6-F(4,3)*p*u*vy-F(4,3)*v*v/rho**2)
 QM=-2*s*uy+j*(vy+4*p*u)/3
 V=p*p/2-6*a*a+6/rho**2
 return dict(rho=rho,p=p,a=a,q=q,u=u,ut=ut,uy=uy,v=v,vy=vy,s=s,j=j,QH=QH,QM=QM,V=V,dy=dy)

def pack(m,n):
 F=m.F
 return [(F(1),F(65,64),True)]+[(F(k,n),F(k+1,n),False) for k in range(n+n//64,2*n-n//64)]+[(F(127,64),F(2),True)]

def enclosure(x):
 return {'exact':x.row(),'display_only':[float(x.lo),float(x.hi)]}

def profiles(m,data,r,tail,Iall,Jall,KH,Kb):
 d=first(m,data,r,tail);rho,p,a,q,dy=[d[k] for k in ('rho','p','a','q','dy')]
 c=2*Iall/(9*Jall);f=q*c;fy=dy(f);ratio=Iall/Jall
 # Integrated-by-parts primitive: exact, with lower endpoint f=f_y=0.
 # K=KH-(I/J)Kb, K_z=-rho^5 QH+(I/J)rho^3 q.
 K0=KH-ratio*Kb;rhs=-rho**5*d['QH']+rho**3*q*ratio
 K=m.J([K0]+[rhs.c[k]/(k+1) for k in range(min(rhs.n+1,3))])
 xi=F(3,2)*(fy-a*f)/p**2-K/(rho**4*p**2)
 xiy=dy(xi);xiyy=dy(xiy);A=f+2*xiy;B=f+2*a*xi
 tau=-d['QM']/(6*a)
 out={k:d[k].val() for k in ('u','ut','uy','v','vy','QH','QM','rho','p','a')}
 out.update(f=f.val(),fy=fy.val(),fyy=dy(fy).val(),fyyy=dy(dy(fy)).val(),xi=xi.val(),xiy=xiy.val(),xiyy=xiyy.val(),A=A.val(),B=B.val(),Ay=dy(A).val(),By=dy(B).val(),tau_angular=tau.val(),KH=KH,Kb=Kb)
 means=m.mean_adapter.initial_means(rho=out['rho'],a=out['a'],u=out['u'],ut=out['ut'],uy=out['uy'],qm=out['QM'],A=out['A'],Ay=out['Ay'],B=out['B'],By=out['By'])
 for vector in ('U0','V0','U0_r'):
  out.update({vector+'_'+name:value for name,value in zip(means['order'],means[vector])})
 return out

def generate(m,data):
 n=1024;targets="3/2,5/4";I=m.I;ELL=m.ELL
 cells=pack(m,n);rows=[];intI=I(0);intJ=I(0)
 for idx,(lo,hi,tail) in enumerate(cells):
  r=I(lo,hi);d=first(m,data,r,tail);w=ELL*(hi-lo)
  i=w*(d['rho']**5*d['QH']).val();j=w*(d['rho']**3*d['q']).val()
  rows.append({'r':[str(lo),str(hi)],'tail':tail,'I':i,'J':j})
  intI=intI+i;intJ=intJ+j
 need(intJ.lo>0,'strict positive J integral')
 need(intI.hi<0,'strict negative I integral verified numerically')
 # r increases toward decreasing y. Suffix sums from r to 2 represent
 # the primitive based at yL=y(2ell). Both cell widths are retained.
 KH=I(0);Kb=I(0)
 for row in reversed(rows):
  lo,hi=map(F,row['r']);insideI=I(min(0,row['I'].lo),max(0,row['I'].hi));insideJ=I(0,max(0,row['J'].hi))
  row['KH_cell']=KH+insideI;row['Kb_cell']=Kb+insideJ
  KH=KH+row['I'];Kb=Kb+row['J']
 points=[]
 for target in map(F,targets.split(',')):
  need(1<target<2,'interior point target')
  # The target must be a partition boundary so its prefix is an exact sum.
  need(any(F(row['r'][0])==target for row in rows),'target must be partition boundary')
  targetI=sum((row['I'] for row in rows if F(row['r'][0])>=target),I(0))
  targetJ=sum((row['J'] for row in rows if F(row['r'][0])>=target),I(0))
  pt=profiles(m,data,I(target),False,intI,intJ,targetI,targetJ)
  points.append({'r':str(target),'profiles':{k:enclosure(v) for k,v in pt.items()}})
 profile_rows=[]
 for row in rows:
  lo,hi=map(F,row['r']);pr=profiles(m,data,I(lo,hi),row['tail'],intI,intJ,row['KH_cell'],row['Kb_cell'])
  profile_rows.append({'r':row['r'],'tail':row['tail'],'integral_I':enclosure(row['I']),'integral_J':enclosure(row['J']), 'profiles':{k:enclosure(v) for k,v in pr.items()}})
 caps={}
 for name,upper in [('packet',F(2)),('evolution_dependency_r_le_5_4',F(5,4))]:
  selected=[row for row in profile_rows if F(row['r'][1])<=upper]
  caps[name]={'r':['0' if upper<2 else '1',str(upper)],'exact_zero_before_packet':True,'components':{k:{'exact_upper':str(max(max(abs(F(row['profiles'][k]['exact'][0])),abs(F(row['profiles'][k]['exact'][1]))) for row in selected)),'display_only':max(max(abs(x) for x in row['profiles'][k]['display_only']) for row in selected)} for k in selected[0]['profiles']}}
 result={'status':'PASS_VALIDATED_EXPLICIT_MEAN_INITIAL_PROFILES','arithmetic':{'outward_grid_denominator':str(m.GRID),'exact_rational':True,'exp':'validated degree28 alternating Taylor after argument halving, repeated outward squaring','bump_tails':'same authenticated monotone s^(2k) exp(4-s), k<=4, s>=4096/63'},'family':'explicit bump Q, modified IFT trial; second parameter derivative convention','domain':{'ell':str(ELL),'r':['1','2'],'dy':'-rho dz','zero_extension':'exact smooth compact support'},'partition':{'N':n,'cells':len(rows),'tail_width':'1/64'},'integrals':{'I':enclosure(intI),'J':enclosure(intJ),'coefficient_2I_over_9J':enclosure(2*intI/(9*intJ)),'metric_moment_2I_over_3':enclosure(2*intI/3)},'targets':points,'global_caps':caps,'profile_cells':profile_rows,'limits':{'background_shell_uncertainty_reduced':False,'evolved_response':False,'shell_energy_or_physical_generator':False,'parameter_radius_effective':False}}
 return result

def upper_sqrt(q):
 GRID=1<<112
 v=q.numerator*GRID*GRID;d=q.denominator;k=math.isqrt(v//d)
 if k*k*d<v:k+=1
 return F(k,GRID)

def enc(v):return {'exact':v.row(),'display_only':[float(v.lo),float(v.hi)]}


def conditioned(m,g,data,st,static_sha256):
 ORDER=('K0','L0','N0','J0','w0');GRAM=(1,2,1,3,1);GRID=1<<112
 py=m.I(*st['y_jets']['p'][1]['exact']);pyy=2*m.I(*st['y_jets']['p'][2]['exact'])
 def row(p):
  vals={k:m.I(*p[k]['exact']) for k in ('rho','a','u','ut','uy','f','fy','xi','xiy','p')};vals['qm']=m.I(*p['QM']['exact']);vals['py']=py
  z=g.conditioned_initial_means(**vals)
  return {vector:{name:enc(v) for name,v in zip(ORDER,z[vector])} for vector in ('U0','V0','U0_r')}
 cells=[{'r':p['r'],'tail':p['tail'],'conditioned':row(p['profiles'])} for p in data['profile_cells']]
 points=[{'r':p['r'],'conditioned':row(p['profiles'])} for p in data['targets']]
 caps={}
 for key,end in (('packet',F(2)),('causal_footprint_r_le_7_4',F(7,4)),('history_r_le_5_4',F(5,4))):
  selected=[p for p in cells if F(p['r'][1])<=end];cc={};norms={}
  for vector in ('U0','V0','U0_r'):
   cc[vector]={name:max(max(abs(F(x)) for x in p['conditioned'][vector][name]['exact']) for p in selected) for name in ORDER}
   norms[vector]=max(upper_sqrt(sum((w*max(abs(F(x)) for x in p['conditioned'][vector][name]['exact'])**2 for w,name in zip(GRAM,ORDER)),F(0))) for p in selected)
  c1=max(sum((upper_sqrt(sum((w*max(abs(F(x)) for x in p['conditioned'][vector][name]['exact'])**2 for w,name in zip(GRAM,ORDER)),F(0))) for vector in ('U0','V0','U0_r')),F(0)) for p in selected)
  caps[key]={'r':['0',str(end)],'component_caps':{vector:{name:{'exact_upper':str(v),'display_only':float(v)} for name,v in values.items()} for vector,values in cc.items()},'Gram_norm_caps':{vector:{'exact_upper':str(v),'display_only':float(v)} for vector,v in norms.items()},'sum_of_Gram_norms_C1_cap':{'exact_upper':str(c1),'display_only':float(c1)}}
 result={'status':'PASS_CONDITIONED_SAME_FROZEN_INITIAL_BOXES','meaning':'Pure second harmonic gauge adapter; same physical family and same retained first tangent. No background/grid/oracle rerun.','original_campaign':{'N':1024,'targets':['3/2','5/4']},'arithmetic':{'exact_rational':True,'outward_grid_denominator':str(GRID),'Gram':list(GRAM),'C1':'sup pointwise [|U0|_Gram+|V0|_Gram+|d_r U0|_Gram]; each vector cap also retained'},'additional_retained_coefficient':{'sha256':static_sha256,'p_y':enc(py),'p_yy':enc(pyy),'scope':'Extra old actual-background coefficient input y_jets p[1], 2 p[2], globally valid on original annulus; not claimed to have appeared as raw profile fields. p_yy retained for possible downstream use, not used here.'},'zero_extension':{'r_less_equal_1':True,'r_greater_equal_2':True,'reason':'Exact analytic compactness and moment identity for the chosen I/J; no subtraction of independent interval I and ratio*J is used to define an exterior value.'},'order':list(ORDER),'targets':points,'global_caps':caps,'profile_cells':cells,'limits':{'same_actual_background_uncertainty_retained':True,'nonlinear_gauge_theorem_retained_conditional':True,'evolution_or_shell_observable_evaluated':False}}
 return result
