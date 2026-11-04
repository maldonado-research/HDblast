#!/usr/bin/env python3
"""Exact polynomial/Fraction identities for a conditional incoming-preparation ledger."""
import argparse,hashlib,json,sys
from fractions import Fraction as F
from pathlib import Path,PurePosixPath
HERE=Path(__file__).resolve().parent
PIN='82bd0eb951cae12c64a0ed56f39f9af3efab8e78c85c56f1a835d4c6af23f530'
DESIGN='f91f854558716fbf5dcf66e66e69b92e5ce427762bd3f5f888d9d2721d5221a7'
class Failure(Exception):pass
def need(v,w):
 if not v:raise Failure(w)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file()and not p.is_symlink()and all(not q.is_symlink()for q in p.absolute().parents)
def read(p):return json.loads(p.read_bytes())
def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+'\n').encode()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
class P:
 def __init__(self,x=0):self.d={():F(x)}if x else {}
 @classmethod
 def fromdict(cls,d):
  p=cls();p.d={k:F(v)for k,v in d.items()if v};return p
 @classmethod
 def var(cls,n):return cls.fromdict({((n,1),):F(1)})
 def __add__(self,b):
  b=b if isinstance(b,P)else P(b);d=dict(self.d)
  for m,c in b.d.items():d[m]=d.get(m,F(0))+c
  return P.fromdict(d)
 __radd__=__add__
 def __neg__(self):return P.fromdict({m:-c for m,c in self.d.items()})
 def __sub__(self,b):return self+-b if isinstance(b,P)else self+P(-b)
 def __rsub__(self,b):return P(b)+-self
 def __mul__(self,b):
  b=b if isinstance(b,P)else P(b);d={}
  for m,c in self.d.items():
   for n,e in b.d.items():
    powers=dict(m)
    for x,k in n:powers[x]=powers.get(x,0)+k
    key=tuple(sorted(powers.items()));d[key]=d.get(key,F(0))+c*e
  return P.fromdict(d)
 __rmul__=__mul__
 def __truediv__(self,b):return self*F(1,b)
 def __pow__(self,n):
  need(type(n)is int and n>=0,'POLYNOMIAL_POWER');p=P(1)
  for j in range(n):p=p*self
  return p
 def coef(self,n,k):
  d={}
  for m,c in self.d.items():
   powers=dict(m)
   if powers.get(n,0)==k:
    powers.pop(n,None);d[tuple(sorted(powers.items()))]=c
  return P.fromdict(d)
 def zero(self):return not self.d
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--design',type=Path,default=HERE/'FROZEN_DESIGN.json');a=ap.parse_args();repo=a.repo_root.absolute();out=a.output.absolute();design=a.design.absolute()
 need(all(not q.is_symlink()for q in [repo,*repo.parents,out,*out.parents,design,*design.parents]),'PATH_SYMLINK');repo=repo.resolve(strict=True);out=out.resolve();design=design.resolve(strict=True)
 need(not out.exists(),'FRESH_OUTPUT_FILE');need(out.parent.is_dir(),'OUTPUT_PARENT_EXISTS');need(all(b!=out and b not in out.parents for b in [repo,HERE]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
 need(regular(HERE/'SOURCE_PINS.json')and sha(HERE/'SOURCE_PINS.json')==PIN,'SOURCE_PIN_FILE_BYTES');need(regular(design)and sha(design)==DESIGN,'FROZEN_DESIGN_BYTES');pins=read(HERE/'SOURCE_PINS.json');need(pins['schema']=='HDBLAST_INCOMING_THEORY_SELECTED_SOURCE_PINS_V1','SOURCE_PIN_SCHEMA')
 for n,row in pins['files'].items():
  rel=PurePosixPath(n);need(not rel.is_absolute()and '..'not in rel.parts and rel.as_posix()==n,'SOURCE_RELATIVE_PATH');p=repo/n
  need(regular(p)and p.stat().st_size==row['bytes']and sha(p)==row['sha256'],'AUTHENTICATED_SOURCE_BYTES')
 # All project source bytes are authenticated before the arithmetic starts.
 d=read(design);par=d['parameters'];checks=[]
 def check(name,v):need(v,name);checks.append(name)
 def identity(name,p):check(name,p.zero())
 ell,T,k,alpha,C,Cf,Cg,expo,q,Nguard,kappaguard=[F(par[n])for n in ['ell','T','tilt_k','boundary_alpha','forced_growth_guard','forced_bulk_square_factor','forced_boundary_square_factor','exponential_guard','retained_q','inverse_N_guard','inverse_kappa_guard']]
 check('declared_dimensions_scales_and_guards',(par['screen_dimension'],ell,T,k,alpha,C,Cf,Cg,expo,q,Nguard,kappaguard)==(3,F(1,10**6),F(5,4),F(1,2),F(1,2),F(301,100),F(1,2),F(5002001,4000000),F(44),F(1201,3200000),F(20000001,20000000),F(47,125000)))
 premises=read(repo/'research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json');check('retained_conditional_operator_assertion',premises['status']=='EXPLICIT_CONDITIONAL_PREMISES_NOT_REPROVED_BY_PORTABLE_CHECKER'and premises['mean_gram']==[1,2,1,3,1]and premises['coefficient_norm_upper']=={'A':'1/10000','B':'1/10000','D':'1/1000000','M':'1/2000'})
 check('future_only_domain_binding',premises['domain']=='0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x'and T<=F(5,4)and T<=F(3,2)and T<=F(7,4))
 A,Ht,B2=[P.var(n)for n in ['A','T0','B2']]
 identity('actual_background_affine_Raychaudhuri',3*(2*A*(Ht-A)+A*A-Ht*Ht-B2/3)+3*(Ht-A)**2+B2)
 ep,th0,th1,th2,p0,p1,p2,sg1,sg2,gm1,gm2,L1,L2=[P.var(n)for n in ['ep','theta0','theta1','theta2','phi0prime','phi1prime','phi2prime','sigma1','sigma2','metric1','metric2','L1','L2']]
 theta=th0+ep*th1+ep**2*th2/2;phipr=p0+ep*p1+ep**2*p2/2;sigma=ep*sg1+ep**2*sg2/2;metric=1+ep*gm1+ep**2*gm2/2
 identity('Raychaudhuri_theta_second_variation',2*(theta*theta/3).coef('ep',2)-(2*th1*th1/3+2*th0*th2/3))
 identity('Raychaudhuri_scalar_second_variation',2*(phipr*phipr).coef('ep',2)-(2*p1*p1+2*p0*p2))
 identity('shear_background_zero_removes_metric_crosses',2*(metric*sigma*sigma).coef('ep',2)-2*sg1*sg1)
 density=1+ep*L1+ep**2*L2/2+ep**2*L1*L1/2
 identity('actual_area_second_variation',2*density.coef('ep',2)-(L2+L1*L1))
 check('area_sign_does_not_follow_from_log_sign_alone',F(-1)<0 and F(-1)+F(2)**2>0)
 R,Rp,th2p=[P.var(n)for n in ['R','Rprime','theta2prime']]
 identity('integrating_factor_R_squared_product_rule',R*(R*R*th2p+2*R*Rp*th2)-(R**3*th2p+2*R*R*Rp*th2))
 V,Q,U,Fc,Vr,Qr,Pc,Mu,G=[P.var(n)for n in ['V','Q','U','F','V_r','Q_r','P','Mu','G']]
 energy=(V*V+Q*Q+2*k*V*Q+U*U)/2;flux=V*Q+k*(V*V+Q*Q)/2
 energy_x=V*(Qr+Fc)+Q*Vr+k*((Qr+Fc)*Q+V*Vr)+U*V;flux_r=Vr*Q+V*Qr+k*(V*Vr+Q*Qr)
 identity('tilted_energy_local_identity',energy_x-flux_r-(V+k*Q)*Fc-U*V)
 identity('moving_outer_characteristic_flux',flux-energy+(1-k)*(V-Q)**2/2+U*U/2)
 identity('sharp_bulk_forcing_square',2*energy-(V+k*Q)**2-((1-k*k)*Q*Q+U*U))
 identity('sharp_bulk_Young_remainder',energy+Fc*Fc/2-(V+k*Q)*Fc-((V+k*Q-Fc)**2/2+(1-k*k)*Q*Q/2+U*U/2))
 v=Pc-Mu-G;radial=Mu+G
 actual=-v*radial-k*(v*v+radial*radial)/2+alpha*U*v
 hom=-k*Pc*Pc/2-(1-k)*Pc*Mu+(1-k)*Mu*Mu+alpha*U*(Pc-Mu)
 source=-(1-k)*Pc*G+2*(1-k)*Mu*G+(1-k)*G*G-alpha*U*G
 identity('actual_nonsymmetric_Robin_source_expansion',actual-hom-source)
 identity('incoming_characteristic_boundary_law',(v-radial)-(Pc-2*(Mu+G)))
 identity('boundary_state_evolution',v-(Pc-Mu-G))
 identity('forced_boundary_P_Young_square',Pc*Pc/8-Pc*G/2+G*G/2-(Pc-2*G)**2/8)
 m=F(premises['coefficient_norm_upper']['M']);beta=alpha+(1-k)*m;Cb=2*(beta*beta/k+(1-k)*m*m+alpha*m)/alpha
 amax=F(1,10000);bmax=amax;dmax=F(1,10**6);Cvbulk=amax+(bmax+k*amax)/2+dmax/2+F(1,2);Cqbulk=k*bmax+(bmax+k*amax)/2+k*dmax/2;Cubulk=dmax/2+k*dmax/2+F(1,2);Cbulk=max(2*Cvbulk/(1-k),2*Cqbulk/(1-k),2*Cubulk)
 check('homogeneous_rate_reconstruction',Cb==F(2003001,1000000)and Cbulk==F(1000351,500000))
 check('forced_growth_rate_guard',max(Cb+1,Cbulk+1)<C);check('actual_boundary_source_square_factor',Cg==1+(m+F(1,2))**2);check('sharp_bulk_source_square_factor',Cf==F(1,2))
 z=C*T;n=32;term=F(1);partial=term
 for j in range(1,n+1):term=term*z/j;partial+=term
 ratio=z/(n+2);upper=partial+(term*z/(n+1))/(1-ratio)
 check('exponential_tail_ratio',0<ratio<1);check('forced_exponential_44_guard',upper<expo)
 r=F(3,2);den=(r-1)*(2-r);denprime=3-2*r
 check('actual_packet_central_derivative_zero',den==F(1,4)and denprime==0 and  4-1/den==0)
 check('actual_packet_both_principal_directions_nonzero',ell>0 and  ell*r>0)
 prior=read(repo/'research/HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/bounds/references/NORMAL.json');con=prior['constants'];rhom=F(con['rho']['exact'][1]);pm=F(con['minus_p']['exact'][1]);am=F(con['a']['exact'][1]);bM=ell*rhom*am;cM=2*ell*rhom*pm/3;delta=T-1;mu=F(20001,20000)
 check('same_initial_causal_constants',rhom==131 and  pm==F(1,7000)and am==F(14,125)and F(con['q']['exact'])==q and F(con['ell']['exact'])==ell)
 check('scalar_parameterized_IVP_matrix_guard',max(1+3*bM,3*ell*ell+bM)<mu)
 zz=mu*delta;tt=F(1);ps=tt
 for j in range(1,7):tt=tt*zz/j;ps+=tt
 eu=ps+(tt*zz/7)/(1-zz/8);check('scalar_parameterized_IVP_exponential',eu<F(4,3))
 Hw=delta*F(4,3)*cM;Hrw=(3*ell*ell+bM)*Hw+cM;Vw=3*ell*ell*Hw+cM;N=max((1+2*Hw+2*T*Vw)/3,1+2*Hrw+2*Vw);kappa=q*Nguard/(1-q)
 check('history_initial_norm_guard',N<Nguard);check('history_resolvent_strict_contraction',0<q<1 and  0<kappa<kappaguard<1 and  delta<=1)
 check('source_bytes_unchanged_after_arithmetic',all(regular(repo/n)and(repo/n).stat().st_size==row['bytes']and sha(repo/n)==row['sha256']for n,row in pins['files'].items()))
 report=enc({'status':'PASS_CONDITIONAL_INCOMING_PREPARATION_IDENTITIES','design_sha256':DESIGN,'source_pins_sha256':PIN,'authenticated_sources':pins['files'],'operator_caps_remain_conditional_future_assertions':True,'background_null_identity':{'screen_dimension':3,'outgoing_generator':'c rho^-2(partial_t+partial_z)','areal_radius':'rho cosh(t)','Ric_kk':'(Phi0prime)^2','static_identity':'a_y=-rho^-2-p^2/3','past_coefficient_enclosures_computed':False},'candidate_null_state_assumptions':d['new_candidate_state_assumptions'],'null_constraint':{'general_second_scalar_cross_term':'-2 Phi0prime Phi2prime retained before setting pulled-back Phi2=0','integrating_factor':'R0^2','optical_second_expansion_nonpositive_under_declared_conditions':True,'log_density_second_nonpositive_under_declared_conditions':True,'actual_area_second_formula':'average(L2+L1^2); no sign from L2 alone','current_shell_expansion_or_H1_bound_established':False,'full_characteristic_solution_constructed':False},'forced_energy':{'tilt_k':k,'boundary_alpha':alpha,'homogeneous_bulk_rate':Cbulk,'homogeneous_boundary_rate':Cb,'forced_growth_guard':C,'bulk_square_factor':Cf,'boundary_square_factor':Cg,'exponential_argument':z,'Taylor_order':n,'Taylor_tail_ratio':ratio,'exponential_rational_upper':upper,'exponential_guard':expo,'scalar_trace_square_factor':2/alpha,'past_source_history_computed':False},'packet_diagnostic':{'central_scaled_radius':r,'Q_center':F(1),'Q_r_center':F(0),'both_scaled_characteristic_derivatives':'ell exp(-3ell/2)>0','not_current_time_symmetric_or_pure_oneway':True,'past_incoming_only_state_excluded':False,'is_physical_energy_flux_measurement':False},'inverse_history':{'norm':'max(sup|value|,sup|first derivative|) on [1,5/4], normalized x for history/r for initial data','cM':cM,'Hw':Hw,'Hrw':Hrw,'Vw':Vw,'N_exact_majorant':N,'N_guard':Nguard,'q_retained':q,'kappa_exact':kappa,'kappa_guard':kappaguard,'observable':'difference of second normalized physical shell scalar means at fixed proper clock','on_admissible_h_only_image':True,'onto_arbitrary_histories':False,'past_state_constructed':False},'verification':{'finite_checks':len(checks),'check_names':checks,'new_PDE_solves':0,'new_background_or_evolved_field_queries':0,'new_bump_queries':0},'limits':d['scope']})
 with out.open('xb')as f:f.write(canonical(report))
 print(json.dumps({'status':'PASS','witness':'INCOMING_PREPARATION_THEORY_IDENTITIES','finite_checks':len(checks),'output_sha256':sha(out)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except (Failure,ArithmeticError,ValueError,OSError,KeyError,TypeError)as e:
  print(json.dumps({'status':'FAIL','witness':str(e)},sort_keys=True));sys.exit(2)
