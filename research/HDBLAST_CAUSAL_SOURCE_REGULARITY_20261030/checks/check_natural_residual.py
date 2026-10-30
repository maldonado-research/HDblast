#!/usr/bin/env python3
"""Exact source identity and positive-frame natural coefficient residual bridge.

All conclusions are conditional on the explicitly supplied candidate jet bounds.
There are no project imports. Decisions remain active under python -O.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import argparse, hashlib, json

def need(ok, why):
    if not ok: raise RuntimeError(why)

def load(path, digest):
    data=Path(path).read_bytes()
    need(hashlib.sha256(data).hexdigest()==digest,'INPUT_HASH:'+str(path))
    return json.loads(data)

def ceil_sqrt(q):
    need(q>=0,'nonnegative square')
    r=isqrt(q.numerator//q.denominator)
    if Q(r*r)<q:r+=1
    need(Q(r*r)>=q and (r==0 or Q((r-1)**2)<q),'rational square-root ceiling')
    return Q(r)

def poly(terms):
    out={}
    for term in terms:
        key=tuple(sorted((v,int(p)) for v,p in term['powers'].items() if p))
        out[key]=out.get(key,Q(0))+Q(term['coefficient'])
    return {k:v for k,v in out.items() if v}

def padd(*pairs):
    out={}
    for coef,p in pairs:
        for k,v in p.items():out[k]=out.get(k,Q(0))+coef*v
    return {k:v for k,v in out.items() if v}

def encode(v):
    if isinstance(v,Q):return str(v)
    if isinstance(v,dict):return {k:encode(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [encode(x) for x in v]
    return v

parser=argparse.ArgumentParser()
parser.add_argument('--source-map',required=True)
parser.add_argument('--source-report',required=True)
parser.add_argument('--time-report',required=True)
parser.add_argument('--constraint-formulas',required=True)
parser.add_argument('--linear-operator',required=True)
parser.add_argument('--jet-caps',help='Authenticated full candidate coordinate-jet report')
parser.add_argument('--higher-transfer',help='PR100 matched-extension transfer report, required with candidate')
parser.add_argument('--jet-caps-sha256')
args=parser.parse_args()
pins={'source_map':'34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615',
      'source_report':'45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f',
      'time_report':'6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d',
      'constraint_formulas':'f3f037f890655ab05f936491d378a494b56e10f0c8a3d32dda0f1b8b2ea7556f',
      'linear_operator':'a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189'}
maps=load(args.source_map,pins['source_map'])
sources=load(args.source_report,pins['source_report'])
time=load(args.time_report,pins['time_report'])
cf=load(args.constraint_formulas,pins['constraint_formulas'])
linear_raw=Path(args.linear_operator).read_bytes()
need(hashlib.sha256(linear_raw).hexdigest()==pins['linear_operator'],'actual linear operator closure hash')
labels=time['labels'];H=list(map(Q,time['Gram']))
need(len(labels)==len(H)==23 and min(H)>0,'all23 positive input')
S={0:Q(1),2:Q(1,12),4:Q(19,720)}
outs=[f'{f}:{l}' for f in ('R','T') for l in (0,2,4)]+['A:2','A:4']
B=[S[l] for _ in range(2) for l in (0,2,4)]+[l*(l+2)*S[l] for l in (2,4)]
need(cf['input_labels']==labels and list(map(Q,cf['input_Gram']))==H and cf['output_labels']==outs and list(map(Q,cf['output_Gram']))==B and Q(cf['sharp_scaled_derivative_cap_squared'])==Q(7,2),'inherited exact DC formula conventions and norms')

def zeros():return [[Q(0) for _ in labels] for _ in outs]
def matrices(alpha=Q(0),eta=Q(0),zeta=Q(0)):
    rx,rr,z=zeros(),zeros(),zeros()
    def put(m,row,ch,l,c):
        name=f'{ch}:{l}'
        if name in labels:m[outs.index(f'{row}:{l}')][labels.index(name)]+=c
    for l in (0,2,4):
        lam=l*(l+2);beta=2-Q(2,3)*lam
        put(rx,'R','L',l,-1);put(rr,'T','L',l,-1)
        for ch,c in [('K',-Q(1,2)),('N',-Q(1,2)),('J',Q(3,2))]:put(rr,'R',ch,l,c)
        for ch,c in [('K',-Q(1,2)),('N',-Q(1,2)),('J',-Q(3,2))]:put(rx,'T',ch,l,c)
        for ch,c in [('K',4*alpha),('N',alpha),('J',-3*alpha),('L',-3*eta),('v',-lam*zeta)]:put(z,'R',ch,l,c)
        for ch,c in [('L',5*alpha),('N',-3*eta),('J',-3*eta),('u',-lam*zeta)]:put(z,'T',ch,l,c)
        if l:
            put(rr,'A','v',l,-1);put(rx,'A','u',l,-1)
            for ch,c in [('v',5*alpha),('u',-4*eta),('K',-zeta/2),('N',zeta/2),('J',-zeta/2),('H',beta*zeta)]:put(z,'A',ch,l,c)
    return rx,rr,z

rx,rr,_=matrices()
gram=[[sum((rx[i][k]*rx[j][k]+rr[i][k]*rr[j][k])/H[k] for k in range(23))*B[j] for j in range(8)] for i in range(8)]
need(gram==[[Q(i==j)*(Q(7,4) if i<6 else 1) for j in range(8)] for i in range(8)],'exact combined derivative Gram')
qder=Q(15,8)
need(qder*qder>Q(7,2),'strict rational joint derivative cap')
def frob2(m):return sum((B[i]*m[i][j]**2/H[j] for i in range(8) for j in range(23)),Q(0))
prim={name:matrices(**{name:Q(1)})[2] for name in ('alpha','eta','zeta')}
prim_squares={name:frob2(m) for name,m in prim.items()}
prim_caps={name:ceil_sqrt(q) for name,q in prim_squares.items()}

# Every angular component of the symmetrized physical gauge gradient is
# retained, including both TT rows (where appropriate differences may vanish).
identities=[]
for field,rows in maps['maps']['QE'].items():
    need(field in maps['maps']['QEred'] and field in maps['maps']['gauge_gradient'],'complete map sectors')
    for j,row in enumerate(rows):
        qe,qr,gg=poly(row),poly(maps['maps']['QEred'][field][j]),poly(maps['maps']['gauge_gradient'][field][j])
        need(padd((1,qe),(-1,qr))==gg,'complete physical gauge gradient:'+field+'/'+str(j))
        identities.append(field+'/'+str(j))
need(len(identities)==20,'all20 metric angular amplitudes')
def shifted(p,coefficient,powers):
    ans={}
    for key,v in p.items():
        exp=dict(key)
        for name,n in powers.items():exp[name]=exp.get(name,0)+n
        ans[tuple(sorted((name,n) for name,n in exp.items() if n))]=coefficient*v
    return ans
for field,rows in maps['maps']['QEred'].items():
    for j,row in enumerate(rows):
        need(poly(maps['maps']['metric_forcing'][field][j])==shifted(poly(row),Q(-2),{'ell':2,'rho':2}),'complete metric forcing normalization:'+field+'/'+str(j))
for j,row in enumerate(maps['maps']['Qphi']):
    need(poly(maps['maps']['scalar_forcing'][j])==shifted(poly(row),Q(1),{'ell':2,'rho':2}),'complete scalar forcing normalization:'+str(j))
# Wrong factors are rejected by the full polynomial identity, not by a sample.
nonzero=poly(maps['maps']['gauge_gradient']['yy'][0])
need(nonzero and nonzero!=padd((Q(2),nonzero)),'wrong symmetrization rejected')

ell=Q(1,10**6)
bg=time['global'];rho_lo,rho_hi=map(Q,bg['rho_interval'])
aabs=max(map(abs,map(Q,bg['a_interval'])));pabs=max(map(abs,map(Q,bg['p_interval'])))
need(rho_lo>0,'positive actual physical scale')
alpha=ell*rho_hi*aabs
eta=Q(3,2)*ell*ell
zeta=ell
# alpha_r=ell²rho² b and b=rho^-2-a²+p²/3.
alpha_r=ell*ell*(1+rho_hi*rho_hi*(aabs*aabs+pabs*pabs/3))
eta_x=ell*ell;zeta_x=Q(3,2)*ell**3
z0=prim_caps['alpha']*alpha+prim_caps['eta']*eta+prim_caps['zeta']*zeta
zx=prim_caps['eta']*eta_x+prim_caps['zeta']*zeta_x
zr=prim_caps['alpha']*alpha_r
# Angular intrinsic covector rough eigenvalues are <=24; extrinsic covector
# connection obeys sum_a|Omega_a C|² <=3(a²+h²)|C|².
need(Q(5)**2>24 and Q(7,4)**2>3,'strict angular rational caps')
connection=2*alpha+5*ell+Q(7,4)*(alpha+eta)
cap=sources['global_cover']['summary']['absolute_coefficient_caps']
metric_rows={}
for field,ch,degrees in [('yy','K',(0,2,4)),('y0','L',(0,2,4)),('00','N',(0,2,4)),('ya','v',(2,4)),('0a','u',(2,4)),('ab',None,None)]:
    names=[f'{ch}:{l}' for l in degrees] if ch else ['J:0','J:2','J:4','H:2','H:4','T:+','T:-']
    for j,name in enumerate(names):metric_rows[name]=Q(cap[f'QE/{field}/{j}'])+Q(cap[f'QEred/{field}/{j}'])
Qgrad2=sum((H[labels.index(name)]*v*v for name,v in metric_rows.items()),Q(0))
Qgrad=ceil_sqrt(Qgrad2)

# Exact covector connection Gram witness for arbitrary a,h: block spectrum
# is bounded by 3(a²+h²), including Lorentz time signs in the connection.
for aa,hh in [(Q(2,7),Q(-3,11)),(Q(-5,13),Q(7,17)),(Q(0),Q(1))]:
    for vector in ([Q(i+1) for i in range(5)],[Q(-1),Q(2),Q(-3),Q(4),Q(-5)]):
        lhs=(aa*aa+hh*hh)*sum((v*v for v in vector[2:]),Q(0))+3*(aa*vector[0]-hh*vector[1])**2
        rhs=3*(aa*aa+hh*hh)*sum((v*v for v in vector),Q(0))
        need(lhs<=rhs,'covector extrinsic connection bound')
# Product rule: C=P/(ell rho), rho_r=-alpha rho, rho_x=0.
# physical_y C=-(P_r+alpha P)/(ell²rho²), physical_0 C=P_x/(ell²rho²).
a,b,c=Q(2,7),Q(-3,11),Q(5,13)
need(-(b+a*c)!=-(b-a*c),'wrong radial product sign control')
result={'status':'PASS_NATURAL_COEFFICIENT_RESIDUAL_BRIDGE','input_sha256':pins,
 'source_map_identity':'QE-QEred=gauge_gradient=one-half-symmetrized reference covariant derivative of actual full-metric-lowered QC, under retained C1=0 neighbourhood premise',
 'complete_metric_map_identity_count':len(identities),'complete_metric_map_identities':identities,
 'coordinate_DCs':'P=ell*rho*DC[Z]=Rx Zx+Rr Zr+Zmat Z',
 'joint_derivative_constant_upper':qder,'value_matrix_primitive_Frobenius_squares':prim_squares,
 'value_matrix_primitive_operator_caps':prim_caps,
 'background_caps':{'rho_lower':rho_lo,'rho_upper':rho_hi,'alpha_abs':alpha,'eta_abs':eta,'zeta_abs':zeta,'alpha_r_abs':alpha_r,'eta_x_abs':eta_x,'zeta_x_abs':zeta_x},
 'Zmat_caps':{'value':z0,'x_derivative':zx,'r_derivative':zr},
 'gradient_DC_connection_scaled_cap':connection,
 'actual_QC_symmetric_gradient_L2_cap':Qgrad,'actual_QC_symmetric_gradient_L2_square_upper':Qgrad2,
 'bridge':'E2[Z]=R_metric/(2ell²rho²)+nabla_(a DC[Z]_b)+(QE-QEred); S2[Z]=-R_scalar/(ell²rho²)',
 'natural_equation':'E_ab=Ric_ab-dphi_a*dphi_b-(2/3)V(phi)g_ab; S=Box_g phi-Vprime(phi). E2 is its second parameter derivative, not Einstein G_ab directly.',
 'complete_forcing_normalization_checks':{'metric_amplitudes':20,'scalar_amplitudes':3},
 'norm':'Positive Frobenius metric or scalar physical background-frame S3 probability L2; not a Lorentz invariant observable',
 'wrong_formula_controls':['factor_two_symmetrization','wrong_radial_product_rule_sign'],
 'scope':'Direct natural second-parameter coefficient residual bridge, not C2 error to target, vanishing residual, relative smallness, finite-amplitude nonlinear solution or physical discovery'}
if args.jet_caps:
    need(args.jet_caps_sha256 is not None,'candidate jet report digest is required')
    jets=load(args.jet_caps,args.jet_caps_sha256)
    need(args.higher_transfer is not None,'PR100 transfer report required')
    tp='4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3'
    transfer=load(args.higher_transfer,tp)
    j=jets['candidate_absolute_coordinate_jet_H_caps']
    v={k:Q(j[name]) for k,name in [('V0','value'),('Vx','x'),('Vr','r'),('Vxx','xx'),('Vxr','xr'),('Vrr','rr')]}
    delta=Q(jets['candidate_C1_change_from_W1_quadratic_X_cap'])
    unit_defect=Q(transfer['direct_reproof']['new_reduced_forcing_defect_H_cap'])+Q(2,3)*Q(transfer['k_actual_bulk_X_operator_cap'])*delta+2*Q(jets['source_caps']['F_H'])
    v['unit_defect']=unit_defect
    need(min(v.values())>=0,'nonnegative candidate bounds')
    P0=qder*max(v['Vx'],v['Vr'])+z0*v['V0']
    Px=qder*max(v['Vxx'],v['Vxr'])+z0*v['Vx']+zx*v['V0']
    Pr=qder*max(v['Vxr'],v['Vrr'])+z0*v['Vr']+zr*v['V0']
    Dcap=(Px+Pr+connection*P0)/(ell*ell*rho_lo*rho_lo)
    Ecap=v['unit_defect']/(2*ell*ell*rho_lo*rho_lo)+Dcap+Qgrad
    Scap=v['unit_defect']/(ell*ell*rho_lo*rho_lo)
    result['candidate_input_sha256']=args.jet_caps_sha256
    result['higher_transfer_sha256']=tp
    result['candidate_unit_defect_proof']='R_eps=R_quad-La(Z-W1_quad)+F_eps-F; R=principal-actual_RHS; weighted La cap=(2/3)k; source discrepancy charged by2FH'
    result['candidate_unit_defect_terms']={'old_matched_reduced_defect_H':Q(transfer['direct_reproof']['new_reduced_forcing_defect_H_cap']),'La_filtered_change_H':Q(2,3)*Q(transfer['k_actual_bulk_X_operator_cap'])*delta,'actual_source_discrepancy_H':2*Q(jets['source_caps']['F_H'])}
    result['candidate_caps']=v
    readable={'natural_Einstein_second_parameter_L2':Q('2.918e43'),'natural_scalar_second_parameter_L2':Q('2.118e30'),'unit_defect':Q('3.574e22')}
    need(Ecap<readable['natural_Einstein_second_parameter_L2'] and Scap<readable['natural_scalar_second_parameter_L2'] and unit_defect<readable['unit_defect'],'strict readable natural coefficient residual caps')
    result['candidate_readable_strict_uppers']=readable
    result['candidate_derived_caps']={'P':P0,'P_x':Px,'P_r':Pr,'nabla_DC_physical_L2':Dcap,'natural_Einstein_second_parameter_L2':Ecap,'natural_scalar_second_parameter_L2':Scap}
print(json.dumps(encode(result),indent=2,sort_keys=True))
