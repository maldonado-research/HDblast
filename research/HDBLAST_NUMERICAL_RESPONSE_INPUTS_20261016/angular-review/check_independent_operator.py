#!/usr/bin/env python3
"""Independent closed-form primitive and formal Laurent operator audit.

Candidate modules are imported only to obtain the matrices being audited.
Expected tables below come from independently derived frame/curvature and
sphere commutation identities documented in ANGULAR_OPERATOR_REVIEW.md.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from pathlib import Path
import json
import sys

sys.dont_write_bytecode=True
N=23
LABELS=[f'{field}:{l}' for field in ('K','L','N','J','w') for l in (0,2,4)]+[f'{field}:{l}' for field in ('v','u') for l in (2,4)]+['H:2','H:4','T:+','T:-']
INDEX={name:i for i,name in enumerate(LABELS)}
VARS=['ell','rho','a','h','c','Ht','p','py','b','d','V','V1','V2','sigma1','sigma2']
NV=len(VARS);ZEXP=(0,)*NV


def need(ok,msg):
    if not ok:raise ValueError(msg)


class L:
    """Sparse rational Laurent polynomial; division only by monomials."""
    def __init__(self,value=0):self.d={e:Q(v) for e,v in (value.d if isinstance(value,L) else value if isinstance(value,dict) else {ZEXP:value}).items() if v}
    def __add__(self,o):
        d=self.d.copy()
        for e,v in L(o).d.items():d[e]=d.get(e,Q(0))+v
        return L(d)
    __radd__=__add__
    def __neg__(self):return L({e:-v for e,v in self.d.items()})
    def __sub__(self,o):return self+-L(o)
    def __rsub__(self,o):return L(o)+-self
    def __mul__(self,o):
        d={}
        for e,v in self.d.items():
            for f,w in L(o).d.items():
                h=tuple(a+b for a,b in zip(e,f));d[h]=d.get(h,Q(0))+v*w
        return L(d)
    __rmul__=__mul__
    def __truediv__(self,o):
        q=L(o);need(len(q.d)==1,'MONOMIAL_DIVISOR')
        e,v=next(iter(q.d.items()));return self*L({tuple(-a for a in e):1/v})
    def __rtruediv__(self,o):return L(o)/self
    def __pow__(self,n):
        if n<0:return L(1)/(self**(-n))
        r=L(1)
        for _ in range(n):r=r*self
        return r
    def __eq__(self,o):return self.d==L(o).d


def zero(n=N,m=None):return [[Q(0) for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n=N):return [[Q(i==j) for j in range(n)] for i in range(n)]
def summat(*ms):return [[sum(m[i][j] for m in ms) for j in range(len(ms[0][0]))] for i in range(len(ms[0]))]
def times(c,m):return [[c*v for v in row] for row in m]
def pair(out,inc,c,m,l):
    a=f'{out}:{l}';b=f'{inc}:{l}'
    if a in INDEX and b in INDEX:m[INDEX[a]][INDEX[b]]+=c
def isotropic(table):
    m=zero()
    for out,inc,c in table:
        for l in (0,2,4):pair(out,inc,Q(c),m,l)
        if out==inc=='H':
            for label in ('T:+','T:-'):m[INDEX[label]][INDEX[label]]=Q(c)
    return m


def independent_primitives(mutation):
    tables={
      'O':[('K','L',-2),('N','L',-2),('L','K',-1),('L','N',-1),('v','u',-1),('u','v',-1)],
      'O2':[('K','K',2),('K','N',2),('N','K',2),('N','N',2),('L','L',4),('v','v',1),('u','u',1)],
      'Sa':[('K','K',-6),('K','J',6),('L','L',-3),('J','K',2),('J','J',-2),('v','v',-6),('u','u',-1),('H','H',-2)],
      'Sh':[('N','N',6),('N','J',6),('L','L',3),('J','N',2),('J','J',2),('v','v',1),('u','u',6),('H','H',2)],
      'Sah':[('K','L',6),('L','K',-3),('L','N',3),('L','J',6),('N','L',-6),('J','L',-4),('v','u',5),('u','v',-5)],
      'Rb':[('K','K',4),('K','N',1),('K','J',-3),('L','L',Q(7,2)),('N','K',1),('N','N',1),('J','K',-1),('J','J',1),('v','v',Q(7,2)),('u','u',1),('H','H',1)],
      'Rd':[('L','L',Q(3,2)),('N','N',3),('N','J',3),('J','N',1),('J','J',1),('v','v',Q(3,2)),('u','u',4),('H','H',4)],
      'Aphi':[('L','w',2),('w','L',1)],
      'Bphi':[('K','w',-4),('w','K',Q(1,2)),('w','N',Q(1,2)),('w','J',-Q(3,2))],
      'Pg':[('K','w',1),('N','w',-1),('J','w',1)],
      'Pww':[('w','w',1)],'W_K':[('w','K',1)],'W_L':[('w','L',1)],
      'By_a':[('N','K',-1),('J','K',1),('K','K',-4),('K','N',-2),('K','J',6)],
      'By_p':[('w','K',Q(1,2))],
      'By_sigma1':[('N','w',-Q(1,3)),('J','w',Q(1,3)),('K','w',Q(4,3))],
      'By_sigma2':[('w','w',-Q(1,2))],
    }
    if mutation=='fibre_time_sign':tables['O'][0]=('K','L',2)
    if mutation=='coordinate_Robin':tables['By_a'].append(('J','J',2))
    m={name:isotropic(table) for name,table in tables.items()}
    m.update({'I':eye(),'Delta':zero(),'Ca':zero(),'Ch':zero(),'Ephi':zero(),'Pmet':zero(),'PD':zero(),'PN':eye()})
    for label,i in INDEX.items():
        if label.startswith('T:'):eigen=6
        else:
            ch,l=label.split(':');l=int(l);eigen=l*(l+2)-(2 if ch in ('v','u') else 6 if ch=='H' else 0)
        m['Delta'][i][i]=-Q(eigen)
        if not label.startswith('w:'):m['Pmet'][i][i]=1
        if label.startswith(('L:','v:')):m['PD'][i][i]=1;m['PN'][i][i]=0
    if mutation=='rough_tensor_shift':m['Delta'][INDEX['H:2']][INDEX['H:2']]-=1
    for l in (0,2,4):
        lam=l*(l+2);beta=2-Q(2,3)*lam
        # Sum_a Omega_a D_a, derived directly by applying the connection
        # to the yy,y0,00,yi,0i,ij components before harmonic projection.
        for o,i,c in [('K','v',2*lam),('L','u',lam),('v','K',1),('v','J',-1),('v','H',-beta),('u','L',1),('J','v',-Q(2,3)*lam),('H','v',2)]:pair(o,i,Q(c),m['Ca'],l)
        for o,i,c in [('L','v',lam),('N','u',2*lam),('v','L',-1),('u','N',-1),('u','J',-1),('u','H',-beta),('J','u',Q(2,3)*lam),('H','u',-2)]:pair(o,i,Q(c),m['Ch'],l)
        pair('v','w',Q(2),m['Ephi'],l);pair('w','v',Q(lam),m['Ephi'],l)
    if mutation=='angular_connection_sign':m['Ca'][INDEX['H:2']][INDEX['v:2']]=-2
    return m


def expected_assembly(v,m,mutation):
    ell,rho,a,h,c,Ht,p,py,b,d,V,V1,V2,sigma1,sigma2=v
    A=summat(times(-3*ell*Ht,m['I']),times(-2*ell*rho*a,m['O']),times(ell*rho*p,m['Aphi']))
    B=summat(times(-3*ell*rho*a,m['I']),times(ell*rho*p,m['Bphi']))
    # Keep independent geometric contributions separated; do not call the
    # candidate's assembly or use its combined scalar row as the derivation.
    connection=summat(times(a*a,m['Sa']),times(h*h,m['Sh']),times(a*h,m['Sah']),times(-a*a,m['O2']),times(-3*a*h,m['O']))
    sign=2 if mutation=='curvature_sign' else -2
    value=summat(connection,times(sign*b,m['Rb']),times(sign*d,m['Rd']),times(Q(4,3)*V,m['Pmet']),times(Q(4,3)*V1,m['Pg']),times(-V2,m['Pww']))
    for l in (0,2,4):
        # -k:Hess(phi0) contributes -py K+ap N-3ap J.
        for inc,coefficient in [('K',-py),('N',a*p),('J',-3*a*p)]:pair('w',inc,coefficient,value,l)
        # -p*dC_y contributes -4ap K-ap N+3ap J+3ph L.
        if mutation!='omit_scalar_DC':
            for inc,coefficient in [('K',-4*a*p),('N',-a*p),('J',3*a*p),('L',3*p*h)]:pair('w',inc,coefficient,value,l)
    D=times(ell*ell,summat(times(1/(c*c),m['Delta']),times(2*rho*a/c,m['Ca']),times(2*rho*h/c,m['Ch']),times(rho*rho,value),times(rho*p/c,m['Ephi'])))
    M=times(-ell*rho,summat(times(a,m['By_a']),times(p,m['By_p']),times(sigma1,m['By_sigma1']),times(sigma2,m['By_sigma2'])))
    return {'A':A,'B':B,'D':D,'M':M,'PD':m['PD'],'PN':m['PN']}


def sphere_expected():
    out={name:zero(rows,cols) for name,rows,cols in [('grad',2,3),('hess',7,3),('div_vector',3,2),('symgrad',7,2),('div_tensor',2,7),('trace',3,7),('lap_scalar',3,3),('lap_vector',2,2),('lap_tensor',7,7)]}
    for i,lam in enumerate((8,24)):
        col=i+1;hrow=i+3
        out['grad'][i][col]=1
        out['hess'][col][col]=-Q(lam,3);out['hess'][hrow][col]=1
        out['div_vector'][col][i]=-lam
        out['symgrad'][col][i]=-Q(2*lam,3);out['symgrad'][hrow][i]=2
        out['div_tensor'][i][col]=1;out['div_tensor'][i][hrow]=2-Q(2*lam,3)
    for i in range(3):out['trace'][i][i]=3
    for key,values in [('lap_scalar',[0,-8,-24]),('lap_vector',[-6,-22]),('lap_tensor',[0,-8,-24,-2,-18,-6,-6])]:
        for i,v in enumerate(values):out[key][i][i]=v
    return out


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def interface_audit(candidate,variables):
    ell,rho=variables[:2];count=0
    fields={'yy':[1,2,3],'y0':[4,5,6],'00':[7,8,9],'ab':[10,11,12,20,21,22,23],
            'scalar':[13,14,15],'ya':[16,17],'0a':[18,19]}
    need(candidate.pack_fields(fields)==list(range(1,24)),'ALL_FIELD_PACKING');count+=23
    metric={k:v for k,v in fields.items() if k!='scalar'}
    forced=candidate.bulk_forcing(ell,rho,metric,fields['scalar'])
    for i,label in enumerate(LABELS):
        expected=(1 if label.startswith('w:') else -2)*ell*ell*rho*(i+1)*rho
        need(L(forced[i])==expected,'FACTORIAL_BULK_SOURCE_PACKING');count+=1
    missing={k:v for k,v in metric.items() if k!='ya'}
    try:candidate.bulk_forcing(ell,rho,missing,[0,0,0])
    except ValueError as error:need(str(error)=='COMPLETE_QERED_REQUIRED','MISSING_QE_DIAGNOSTIC')
    else:raise ValueError('MISSING_QE_ACCEPTED')
    source={'GN':[1,2,3],'QK00':[7,8,9],'QKab':[10,11,12,20,21,22,23],
            'Qp':[13,14,15],'QK0a':[18,19]}
    result=candidate.boundary_forcing(ell,rho,source)
    for i,label in enumerate(LABELS):
        factor=0 if label.startswith(('L:','v:')) else -1 if label.startswith('K:') else 1 if label.startswith('w:') else 2
        need(L(result[i])==factor*ell*rho*(i+1),'ALL_BOUNDARY_SOURCE_PACKING');count+=1
    initial=dict(fields);del initial['ab']
    initial.update({'ab_trace_amplitude':[10,11,12],'ab_HessTF':[20,21],'ab_TT':[22,23]})
    u,v=candidate.pack_initial({'values':initial,'t_velocities':initial},ell)
    for i in range(23):
        need(u[i]==i+1 and L(v[i])==ell*(i+1),'INITIAL_SCALED_VELOCITY_PACKING');count+=2
    return count


def run(mutation):
    base=Path(__file__).resolve().parent;pins=json.loads((base/'INPUT_PINS.json').read_text())
    for key,item in pins.items():need(sha256((base/item['path']).read_bytes()).hexdigest()==item['sha256'],'INPUT_HASH_'+key)
    candidate=load(base/'inputs/linear_operator.py','candidate_linear_operator')
    need(candidate.LABELS==LABELS,'BASIS_ORDER')
    expected=independent_primitives(mutation);actual=candidate.matrices()
    need(set(expected)==set(actual),'PRIMITIVE_ROSTER')
    comparisons=0
    for name in expected:
        need(expected[name]==actual[name],'PRIMITIVE_'+name);comparisons+=N*N
    norms={0:Q(1),2:Q(1,12),4:Q(19,720)};g=[]
    for name in LABELS:
        if name.startswith('T:'):q=Q(1,6)
        else:
            ch,l=name.split(':');l=int(l);q=norms[l]
            q*= {'K':1,'L':1 if mutation=='mixed_weight' else 2,'N':1,'J':3,'w':1,'v':2*l*(l+2),'u':2*l*(l+2),'H':Q(2,3)*l*(l+2)*(l*(l+2)-3)}[ch]
        g.append(q)
    need(g==candidate.gram(),'PHYSICAL_GRAM')
    variables=[]
    for i in range(NV):
        e=[0]*NV;e[i]=1;variables.append(L({tuple(e):1}))
    direct=expected_assembly(variables,expected,mutation);assembled=candidate.assemble(*variables)
    coefficients=0
    for name in direct:
        for i in range(N):
            for j in range(N):
                need(L(direct[name][i][j])==L(assembled[name][i][j]),'FORMAL_ASSEMBLY_'+name)
                coefficients+=len(L(direct[name][i][j]).d)
    interface_checks=interface_audit(candidate,variables)
    angular=load(base/'inputs/angular_exact.py','candidate_angular_exact')
    primitive=angular.primitive_matrices();sphere=sphere_expected()
    for key in sphere:need(primitive[key]==sphere[key],'INTRINSIC_PRIMITIVE_'+key)
    return {'status':'PASS_INDEPENDENT_23_OPERATOR_AND_PRIMITIVES',
      'primitive_matrix_entry_comparisons':comparisons,'formal_full_matrix_entry_comparisons':6*N*N,
      'formal_nonzero_Laurent_coefficients':coefficients,'formal_variables':VARS,
      'linear_interface_entry_checks':interface_checks,'incomplete_Einstein_source_rejected':True,
      'sphere_primitives':list(sphere),'physical_Gram':list(map(str,g)),
      'method':'Candidate exact polynomial reconstruction is compared against independently derived sphere, frame and curvature tables. Full assembly is compared as formal Laurent polynomials, not sampled backgrounds.',
      'scope':'Homogeneous 23-amplitude bulk/boundary operator only; actual source projection, initial adapter and numerical evolution require separate acceptance.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--inject',choices=['fibre_time_sign','rough_tensor_shift','coordinate_Robin','angular_connection_sign','mixed_weight','curvature_sign','omit_scalar_DC']);a=p.parse_args()
    need(not a.output.exists(),'OUTPUT_EXISTS')
    need(not a.output.resolve().is_relative_to(Path(__file__).resolve().parent),'OUTPUT_MUST_BE_EXTERNAL')
    result=run(a.inject);need(a.inject is None,'CONTROL_NOT_REJECTED')
    result['checker_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(result['status'])


if __name__=='__main__':main()
