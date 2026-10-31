#!/usr/bin/env python3
"""Private exact TT block/current-data feasibility; no PDE or curvature replay."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PINS_SHA='92d315be572092b9c4e39ac2c6ed7ce2cccff81edeecb0a03bb3ff5e4ebee7e2'
G='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/'
A='research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/'
E='research/HDBLAST_EXACT_PACKET_INITIAL_20261031/'

def need(ok,why):
    if not ok:raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_bytes())
def authenticate(root,velocity):
    need(sha(HERE/'SOURCE_PINS.json')==PINS_SHA,'pinned feasibility inputs');pins=read(HERE/'SOURCE_PINS.json');count=0
    for name,pin in pins['source_manifest_sha256'].items():
        cap=root/'research'/name;need(sha(cap/'MANIFEST.json')==pin,'manifest pin '+name);seen={'MANIFEST.json'}
        for row in read(cap/'MANIFEST.json')['files']:
            rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'safe manifest member');seen.add(str(rel));p=cap/rel
            need(p.is_file() and not p.is_symlink() and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'authenticated source member '+str(rel));count+=1
        need(all(not p.is_symlink() and (not p.is_file() or str(p.relative_to(cap)) in seen) for p in cap.rglob('*')),'closed source capsule '+name)
    for rel,pin in pins['source_sha256'].items():need(sha(root/rel)==pin,'selected original source '+rel)
    for rel,pin in pins['velocity_sha256'].items():need(sha(velocity/rel)==pin,'selected current Cauchy source '+rel)
    return {'source_pin_map_sha256':PINS_SHA,'closed_capsule_members':count,'manifest_count':len(pins['source_manifest_sha256'])}
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

class Expr:
    """Expression DAG with only exact zero/constant rules; never numerical sampling."""
    def __init__(self,value=None,op='symbol',args=()):self.value=None if value is None else Q(value);self.op=op;self.args=args
    @staticmethod
    def cv(x):return x if isinstance(x,Expr) else Expr(x)
    def zero(self):return self.value==0
    def __add__(self,y):
        y=Expr.cv(y)
        if self.zero():return y
        if y.zero():return self
        if self.value is not None and y.value is not None:return Expr(self.value+y.value)
        return Expr(op='+',args=(self,y))
    __radd__=__add__
    def __neg__(self):return Expr(-self.value) if self.value is not None else Expr(op='-',args=(self,))
    def __sub__(self,y):return self+-Expr.cv(y)
    def __rsub__(self,y):return Expr.cv(y)+-self
    def __mul__(self,y):
        y=Expr.cv(y)
        if self.zero() or y.zero():return Expr(0)
        if self.value is not None and y.value is not None:return Expr(self.value*y.value)
        return Expr(op='*',args=(self,y))
    __rmul__=__mul__
    def __truediv__(self,y):
        y=Expr.cv(y);need(not y.zero(),'formal denominator nonzero')
        if self.zero():return Expr(0)
        if self.value is not None and y.value is not None:return Expr(self.value/y.value)
        return Expr(op='/',args=(self,y))
    def __rtruediv__(self,y):return Expr.cv(y)/self
    def __pow__(self,n):
        if self.value is not None:return Expr(self.value**n)
        return Expr(1) if n==0 else Expr(op='power',args=(self,n))
    def dz(self):return Expr(0) if self.value is not None else Expr(op='formal_derivative',args=(self,))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);ap.add_argument('--velocity-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--fault',choices=('bulk_cross_block','boundary_cross_block','registered_TT_initial','source_orientation'))
    x=ap.parse_args();need(not x.output.exists(),'fresh output');before=authenticate(x.source_root,x.velocity_root)
    # All selected source bytes and closed capsules are authenticated before imports.
    op=load('tt_actual_registered_operator',x.source_root/G/'projection/inputs/linear_operator.py')
    adapter=load('tt_actual_registered_adapters',x.source_root/G/'projection/inputs/source_adapters.py')
    algebra=load('tt_retained_exact_Laurent',x.source_root/A/'structure/exact_angular_structure.py')
    initial=load('tt_completed_initial_expression',x.source_root/A/'initial/full23_initial.py')
    tt=[op.IDX['T:+'],op.IDX['T:-']];need(tt==[21,22] and len(op.LABELS)==23,'actual registered23 ordering')
    primitives=op.matrices()
    if x.fault=='bulk_cross_block':primitives['Aphi'][21][12]=Q(1)
    for name,m in primitives.items():
        need(all(m[i][j]==0 and m[j][i]==0 for i in tt for j in range(23) if j!=i),'TT primitive row/column closure:'+name)
        need(m[21][21]==m[22][22],'same TT orientation operator:'+name)
    z=algebra.symbols('ell rho a h c Ht p py b d V V1 V2 sigma1 sigma2');mat=op.assemble(**z)
    if x.fault=='boundary_cross_block':mat['M'][21][6]=algebra.L(1)
    for name,m in mat.items():need(all(not m[i][j] and not m[j][i] for i in tt for j in range(23) if j!=i),'assembled TT block:'+name)
    need(all(not mat['M'][i][j] and not mat['M'][j][i] for i in tt for j in range(23)),'actual TT M row/column exactly zero')
    expected={'A':-3*z['ell']*z['Ht'],'B':-3*z['ell']*z['rho']*z['a'],'D':z['ell']*z['ell']*(-6/(z['c']*z['c'])+z['rho']*z['rho']*(-2*z['a']*z['a']+2*z['h']*z['h']-2*z['b']-8*z['d']+Q(4,3)*z['V'])),'M':algebra.L(0)}
    for name,value in expected.items():need(mat[name][21][21]==mat[name][22][22]==value,'actual TT symbolic coefficient '+name)
    need(mat['PD'][21][21]==mat['PD'][22][22]==0 and mat['PN'][21][21]==mat['PN'][22][22]==1,'Neumann TT projectors')
    stored=read(x.source_root/A/'compute/response-normal/RESPONSE_REPORT.json')
    for key,m in stored['constant_matrices'].items():
        need(all(Q(m[i][j])==0 and Q(m[j][i])==0 for i in tt for j in range(23) if j!=i),'stored frozen L0 closure '+key)
        if key=='M':need(all(Q(m[i][j])==0 and Q(m[j][i])==0 for i in tt for j in range(23)),'stored frozen M exact zero TT')
    for key,m in stored['actual_coefficient_boxes'].items():
        need(all(all(Q(t)==0 for t in m[i][j]+m[j][i]) for i in tt for j in range(23) if j!=i),'actual coefficient interval closure '+key)
        if key=='M':need(all(all(Q(t)==0 for t in m[i][j]+m[j][i]) for i in tt for j in range(23)),'actual shell M interval zero TT')
    maps=read(x.source_root/G/'projection/GEOMETRIC_COEFFICIENT_MAPS.json')['maps']
    if x.fault=='source_orientation':maps['metric_forcing']['ab'][5]=[]
    equalmaps={}
    for name,fields in maps.items():
        if isinstance(fields,dict) and 'ab' in fields:
            need(len(fields['ab'])==7 and fields['ab'][5]==fields['ab'][6],'complete geometric TT source orientation '+name);equalmaps[name]=len(fields['ab'][5])
    need('metric_forcing' in equalmaps and 'QEred' in equalmaps,'complete geometric source maps')
    sh=algebra.symbols('a p sigmap sigmapp sigmappp h X V L d v');shell=op.boundary_forcing(z['ell'],z['rho'],adapter.shell_sources(**sh))
    gt=-Q(8,3)*z['ell']*z['rho']*(sh['a']*sh['h']+sh['sigmap']*sh['v']/3)*sh['d']
    need(shell[21]==shell[22]==gt and bool(gt),'equal generally nonzero actual TT shell forcing')
    # Evaluate the complete registered expression on generic expression DAGs.
    # The only simplifications used to certify zero are exact algebraic identities.
    mean_ids=[op.IDX[k+':0'] for k in ('K','L','N','J','w')];need(mean_ids==[0,3,6,9,12],'only five scalar mean overrides')
    dummy=SimpleNamespace(init=SimpleNamespace(UP={1:'U1',2:'U2'},poly=lambda coefficients,eta:Expr(op='fixed_polynomial',args=(coefficients,eta))),op=op,adapter=adapter,ell=Q(1,10**6),mean_ids=mean_ids)
    displacement,velocity=initial.Context._registered(dummy,{k:Expr(op=k) for k in ('rho','eta','p','a','alpha')},Expr(op='exact_E'),Expr(op='exact_q'))
    if x.fault=='registered_TT_initial':velocity[21]=Expr(1)
    need(all(displacement[i].zero() and velocity[i].zero() for i in tt),'completed registered TT initial value/velocity identically zero')
    profile=read(x.source_root/A/'initial/results/INITIAL_PROFILES.json');need(profile['labels']==op.LABELS,'same old P1 all23 labels')
    need(all(Q(row[i])==0 for rows in profile['profiles'].values() for row in rows for i in tt),'old all23 stored P1 TT channels zero')
    old=read(x.source_root/E/'reference/INITIAL.json')['campaigns'][0];new=read(x.velocity_root/'INITIAL_NORMAL.json')['campaigns'][0]
    for label,report in [('old_exact_packet',old),('new_exact_Cauchy',new)]:
        need(report['labels']==op.LABELS,'current candidate label identity')
        need(all(Q(caps[i])==0 for leaf in report['leaves'] for caps in leaf['absolute_component_caps'] for i in tt),'completed candidate all-leaf TT zero '+label)
    geometry=read(x.source_root/A/'structure/ANGULAR_STRUCTURE_NORMAL.json')['TT_geometry']
    need(geometry['trace']==geometry['divergence']==0,'retained exact TT trace/divergence identities')
    need(authenticate(x.source_root,x.velocity_root)==before,'sources unchanged after exact queries')
    result={'status':'PASS_PRIVATE_TT_BLOCK_AND_CURRENT_DATA_FEASIBILITY','authentication':before,'labels':op.LABELS,'TT_indices':tt,
      'primitive_matrices_checked':len(primitives),'actual_and_frozen_TT_rows_and_columns_close':True,'actual_and_frozen_M_TT_rows_columns_zero':True,
      'symbolic_TT_coefficients':{key:value.json() for key,value in expected.items()},'TT_shell_forcing':gt.json(),'TT_shell_forcing_identically_zero':False,
      'complete_geometric_TT_orientation_map_term_counts':equalmaps,'completed_registered_old_new_TT_initial_functions_zero':True,'mean_override_indices':mean_ids,
      'stored_old_P1_channels_TT_zero':True,'derivative_extension_convolution_and_parity_preserve_TT_zero':'linear componentwise operations preserve the identically zero completed functions and all their jets',
      'constraint_scope':'Retained trace/divergence-free TT identities eliminate its linear scalar/vector subsidiary projections; physical exact-target completion remains conditional. Nominal full candidates are not newly certified to satisfy constraints.',
      'conditional_finite_response_identity':{'P_T':'identity because complete TT M row is zero','W0_T':'L G_T','W1_T':'L G_T + K(La_T L G_T + F_T)','old_vs_new':'Same operators and original F/G; exact TT finite functions equal when only nonTT initial data change','protected_cone':'0<=x<=1/4,1+x<=r<=3/2-x','protected_W0_T':'0','protected_W1_T':'K F_T'},
      'scope':'Private feasibility. No new PDE existence/uniqueness proof, screen-curvature certificate, point table, nonlinear shell-contact or cosmological claim.'}
    x.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':result['status'],'primitive_matrices':len(primitives)}))
if __name__=='__main__':main()
