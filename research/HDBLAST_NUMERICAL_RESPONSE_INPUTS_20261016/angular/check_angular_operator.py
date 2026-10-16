#!/usr/bin/env python3
"""Independent all-point angular reconstruction and symbolic PDE checks."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as Q
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from exact_laurent import variables
import angular_exact as a
import linear_operator as l

CONTROLS={'rough_tensor_shift':'PRIMITIVE_LAPLACIAN','cross_sign':'CROSS_RECONSTRUCTION','drop_TT':'TT_RETAINED','omit_scalar_DC':'FULL_SYMBOLIC_BULK','curvature_mixed':'CURVATURE','boundary_frame':'FULL_SYMBOLIC_BOUNDARY','dirichlet_mixed':'BOUNDARY_SPLIT'}
def need(ok,label):
    if not ok:raise ValueError(label)
def trans(t):return {(i,j):t[j,i] for i,j in t}
def zero_section():return {k:a.zero_tensor(r) for k,r in [('K',0),('L',0),('N',0),('w',0),('v',1),('u',1),('m',2)]}
def basis_section(name):
    out=zero_section();field,degree=name.split(':')
    if field in ('K','L','N','w'):out[field]=a.SCALARS[(0,2,4).index(int(degree))]
    elif field=='J':out['m']=a.TENSORS[(0,2,4).index(int(degree))]
    elif field in ('v','u'):out[field]=a.GRADS[(2,4).index(int(degree))]
    elif field=='H':out['m']=a.HESS[(2,4).index(int(degree))]
    else:out['m']=a.TT[(0 if degree=='+' else 1)]
    return out
def sum_section(coeff):
    out=zero_section()
    for c,b in zip(coeff,BASE):
        if c:
            for k in out:out[k]=a.tadd(out[k],a.tscale(c,b[k]))
    return out
BASE=[basis_section(name) for name in l.LABELS]
def cross(v,which):
    o=zero_section();dv=a.div(v['v']);du=a.div(v['u']);dm=a.div(v['m'])
    if which=='a':
        o['K']=a.tscale(-2,dv);o['L']=a.tscale(-1,du)
        o['v']=a.tadd(a.covd(v['K']),a.tscale(-1,dm));o['u']=a.covd(v['L'])
        d=a.covd(v['v']);o['m']=a.tadd(d,trans(d))
    else:
        o['L']=a.tscale(-1,dv);o['N']=a.tscale(-2,du)
        o['v']=a.tscale(-1,a.covd(v['L']))
        o['u']=a.tscale(-1,a.tadd(a.covd(v['N']),dm))
        d=a.covd(v['u']);o['m']=a.tscale(-1,a.tadd(d,trans(d)))
    return o
def expected_channel(entries):
    m=l.zeros()
    for out,inc,c in entries:
        for degree in (0,2,4):l.put(m,out,inc,c,degree)
        if out==inc=='H':
            for name in ['T:+','T:-']:m[l.IDX[name]][l.IDX[name]]=Q(c)
    return m

def symbolic_bulk(m,mutation):
    ell,rho,aa,h,c,Ht,p,py,b,d,V,V1,V2,sg1,sg2=variables()
    assembled=l.assemble(ell,rho,aa,h,c,Ht,p,py,b,d,V,V1,V2,sg1,sg2)
    if mutation=='omit_scalar_DC':
        for degree in (0,2,4):assembled['A'][l.IDX[f'w:{degree}']][l.IDX[f'L:{degree}']]=0
    if mutation=='boundary_frame':
        for degree in (0,2,4):assembled['M'][l.IDX[f'J:{degree}']][l.IDX[f'J:{degree}']]=-2*ell*rho*aa
    add,scale,mm=l.add,l.scale,l.mm
    Om=scale(aa,m['O'])
    S2=add(scale(aa**2,m['Sa']),scale(aa*h,m['Sah']),scale(h**2,m['Sh']))
    wave_x=scale(-1/ell,add(scale(3*Ht,m['I']),scale(2*rho,Om)))
    wave_r=scale(-3*aa*rho/ell,m['I'])
    wave_0=add(scale(1/c**2,m['Delta']),scale(2*rho/c,add(scale(aa,m['Ca']),scale(h,m['Ch']))),scale(rho**2,add(scale(-1,mm(Om,Om)),scale(-3*h,Om),S2)))
    rhs_0=add(scale(2,add(scale(b,m['Rb']),scale(d,m['Rd']))),scale(-Q(4,3)*V,m['Pmet']),scale(-Q(4,3)*V1,m['Pg']))
    rhs_x=l.zeros();rhs_r=l.zeros();boundary=l.zeros()
    for degree in (0,2,4):
        idx=lambda ch:l.IDX[f'{ch}:{degree}']
        # Uncombined Hess(phi0): py*K + a*p*(-N+3J).
        rhs_0[idx('w')][idx('w')]+=V2
        rhs_0[idx('w')][idx('K')]+=py
        rhs_0[idx('w')][idx('N')]+=-aa*p
        rhs_0[idx('w')][idx('J')]+=3*aa*p
        # Independent divergence trace formula, with all derivative slots.
        rhs_x[idx('w')][idx('L')]+=-p/(ell*rho)
        for ch,factor in [('K',-Q(1,2)),('N',-Q(1,2)),('J',Q(3,2))]:rhs_r[idx('w')][idx(ch)]+=factor*p/(ell*rho)
        for ch,factor in [('K',4),('N',1),('J',-3)]:rhs_0[idx('w')][idx(ch)]+=factor*aa*p
        rhs_0[idx('w')][idx('L')]+=-3*h*p
        if degree:rhs_0[idx('w')][idx('v')]+=-p*degree*(degree+2)/(rho*c)
        # -2 (dphi tensor dw+dw tensor dphi), before solving for U_xx.
        rhs_r[idx('K')][idx('w')]+=4*p/(ell*rho)
        rhs_x[idx('L')][idx('w')]+=-2*p/(ell*rho)
        if degree:rhs_0[idx('v')][idx('w')]+=-2*p/(rho*c)
        boundary[idx('N')][idx('K')]=-aa;boundary[idx('N')][idx('w')]=-sg1/3
        boundary[idx('J')][idx('K')]=aa;boundary[idx('J')][idx('w')]=sg1/3
        boundary[idx('w')][idx('K')]=p/2;boundary[idx('w')][idx('w')]=-sg2/2
        boundary[idx('K')][idx('K')]=-4*aa;boundary[idx('K')][idx('N')]=-2*aa
        boundary[idx('K')][idx('J')]=6*aa;boundary[idx('K')][idx('w')]=4*sg1/3
    for key,wave,rhs in [('A',wave_x,rhs_x),('B',wave_r,rhs_r),('D',wave_0,rhs_0)]:
        wanted=scale(ell**2,add(wave,scale(-rho**2,rhs)))
        need(wanted==assembled[key],'FULL_SYMBOLIC_BULK')
    need(scale(-ell*rho,boundary)==assembled['M'],'FULL_SYMBOLIC_BOUNDARY')
    return {'formal_coefficient_rows':23,'independent_jet_slots':69,'background_symbols':15,'bulk_identity':True,'all_boundary_rows':True,'arithmetic':'Exact rational Laurent polynomials; no parameter sampling'}

def run(mutation=None):
    prim=a.primitive_matrices(mutation)
    need(prim['lap_tensor']==[[Q(-v) if i==j else Q(0) for j,v in enumerate([0,8,24,2,18,6,6])] for i in range(7)],'PRIMITIVE_LAPLACIAN')
    gramsets=[(a.SCALARS,[Q(1),Q(1,12),Q(19,720)]),(a.GRADS,[Q(2,3),Q(19,30)]),(a.TENSORS,[Q(3),Q(1,4),Q(19,240),Q(20,9),Q(133,15),Q(1,6),Q(1,6)])]
    for basis,diag in gramsets:
        for i,t in enumerate(basis):
            for j,u in enumerate(basis):need(a.inner(t,u)==(diag[i] if i==j else 0),'EXACT_GRAM')
    m=l.matrices()
    expected_o=expected_channel([('K','L',-2),('L','K',-1),('L','N',-1),('N','L',-2),('v','u',-1),('u','v',-1)])
    need(m['O']==expected_o,'TIME_CONNECTION')
    for i,v in enumerate(BASE):
        g=sum((a.inner(v[ch],v[ch])*(2 if ch in ('L','v','u') else 1) for ch in v),Q(0))
        need(g==l.gram()[i],'FULL_FROBENIUS_GRAM')
    if mutation=='cross_sign':m['Ca'][l.IDX['v:2']][l.IDX['H:2']]*=-1
    if mutation=='drop_TT':m['Delta'][l.IDX['T:+']][l.IDX['T:+']]=0
    if mutation=='curvature_mixed':m['Rb'][l.IDX['L:2']][l.IDX['L:2']]=Q(5,2)
    if mutation=='dirichlet_mixed':m['PD'][l.IDX['v:2']][l.IDX['v:2']]=0
    for which in ('a','h'):
        for j,v in enumerate(BASE):
            direct=cross(v,which);recon=sum_section([row[j] for row in m['C'+which]])
            for k in direct:a.expect(direct[k],recon[k],'CROSS_RECONSTRUCTION')
    need(m['Delta'][21][21]==m['Delta'][22][22]==-6,'TT_RETAINED')
    expect_b=expected_channel([('K','K',4),('K','N',1),('K','J',-3),('L','L',Q(7,2)),('v','v',Q(7,2)),('N','K',1),('N','N',1),('J','K',-1),('J','J',1),('u','u',1),('H','H',1)])
    expect_d=expected_channel([('L','L',Q(3,2)),('v','v',Q(3,2)),('N','N',3),('N','J',3),('J','N',1),('J','J',1),('u','u',4),('H','H',4)])
    need(m['Rb']==expect_b and m['Rd']==expect_d,'CURVATURE')
    need(sum(m['PD'][i][i] for i in range(23))==5 and sum(m['PN'][i][i] for i in range(23))==18 and l.add(m['PD'],m['PN'])==l.ident(),'BOUNDARY_SPLIT')
    # Compare fibre-generated connection squares to an independent channel table.
    sa=expected_channel([('K','K',-6),('K','J',6),('L','L',-3),('v','v',-6),('u','u',-1),('J','K',2),('J','J',-2),('H','H',-2)])
    sh=expected_channel([('L','L',3),('N','N',6),('N','J',6),('v','v',1),('u','u',6),('J','N',2),('J','J',2),('H','H',2)])
    sah=expected_channel([('K','L',6),('L','K',-3),('L','N',3),('L','J',6),('N','L',-6),('J','L',-4),('v','u',5),('u','v',-5)])
    need(m['Sa']==sa and m['Sh']==sh and m['Sah']==sah,'CONNECTION_SQUARES')
    symbolic=symbolic_bulk(m,mutation)
    return {'status':'PASS_EXACT_23_AMPLITUDE_LINEAR_OPERATOR','scope':'Complete homogeneous bulk and boundary operator. Full tensor Einstein quadratic source adapter remains separate. No evaluated response.','basis':l.LABELS,'Gram_diagonal':list(map(str,l.gram())),'primitive_matrices':{k:[[str(v) for v in row] for row in mat] for k,mat in prim.items()},'matrix_nonzero_entries':{key:[[i,j,str(v)] for i,row in enumerate(mat) for j,v in enumerate(row) if v] for key,mat in m.items()},'symbolic_covariant_reconstruction':symbolic,'angular_reconstruction':'Exact zero polynomials modulo |n|²−1 for every primitive and all 46 angular connection columns','boundary_counts':{'Dirichlet':5,'Robin':18},'pointwise_factor':17}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--inject',choices=sorted(CONTROLS));args=p.parse_args()
    source_root=Path(__file__).resolve().parent
    need(source_root not in args.output.resolve().parents and args.output.resolve()!=source_root,'OUTPUT_MUST_BE_EXTERNAL')
    need(not args.output.exists(),'OUTPUT_EXISTS');r=run(args.inject);need(args.inject is None,'CONTROL_NOT_REJECTED')
    r['pins']={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['angular_exact.py','linear_operator.py','exact_laurent.py','check_angular_operator.py']}
    args.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(r['status'])
if __name__=='__main__':main()
