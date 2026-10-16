#!/usr/bin/env python3
"""Independent exact linearization witnesses and full-response norm audit.

No producer code is imported. Curved, nonnormal coordinate jets are expanded
directly and compared with the covariant linear operator. Finite witnesses
supplement, rather than replace, the all-point proof in EVOLUTION_REVIEW.md.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

N=5
I=range(N)
P2=list(product(I,repeat=2))
P3=list(product(I,repeat=3))
P4=list(product(I,repeat=4))
CONTROLS={
    'curvature_sign':'CURVED_RICCI_LINEARIZATION',
    'derivative_slot':'CURVED_RICCI_LINEARIZATION',
    'symmetrization':'CURVED_RICCI_LINEARIZATION',
    'scalar_DC':'SCALAR_LINEARIZATION',
    'omit_TT':'ANGULAR_MULTIPLICITY',
    'scalar_component_cutoff':'ANGULAR_MULTIPLICITY',
    'evaluation_factor':'EVALUATION_TRACE',
    'missing_boundary_weight':'BOUNDARY_TRACE_WEIGHT',
    'wrong_initial_frame_sign':'INITIAL_FRAME_CONVERSION',
}


def need(ok,message):
    if not ok: raise ValueError(message)


def zero(): return [[Q(0) for _ in I] for _ in I]
def mul(a,b): return [[sum(a[i][k]*b[k][j] for k in I) for j in I] for i in I]
def trans(a): return [list(row) for row in zip(*a)]
def neg(a): return [[-v for v in row] for row in a]


def inv(a):
    x=[list(a[i])+[Q(i==j) for j in I] for i in I]
    for j in I:
        k=next(k for k in range(j,N) if x[k][j])
        x[k],x[j]=x[j],x[k]
        v=x[j][j]; x[j]=[z/v for z in x[j]]
        for i in I:
            if i!=j:
                v=x[i][j]; x[i]=[u-v*w for u,w in zip(x[i],x[j])]
    return [row[N:] for row in x]


def jets(seed):
    gen=random.Random(seed)
    def get(): return Q(gen.randrange(-3,4),gen.randrange(1,4))
    def symmetric():
        a=zero()
        for i in I:
            for j in range(i,N): a[i][j]=a[j][i]=get()
        return a
    def second():
        a=[[None for _ in I] for _ in I]
        for i in I:
            for j in range(i,N): a[i][j]=a[j][i]=symmetric()
        return a
    l=zero(); eta=zero()
    for i in I:
        eta[i][i]=Q(-1 if i==1 else 1); l[i][i]=Q(1)
        for j in range(i+1,N): l[i][j]=get()/5
    g=mul(trans(l),mul(eta,l))
    return (g,[symmetric() for _ in I],second(),symmetric(),
            [symmetric() for _ in I],second(),[get() for _ in I],
            symmetric(),get(),[get() for _ in I],symmetric(),get(),get(),get())


def curved_audit(seed,mutation):
    g,dg,ddg,k,dk,ddk,dp,ddp,w,dw,ddw,V,Vp,Vpp=jets(seed)
    gi=inv(g); iu=neg(mul(mul(gi,k),gi))
    # Direct first coefficients of inverse, its spatial derivative, connection.
    ims=[gi,iu]; gs=[g,k]; dgs=[dg,dk]; ddgs=[ddg,ddk]
    dis=[]
    for n in (0,1):
        dis.append([[[ -sum(ims[j][a][c]*dgs[l][e][c][d]*ims[n-j-l][d][b]
                 for j in range(n+1) for l in range(n-j+1) for c,d in P2)
                 for b in I] for a in I] for e in I])
    con=[]; dcon=[]
    for n in (0,1):
        con.append({(a,b,c):sum(ims[j][a][d]*(dgs[n-j][b][c][d]+dgs[n-j][c][b][d]-dgs[n-j][d][b][c])
                    for j in range(n+1) for d in I)/2 for a,b,c in P3})
        dcon.append({(e,a,b,c):sum(
             dis[j][e][a][d]*(dgs[n-j][b][c][d]+dgs[n-j][c][b][d]-dgs[n-j][d][b][c])+
             ims[j][a][d]*(ddgs[n-j][e][b][c][d]+ddgs[n-j][e][c][b][d]-ddgs[n-j][e][d][b][c])
             for j in range(n+1) for d in I)/2 for e,a,b,c in P4})
    G=con[0]; DG=dcon[0]
    direct={(a,b):sum(dcon[1][c,c,a,b]-dcon[1][b,c,a,c] for c in I)+
            sum(con[j][c,c,d]*con[1-j][d,a,b]-con[j][c,b,d]*con[1-j][d,a,c]
                for j in (0,1) for c,d in P2) for a,b in P2}
    # Recipe side: both covariant derivative slots, including the first slot.
    H={(a,b,c):dk[a][b][c]-sum(G[d,a,b]*k[d][c]+G[d,a,c]*k[b][d] for d in I) for a,b,c in P3}
    partialH={(e,a,b,c):ddk[e][a][b][c]-sum(
        DG[e,d,a,b]*k[d][c]+G[d,a,b]*dk[e][d][c]+
        DG[e,d,a,c]*k[b][d]+G[d,a,c]*dk[e][b][d] for d in I) for e,a,b,c in P4}
    slot=0 if mutation=='derivative_slot' else 1
    HH={(e,a,b,c):partialH[e,a,b,c]-sum(slot*G[d,e,a]*H[d,b,c]+
        G[d,e,b]*H[a,d,c]+G[d,e,c]*H[a,b,d] for d in I) for e,a,b,c in P4}
    box={(a,b):sum(gi[c][d]*HH[c,d,a,b] for c,d in P2) for a,b in P2}
    DC={a:sum(gi[b][c]*(H[b,a,c]-H[a,b,c]/2) for b,c in P2) for a in I}
    NDC={(e,a):sum(gi[b][c]*(HH[e,b,a,c]-HH[e,a,b,c]/2) for b,c in P2) for e,a in P2}
    R={(a,b,c,d):DG[c,a,d,b]-DG[d,a,c,b]+sum(G[a,c,e]*G[e,d,b]-G[a,d,e]*G[e,c,b] for e in I)
       for a,b,c,d in P4}
    ric={(a,b):sum(R[c,a,c,b] for c in I) for a,b in P2}
    ku=mul(mul(gi,k),gi)
    lr={(a,b):sum((ric[a,c]*gi[c][d]*k[d][b]+ric[b,c]*gi[c][d]*k[d][a])/2 for c,d in P2)-
          sum(g[a][e]*R[e,c,b,d]*ku[c][d] for c,d,e in P3) for a,b in P2}
    factor=-1 if mutation=='curvature_sign' else 1
    sym=Q(1) if mutation=='symmetrization' else Q(1,2)
    need(any(lr.values()) and any(DC.values()) and any(G.values()),'NONTRIVIAL_CURVED_WITNESS')
    count=0
    for a,b in P2:
        recipe=-box[a,b]/2+sym*(NDC[a,b]+NDC[b,a])+factor*lr[a,b]
        need(direct[a,b]==recipe,'CURVED_RICCI_LINEARIZATION');count+=1
        reduced=direct[a,b]-(NDC[a,b]+NDC[b,a])/2
        need(reduced==-box[a,b]/2+lr[a,b],'REDUCED_RICCI_LINEARIZATION');count+=1
        Edirect=reduced-dp[a]*dw[b]-dw[a]*dp[b]-Q(2,3)*(Vp*w*g[a][b]+V*k[a][b])
        Erecipe=-box[a,b]/2+lr[a,b]-dp[a]*dw[b]-dw[a]*dp[b]-Q(2,3)*(Vp*w*g[a][b]+V*k[a][b])
        need(Edirect==Erecipe,'REDUCED_EINSTEIN_LINEARIZATION');count+=1
    hphi={(a,b):ddp[a][b]-sum(G[c,a,b]*dp[c] for c in I) for a,b in P2}
    hw={(a,b):ddw[a][b]-sum(G[c,a,b]*dw[c] for c in I) for a,b in P2}
    Sdirect=sum(iu[a][b]*hphi[a,b]+gi[a][b]*(hw[a,b]-sum(con[1][c,a,b]*dp[c] for c in I)) for a,b in P2)-Vpp*w
    scalarfactor=0 if mutation=='scalar_DC' else 1
    Srecipe=sum(gi[a][b]*hw[a,b]-ku[a][b]*hphi[a,b] for a,b in P2)-scalarfactor*sum(gi[a][b]*DC[a]*dp[b] for a,b in P2)-Vpp*w
    need(Sdirect==Srecipe,'SCALAR_LINEARIZATION');count+=1
    # Independent direct relative-vector contraction checks lowered dC.
    for a in I:
        need(DC[a]==sum(g[a][d]*gi[b][c]*con[1][d,b,c] for b,c,d in P3),'RELATIVE_VECTOR_LINEARIZATION');count+=1
    return {'seed':seed,'exact_comparisons':count,'nonnormal_curved_Lorentzian_metric':True}


def angular_audit(mutation):
    irreps=[(0,0),(1,1),(2,2),(2,0),(0,2)]
    if mutation=='omit_TT': irreps=irreps[:3]
    rows=[]; spin_dims={}
    for s,multiplicity in [(0,5),(1,2),(2,1)]:
        selected=[]
        for l,r in irreps:
            if abs(l-r)<=s<=l+r:
                d=(2*l+1)*(2*r+1)
                eig=2*(l*(l+1)+r*(r+1))-s*(s+1)
                selected.append({'irrep':[l,r],'dimension':d,'rough_eigenvalue':eig})
        spin_dims[s]=sum(item['dimension'] for item in selected)
        rows.append({'spin':s,'fibre_multiplicity':multiplicity,'sectors':selected})
    total=sum(row['fibre_multiplicity']*sum(q['dimension'] for q in row['sectors']) for row in rows)
    if mutation=='scalar_component_cutoff': total=16*35
    need(spin_dims=={0:35,1:34,2:44} and total==287,'ANGULAR_MULTIPLICITY')
    need(max(q['rough_eigenvalue'] for row in rows for q in row['sectors'])==24,'ANGULAR_SPECTRUM')
    evaluation=16 if mutation=='evaluation_factor' else 17
    need(evaluation**2>total,'EVALUATION_TRACE')
    return {'branching':rows,'dimension':total,'Dirichlet':spin_dims[0]+spin_dims[1],
            'Robin':total-spin_dims[0]-spin_dims[1],'rough_operator_cap':24,'evaluation_cap':evaluation}


def rational_audit(base,mutation):
    pins=json.loads((base/'INPUT_PINS.json').read_text());data={}
    for key,item in pins.items():
        raw=(base/item['path']).read_bytes()
        need(sha256(raw).hexdigest()==item['sha256'],'HASH_'+key)
        data[key]=json.loads(raw)
    bg=data['background']; bulk=data['bulk']; shell=data['shell']; producer=data['producer']
    need(bg['domain']['ell']=='1/1000000' and bg['domain']['scaled_x']==[0,3] and bg['domain']['scaled_r']==[0,8],'WHOLE_COLLAR_DOMAIN')
    ell=Q(1,10**6); rho=Q(131); alpha=Q(14,125); p=Q(7,50000); ht=3*ell; h=ht/129
    for key,cap in [('rho',rho),('a',alpha),('p',p),('a_y',Q(1,10000))]:
        need(Q(bg['static_actual_derivative_caps'][key][0])<cap,'BACKGROUND_'+key)
    lo,hi=map(Q,bg['first_exit']['rho']['outer']['exact'])
    need(lo>=129 and hi<=rho,'RHO_LOWER')
    lo,hi=map(Q,bg['first_exit']['eta']['outer']['exact'])
    need(-Q(9,10**5)<=lo<=hi<=-Q(8,10**5),'SIGNED_PHI_TUBE')
    need(0<1+lo<=1+hi<1,'POSITIVE_SUBUNIT_PHI')
    # Differentiate polynomial coefficient lists, then apply a positive norm.
    poly=[-Q(2,27),Q(0),Q(14,9),Q(50,27),-Q(1,6),-Q(4,9),-Q(2,27)]
    pot=[]
    for _ in range(3):
        pot.append(sum(abs(v)*Q(9,10**5)**j for j,v in enumerate(poly)))
        poly=[(j+1)*poly[j+1] for j in range(len(poly)-1)]
    need(pot[0]<1 and pot[1]<Q(1,1000) and pot[2]<4,'POTENTIAL_CAPS')
    py=pot[1]+4*alpha*p
    need(py<Q(1,500) and py+2*alpha*p<1,'SCALAR_HESSIAN')
    curvature=alpha*alpha+Q(1,10000)+Q(1,129**2)
    need(curvature<Q(13,1000) and 29*curvature<1,'CURVATURE_OPERATOR')
    # Sum independently identified blocks, including both off-diagonal systems.
    blocks_A=[3*ht,4*rho*alpha,4*rho*p,3*rho*p]
    blocks_B=[3*rho*alpha,4*rho*p,3*rho*p]
    z=(2*alpha)**2+3*h*2*alpha+3*(2*(alpha+h))**2
    blocks_D={'angular':Q(24),'angular_connection':40*rho*(alpha+h),
              'connection_square':rho*rho*z,'curvature':2*rho*rho,
              'metric_potential':Q(4,3)*rho*rho,'metric_scalar_potential':4*rho*rho,
              'scalar_potential':4*rho*rho,'scalar_Hessian':rho*rho,
              'scalar_DC_connection':rho*rho*p*(8*alpha+6*h),'angular_derivative_couplings':30*rho*p}
    A=ell*sum(blocks_A);B=ell*sum(blocks_B);D=ell*ell*sum(blocks_D.values())
    M=ell*rho*(10*alpha+4*p+p/2+2)
    coefs={'A':A,'B':B,'D':D,'M':M}; caps={'A':Q(1,10000),'B':Q(1,10000),'D':ell,'M':Q(1,2000)}
    for key in coefs:
        need(coefs[key]<caps[key],'OPERATOR_'+key)
        need(coefs[key]==Q(producer['coefficient_operator_upper'][key]),'PRODUCER_OPERATOR_'+key)
    trace=1 if mutation=='missing_boundary_weight' else 3
    need(trace==3,'BOUNDARY_TRACE_WEIGHT')
    q=3*(caps['A']+caps['B']+3*caps['D'])+trace*caps['M']
    need(q==Q(2109,10**6),'CONTRACTION')
    metric=Q(bulk['actual_source_bounds']['QEred']['strict_simple_upper'])
    scalar=Q(bulk['actual_source_bounds']['Qphi']['strict_simple_upper'])
    gauge=Q(bulk['actual_source_bounds']['QC']['strict_simple_upper'])
    source=shell['bounds']['one_copy_outward_residual_hessian']
    qr=max(Q(v['exact_upper']) for v in source['metric']['value'].values())
    qp=Q(source['scalar']['value']['exact_upper'])
    gn=8*qr+2*gauge
    need(metric<10**21 and scalar<4*10**15 and gauge<6000,'ACTUAL_BULK_JOIN')
    need(qr<15*10**8 and qp<44*10**8 and gn<13*10**9,'ACTUAL_SHELL_JOIN')
    f=ell**2*rho**2*(10*10**21+4*10**15)
    g=ell*rho*(8*15*10**8+44*10**8+13*10**9)
    need(f<2*10**14 and g<4*10**6,'SOURCE_HILBERT_JOIN')
    init=data['initial']['convenient_strict_caps']
    k=Q(init['metric2_components']);dk=Q(init['metric2_full_initial_covariant_derivative_components'])
    w=Q(init['scalar2']);dw=Q(init['scalar2_full_initial_covariant_derivative_components'])
    u0=5*k+w;ur=ell*rho*(5*dk+dw);ux=ell*rho*(5*dk+dw+2*alpha*5*k)
    # For positive majorants the subtracted frame connection is added in norm.
    checkux=ell*rho*(5*dk+dw+(-1 if mutation=='wrong_initial_frame_sign' else 1)*2*alpha*5*k)
    need(checkux==ux,'INITIAL_FRAME_CONVERSION')
    free=max(u0/3+ux,ur+ux)
    need(u0<6*10**33 and ur<5*10**35 and ux<5*10**35 and free<10**36,'INITIAL_HILBERT_JOIN')
    response=(10**36+3*2*10**14+4*10**6)/(1-q)
    need(response==Q(producer['response']['X_exact_upper']) and response<Q('1.02e36'),'RESPONSE_CAP')
    need(51*response<6*10**37 and 17*response<2*10**37,'POINTWISE_CAP')
    return {'input_sha256':{key:item['sha256'] for key,item in pins.items()},
            'potential_caps':list(map(str,pot)),'curvature_component_cap':str(curvature),
            'scalar_Hessian_cap':str(py+2*alpha*p),'sigma_second_half_cap':'2',
            'sigma_reason':'authenticated 0<phi0<1, sigma_second/2=2phi0',
            'zeroth_blocks_before_ell_squared':{k:str(v) for k,v in blocks_D.items()},
            'operator_majorants':{k:str(v) for k,v in coefs.items()},'contraction':str(q),
            'resolvent':str(1/(1-q)),'initial_value_H':str(u0),'initial_radial_H':str(ur),
            'initial_time_H':str(ux),'free_X':str(free),'F_H':str(f),'G_H':str(g),
            'response_X':str(response),'pointwise_value_cap':'6e37','pointwise_scaled_derivative_cap':'2e37'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inject',choices=sorted(CONTROLS));args=parser.parse_args()
    need(not args.output.exists(),'OUTPUT_EXISTS')
    base=Path(__file__).resolve().parent
    curved=[curved_audit(seed,args.inject) for seed in (161640,161641)]
    angular=angular_audit(args.inject); quantitative=rational_audit(base,args.inject)
    need(args.inject is None,'MUTATION_NOT_REJECTED')
    result={'status':'PASS_INDEPENDENT_EVOLUTION_OPERATOR_AND_NORMS','curved_finite_jet_checks':curved,
            'angular':angular,'quantitative':quantitative,
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'All-point operator and norm arguments are in EVOLUTION_REVIEW.md; finite jets test algebra, not PDE existence. Qualitative regularity/subsidiary acceptance is a separate review.'}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'])


if __name__=='__main__': main()
