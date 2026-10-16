#!/usr/bin/env python3
"""Exact channelwise continuous C3 raw-gauge bounds from frozen coefficients.

Only standard-library Fraction arithmetic. No old scientific module is
imported and no differentiation of a value/error bound is performed.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
PINS={
 'GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json':'4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53',
 'COEFFICIENT_VERIFICATION_RECEIPT.json':'30b970991fe5e189892cbf007d9063ad78deb367f09cb8a5a4a6ce9988c55197',
 'MIXED_WAVE_REGULARITY.md':'00f654647d3cd5578e3f9afcb00dc3ec81739036f40adbe43f53e84b87506e98',
}
X=F(3)
ELL=F(1,10**6)
THETA=X/2


def need(ok,message):
    if not ok:raise ValueError(message)


def vadd(*vectors):return [sum(v[i] for v in vectors) for i in range(3)]
def vscale(c,v):return [c*x for x in v]
def mv(a,v):return [sum((a[i][j]*v[j] for j in range(3)),F(0)) for i in range(3)]
def madd(*matrices):return [[sum(a[i][j] for a in matrices) for j in range(3)] for i in range(3)]
def mscale(c,a):return [[c*x for x in row] for row in a]
def ident():return [[F(int(i==j)) for j in range(3)] for i in range(3)]


def inverse3(a):
    aug=[list(row)+ident()[i] for i,row in enumerate(a)]
    for k in range(3):
        need(aug[k][k]!=0,'exact inverse pivot')
        pivot=aug[k][k];aug[k]=[x/pivot for x in aug[k]]
        for i in range(3):
            if i!=k:
                q=aug[i][k];aug[i]=[x-q*y for x,y in zip(aug[i],aug[k])]
    out=[row[3:] for row in aug]
    for i in range(3):
        for j in range(3):
            need(sum(a[i][k]*out[k][j] for k in range(3))==int(i==j),'exact inverse product')
    return out


def pclean(p):return {k:F(v) for k,v in p.items() if v}
def padd(*ps):
    out={}
    for p in ps:
        for k,c in p.items():out[k]=out.get(k,F(0))+c
    return pclean(out)
def pscale(a,p):return pclean({k:a*v for k,v in p.items()})
def pmul(a,b):
    out={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():out[i+k,j+l]=out.get((i+k,j+l),F(0))+x*y
    return pclean(out)
def pd(p,ix=0,ir=0):
    out=p
    for axis,times in ((0,ix),(1,ir)):
        for _ in range(times):
            nxt={}
            for key,c in out.items():
                if key[axis]:
                    kk=list(key);kk[axis]-=1;nxt[tuple(kk)]=c*key[axis]
            out=nxt
    return pclean(out)
def pzero(p):return pclean({(0,j):c for (i,j),c in p.items() if i==0})


def differential_identity_checks(injection):
    # An exact polynomial with zero initial displacement and velocity.
    # The source is defined from its full variable-coefficient wave residual.
    # This checks product-rule and recovery identities, not numerical samples.
    u={(2,0):F(1,2),(2,1):F(1,2),(2,3):F(1,2),
       (3,0):F(1,3),(3,2):F(1,6),(4,1):F(1,24),(5,0):F(1,120)}
    A={(0,0):F(1),(1,0):F(1),(0,1):F(2),(1,1):F(1)}
    B={(0,0):F(2),(1,0):F(-1),(0,2):F(1)}
    D={(0,0):F(3),(2,0):F(1),(0,1):F(-1)}
    ux,ur=pd(u,1),pd(u,0,1)
    forcing=padd(pmul(A,ux),pmul(B,ur),pmul(D,u))
    source=padd(pd(u,2),pscale(-1,pd(u,0,2)),pscale(-1,forcing))
    v=ux;w=pd(u,2)
    homogeneous=lambda z:padd(pd(z,2),pscale(-1,pd(z,0,2)),
                              pscale(-1,pmul(A,pd(z,1))),pscale(-1,pmul(B,pd(z,0,1))),
                              pscale(-1,pmul(D,z)))
    r1=padd(pmul(pd(A,1),ux),pmul(pd(B,1),ur),pmul(pd(D,1),u),pd(source,1))
    need(not padd(homogeneous(v),pscale(-1,r1)),'first-time commutator identity')
    factor=1 if injection=='drop_second_commutator_two' else 2
    r2=padd(pscale(factor,pmul(pd(A,1),pd(v,1))),
            pscale(factor,pmul(pd(B,1),pd(v,0,1))),pscale(factor,pmul(pd(D,1),v)),
            pmul(pd(A,2),ux),pmul(pd(B,2),ur),pmul(pd(D,2),u),pd(source,2))
    need(not padd(homogeneous(w),pscale(-1,r2)),'second-time commutator factor two')
    s0=pzero(source)
    init_A=pzero(pd(A,1)) if injection=='use_Ax0_initial_velocity' else pzero(A)
    expected=padd(pmul(init_A,s0),pzero(pd(source,1)))
    need(not padd(pzero(pd(u,3)),pscale(-1,expected)),'second-time initial velocity uses A(0), not A_x(0)')
    need(pzero(w)==s0 and not pzero(u) and not pzero(v),'differentiated initial displacement identities')
    Fr=padd(pmul(pd(A,0,1),ux),pmul(A,pd(u,1,1)),
             pmul(pd(B,0,1),ur),pmul(B,pd(u,0,2)),
             pmul(pd(D,0,1),u),pmul(D,ur))
    if injection!='omit_source_radial_derivative':Fr=padd(Fr,pd(source,0,1))
    Fx=pd(padd(forcing,source),1)
    need(not padd(pd(u,1,2),pscale(-1,pd(u,3)),Fx),'mixed third radial recovery')
    need(not padd(pd(u,0,3),pscale(-1,pd(u,2,1)),Fr),'third radial recovery retains source radial derivative')
    return {'first_time_equation':'PASS_EXACT_POLYNOMIAL_IDENTITY',
            'second_time_equation':'PASS_EXACT_POLYNOMIAL_IDENTITY',
            'initial_data':'PASS_EXACT_POLYNOMIAL_IDENTITIES',
            'third_radial_recoveries':'PASS_EXACT_POLYNOMIAL_IDENTITIES'}


def run(injection='none'):
    for name,pin in PINS.items():
        need(sha256((ROOT/'inputs'/name).read_bytes()).hexdigest()==pin,'authenticated input: '+name)
    coef=json.loads((ROOT/'inputs/GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json').read_text())
    receipt=json.loads((ROOT/'inputs/COEFFICIENT_VERIFICATION_RECEIPT.json').read_text())
    need(receipt['report_sha256']==PINS['GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json'],
         'coefficient receipt matches exact report')
    need(coef['status']=='PASS_WHOLE_COLLAR_GAUGE_COEFFICIENT_AND_SOURCE_CERTIFICATE','coefficient status')
    need(F(coef['domain']['ell'])==ELL and coef['domain']['scaled_x']==[0,3], 'scaled coordinate domain')
    checks=differential_identity_checks(injection)
    def matrix(name,i=0,j=0):
        return [[F(entry[f'{i},{j}']['abs_upper']) for entry in row]
                for row in coef['matrix_mixed_derivative_intervals'][name]]
    a,b,d=[matrix(n) for n in ('A','B','D')]
    ax,bx,dx=[matrix(n,1,0) for n in ('A','B','D')]
    axx,bxx,dxx=[matrix(n,2,0) for n in ('A','B','D')]
    ar,br,dr=[matrix(n,0,1) for n in ('A','B','D')]
    theta=F(1) if injection=='drop_field_value_factor' else THETA
    L0=madd(a,b,mscale(theta,d));M=mscale(X,L0)
    need(M==[[F(x) for x in row] for row in coef['weighted_C1_matrix']],
         'weighted field value factor is X/2=3/2')
    need(max(map(sum,M))<1,'same matrix contraction')
    inv=inverse3(madd(ident(),mscale(-1,M)))
    need(all(v>=0 for row in inv for v in row),'nonnegative same resolvent')
    source=lambda i,j:[F(coef['source_strict_absolute_mixed_derivative_caps'][f'{i},{j}']),F(0),F(0)]
    s0,sx,sr,sxx=source(0,0),source(1,0),source(0,1),source(2,0)
    n0=mv(inv,vscale(X,s0))
    need(n0==list(map(F,coef['conditional_weighted_C1_channel_bounds'])),'recomputed original channel bound')
    Lx=madd(ax,bx,mscale(theta,dx))
    Lxx=madd(axx,bxx,mscale(theta,dxx))
    r1=vadd(mv(Lx,n0),sx)
    free_v=vscale(2,s0)
    n1=mv(inv,vadd(free_v,vscale(X,r1)))
    # Conservative global majorants for the initial data are legitimate.
    # A(0)S0 is bounded by the global elementwise A bound, never by A_x.
    h0=s0;h0r=sr;h1=vadd(mv(a,s0),sx)
    free_w=[max(2*h0[i]/X+2*h1[i],h0r[i]+h1[i]) for i in range(3)]
    r2=vadd(vscale(2,mv(Lx,n1)),mv(Lxx,n0),sxx)
    n2=mv(inv,vadd(free_w,vscale(X,r2)))
    for n,rhs in [(n0,vscale(X,s0)),(n1,vadd(free_v,vscale(X,r1))),
                  (n2,vadd(free_w,vscale(X,r2)))]:
        need(vadd(n,vscale(-1,mv(M,n)))==rhs,'exact nonnegative resolvent majorant equation')
    bounds={(0,0):vscale(theta,n0),(1,0):n0,(0,1):n0,
            (2,0):[min(n1[i],theta*n2[i]) for i in range(3)],(1,1):n1,
            (3,0):n2,(2,1):n2}
    bounds[0,2]=vadd(bounds[2,0],mv(a,bounds[1,0]),mv(b,bounds[0,1]),
                     mv(d,bounds[0,0]),s0)
    fx=vadd(mv(ax,bounds[1,0]),mv(a,bounds[2,0]),
            mv(bx,bounds[0,1]),mv(b,bounds[1,1]),
            mv(dx,bounds[0,0]),mv(d,bounds[1,0]),sx)
    fr=vadd(mv(ar,bounds[1,0]),mv(a,bounds[1,1]),
            mv(br,bounds[0,1]),mv(b,bounds[0,2]),
            mv(dr,bounds[0,0]),mv(d,bounds[0,1]),sr)
    bounds[1,2]=vadd(bounds[3,0],fx)
    bounds[0,3]=vadd(bounds[2,1],fr)
    unscaled={}
    for key,values in bounds.items():
        order=sum(key)
        factor=ELL**order if injection=='wrong_physical_scaling' else ELL**(-order)
        unscaled[key]=vscale(factor,values)
        need(vscale(ELL**order,unscaled[key])==values,'unscaled derivatives require ell to negative total order')
    channels=('f','b','d')
    ordercaps={name:[max(v[channel] for key,v in bounds.items() if sum(key)==order)
                     for order in range(4)] for channel,name in enumerate(channels)}
    # Display caps use exact decimal ceiling with six significant figures.
    def round_up_decimal(value):
        need(value>0,'positive readable cap')
        exponent=0;scaled=value
        while scaled>=10:scaled/=10;exponent+=1
        while scaled<1:scaled*=10;exponent-=1
        step=F(10)**(exponent-5)
        return ((value/step).__floor__()+1)*step
    readable={name:list(map(round_up_decimal,values)) for name,values in ordercaps.items()}
    need(all(x<y for name in channels for x,y in zip(ordercaps[name],readable[name])),
         'strict readable total-order caps')
    records=lambda data:{f'{i},{j}':list(map(str,v)) for (i,j),v in sorted(data.items(),key=lambda z:(sum(z[0]),z[0]))}
    all_order=[max(ordercaps[name][order] for name in channels) for order in range(4)]
    all_readable=list(map(round_up_decimal,all_order))
    return {
      'status':'PASS_CONTINUOUS_RAW_GAUGE_ALL_MIXED_C3_CAPS',
      'source_pins':PINS,'inherited_actual_source_pins':coef['input_sha256'],
      'domain':{'ell':str(ELL),'scaled':'0<=x<=3, r>=0, x+r<=8',
                'unscaled':'0<=t<=3ell, z>=0, t+z<=8ell','boundary':'physical shell r=0; one-sided classical derivatives'},
      'channels':list(channels),
      'variables':'raw contravariant amplitudes xi^y=fY,xi^t=bY,xi^i=d D^i_kappa Y; no rho rescaling',
      'norm':'N_i(V)=max((2/3)||V_i||,||partial_x V_i||,||partial_r V_i||)',
      'weighted_C1_matrix':[[str(x) for x in row] for row in M],
      'same_resolvent':[[str(x) for x in row] for row in inv],
      'weighted_C1_bounds':{'U':list(map(str,n0)),'U_x':list(map(str,n1)),'U_xx':list(map(str,n2))},
      'initial_majorants':{'v_initial_displacement':['0']*3,'v_initial_velocity':list(map(str,s0)),
                           'v_free_weighted_C1':list(map(str,free_v)),
                           'w_initial_displacement':list(map(str,h0)),
                           'w_initial_displacement_r':list(map(str,h0r)),
                           'w_initial_velocity_A_at_zero_S0_plus_Sx0':list(map(str,h1)),
                           'w_free_weighted_C1':list(map(str,free_w))},
      'commutator_forcing_majorants':{'for_Ux':list(map(str,r1)),'for_Uxx':list(map(str,r2))},
      'scaled_mixed_partial_caps':records(bounds),
      'unscaled_coordinate_mixed_partial_caps':records(unscaled),
      'scaled_total_order_caps_by_channel':{k:list(map(str,v)) for k,v in ordercaps.items()},
      'strict_readable_scaled_total_order_caps_by_channel':{k:list(map(str,v)) for k,v in readable.items()},
      'unscaled_total_order_caps_by_channel':{k:[str(value*ELL**(-order)) for order,value in enumerate(v)] for k,v in ordercaps.items()},
      'scaled_total_order_maximum_all_channels':list(map(str,all_order)),
      'strict_readable_scaled_total_order_maximum_all_channels':list(map(str,all_readable)),
      'strict_readable_unscaled_total_order_maximum_all_channels':
          [str(value*ELL**(-order)) for order,value in enumerate(all_readable)],
      'recovery_majorants':{'F_x':list(map(str,fx)),'F_r':list(map(str,fr))},
      'initial_velocity_caution':'Uxx initial velocity is A(0,r)S(0,r)+S_x(0,r), not A_x(0,r)S(0,r)+S_x(0,r).',
      'regularity_basis':'Pinned compatible C2-coefficient/source mixed-wave theorem constructs U,U_x,U_xx with C1 convergence; wave equation recovers continuous radial derivatives. Source flatness on x+r<1 supplies the corner compatibility.',
      'identity_checks':checks,
      'arithmetic':'exact fractions.Fraction only; no floating-point computation or scientific module imports',
      'scope':'Pointwise continuous raw-amplitude coordinate derivative caps through total order3, conditional on retained actual branch/master and the pinned gauge-wave construction. Not frame/covariant/angular derivative norms, energy norms, nonlinear amplitude/remainder bounds, or a full second-order contact solution.'
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--inject-error',choices=['none','use_Ax0_initial_velocity','drop_second_commutator_two',
                         'drop_field_value_factor','omit_source_radial_derivative','wrong_physical_scaling'],default='none')
    args=parser.parse_args();destination=args.output.resolve()
    need(destination!=ROOT and ROOT not in destination.parents,'output must be outside higher-norm input tree')
    need(not destination.exists(),'refuse to overwrite output')
    report=run(args.inject_error)
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'strict_readable_scaled_total_order_caps_by_channel':
                     report['strict_readable_scaled_total_order_caps_by_channel'],'output':str(destination)},indent=2))


if __name__=='__main__':main()
