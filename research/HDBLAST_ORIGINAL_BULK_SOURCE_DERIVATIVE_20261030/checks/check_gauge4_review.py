#!/usr/bin/env python3
"""Independent ordinary-derivative algebra and actual gauge fourth-jet ledger."""
from fractions import Fraction as Q
from pathlib import Path
from math import comb
import argparse, hashlib, json

checks=[]
controls=[]
def need(ok, name):
    if not ok: raise RuntimeError(name)
    checks.append(name)
def reject(ok, name):
    need(not ok, 'reject '+name)
    controls.append(name)
def load(path, pin):
    b=path.read_bytes()
    need(hashlib.sha256(b).hexdigest()==pin,'authenticate input '+pin[:16])
    return json.loads(b)
def encode(z):
    if isinstance(z,Q):return str(z)
    if isinstance(z,dict):return {k:encode(v) for k,v in z.items()}
    if isinstance(z,(list,tuple)):return [encode(v) for v in z]
    return z

# Polynomials in H with exact ordinary time derivative dH/dt=1-H^2.
def padd(*terms):
    d={}
    for c,p in terms:
        for k,v in p.items():d[k]=d.get(k,Q(0))+c*v
    return {k:v for k,v in d.items() if v}
def pmul(p,q):
    d={}
    for k,v in p.items():
        for j,w in q.items():d[k+j]=d.get(k+j,Q(0))+v*w
    return {k:v for k,v in d.items() if v}
def pd(p):return {k-1:k*v for k,v in p.items() if k}
def dt(p):return pmul({0:Q(1),2:Q(-1)},pd(p))
def dq(p):return padd((1,dt(p)),(-Q(3,2),pmul({1:Q(1)},p)))
H={1:Q(1)};J={0:Q(1),2:Q(-1)}
H3=dt(dt(dt(H)));J3=dt(dt(dt(J)));HJ3=dt(dt(dt(pmul(H,J))))
need(H3=={0:Q(-2),2:Q(8),4:Q(-6)},'ordinary third derivative of tanh')
need(J3=={1:Q(16),3:Q(-40),5:Q(24)},'ordinary third derivative of sech squared')
need(HJ3=={0:Q(-8),2:Q(68),4:Q(-120),6:Q(60)},'ordinary third derivative of HJ')
qpoly={0:Q(1)}
qpolys=[qpoly]
for _ in range(3):qpolys.append(dq(qpolys[-1]))
need(qpolys[3]=={1:Q(39,4),3:Q(-105,8)},'ordinary third derivative of cosh to minus three halves')
Hsq3=dt(dt(dt(pmul(H,H))))
need(padd((1,Hsq3),(1,J3))=={},'H squared third derivative equals minus J third')
reject(padd((Q(1,6),qpolys[3]))==qpolys[3],'normalized third Taylor coefficient used as derivative')

# Noncommutative matrix/vector polynomial algebra. Words preserve matrix
# multiplication order; the final factor is a source vector derivative.
def fa(*terms):
    out={}
    for c,p in terms:
        for word,v in p.items():out[word]=out.get(word,Q(0))+c*v
    return {w:c for w,c in out.items() if c}
def fm(p,q):
    out={}
    for w,c in p.items():
        for v,d in q.items():out[w+v]=out.get(w+v,Q(0))+c*d
    return {w:c for w,c in out.items() if c}
def atom(name,x=0,r=0):return {((name,x,r),):Q(1)}
def fd(p,axis):
    out={}
    for word,c in p.items():
        for i,factor in enumerate(word):
            f=list(factor);f[axis]+=1
            w=word[:i]+(tuple(f),)+word[i+1:]
            out[w]=out.get(w,Q(0))+c
    return out
cache={}
def initial(a,b):
    if a<2:return {}
    if (a,b) in cache:return cache[a,b]
    # Differentiate the original PDE, not an asserted fourth-jet formula.
    x=a-2;r=b
    result=fa((1,initial(x,r+2)),(1,atom('s',x,r)))
    for k in range(x+1):
        for j in range(r+1):
            c=Q(comb(x,k)*comb(r,j))
            result=fa((1,result),
                (c,fm(atom('P',x-k,r-j),initial(k+1,j))),
                (c,fm(atom('Q',x-k,r-j),initial(k,j+1))),
                (c,fm(atom('R',x-k,r-j),initial(k,j))))
    cache[a,b]=result
    return result
s=atom('s');P=atom('P')
d0=fa((1,fm(P,s)),(1,atom('s',1)))
v0=fa((1,atom('s',0,2)),(1,fm(P,d0)),(2,fm(atom('P',1),s)),
      (1,fm(atom('Q'),atom('s',0,1))),(1,fm(atom('R'),s)),(1,atom('s',2)))
d0r=fd(d0,2)
need(initial(2,0)==s,'zero Cauchy data retain initial acceleration source')
need(initial(3,0)==d0,'independent third time initial displacement')
need(initial(4,0)==v0,'noncommutative fourth time initial velocity')
need(initial(3,1)==d0r,'independent radial derivative of third initial displacement')
reject(fa((1,v0),(-1,fm(P,fm(P,s))))==v0,'drop P squared times source')
reject(fa((1,v0),(-1,fm(atom('P',1),s)))==v0,'drop factor two in P_x times source')
reject(fa((1,d0r),(-1,fm(atom('P',0,1),s)))==d0r,'drop P_r times source')
full=fa((1,fm(atom('P'),atom('v',1))),
        (1,fm(atom('Q'),atom('v',0,1))),
        (1,fm(atom('R'),atom('v'))),(1,atom('s')))
third=fd(fd(fd(full,1),1),1)
principal=fa((1,fm(atom('P'),atom('v',4))),
             (1,fm(atom('Q'),atom('v',3,1))),
             (1,fm(atom('R'),atom('v',3))))
comm3=fa((3,fm(atom('P',1),atom('v',3))),
         (3,fm(atom('Q',1),atom('v',2,1))),
         (3,fm(atom('R',1),atom('v',2))),
         (3,fm(atom('P',2),atom('v',2))),
         (3,fm(atom('Q',2),atom('v',1,1))),
         (3,fm(atom('R',2),atom('v',1))),
         (1,fm(atom('P',3),atom('v',1))),
         (1,fm(atom('Q',3),atom('v',0,1))),
         (1,fm(atom('R',3),atom('v'))),(1,atom('s',3)))
need(fa((1,third),(-1,principal))==comm3,'independent noncommutative third time commutator')
reject(fa((1,comm3),(-2,fm(atom('P',1),atom('v',3))))==comm3,
       'drop binomial factor3 from highest field commutator')

ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--master4-report',type=Path)
ap.add_argument('--producer-report',type=Path,required=True)
a=ap.parse_args()
raw=a.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs'
pincoeff='4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53'
pinmaster='d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e'
coeff=load(raw/'GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json',pincoeff)
master=load(a.master4_report or a.repo_root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json',pinmaster)
need(coeff['status']=='PASS_WHOLE_COLLAR_GAUGE_COEFFICIENT_AND_SOURCE_CERTIFICATE','actual coefficient certificate status')
need(master['status']=='PASS_CONDITIONAL_ACTUAL_MASTER_CONTINUOUS_FOURTH_JETS','actual master fourth certificate status')
need(coeff['domain']['ell']=='1/1000000' and coeff['domain']['scaled_x']==[0,3],'same scaling and triangle time')
ell=Q(1,10**6);X=Q(3);theta=X/2;h=X*ell;lam=Q(8)
static=coeff['static_actual_derivative_caps']
rho,aa,C,Cz=(Q(static[k][0]) for k in ('rho','a','C','C_z'))
Cs=6*ell*rho*C
C0=6*ell**2*rho*(Cz+4*rho*aa*C)
need(h*h<Q(1,3),'small strip permits tanh third derivative bound2')
capH3=Q(2)
capJ3=16*h+40*h**3+24*h**5
capHJ3=8+68*h*h+120*h**4+60*h**6
qder=[Q(1),Q(3,2)*h,Q(3,2)+Q(15,4)*h*h,Q(39,4)*h+Q(105,8)*h**3]
M=list(map(Q,coeff['master_scaled_order_strict_caps']))+[max(map(Q,master['exact']['all_scaled_order_four_mixed_caps'].values()))]
S3=sum(Q(comb(3,j))*ell**j*qder[j]*(Cs*M[4-j]+C0*M[3-j]) for j in range(4))
mat={k:{j:Q(v) for j,v in rows.items()} for k,rows in coeff['matrix_maximum_absolute_row_sum_by_mixed_derivative'].items()}
mat['A']['3,0']=5*ell**4*capH3;mat['B']['3,0']=Q(0)
mat['D']['3,0']=ell**5*max(lam*capJ3+6*rho*rho*aa*capH3,
    (3+lam)*capJ3+2*lam*capH3,2*aa*capJ3+2*capHJ3+abs(4-lam)*capJ3)
src={k:Q(v) for k,v in coeff['source_strict_absolute_mixed_derivative_caps'].items()}
src['3,0']=S3
L1=lambda k:mat['A'][k]+mat['B'][k]
L0=lambda k:mat['D'][k]
alpha=lambda k:L1(k)+theta*L0(k)
q=X*alpha('0,0')
need(0<q<Q(3,250),'same actual scalar contraction')
N0=X*src['0,0']/(1-q)
N1=(2*src['0,0']+X*(alpha('1,0')*N0+src['1,0']))/(1-q)
v2=mat['A']['0,0']*src['0,0']+src['1,0']
free2=max(2*src['0,0']/X+2*v2,src['0,1']+v2)
N2=(free2+X*(2*alpha('1,0')*N1+alpha('2,0')*N0+src['2,0']))/(1-q)
B2=N1+alpha('0,0')*N0+src['0,0']
Tx=L1('0,0')*N1+L0('0,0')*N0+alpha('1,0')*N0+src['1,0']
Tr=L1('0,0')*B2+L0('0,0')*N0+alpha('0,1')*N0+src['0,1']
B3=N2+max(Tx,Tr)
v3=src['0,2']+mat['A']['0,0']*v2+2*mat['A']['1,0']*src['0,0']+mat['B']['0,0']*src['0,1']+mat['D']['0,0']*src['0,0']+src['2,0']
dr3=mat['A']['0,1']*src['0,0']+mat['A']['0,0']*src['0,1']+src['1,1']
free3=max(2*v2/X+2*v3,dr3+v3)
N3=(free3+X*(3*alpha('1,0')*N2+3*alpha('2,0')*N1+alpha('3,0')*N0+S3))/(1-q)
Txx=L1('0,0')*N2+L0('0,0')*N1+2*L1('1,0')*N1+2*L0('1,0')*N0+alpha('2,0')*N0+src['2,0']
Txr=L1('0,0')*B3+L0('0,0')*B2+(L1('1,0')+L1('0,1'))*B2+(L0('1,0')+L0('0,1'))*N0+alpha('1,1')*N0+src['1,1']
Trr=L1('0,0')*B3+L0('0,0')*B2+2*L1('0,1')*B2+2*L0('0,1')*N0+alpha('0,2')*N0+src['0,2']
four={'4,0':N3,'3,1':N3,'2,2':N3+Txx,'1,3':N3+Txr,'0,4':N3+Txx+Trr}
need(all(v>0 for v in four.values()),'positive allfive recovered fourth jet caps')
reject(N3==(free3+X*(alpha('1,0')*N2+3*alpha('2,0')*N1+alpha('3,0')*N0+S3))/(1-q),'omit binomial3 in highest commutator')
reject(N3==free3+X*(3*alpha('1,0')*N2+3*alpha('2,0')*N1+alpha('3,0')*N0+S3),'omit same operator inverse')
reject(S3==sum(Q(comb(3,j))*ell**j*qder[j]*(Cs*M[3-j]+C0*M[3-j]) for j in range(4)),'wrong master order in source radial term')
result={'status':'PASS_INDEPENDENT_ACTUAL_GAUGE_FOURTH_REVIEW','input_sha256':{'coefficients':pincoeff,'master4':pinmaster},
 'domain':'0<=x<=3,r>=0,x+r<=8; physical side and one-sided traces',
 'source_constants':{'Cs':Cs,'C0':C0,'ordinary_temporal_q_derivative_caps':qder,'source_xxx':S3},
 'coefficient_third_time_caps':{k:mat[k]['3,0'] for k in ('A','B','D')},
 'q':q,'N0':N0,'N1':N1,'N2':N2,'N3':N3,'B2':B2,'B3':B3,
 'initial_third_time_caps':{'d0':v2,'v0':v3,'d0r':dr3,'free':free3},
 'differentiated_forcing_caps':{'Txx':Txx,'Txr':Txr,'Trr':Trr},
 'fourth_mixed_caps':four,'named_checks':checks,'wrong_formula_controls_rejected':controls,
 'scope':'conditional actual first-gauge coordinate fourth-jet ledger; no actual full bulk F_x, second-target C2 error, nonlinear or cosmological inference'}
producerpin='94579087d39d9c849e9264d40bd1d71d64892a668fe3e2e08aade9ca4296b873'
p=load(a.producer_report,producerpin)
need(p['status']=='PASS_CONDITIONAL_ACTUAL_FIRST_GAUGE_CONTINUOUS_FOURTH_JETS'
     and p['actual_source_unchanged'] is True,'actual producer scope and unchanged source')
need(Q(p['q0'])==q and Q(p['resolvent'])==1/(1-q),'independent same contraction and inverse')
need(list(map(Q,p['temporal_ordinary_derivative_caps']['cosh_minus_three_halves_orders_0_to_3']))==qder,
     'ordinary temporal q derivative caps agree independently')
need({k:Q(p['static_source_caps'][k]) for k in ('Cs','C0')}=={'Cs':Cs,'C0':C0},'independent actual static source caps')
for k in ('A','B','D'):
    need(Q(p['ordinary_matrix_row_caps'][k]['3,0'])==mat[k]['3,0'],'independent matrix third-time row cap '+k)
need(Q(p['ordinary_source_mixed_caps']['3,0'])==S3,'independent source third-time cap')
for k,v in enumerate((N0,N1,N2,N3)):
    need(Q(p['N_time_derivative_caps'][str(k)])==v,'independent actual N'+str(k))
need(Q(p['free_second_time_field_N_cap'])==free2 and Q(p['free_third_time_field_N_cap'])==free3,'independent stabilized free-wave ledgers')
need({k:Q(p['initial_third_time_field_caps'][k]) for k in ('displacement','velocity','radial_displacement')}==
     {'displacement':v2,'velocity':v3,'radial_displacement':dr3},'independent nonzero initial gauge third-time data')
need({k:Q(p['gauge_forcing_second_derivative_caps'][k]) for k in ('2,0','1,1','0,2')}==
     {'2,0':Txx,'1,1':Txr,'0,2':Trr},'independent second differentiated forcing caps')
need({k:Q(v) for k,v in p['actual_gauge_scaled_mixed_order_four_caps'].items()}==four,'all five fourth mixed derivative caps agree')
mix={'0,0':theta*N0,'1,0':N0,'0,1':N0,'2,0':N1,'1,1':N1,'0,2':B2,
     '3,0':N2,'2,1':N2,'1,2':N2+Tx,'0,3':N2+Tr,**four}
need(set(mix)==set(p['actual_mixed_scaled_derivative_caps']),'all15 mixed derivative slots present')
for k,v in mix.items():
    need(list(map(Q,p['actual_mixed_scaled_derivative_caps'][k]))==[v,v,v], 'all three amplitude rows '+k)
need(list(map(Q,p['actual_gauge_uniform_scaled_order_caps_0_to_4']))==[theta*N0,N0,B2,B3,max(four.values())],
     'independent uniform actual gauge derivative caps through order4')
readables=list(map(Q,p['readable_strict_uniform_scaled_order_caps_0_to_4']))
need(all(v<c for v,c in zip([theta*N0,N0,B2,B3,max(four.values())],readables)),'strict readable actual gauge cap displays')
result['producer_report_sha256']=producerpin
result['all15_mixed_slots_and_three_amplitudes_match']=True
print(json.dumps(encode(result),sort_keys=True,indent=2))
