"""Finite rational coefficient-jet majorants and exact reduced-row test jets.

Analytic coefficient extraction and the background/physical patch are
premises of COEFFICIENTS.md; the finite arithmetic is recomputed here.
"""
from fractions import Fraction as F
from math import factorial, comb
import json


def evaluate(root, data, ck):
    actual=json.loads((root/'inputs'/'ACTUAL_POTENTIAL_JET_BOUNDS.json').read_text())
    def mul(a,b,N):
        c=[F(0)]*(N+1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                if i+j<=N:c[i+j]+=x*y
        return c
    def add(*aa):return [sum(v,F(0)) for v in zip(*aa)]
    def scale(a,c):return [c*x for x in a]
    def power(a,n,N):
        o=[F(1)]+[F(0)]*N
        for _ in range(n):o=mul(o,a,N)
        return o
    N=data['background_order']
    U=[F(-2,27),F(0),F(14,9),F(50,27),F(-1,6),F(-4,9),F(-2,27)]
    ck('source actual potential coefficients',U==[F(x) for x in actual['actual_model']['U_coefficients']])
    W=[F(1,3),F(0),F(1),F(1,3)]+[F(0)]*3
    Wp=[(i+1)*W[i+1] for i in range(6)]+[F(0)]
    ck('actual potential from superpotential',add(scale(mul(Wp,Wp,6),F(1,2)),scale(mul(W,W,6),F(-2,3)))==U)
    a0=dict(rho=F(131),b=F(1,129),a=F(14,125),p=F(7,50000),eta=F(9,100000))
    for k in ['rho','a','p','eta']:
        ck('source all-point box '+k,a0[k]==max(abs(F(v)) for v in actual['background_intervals'][k]['exact']))
    ck('source reciprocal rho box',a0['b']==1/F(actual['background_intervals']['rho']['exact'][0]))
    y={k:[v]+[F(0)]*N for k,v in a0.items()}
    for n in range(N):
        rho,b,a,p,eta=[y[k] for k in ['rho','b','a','p','eta']]
        up=[F(0)]*(N+1)
        for j in range(1,7):up=add(up,scale(power(eta,j-1,N),abs(j*U[j])))
        rhs=dict(rho=mul(a,mul(rho,rho,N),N),b=a,a=add(b,scale(mul(rho,mul(p,p,N),N),F(1,3))),p=add(mul(rho,up,N),scale(mul(rho,mul(a,p,N),N),4)),eta=mul(rho,p,N))
        for k in y:
            raw=rhs[k][n]/(n+1)
            upper=F(-(-raw.numerator//raw.denominator))
            ck(f'upward background {k} order {n+1}',upper>=raw)
            y[k][n+1]=upper
    for k in y:
        for n in range(N+1):ck(f'background derivative {k} order {n} below 1e100',factorial(n)*y[k][n]<10**data['background_derivative_decade'])
    ck('complex chart denominator margin',6*F(1,100)+3*F(1,10000)<F(1,10))
    ck('retained mixed flow jets through eighteen', all(8*factorial(m)*factorial(k)*10**(6*m)*16000**k < 10**data['flow_derivative_decade'] for m in range(19) for k in range(19-m)))
    ck('frame derivative cap order19',8*factorial(19)*100**19<10**data['frame_derivative_decade'])
    ck('reconstruction raw derivative cap',6**18*10**(2*data['background_derivative_decade']+2*data['flow_derivative_decade']+2*data['frame_derivative_decade'])<10**data['reconstruction_derivative_decade'])
    ck('reference raw metric derivative cap',4**19*10**(2*data['background_derivative_decade']+2*data['frame_derivative_decade'])<10**data['reference_metric_derivative_decade'])
    ck('reference connection raw derivative cap',125*3*2**17*10**data['reconstruction_derivative_decade']<10**data['reference_connection_derivative_decade'])
    ck('all multiindex normalization count',comb(21+18,18)<10**12)
    ck('raw to aggregate input budget',10**6*10**12*10**data['reference_connection_derivative_decade']<10**data['aggregate_input_decade'])
    ck('inverse Jacobian entry sum',16*4*F(101,100)<65)
    ck('inverse jet bound',65*(1+65*10**data['aggregate_input_decade'])**16<10**data['inverse_jet_decade'])
    ck('final derivative reserve',factorial(16)*10**30*10**(8*data['inverse_jet_decade'])<10**data['final_coefficient_jet_decade'])
    ck('coefficient derivative order within derived range', N >= 19 and 0 <= data['coefficient_jet_order'] <= 16)
    # Exact full first/second jet reduction, with nonzero bulk reference data.
    d=5
    G=[[F(0) for _ in range(d)] for _ in range(d)]
    for i in range(d):G[i][i]=F(-2 if i==1 else i+2)
    G[0][2]=G[2][0]=F(1,7)
    def inv(A):
        a=[row[:]+[F(i==j) for j in range(len(A))] for i,row in enumerate(A)]
        for i in range(len(A)):
            p=next(j for j in range(i,len(A)) if a[j][i]);a[i],a[p]=a[p],a[i]
            z=a[i][i];a[i]=[x/z for x in a[i]]
            for j in range(len(A)):
                if i!=j:
                    z=a[j][i];a[j]=[x-z*y for x,y in zip(a[j],a[i])]
        return [r[len(A):] for r in a]
    q=inv(G)
    k=[[[F((i+2)*(a+b+1)%11-5,13) for b in range(d)] for a in range(d)] for i in range(d)]
    h=[[[[F((i+j+1)*(a+b+2)%13-6,17) for b in range(d)] for a in range(d)] for j in range(d)] for i in range(d)]
    K=[[[F((a+1)*(b+c+2)%7-3,19) for c in range(d)] for b in range(d)] for a in range(d)]
    Kx=[[[[F((i+1)*(a+2)*(b+c+1)%11-5,23) for c in range(d)] for b in range(d)] for a in range(d)] for i in range(d)]
    zero=[[[[F(0) for _ in range(d)] for _ in range(d)] for _ in range(d)] for _ in range(d)]
    def reduced(hh,kkx=Kx):
        dq=[[[-sum(q[a][r]*k[i][r][s]*q[s][b] for r in range(d) for s in range(d)) for b in range(d)] for a in range(d)] for i in range(d)]
        gam=[[[sum(q[a][r]*(k[b][r][c]+k[c][r][b]-k[r][b][c])/2 for r in range(d)) for c in range(d)] for b in range(d)] for a in range(d)]
        dg=[[[[sum(dq[i][a][r]*(k[b][r][c]+k[c][r][b]-k[r][b][c])/2+q[a][r]*(hh[i][b][r][c]+hh[i][c][r][b]-hh[i][r][b][c])/2 for r in range(d)) for c in range(d)] for b in range(d)] for a in range(d)] for i in range(d)]
        C=[sum(q[b][c]*(gam[a][b][c]-K[a][b][c]) for b in range(d) for c in range(d)) for a in range(d)]
        dc=[[sum(dq[i][b][c]*(gam[a][b][c]-K[a][b][c])+q[b][c]*(dg[i][a][b][c]-kkx[i][a][b][c]) for b in range(d) for c in range(d)) for a in range(d)] for i in range(d)]
        low=[sum(G[b][a]*C[a] for a in range(d)) for b in range(d)]
        nc=[[sum(k[i][b][a]*C[a]+G[b][a]*dc[i][a]-gam[a][i][b]*low[a] for a in range(d)) for b in range(d)] for i in range(d)]
        ric=[[sum(dg[c][c][a][b]-dg[b][c][a][c] for c in range(d))+sum(gam[c][c][r]*gam[r][a][b]-gam[c][b][r]*gam[r][a][c] for c in range(d) for r in range(d)) for b in range(d)] for a in range(d)]
        return [[ric[a][b]-(nc[a][b]+nc[b][a])/2 for b in range(d)] for a in range(d)]
    full,lo=reduced(h),reduced(zero)
    for a in range(d):
        for b in range(a,d):
            principal=-sum(q[i][j]*h[i][j][a][b] for i in range(d) for j in range(d))/2
            ck(f'reduced principal row {a}{b}',full[a][b]-lo[a][b]==principal)
    ck('reject wrong principal sign',full[0][0]-lo[0][0]!=sum(q[i][j]*h[i][j][0][0] for i in range(d) for j in range(d))/2)
    frozen=reduced(zero,zero)
    ck('reject frozen bulk reference derivative',any(frozen[a][b]!=lo[a][b] for a in range(d) for b in range(d)))
    # A rational moving point-map jet checks all 16 chain-rule rows.
    P=[[F((i+1)*(a+2)%13-6,7) for a in range(16)] for i in range(d)]
    J=[[F(i==j)+F((i+1)*(j+1),10000) for j in range(16)] for i in range(16)]
    Ji=inv(J)
    pxx=[[[F((i+j+2)*(a+1)%7-3,11) for a in range(16)] for j in range(d)] for i in range(d)]
    pxv=[[[F((i+1)*(a+b+2)%7-3,17) for b in range(16)] for a in range(16)] for i in range(d)]
    pvv=[[[F((a+1)*(b+c+1)%7-3,19) for c in range(16)] for b in range(16)] for a in range(16)]
    H=[[[F((i+j+1)*(a+1)%11-5,23) for a in range(16)] for j in range(d)] for i in range(d)]
    Q=[sum(q[i][j]*(sum(J[a][b]*H[i][j][b] for b in range(16))+pxx[i][j][a]+sum(pxv[i][a][b]*P[j][b]+pxv[j][a][b]*P[i][b] for b in range(16))+sum(pvv[a][b][c]*P[i][b]*P[j][c] for b in range(16) for c in range(16))) for i in range(d) for j in range(d)) for a in range(16)]
    def sub(factor):return [Q[a]-sum(q[i][j]*(pxx[i][j][a]+factor*sum(pxv[i][a][b]*P[j][b] for b in range(16))+sum(pvv[a][b][c]*P[i][b]*P[j][c] for b in range(16) for c in range(16))) for i in range(d) for j in range(d)) for a in range(16)]
    res=sub(data['point_map_mixed_factor']);wrong=sub(1)
    for a in range(16):ck(f'point map chain rule row {a}',sum(Ji[a][b]*res[b] for b in range(16))==sum(q[i][j]*H[i][j][a] for i in range(d) for j in range(d)))
    ck('reject missing chain rule factor two',any(sum(Ji[a][b]*wrong[b] for b in range(16))!=sum(q[i][j]*H[i][j][a] for i in range(d) for j in range(d)) for a in range(16)))
    return {'background_raw_derivative_maxima':{k:str(max(factorial(n)*v for n,v in enumerate(y[k]))) for k in y}, 'background_order':N, 'coefficient_order':data['coefficient_jet_order'], 'coefficient_jet_ceiling':'10^'+str(data['final_coefficient_jet_decade'])}
