"""Exact independent indexed contractions of source-adapter formulas."""
from fractions import Fraction as Q
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from source_adapters import connection_contractions,gauge_source,scalar_bulk_source,shell_sources,Y2
sgn=(1,-1,1,1,1)

def check(seed, symbolic=False, point=None):
    keys=('A','B','C','D','E','L','M','Ay','At','By','Bt','Cy','Ct','Dy','Dt','Ey','Et','Ly','Lt','My','Mt')
    m={k:Q((i+1)*(seed+2)%17-8, i%4+2) for i,k in enumerate(keys)}
    b={'rho':Q(13,2),'c':Q(5,4),'a':Q(-2,7),'theta':Q(3,5),'p':Q(-3,11),'py':Q(2,13),'Vppp':Q(-7,9)}
    s={k:Q((i+3)*(seed+1)%13-6,i+2) for i,k in enumerate(('s','sy','st','syy','syt','stt'))}
    if symbolic:
        import source_laurent as ring
        symbols=dict(zip(ring.NAMES,ring.variables()))
        m={k:symbols[k] for k in m}
        b={k:symbols[k] for k in b}
        s={k:symbols[k] for k in s}
    n=point or (Q(3,5),Q(0),Q(4,5),Q(0))
    if n[3] != 0:raise RuntimeError('This checker uses the n4=0 Householder section.')
    Es=tuple(tuple(Q(i==j)-n[i]*(n[j]-Q(j==3)) for j in range(4)) for i in range(3))
    AA=(1,-1,0,0)
    Y=sum(AA[i]*n[i]*n[i] for i in range(4));v=[2*sum(e[i]*AA[i]*n[i] for i in range(4)) for e in Es]
    HH=[[2*sum(Es[i][k]*AA[k]*Es[j][k] for k in range(4))-2*Y*(i==j) for j in range(3)] for i in range(3)]
    g2=sum(x*x for x in v);h2=sum(x*x for row in HH for x in row)
    rho,c,a,h=b['rho'],b['c'],b['a'],b['theta']/b['rho']
    def tensor(suffix=''):
        A,B,C,D,E,L,M=(m[k+suffix] for k in ('A','B','C','D','E','L','M'))
        q=[[Q(0) for _ in range(5)] for _ in range(5)]
        q[0][0]=A*Y;q[0][1]=q[1][0]=B*Y;q[1][1]=D*Y
        for i in range(3):
            q[0][i+2]=q[i+2][0]=C*v[i];q[1][i+2]=q[i+2][1]=E*v[i]
            for j in range(3):q[i+2][j+2]=L*Y*(i==j)+M*HH[i][j]
        return q
    q=tensor();raw=[tensor('y'),[[x/rho for x in row] for row in tensor('t')]]
    for k in range(3):
        qq=[[Q(0) for _ in range(5)] for _ in range(5)]
        qq[0][0]=m['A']*v[k];qq[0][1]=qq[1][0]=m['B']*v[k];qq[1][1]=m['D']*v[k]
        for i in range(3):
            qq[0][i+2]=qq[i+2][0]=m['C']*HH[k][i];qq[1][i+2]=qq[i+2][1]=m['E']*HH[k][i]
            for j in range(3):qq[i+2][j+2]=m['L']*v[k]*(i==j)-m['M']*((k==i)*v[j]+(k==j)*v[i]+2*(i==j)*v[k])
        raw.append([[x/(rho*c) for x in row] for row in qq])
    G=[[[Q(0) for _ in range(5)] for _ in range(5)] for _ in range(5)]
    G[1][0][1]=a;G[1][1][0]=a
    for k in range(3):
        ii=k+2;G[ii][0][ii]=a;G[ii][1][ii]=h;G[ii][ii][1]=h;G[ii][ii][0]=-a
    dq=[[[raw[k][i][j]-sum(G[k][i][l]*q[l][j]+G[k][j][l]*q[i][l] for l in range(5)) for j in range(5)] for i in range(5)] for k in range(5)]
    qP=[sum(sgn[i]*sgn[j]*q[i][j]*(dq[i][j][k]+dq[j][i][k]-dq[k][i][j])/2 for i in range(5) for j in range(5)) for k in range(5)]
    co=connection_contractions(m,b)
    predict=[co[key][0]*Y*Y+co[key][1]*g2+co[key][2]*h2 for key in ('y','0')]
    predict.extend(co['a'][0]*Y*v[i]+co['a'][1]*sum(HH[i][j]*v[j] for j in range(3)) for i in range(3))
    if symbolic:
        if any(x-y!=0 for x,y in zip(qP,predict)):raise ValueError('SOURCE_QP_IDENTITY')
    elif qP!=predict:raise ValueError('SOURCE_QP_IDENTITY')
    # Independent scalar Hessian via e_i e_j eta and background connection.
    z,zy,zt,zyy,zyt,ztt=(s[k] for k in ('s','sy','st','syy','syt','stt'))
    grad=[zy*Y,zt/rho*Y]+[z/(rho*c)*x for x in v]
    eh=[[Q(0) for _ in range(5)] for _ in range(5)]
    eh[0][0]=zyy*Y;eh[0][1]=(zyt-a*zt)/rho*Y;eh[1][0]=zyt/rho*Y;eh[1][1]=ztt/(rho*rho)*Y
    for i in range(3):
        eh[0][i+2]=(zy-a*z)/(rho*c)*v[i];eh[i+2][0]=zy/(rho*c)*v[i]
        eh[1][i+2]=(zt-b['theta']*z)/(rho*rho*c)*v[i];eh[i+2][1]=zt/(rho*rho*c)*v[i]
        for j in range(3):eh[i+2][j+2]=z/(rho*rho*c*c)*HH[i][j]
    eh=[[eh[i][j]-sum(G[i][j][k]*grad[k] for k in range(5)) for j in range(5)] for i in range(5)]
    phiH=(b['py'],-a*b['p'],a*b['p'],a*b['p'],a*b['p'])
    sq=sum(phiH[i]*sgn[j]*q[i][j]**2 for i in range(5) for j in range(5))
    qh=sum(sgn[i]*sgn[j]*q[i][j]*eh[i][j] for i in range(5) for j in range(5))
    want=2*(sq-qh+b['p']*qP[0])-b['Vppp']*z*z*Y*Y
    coef=scalar_bulk_source(m,s,b);S=n[0]**2+n[1]**2;Z=S-Q(1,2);F=Y*Y-S/2+Q(1,12)
    actual=sum(x*y for x,y in zip(coef,(1,Z,F)))
    if symbolic:
        if want-actual!=0:raise ValueError('SOURCE_QPHI_IDENTITY')
    elif want!=actual:raise ValueError('SOURCE_QPHI_IDENTITY')


def rank(rows):
    rows=[list(row) for row in rows];r=0
    for col in range(len(rows[0])):
        pivot=next((i for i in range(r,len(rows)) if rows[i][col]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        divisor=rows[r][col];rows[r]=[x/divisor for x in rows[r]]
        for i in range(r+1,len(rows)):
            factor=rows[i][col];rows[i]=[x-factor*y for x,y in zip(rows[i],rows[r])]
        r+=1
    return r
def validate():
    for seed in range(8):check(seed)
    print('PASS: 8 exact indexed metric-jet witnesses verify all five qP components and Qphi.')
    
    # The invariant angular primitive spaces are proved by the derivative/product
    # identities in check_source_identities.py. These exact angular witnesses span
    # each primitive space, and all radial/time coefficients and backgrounds here
    # are independent symbolic variables, not numerical sample amplitudes.

    points=((Q(1),Q(0),Q(0),Q(0)),(Q(0),Q(0),Q(1),Q(0)),(Q(3,5),Q(0),Q(4,5),Q(0)),(Q(1,3),Q(2,3),Q(2,3),Q(0)))
    scalar_rows=[]
    vector_rows=[]
    for point in points:
        Y=point[0]**2-point[1]**2;S=point[0]**2+point[1]**2
        scalar_rows.append((Y*Y,4*(S-Y*Y),8-8*S+24*Y*Y))
        Es=tuple(tuple(Q(i==j)-point[i]*(point[j]-Q(j==3)) for j in range(4)) for i in range(3))
        DY=tuple(2*sum(E[k]*(1,-1,0,0)[k]*point[k] for k in range(4)) for E in Es)
        H=tuple(tuple(2*sum(Es[i][k]*(1,-1,0,0)[k]*Es[j][k] for k in range(4))-2*Y*(i==j) for j in range(3)) for i in range(3))
        vector_rows.extend((Y*DY[i],sum(H[i][j]*DY[j] for j in range(3))) for i in range(3))
        check(0,symbolic=True,point=point)
    if rank(scalar_rows)!=3 or rank(vector_rows)!=2:raise RuntimeError('Angular witnesses do not span primitive spaces.')
    print('PASS: universal symbolic radial/time/background identities on rank-complete angular witnesses (scalar rank 3, vector rank 2).')

if __name__ == '__main__':validate()
