"""Exact SO(4) scalar-module, natural-tensor and temporal-jet identities.

Finite algebra supports, rather than proves, the analytic equivariance and
uniqueness arguments. All sphere integrals below use probability measure;
the full geometric volume cancels in every checked identity.
"""
from fractions import Fraction as F
from math import factorial

DIM = 12
ZERO = (0,) * DIM


class P:
    def __init__(self, terms=None):
        self.d = {k:F(v) for k,v in (terms or {}).items() if v}
    @staticmethod
    def cast(v):
        return v if isinstance(v,P) else P({ZERO:F(v)})
    @staticmethod
    def var(i):
        k=list(ZERO);k[i]=1
        return P({tuple(k):1})
    def __add__(self,v):
        d=dict(self.d)
        for k,c in P.cast(v).d.items():d[k]=d.get(k,F(0))+c
        return P(d)
    __radd__=__add__
    def __neg__(self):return P({k:-c for k,c in self.d.items()})
    def __sub__(self,v):return self+-P.cast(v)
    def __rsub__(self,v):return P.cast(v)+-self
    def __mul__(self,v):
        d={}
        for a,b in self.d.items():
            for c,e in P.cast(v).d.items():
                k=tuple(x+y for x,y in zip(a,c));d[k]=d.get(k,F(0))+b*e
        return P(d)
    __rmul__=__mul__
    def __eq__(self,v):return self.d==P.cast(v).d
    def diff(self,i):
        d={}
        for a,c in self.d.items():
            if a[i]:
                k=list(a);k[i]-=1;d[tuple(k)]=c*a[i]
        return P(d)
    def at_zero(self,i):return P({k:c for k,c in self.d.items() if not k[i]})


def integral(p):
    out=F(0)
    for k,c in p.d.items():
        if any(k[4:]):raise ValueError('Non-sphere variable in sphere integral')
        if any(x%2 for x in k[:4]):continue
        numerator=1
        for x in k[:4]:
            for j in range(1,x,2):numerator*=j
        denominator=1
        for j in range(sum(k[:4])//2):denominator*=4+2*j
        out+=c*F(numerator,denominator)
    return out


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def tr(a):return [list(x) for x in zip(*a)]


def evaluate(data,check):
    n=[P.var(i) for i in range(4)]
    pairs=[(a,b) for a in range(4) for b in range(a+1,4)]
    r2=sum((x*x for x in n),P())
    def rotate(p,a,b):return n[a]*p.diff(b)-n[b]*p.diff(a)
    Y=n[0]*n[0]-n[1]*n[1]
    casimir=-sum((rotate(rotate(Y,a,b),a,b) for a,b in pairs),P())
    check('selected_quadrupole_casimir',casimir==data['casimir']*Y)
    basis=[n[i]*n[j] for i,j in pairs]+[n[i]*n[i]-n[3]*n[3] for i in range(3)]
    matrices=[]
    for a,b in pairs:
        cols=[]
        for q in basis:
            v=rotate(q,a,b);coeff=[]
            for i,j in pairs:
                k=list(ZERO);k[i]+=1;k[j]+=1;coeff.append(v.d.get(tuple(k),F(0)))
            for i in range(3):
                k=list(ZERO);k[i]=2;coeff.append(v.d.get(tuple(k),F(0)))
            check('module_closed:'+str(a)+str(b)+':'+str(len(cols)),
                  v==sum((c*q for c,q in zip(coeff,basis)),P()))
            cols.append(coeff)
        matrices.append(tr(cols))
    gram=[[integral(x*y) for y in basis] for x in basis]
    check('quadrupole_probability_norm',integral(Y*Y)==F(1,6))
    for k,R in enumerate(matrices):
        left=mm(tr(R),gram);right=mm(gram,R)
        check('skew_adjoint_module:'+str(k),all(left[i][j]+right[i][j]==0
               for i in range(9) for j in range(9)))
    word=gram
    for order in range(15):
        check('all_ordered_word_Gram:'+str(order),word==[[data['casimir']**order*x for x in row] for row in gram])
        children=[mm(mm(tr(R),word),R) for R in matrices]
        word=[[sum(c[i][j] for c in children) for j in range(9)] for i in range(9)]
    K=[]
    for a,b in pairs:
        v=[P() for _ in range(4)];v[a]=-n[b];v[b]=n[a];K.append(v)
    for i in range(4):
        check('Killing_tight_frame_row:'+str(i),all(sum((v[i]*v[j] for v in K),P())==
              ((r2 if i==j else P())-n[i]*n[j]) for j in range(4)))
    metric=[[(1 if i==j else 0)-n[i]*n[j] for j in range(4)] for i in range(4)]
    grad=[sum((metric[i][j]*Y.diff(j) for j in range(4)),P()) for i in range(4)]
    for k,(a,b) in enumerate(pairs):
        v=K[k]
        lie_grad=[sum((v[j]*grad[i].diff(j)+v[j].diff(i)*grad[j] for j in range(4)),P()) for i in range(4)]
        target=[sum((metric[i][j]*rotate(Y,a,b).diff(j) for j in range(4)),P()) for i in range(4)]
        check('natural_projected_gradient_intertwiner:'+str(k),lie_grad==target)
        lie_metric=[[sum((v[h]*metric[i][j].diff(h)+v[h].diff(i)*metric[h][j]+v[h].diff(j)*metric[i][h]
                          for h in range(4)),P()) for j in range(4)] for i in range(4)]
        check('natural_sphere_metric_invariance:'+str(k),all(q==0 for row in lie_metric for q in row))
    # At the north pole translating generators have zero tangential spin,
    # while stabilizer generators have zero vector value. Thus cross terms
    # vanish for every tensor representation, not just scalar components.
    north=[0,0,0,1];Pi=[[F(int(i==j and i<3)) for j in range(4)] for i in range(4)]
    cross=[[[F(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for a,b in pairs:
        M=[[F(0) for _ in range(4)] for _ in range(4)];M[a][b]=-1;M[b][a]=1
        value=[sum(M[i][j]*north[j] for j in range(4)) for i in range(4)]
        spin=mm(mm(Pi,M),Pi)
        for i in range(4):
            for j in range(4):
                for h in range(4):cross[i][j][h]+=value[i]*spin[j][h]
    check('north_pole_translation_spin_cross_zero',all(x==0 for row in cross for col in row for x in col))
    t=P.var(4);a0,a1,a2=[P.var(i) for i in (5,6,7)];u0,u1,u2,u3=[P.var(i) for i in (8,9,10,11)]
    coeff=a0+a1*t+F(1,2)*a2*t*t
    u=u0+u1*t+F(1,2)*u2*t*t+F(1,6)*u3*t*t*t
    actual=(coeff*u.diff(4)).diff(4).diff(4).at_zero(4)
    check('temporal_Leibniz_cross_factor',actual==a0*u3+data['time_cross_factor']*a1*u2+a2*u1)
    return {'scalar_module_dimension':9,'casimir':data['casimir'],'ordered_word_orders_checked':list(range(15)),
            'sphere_measure':'probability; common full volume cancels','scope':'Finite SO(4) and temporal-jet algebra. Analytic bundle equivariance, uniqueness and PDE estimates are written proof inputs.'}
