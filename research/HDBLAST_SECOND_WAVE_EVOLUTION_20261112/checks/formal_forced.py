"""Exact forced temporal jets and mixed product identities.

The time recurrence uses independent formal source and spatial jets in a
dense two-fiber coefficient example. It is checked against independently
multiplied normalized Taylor series. Mixed differential products are exact
bivariate polynomial identities, including both u*u_tt and u*u_zz.
"""
from fractions import Fraction as F
from math import comb,factorial


class Linear:
    """Rational linear form in named, independent jet indeterminates."""
    def __init__(self,terms=None):
        self.terms={k:F(v) for k,v in (terms or {}).items() if v}
    def __add__(self,other):
        if other==0:return self
        out=dict(self.terms)
        for key,value in other.terms.items():out[key]=out.get(key,F(0))+value
        return Linear(out)
    __radd__=__add__
    def __neg__(self):return Linear({k:-v for k,v in self.terms.items()})
    def __sub__(self,other):return self+-other
    def __mul__(self,scalar):return Linear({k:v*F(scalar) for k,v in self.terms.items()})
    __rmul__=__mul__
    def __truediv__(self,scalar):return self* (1/F(scalar))
    def __eq__(self,other):
        if isinstance(other,(int,F)):return other==0 and not self.terms
        return isinstance(other,Linear) and self.terms==other.terms


def matvec(a,v):
    return [sum((a[i][j]*v[j] for j in range(2)),Linear()) for i in range(2)]


def temporal(data,check):
    last=11
    def jet(name,j):return [Linear({name+':'+str(j)+':'+str(a):1}) for a in range(2)]
    # These arrays are normalized Taylor coefficients, not derivatives.
    zz=[jet('zz',j) for j in range(last+1)]
    lap=[jet('lap',j) for j in range(last+1)]
    z=[jet('z',j) for j in range(last+1)]
    ang=[jet('ang',j) for j in range(last+1)]
    force=[jet('F',j) for j in range(last+1)]
    s=[F((-1)**j*(j+2),j+3) for j in range(last+1)]
    maps={name:[[[F((j+1)*(a+2)+b+offset,17+offset) for b in range(2)]
                 for a in range(2)] for j in range(last+1)]
          for name,offset in [('t',1),('z',3),('ang',5),('zero',7)]}
    rho2=F(16900)
    u=[jet('U',0),jet('U',1)]
    # Compute u coefficients from literal Taylor convolution, with no
    # binomial coefficients. The source coefficient is always present.
    for n in range(last+1):
        row=[zz[n][a]+rho2*force[n][a] for a in range(2)]
        for k in range(n+1):
            dt=matvec(maps['t'][k],[(n-k+1)*v for v in u[n-k+1]])
            dz=matvec(maps['z'][k],z[n-k])
            da=matvec(maps['ang'][k],ang[n-k])
            zero=matvec(maps['zero'][k],u[n-k])
            row=[row[a]+s[k]*lap[n-k][a]-dt[a]-dz[a]-da[a]-zero[a] for a in range(2)]
        u.append([v/F((n+1)*(n+2)) for v in row])
    for j in range(last+1):
        rhs=[factorial(j)*zz[j][a]+data['forcing_initial_factor']*rho2*factorial(j)*force[j][a] for a in range(2)]
        for k in range(j+1):
            dt=matvec([[factorial(k)*v for v in row] for row in maps['t'][k]],
                      [factorial(j-k+1)*v for v in u[j-k+1]])
            dz=matvec([[factorial(k)*v for v in row] for row in maps['z'][k]],
                      [factorial(j-k)*v for v in z[j-k]])
            da=matvec([[factorial(k)*v for v in row] for row in maps['ang'][k]],
                      [factorial(j-k)*v for v in ang[j-k]])
            zero=matvec([[factorial(k)*v for v in row] for row in maps['zero'][k]],
                        [factorial(j-k)*v for v in u[j-k]])
            rhs=[rhs[a]+comb(j,k)*(factorial(k)*s[k]*factorial(j-k)*lap[j-k][a]
                                  -dt[a]-dz[a]-da[a]-zero[a]) for a in range(2)]
        for a in range(2):
            check('forced_initial_recurrence:'+str(j)+':'+str(a),factorial(j+2)*u[j+2][a]==rhs[a])
            other=rhs[a]-factorial(j)*zz[j][a]-data['forcing_initial_factor']*rho2*factorial(j)*force[j][a]
            normal_rhs=factorial(j+2)*u[j+2][a]-other+data['normal_forcing_sign']*rho2*factorial(j)*force[j][a]
            check('forced_normal_source_sign:'+str(j)+':'+str(a),factorial(j)*zz[j][a]==normal_rhs)
        check('nonzero_source_jet:'+str(j),all(force[j][a]!=0 for a in range(2)))
    return {'highest_initial_time_jet':13,'recurrence_steps':12,
            'formal_two_fiber_source_jet_count':24,
            'largest_temporal_linear_form_terms':max(len(x.terms) for row in u for x in row)}


class P:
    def __init__(self,terms=None):self.terms={k:F(v) for k,v in (terms or {}).items() if v}
    @staticmethod
    def cast(x):return x if isinstance(x,P) else P({(0,0):F(x)})
    def __add__(self,other):
        out=dict(self.terms)
        for key,value in P.cast(other).terms.items():out[key]=out.get(key,F(0))+value
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.terms.items()})
    def __sub__(self,other):return self+-P.cast(other)
    def __mul__(self,other):
        out={}
        for (a,b),v in self.terms.items():
            for (c,d),w in P.cast(other).terms.items():out[a+c,b+d]=out.get((a+c,b+d),F(0))+v*w
        return P(out)
    __rmul__=__mul__
    def derivative(self,t,z):
        return P({(a-t,b-z):v*F(factorial(a),factorial(a-t))*F(factorial(b),factorial(b-z))
                  for (a,b),v in self.terms.items() if a>=t and b>=z})
    def __eq__(self,other):return self.terms==P.cast(other).terms


def mixed(data,check):
    def sample(seed,degree):
        return P({(t,z):F((seed+t+1)*(z+2)+t*z,seed+11)
                  for t in range(degree+1) for z in range(degree+1-t)})
    u=[sample(2,5),sample(5,6)]
    coeff=[[[sample(7+3*a+b+2*c,3) for c in range(2)] for b in range(2)] for a in range(2)]
    cases={
        'value_value':(u,u),
        'value_first':(u,[v.derivative(1,0) for v in u]),
        'first_first':([v.derivative(1,0) for v in u],[v.derivative(0,1) for v in u]),
        'principal_tt':(u,[v.derivative(2,0) for v in u]),
        'principal_zz':(u,[v.derivative(0,2) for v in u]),
    }
    count=0
    for label,(left,right) in cases.items():
        direct=[sum((coeff[a][b][c]*left[b]*right[c] for b in range(2) for c in range(2)),P()) for a in range(2)]
        for t in range(5):
            for z in range(5-t):
                candidate=[P(),P()]
                for at in range(t+1):
                    for bt in range(t-at+1):
                        ct=t-at-bt
                        for az in range(z+1):
                            for bz in range(z-az+1):
                                cz=z-az-bz
                                multiplicity=comb(t,at)*comb(t-at,bt)*comb(z,az)*comb(z-az,bz)
                                # Setting the factor to zero intentionally
                                # replaces every multinomial weight by one.
                                weight=1+data['mixed_product_binomial_factor']*(multiplicity-1)
                                for a in range(2):
                                    candidate[a]+=weight*sum((coeff[a][b][c].derivative(at,az)*
                                        left[b].derivative(bt,bz)*right[c].derivative(ct,cz)
                                        for b in range(2) for c in range(2)),P())
                for a in range(2):
                    check('mixed_Leibniz:'+label+':'+str(t)+':'+str(z)+':'+str(a),
                          direct[a].derivative(t,z)==candidate[a])
                    count+=1
    # Independently multiply amplitude Taylor arrays of q(V) and H(V).
    q0=P.cast(F(3,7));q1=u[0]+2*u[1]
    v=[sample(3,4),sample(9,4)]
    zero=[P(),P()]
    ve=[zero,u,[F(1,2)*x for x in v]]
    du=[x.derivative(1,0)+x.derivative(0,1) for x in u]
    dv=[x.derivative(1,0)+x.derivative(0,1) for x in v]
    dve=[zero,du,[F(1,2)*x for x in dv]]
    def bilinear(left,right):
        return [sum((coeff[a][b][c]*left[b]*right[c] for b in range(2) for c in range(2)),P()) for a in range(2)]
    c2uu,b1udu,d2du=bilinear(u,u),bilinear(u,du),bilinear(du,du)
    for direction in ((2,0),(0,2)):
        for a in range(2):
            hu=u[a].derivative(*direction);hv=v[a].derivative(*direction)
            exact_second=2*(F(1,2)*q0*hv+q1*hu)
            candidate=q0*hv+data['principal_source_factor']*q1*hu
            check('principal_second_variation:'+str(direction[0])+':'+str(direction[1])+':'+str(a),exact_second==candidate)
            # Build the second coefficient from amplitude Taylor products
            # q(V)H(V)-c(V)-b(V)DV-d(DV,DV), with generic dense bilinear
            # maps and c(V) quadratic coefficient one half.
            c_linear=sum((F((a+1)*(b+2),5)*ve[2][b] for b in range(2)),P())
            b_linear=sum((F((a+2)*(b+1),7)*dve[2][b] for b in range(2)),P())
            quadratic=P()
            for k in range(3):
                quadratic+=F(1,2)*bilinear(ve[k],ve[2-k])[a]
                quadratic+=bilinear(ve[k],dve[2-k])[a]+bilinear(dve[k],dve[2-k])[a]
            direct_residual=2*(F(1,2)*q0*hv+q1*hu-c_linear-b_linear-quadratic)
            linear_plus_v=q0*hv-2*c_linear-2*b_linear
            energy_force=data['principal_source_factor']*q1*hu-c2uu[a]-data['mixed_source_factor']*b1udu[a]-data['first_jet_source_factor']*d2du[a]
            check('complete_energy_sign_second_source:'+str(direction[0])+':'+str(direction[1])+':'+str(a),
                  direct_residual==linear_plus_v+energy_force)
    return {'mixed_polynomial_rows':count,'mixed_total_differential_order_tested':4,
            'principal_second_derivatives_tested':['tt','zz'],
            'scope':'Exact polynomial identities; full order-twelve derivative allocation is checked separately.'}


def evaluate(data,check):
    return {'temporal':temporal(data,check),'mixed':mixed(data,check)}
