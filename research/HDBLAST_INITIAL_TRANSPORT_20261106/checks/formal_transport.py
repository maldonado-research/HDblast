"""Exact truncated-polynomial checks of initial chart transport.

The independent composition uses sqrt(1+x)-1 through degree two and a
nontrivial three-state quadratic map. Coefficients of the Cauchy data are
indeterminates, not sampled data. This supplements the general chain rule
in the proof; it does not certify a continuous flow or Sobolev theorem.
"""
from fractions import Fraction as F


# e,t followed by a_i,b_i,c_i,d_i, i=0,1,2.
DIM = 14
ZERO = (0,)*DIM


class P:
    def __init__(self, terms=None):
        self.terms = {k:F(v) for k,v in (terms or {}).items() if v}

    @staticmethod
    def cast(value):
        return value if isinstance(value,P) else P({ZERO:F(value)})

    @staticmethod
    def var(index):
        power = list(ZERO)
        power[index] = 1
        return P({tuple(power):1})

    def __add__(self, other):
        other=P.cast(other)
        result=dict(self.terms)
        for key,value in other.terms.items():
            result[key]=result.get(key,F(0))+value
        return P(result)

    __radd__=__add__

    def __neg__(self):
        return P({k:-v for k,v in self.terms.items()})

    def __sub__(self, other):
        return self+-P.cast(other)

    def __mul__(self, other):
        other=P.cast(other)
        result={}
        for ka,va in self.terms.items():
            for kb,vb in other.terms.items():
                key=tuple(a+b for a,b in zip(ka,kb))
                if key[0]>2 or key[1]>1:
                    continue
                result[key]=result.get(key,F(0))+va*vb
        return P(result)

    __rmul__=__mul__

    def extract(self, e_order, t_order):
        return P({(0,0)+key[2:]:value for key,value in self.terms.items()
                  if key[:2]==(e_order,t_order)})

    def zero_variable(self, index):
        return P({key:value for key,value in self.terms.items() if key[index]==0})

    def __eq__(self, other):
        return self.terms==P.cast(other).terms


def evaluate(data, check):
    e,t=P.var(0),P.var(1)
    a=[P.var(2+i) for i in range(3)]
    b=[P.var(5+i) for i in range(3)]
    c=[P.var(8+i) for i in range(3)]
    d=[P.var(11+i) for i in range(3)]
    y=[e*a[i]+F(1,2)*e*e*b[i]+t*(e*c[i]+F(1,2)*e*e*d[i])
       for i in range(3)]
    # These are independently derived sqrt Taylor coefficients, not DATA.
    r=F(1,2)*y[0]-F(1,8)*y[0]*y[0]
    check('root_square_identity_through_e2_t1',(1+r)*(1+r)==1+y[0])
    check('root_first_derivative',r.extract(1,0)==F(data['root_first_derivative'])*a[0])
    check('root_second_derivative',2*r.extract(2,0)==
          F(data['root_first_derivative'])*b[0]+F(data['root_second_derivative'])*a[0]*a[0])
    check('root_second_velocity',2*r.extract(2,1)==
          F(data['root_first_derivative'])*d[0]+
          F(data['velocity_cross_factor'])*F(data['root_second_derivative'])*a[0]*c[0])
    u=[r,y[1],y[2]]

    # Non-diagonal first differential and symmetric second differential,
    # with independent explicit time dependence to exercise A_t and B_t.
    l0=[[F((i+1)*(j+2)+1,7) for j in range(3)] for i in range(3)]
    lt=[[F((i+2)*(j+1)-4,11) for j in range(3)] for i in range(3)]
    h0=[[[F((i+1)*(j+k+3)+j*k+1,13) for k in range(3)]
          for j in range(3)] for i in range(3)]
    ht=[[[F((i+2)*(j+k+2)+j*k-3,17) for k in range(3)]
          for j in range(3)] for i in range(3)]

    def linear(matrix, values):
        return [sum((matrix[i][j]*values[j] for j in range(3)),P()) for i in range(3)]

    def bilinear(tensor,left,right):
        return [sum((tensor[i][j][k]*left[j]*right[k] for j in range(3) for k in range(3)),P())
                for i in range(3)]

    first=F(data['root_first_derivative'])
    second=F(data['root_second_derivative'])
    def root_linear(values):
        return [first*values[0],values[1],values[2]]

    def root_bilinear(left,right):
        return [second*left[0]*right[0],P(),P()]

    def composed_a(matrix,values):
        return linear(matrix,root_linear(values))

    def composed_b(matrix,tensor,left,right):
        one=linear(matrix,root_bilinear(left,right))
        two=bilinear(tensor,root_linear(left),root_linear(right))
        return [one[i]+two[i] for i in range(3)]

    term_counts={}
    for label,linear_time,hessian_time in [('generic_time',lt,ht),
        ('even_frame',[[F(0)]*3 for _ in range(3)],[[[F(0)]*3 for _ in range(3)] for _ in range(3)])]:
        L=[[l0[i][j]+t*linear_time[i][j] for j in range(3)] for i in range(3)]
        H=[[[h0[i][j][k]+t*hessian_time[i][j][k] for k in range(3)] for j in range(3)] for i in range(3)]
        raw_linear=linear(L,u)
        raw_quadratic=bilinear(H,u,u)
        actual=[raw_linear[i]+F(1,2)*raw_quadratic[i] for i in range(3)]
        aa,ab,ac,ad=(composed_a(l0,v) for v in (a,b,c,d))
        ata,atb=(composed_a(linear_time,v) for v in (a,b))
        baa=composed_b(l0,h0,a,a)
        bac=composed_b(l0,h0,a,c)
        btaa=composed_b(linear_time,hessian_time,a,a)
        for i in range(3):
            check(label+':first_value:'+str(i),actual[i].extract(1,0)==aa[i])
            check(label+':second_value:'+str(i),2*actual[i].extract(2,0)==ab[i]+baa[i])
            check(label+':first_velocity:'+str(i),actual[i].extract(1,1)==ac[i]+ata[i])
            check(label+':second_velocity:'+str(i),2*actual[i].extract(2,1)==
                  ad[i]+atb[i]+F(data['velocity_cross_factor'])*bac[i]+btaa[i])
            term_counts[label+':'+str(i)]=len(actual[i].terms)
        # A false formula with cross factor one differs as a polynomial.
        check(label+':wrong_velocity_factor_detected',
              2*actual[0].extract(2,1)!=ad[0]+atb[0]+bac[0]+btaa[0])

    # Temporal orthonormal angular factors are even: sech t has no t term.
    # Check their exact quadratic inverses separately without discarding t^2.
    for angular_slots in (0,1,2):
        cosh_coefficient=F(angular_slots,2)
        sech_coefficient=-F(angular_slots,2)
        check('even_frame_quadratic_inverse_slots_'+str(angular_slots),
              cosh_coefficient+sech_coefficient==0)
    # c_1 represents a nonzero mixed temporal-radial first velocity. The
    # non-diagonal Hessian makes dropping it change second transport.
    second_velocity=2*actual[0].extract(2,1)
    check('mixed_temporal_radial_velocity_omission_detected',
          second_velocity!=second_velocity.zero_variable(9))
    return {'formal_state_dimension':3,'formal_indeterminate_count':DIM,
            'root_first_derivative':str(first),'root_second_derivative':str(second),
            'second_velocity_cross_factor':data['velocity_cross_factor'],
            'composition_monomial_counts':term_counts,
            'scope':'Exact symbolic coefficient identities for non-diagonal quadratic map jets; no numerical flow sampling.'}
