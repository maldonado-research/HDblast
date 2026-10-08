"""Generic certified finite-source target engine; performs no I/O or source callbacks.

Exact rational real polynomial coefficients are supplied by the guarded caller.
This module does not authenticate a public freeze and does not itself authorize
construction or decoding. All actual entry work belongs behind that outer guard.
"""
from fractions import Fraction as Q
from math import factorial
from flint import arb, acb, arb_poly, acb_poly, fmpq, ctx

BITS=512
SOURCE_DEGREE=112
MOMENTUM_DEGREE=832
MAX_K=Q(256)
EXPORT_BITS=96
EXPORT_SCALE=2**EXPORT_BITS
EXPORTED_L1_GATE=Q(1,10**18)
DELTA=Q(1,128)
SOURCES=('positive_B','signed_uB')

def require(ok,message):
    if not ok: raise ValueError(message)

def av(q):
    require(type(q) is Q,'Exact Fraction required')
    return arb(fmpq(q.numerator,q.denominator))

def geometry():
    out=[];left=DELTA
    while left<Q(1,2):
        right=min(3*left/2,Q(1,2))
        out.append({'left':left-5,'right':right-5,
                    'center':(left+right)/2-5,'half':(right-left)/2})
        left=right
    require(len(out)==11,'Eleven source parents required')
    return out

def exponential_tail(x,n):
    """Upper bound for sum_{j=n+1}^infinity x^j/j!, all exact rational."""
    require(type(x) is Q and x>=0 and type(n) is int and n>=0 and x<n+2,
            'Invalid positive factorial-tail premises')
    return x**(n+1)/factorial(n+1)/(1-x/Q(n+2))

def kernel_tails(source_l1_bound):
    require(type(source_l1_bound) is Q and source_l1_bound>=0,'Nonnegative exact source L1 bound required')
    n=MOMENTUM_DEGREE;x=MAX_K
    return {'W':source_l1_bound*exponential_tail(x,n),
            'U':source_l1_bound*x**(n+1)/(2*factorial(n+2))/(1-x/Q(n+3))}

def cap_bounds(source,k):
    require(source in SOURCES and type(k) is Q and 0<=k<=MAX_K,'Registered source and exact band node required')
    d=DELTA*(2-DELTA);ell=1/(5-DELTA);h=Q(3,8)**63
    pd=2*(1-DELTA)/d**2;rho=SOURCES.index(source)
    return {'W':h*(pd+rho+2*ell+2*k+DELTA*(4*k*k+4*k*ell+6*ell*ell)),
            'U':h*((pd+rho+2*ell)/2+1+DELTA*(3*ell*ell+2*ell+2*k))}

class EntireTarget:
    """Finite real-polynomial Duhamel response plus caller-certified target error.

The 11 rows are ordered dictionaries with exact Fraction 'center', 'half',
'coefficients' (113 ordinary powers in s-center), and 'uniform_error'.
Authentication, coefficient SHA256 comparison, and the certificate establishing
uniform_error are mandatory caller responsibilities, not inferred from values.
"""
    def __init__(self,rows):
        ctx.prec=BITS
        require(type(rows) is list and len(rows)==11,'All eleven ordered parents required')
        geom=geometry();self.source_l1=Q(0);self.source_error_l1=Q(0)
        moments=[arb(0) for _ in range(MOMENTUM_DEGREE+2)]
        for index,(row,g) in enumerate(zip(rows,geom)):
            require(row['center']==g['center'] and row['half']==g['half'],'Changed parent geometry')
            co=row['coefficients'];error=row['uniform_error']
            require(type(co) in (list,tuple) and len(co)==SOURCE_DEGREE+1 and all(type(c) is Q for c in co),
                    'Exactly113 exact real coefficients required')
            require(type(error) is Q and error>=0,'Nonnegative certified source error required')
            # This independently computed polynomial L1 majorant does not rely
            # on source samples, disk estimates, or a small coefficient premise.
            normalized=[c*g['half']**j for j,c in enumerate(co)]
            self.source_l1+=2*g['half']*sum((abs(c) for c in normalized),Q(0))
            self.source_error_l1+=2*g['half']*error
            p=arb_poly([av(c) for c in normalized]);r=av(g['half'])
            affine=arb_poly([av(-9-2*g['center']),-2*r])
            for n in range(MOMENTUM_DEGREE+2):
                integral=p.integral()
                moments[n]+=r*(integral(1)-integral(-1))
                if n<MOMENTUM_DEGREE+1:p=p*affine
        u=[];w=[];fact=1
        for n in range(MOMENTUM_DEGREE+1):
            if n:fact*=n
            uv=-moments[n+1]/(2*(n+1)*fact);wv=-moments[n]/fact
            def phase(v):
                if n%4==0:return acb(v,0)
                if n%4==1:return acb(0,v)
                if n%4==2:return acb(-v,0)
                return acb(0,-v)
            u.append(phase(uv));w.append(phase(wv))
        self.U=acb_poly(u);self.W=acb_poly(w)
        self.tails=kernel_tails(self.source_l1)

    def evaluate(self,k,source=None):
        """Return exact outward dyadic96 rectangles at one exact rational node.

source=None means the supplied finite-polynomial model (fabricated testing).
source=registered label adds certified source errors and the prescribed cap.
There is no source callback, data read, or adaptive precision/degree escalation.
"""
        require(ctx.prec==BITS,'Arithmetic precision changed')
        require(type(k) is Q and 0<=k<=MAX_K,'Exact node outside registered band')
        require(source is None or source in SOURCES,'Unknown source')
        extra={'U':self.tails['U'],'W':self.tails['W']}
        if source is not None:
            cap=cap_bounds(source,k)
            extra['U']+=cap['U']+self.source_error_l1/2
            extra['W']+=cap['W']+self.source_error_l1
        results={};x=av(k)
        for name,p in (('U',self.U),('W',self.W)):
            z=p(x)
            results[name]=export_rectangle(z,extra[name],real=(k==0))
        return results

def exact_q(ball):
    q=ball.fmpq()
    return Q(int(q.numerator),int(q.denominator))

def floor_grid(q):return (q.numerator*EXPORT_SCALE)//q.denominator

def ceil_grid(q):return -floor_grid(-q)

def export_rectangle(z,extra,real=False):
    require(type(extra) is Q and extra>=0,'Nonnegative exact analytic error required')
    out={}
    for label,part in (('real',z.real),('imag',z.imag)):
        if real and label=='imag':
            require(part.is_zero(),'k=0 imaginary arithmetic must be exact zero')
            out[label]=[0,0]
        else:
            out[label]=[floor_grid(exact_q(part.lower())-extra),
                        ceil_grid(exact_q(part.upper())+extra)]
        require(out[label][0]<=out[label][1],'Reversed exported interval')
    rad=Q(sum(out[label][1]-out[label][0] for label in ('real','imag')),2*EXPORT_SCALE)
    require(rad<=EXPORTED_L1_GATE,'Complete exported target L1 radius exceeds fixed gate')
    # Grid exponent is package-level metadata. No rounded decimal values appear.
    return out

def stored_state_coordinates(k,u,w,epsilon):
    """Exact normalization defect and interval-ready stable scaled amplitude.

u,w are pairs of exact rationals decoded from original represented components.
This function neither changes the stored inputs nor projects normalization.
"""
    require(type(k) is Q and k>0,'Positive exact stored momentum required')
    require(type(epsilon) is Q and epsilon!=0,'Exact represented epsilon required')
    require(all(type(v) is Q for v in (*u,*w)),'Exact decoded components required')
    U=tuple(v/epsilon for v in u);W=tuple(v/epsilon for v in w)
    c=U[0]-W[1]/(2*k)
    return {'U_saved':U,'W_saved':W,'c':c,
            'first_order_wronskian_defect':2*epsilon*c,
            'kA_saved':(W[1]/2,-W[0]/2)}
