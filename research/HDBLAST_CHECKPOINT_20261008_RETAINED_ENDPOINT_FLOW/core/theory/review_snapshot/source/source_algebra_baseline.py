"""Independent exact-rational analytic models; definitions only before freeze.

No module-level registered source evaluation occurs. The registered builder
requires a root-authenticated prospective freeze receipt; this module does not
fetch a registration or substitute for the root entry-point integrity guard.
"""
from fractions import Fraction as Q
from math import factorial
import re


def exp_rational_interval(q,terms=200):
    """Exact alternating enclosure for exp(q), -1≤q≤0."""
    if not isinstance(q,Q) or not -1<=q<=0 or terms<2 or terms%2:
        raise ValueError('Exact q in [-1,0] and positive even degree required')
    value=Q(1);term=Q(1)
    for n in range(1,terms+1):
        term=term*q/n;value+=term
    upper=value
    lower=value+term*q/(terms+1)
    if lower<=0 or upper<lower:
        raise ValueError('Invalid exact exponential enclosure')
    return lower,upper


def normalized_exp_coefficients(f,count):
    """Coefficients of exp(f(x)-f(0)), exact rational formal algebra."""
    if len(f)<count+1:
        raise ValueError('Full formal exponent coefficient inventory required')
    beta=[Q(1)]
    for n in range(1,count+1):
        beta.append(sum((j*f[j]*beta[n-j] for j in range(1,n+1)),Q(0))/n)
    return tuple(beta)


def reciprocal_coefficients(d,count):
    if d[0]==0:
        raise ValueError('Nonzero formal denominator required')
    out=[Q(1)/d[0]]
    for n in range(1,count+1):
        out.append(-sum((d[j]*out[n-j] for j in range(1,min(n,len(d)-1)+1)),Q(0))/d[0])
    return tuple(out)


def dyadic_point(q,bits=512):
    return Q((q.numerator*(1<<bits))//q.denominator,1<<bits)


def _authorize(receipt):
    if not isinstance(receipt,dict) or receipt.get('status')!='PASS_REMOTE_REGISTERED_SOURCE_GO':
        raise ValueError('New independently read-back remote freeze required before source evaluation')
    if not re.fullmatch(r'[0-9a-f]{40}',receipt.get('freeze_commit','')):
        raise ValueError('Exact freeze commit required')
    if not re.fullmatch(r'[0-9a-f]{64}',receipt.get('registration_sha256','')):
        raise ValueError('Exact registration hash required')
    if receipt.get('input_scope')!='NO_RETAINED_ARRAYS_SOURCE_OPERATOR_ONLY':
        raise ValueError('Narrow source-only input scope required')


def registered_source_model(center,source,authorization,degree=24,bits=512):
    """Fresh own source coefficients and enclosure; only call after new freeze.

    Exp is enclosed by an exact rational alternating series. All other
    normalized coefficient algebra is exact. The only approximation in the
    chosen polynomial is an explicit dyadic point/radius operation.
    """
    _authorize(authorization)
    if degree!=24 or bits!=512:
        raise ValueError('Candidate degree24 / coefficient512bits required')
    if source not in ('positive_B','signed_uB'):
        raise ValueError('Registered real source identifier required')
    if not isinstance(center,Q) or not Q(-9,2)<=center<=Q(-7,2):
        raise ValueError('Exact audited negative source center required')
    z=center+4
    d=(1-z*z,-2*z,Q(-1))
    inv=reciprocal_coefficients(d,degree+2)
    exponent=(1-inv[0],)+tuple(-a for a in inv[1:])
    b0lo,b0hi=exp_rational_interval(exponent[0])
    beta=normalized_exp_coefficients(exponent,degree+2)
    if source=='positive_B':
        h=beta
    else:
        h=tuple(z*beta[n]+(beta[n-1] if n else 0) for n in range(degree+3))
    # L=-1/(center+x), not a frozen geometry value.
    ls=tuple(-(-Q(1))**n/center**(n+1) for n in range(degree+1))
    l2=tuple(sum((ls[j]*ls[n-j] for j in range(n+1)),Q(0)) for n in range(degree+1))
    gamma=[]
    for n in range(degree+1):
        gamma.append(4*sum((l2[j]*h[n-j] for j in range(n+1)),Q(0))
                     -2*sum((ls[j]*(n-j+1)*h[n-j+1] for j in range(n+1)),Q(0))
                     -(n+1)*(n+2)*h[n+2])
    work_gamma=tuple(sum((ls[j]*gamma[n-j] for j in range(n+1)),Q(0))
                     for n in range(degree+1))
    half=Q(1,128)
    def choose(multiplier,analytic_remainder):
        point=[];radii=[]
        for a in multiplier:
            lo,hi=(a*b0lo,a*b0hi) if a>=0 else (a*b0hi,a*b0lo)
            chosen=dyadic_point((lo+hi)/2,bits)
            point.append(chosen);radii.append(max(chosen-lo,hi-chosen))
        error=analytic_remainder+sum((r*half**n for n,r in enumerate(radii)),Q(0))
        return {'center_coefficients':tuple(point),
                'coefficient_radii':tuple(radii),
                'uniform_source_error':error,
                'analytic_tail':analytic_remainder}
    return {'forcing':choose(gamma,Q(1,15*2**90)),
            'Lg':choose(work_gamma,Q(1,15*2**91)),
            'center':center,'source':source,'half_width':half,
            'degree':degree,'coefficient_bits':bits,
            'coefficient_algorithm':'own exact rational inverse/exp/geometry recurrences',
            'exp_algorithm':'own exact alternating rational series degree200'}
