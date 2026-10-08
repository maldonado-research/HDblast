"""New-scope authenticated later-source construction; no source at import.

The byte-identical upstream helper supplies only exact algebra/exp primitives.
Its old physical entry and old scope-specific GO are never called here.
"""
from fractions import Fraction as Q
import hashlib
import json
from source_algebra_baseline import reciprocal_coefficients, normalized_exp_coefficients, exp_rational_interval

SOURCES = ('positive_B','signed_uB')
DEGREE = 24
BITS = 512
HALF = Q(1,128)
RADIUS = Q(1,8)
CAUCHY_MAJORANT = Q(64)
TAIL = CAUCHY_MAJORANT/Q(15*16**24)

def require(ok, text):
    if not ok:
        raise ValueError(text)

def canonical(q):
    require(type(q) is Q, 'exact Fraction required')
    return f'{q.numerator}/{q.denominator}'

def coefficient_sha(co):
    return hashlib.sha256(json.dumps([canonical(v) for v in co],separators=(',',':')).encode()).hexdigest()

def outward_error(value):
    require(type(value) is Q and value>=0,'nonnegative exact source error required')
    scale=1<<BITS
    return Q((value.numerator*scale+value.denominator-1)//value.denominator,scale)

def source_cell(center, source, authorization):
    require(callable(authorization), 'authenticated new-scope source callback required')
    require(source in SOURCES and type(center) is Q and Q(-9,2)<center<Q(-7,2),
            'registered exact source/center required')
    authorization('later_source_taylor', {'source':source, 'center':canonical(center),
                  'scope':'LATER_ENDPOINT_SOURCE64CELLS_DEGREE24'})
    z = center+4
    inv = reciprocal_coefficients((1-z*z,-2*z,Q(-1)),DEGREE+2)
    exponent = (1-inv[0],)+tuple(-a for a in inv[1:])
    lo, hi = exp_rational_interval(exponent[0],200)
    beta = normalized_exp_coefficients(exponent,DEGREE+2)
    h = beta if source == SOURCES[0] else tuple(z*beta[n]+(beta[n-1] if n else 0)
                                                for n in range(DEGREE+3))
    ls = tuple(-(-Q(1))**n/center**(n+1) for n in range(DEGREE+1))
    l2 = tuple(sum((ls[j]*ls[n-j] for j in range(n+1)),Q(0)) for n in range(DEGREE+1))
    gamma = tuple(4*sum((l2[j]*h[n-j] for j in range(n+1)),Q(0))
                  -2*sum((ls[j]*(n-j+1)*h[n-j+1] for j in range(n+1)),Q(0))
                  -(n+1)*(n+2)*h[n+2] for n in range(DEGREE+1))
    chosen = []; discrepancy = Q(0); scale = 1<<BITS
    for n, v in enumerate(gamma):
        a,b = (v*lo,v*hi) if v >= 0 else (v*hi,v*lo)
        middle = (a+b)/2
        c = Q(middle.numerator*scale//middle.denominator,scale)
        chosen.append(c)
        discrepancy += max(c-a,b-c)*HALF**n
    raw_error=TAIL+discrepancy;error=outward_error(raw_error)
    return tuple(chosen),error, {'analytic_tail':canonical(TAIL),
               'coefficient_error':canonical(discrepancy),'source_error_rounding':canonical(error-raw_error),
               'source_error':canonical(error),
               'coefficient_sha256':coefficient_sha(chosen),
               'exp_degree':200,'coefficient_bits':BITS}

def manufactured_cell(center, source):
    # Deliberately all25 coefficients nonzero; bounded normalized coefficients
    # exercise the same degree/conditioning without constructing a real source.
    sid = SOURCES.index(source)
    co = tuple(Q((-1)**(n+sid),n+1)/(2*HALF)**n for n in range(DEGREE+1))
    return co,Q(0),{'analytic_tail':'0/1','coefficient_error':'0/1','source_error_rounding':'0/1','source_error':'0/1',
                   'coefficient_sha256':coefficient_sha(co),'manufactured':True}

def build_models(authorization=None, manufactured=False, check=lambda:None):
    require(type(manufactured) is bool, 'explicit source mode required')
    require(manufactured or callable(authorization), 'new source authorization required')
    models = {}; proof = []
    for source in SOURCES:
        rows = []
        for j in range(64):
            check()
            center = Q(-9,2)+Q(2*j+1,128)
            co, error, evidence = manufactured_cell(center,source) if manufactured else source_cell(center,source,authorization)
            row = {'center':center,'half':HALF,'coefficients':co,'uniform_error':error}
            rows.append(row)
            proof.append({'source':source,'cell':j,'center':canonical(center),'half':canonical(HALF),
                          'radius':canonical(RADIUS),'cauchy_majorant':canonical(CAUCHY_MAJORANT),
                          'coefficients':[canonical(v) for v in co],**evidence})
        models[source] = rows
    return models, {'schema_version':1,'manufactured':manufactured,'source_count':2,
                     'cell_count_per_source':64,'degree':DEGREE,'rows':proof}
