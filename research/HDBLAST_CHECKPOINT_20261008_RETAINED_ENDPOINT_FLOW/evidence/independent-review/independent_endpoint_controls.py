"""Exact fabricated IVPs; no retained file or physical source construction.

Prescribe a polynomial solution, derive its real source by differentiation,
then check the production propagator against exact rational endpoint values.
This oracle uses neither exponential/moment expansion nor numerical integration.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import hashlib
import json
import sys
import time

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parent
ROOT = BASE.parent/'implementation'
sys.path[:0] = [str(ROOT/'engine')]
from endpoint_engine import EndpointFlow, geometry, SCALE, ENDPOINTS, ANCHOR
import endpoint_aggregate as aggregate

def require(ok, text):
    if not ok: raise RuntimeError(text)
def evaluate(co,x): return sum((v*x**j for j,v in enumerate(co)),Q(0))
def derivative(co): return tuple((j+1)*v for j,v in enumerate(co[1:]))
def primitive(co): return (Q(0),)+tuple(v/(j+1) for j,v in enumerate(co))
def to_local(co,offset):
    return tuple(sum((co[n]*comb(n,j)*offset**(n-j) for n in range(j,len(co))),Q(0))
                 for j in range(25))
def contained(rect,z):
    return all(Q(rect[a][0],SCALE)<=z[j]<=Q(rect[a][1],SCALE)
               for j,a in enumerate(('real','imag')))
def rows_for(co):
    return [dict(g,coefficients=to_local(co,g['center']-ANCHOR),uniform_error=Q(0))
            for g in geometry()]
def exact_rectangle(z):
    require(all((v*SCALE).denominator==1 for v in z),'dyadic exact rectangle required')
    return {axis:[int(v*SCALE),int(v*SCALE)] for axis,v in zip(('real','imag'),z)}

def run():
    started=time.monotonic();checks=[]
    # P contains every power through24, so signs, top degree, factorial order,
    # lag translation, integration Jacobian and incoming transport all matter.
    P=tuple(Q((-1)**j,(j+1)*2**j) for j in range(25))
    P1=derivative(P);P2=derivative(P1);I=primitive(P)
    C=Q(3,11);D=Q(-5,13)
    for k in (Q(0),Q(3,7),Q(256)):
        source=tuple(-4*k*k*P[j]-(P2[j] if j<len(P2) else 0)-(2*k*C if j==0 else 0)
                     for j in range(25))
        incoming_u=(P[0],D);incoming_w=(P1[0],2*k*P[0]+C)
        for endpoint in ENDPOINTS:
            flow=EndpointFlow(rows_for(source),endpoint)
            x=endpoint-ANCHOR
            expected_u=(evaluate(P,x),2*k*evaluate(I,x)+C*x+D)
            expected_w=(evaluate(P1,x),2*k*evaluate(P,x)+C)
            result=flow.evaluate(k,incoming_u,incoming_w)
            require(contained(result['U'],expected_u),'prescribed polynomial U escaped')
            require(contained(result['W'],expected_w),'prescribed polynomial W escaped')
            # A sign-inverted forcing gives a distinct endpoint for this IVP;
            # it cannot also lie in the tiny reported certified rectangle.
            if k==0:
                wrong_u=(2*(incoming_u[0]+x*incoming_w[0])-expected_u[0],expected_u[1])
                require(not contained(result['U'],wrong_u),'forcing-sign counterexample not separated')
            checks.append({'case':'prescribed_degree24_real_forcing','k':str(k),'endpoint':str(endpoint)})
    # A jump at every cell boundary and independently varying high coefficients
    # expose errors hidden by a single smooth global forcing at k=0.
    cells=[]
    for j,g in enumerate(geometry()):
        co=tuple(Q((-1)**(j+n)*(j+1),n+3)*128**n for n in range(25))
        cells.append(dict(g,coefficients=co,uniform_error=Q(0)))
    for endpoint in ENDPOINTS:
        u=w=Q(0)
        for row in cells:
            if row['right']>endpoint:continue
            h=row['half'];lag=endpoint-row['center']
            integral=sum((a*(h**(n+1)-(-h)**(n+1))/(n+1) for n,a in enumerate(row['coefficients'])),Q(0))
            first=sum((a*(h**(n+2)-(-h)**(n+2))/(n+2) for n,a in enumerate(row['coefficients'])),Q(0))
            w-=integral;u-=lag*integral-first
        result=EndpointFlow(cells,endpoint).evaluate(Q(0),(Q(0),Q(0)),(Q(0),Q(0)))
        require(contained(result['U'],(u,Q(0))) and contained(result['W'],(w,Q(0))),
                'piecewise high-degree exact integral escaped')
        checks.append({'case':'discontinuous_degree24_zero_k','endpoint':str(endpoint)})
    # Uniform uncertainty at zero k is attained by a same-sign piecewise real
    # source. Check its full mass and first moment, including the prefix split.
    uncertain=[dict(g,coefficients=(Q(0),)*25,uniform_error=Q(j+1,2**100))
               for j,g in enumerate(geometry())]
    for endpoint in ENDPOINTS:
        flow=EndpointFlow(uncertain,endpoint)
        ew=sum((2*r['half']*r['uniform_error'] for r in uncertain if r['right']<=endpoint),Q(0))
        eu=sum((2*r['half']*(endpoint-r['center'])*r['uniform_error'] for r in uncertain if r['right']<=endpoint),Q(0))
        require(flow.source_error_U==eu and flow.source_error_W==ew,'source uncertainty mass/first moment differs')
        result=flow.evaluate(Q(0),(Q(0),Q(0)),(Q(0),Q(0)))
        for sign in (-1,1):
            require(contained(result['U'],(sign*eu,Q(0))) and contained(result['W'],(sign*ew,Q(0))),
                    'attainable source uncertainty escaped')
        checks.append({'case':'source_uncertainty_attainment','endpoint':str(endpoint)})
    broad=[dict(g,coefficients=(Q(0),)*25,uniform_error=Q(1,10**6)) for g in geometry()]
    try:EndpointFlow(broad,ENDPOINTS[0]).evaluate(Q(0),(Q(0),Q(0)),(Q(0),Q(0)))
    except ValueError as error:require('radius exceeds' in str(error),'unexpected gate rejection')
    else:raise RuntimeError('overwide uncertainty silently passed fixed gate')
    checks.append({'case':'overwide_complete_export_rejected'})
    # Independent direct stress function checks correlated stress algebra and
    # proves nonzero canonical drift is retained rather than projected away.
    for k in (Q(1,8),Q(2),Q(255)):
        for endpoint in ENDPOINTS:
            target_u=(Q(3,8),Q(5,16));target_w=(Q(7,32),Q(1,4))
            incoming_u=(target_u[0],Q(0));incoming_w=(Q(0),target_w[1])
            saved_u=(Q(-1,4),Q(3,8));saved_w=(Q(-5,16),Q(7,8));weight=Q(3,64)
            result=aggregate.node_terms(k,weight,endpoint,saved_u,saved_w,
                {'U':exact_rectangle(target_u),'W':exact_rectangle(target_w)},incoming_u,incoming_w)
            L=-1/endpoint;mu=weight*k*k/(2*aggregate.PI*aggregate.PI)
            def stress(u,w):
                return ((k+3*L*L/(2*k))*u[0]-L*w[0]/(2*k)-w[1]/2,
                        (k/3-L*L/(2*k))*u[0]-L*w[0]/(2*k)-w[1]/2)
            expected=tuple(mu*(a-b) for a,b in zip(stress(saved_u,saved_w),stress(target_u,target_w)))
            require(result['R_error']==(expected[0],)*2 and result['P_error']==(expected[1],)*2,
                    'direct and canonical stress disagree')
            require(result['canonical_drift']!=0,'manufactured nonzero drift vanished')
            checks.append({'case':'exact_stress_and_nonzero_canonical_drift','k':str(k),'endpoint':str(endpoint)})
    bad={'U':exact_rectangle((Q(0),Q(0))),'W':exact_rectangle((Q(0),Q(0)))}
    try:aggregate.node_terms(Q(1),Q(1),ENDPOINTS[0],(Q(0),Q(0)),(Q(0),Q(0)),bad,(Q(1),Q(0)),(Q(0),Q(0)))
    except ValueError as error:require('disjoint' in str(error),'unexpected canonical rejection')
    else:raise RuntimeError('impossible canonical target accepted')
    checks.append({'case':'inconsistent_canonical_target_rejected'})
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
          for p in (ROOT/'engine/endpoint_engine.py',ROOT/'engine/endpoint_aggregate.py')}
    result={'status':'PASS_INDEPENDENT_EXACT_FABRICATED_CONTROLS','checks':len(checks),'cases':checks,
            'source_pins':pins,'wall_seconds':time.monotonic()-started,'retained_decodes':0,'physical_source_calls':0}
    output=BASE/('EXACT_CONTROLS_OPTIMIZED.json' if not __debug__ else 'EXACT_CONTROLS_NORMAL.json')
    with output.open('x') as out:json.dump(result,out,sort_keys=True,indent=2);out.write('\n')
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':run()
