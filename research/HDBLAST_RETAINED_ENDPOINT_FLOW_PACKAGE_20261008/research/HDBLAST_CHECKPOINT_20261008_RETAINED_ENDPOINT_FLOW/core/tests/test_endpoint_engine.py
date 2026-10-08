"""Manufactured exact functions and mutation checks; never a physical source."""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'engine'))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'source'))
from endpoint_engine import EndpointFlow,geometry,SCALE,av,cq,rectangle,ENDPOINTS,ANCHOR,BITS
from later_source import outward_error
from flint import arb,acb,ctx

def rows(coeff):
    return [dict(g,coefficients=tuple(coeff)+tuple(Q(0) for _ in range(25-len(coeff))),uniform_error=Q(0)) for g in geometry()]

def polynomial_rows(k):
    out=[]
    # Exact solution W=x+i*k*x², U=x²/2+i*k*x³/3 with real
    # g=-1-2*k²*x². Varying cell centers test moment translation and signs.
    for g in geometry():
        a=g['center']-ANCHOR
        co=(Q(-1)-2*k*k*a*a,-4*k*k*a,-2*k*k)
        out.append(dict(g,coefficients=co+tuple(Q(0) for _ in range(22)),uniform_error=Q(0)))
    return out

def contained(rect,z):
    for axis, a in (('real',z.real),('imag',z.imag)):
        lo,hi=rect[axis]
        if not (av(Q(lo,SCALE)) <= a.lower() and a.upper() <= av(Q(hi,SCALE))):
            raise RuntimeError('independent closed-form outside rectangle')

def run():
    passed=[]
    # An exact constant forcing has an independent elementary flow. Including
    # tiny and zero k exercises the entire kernel without division by k.
    for endpoint in ENDPOINTS:
        target=EndpointFlow(rows([Q(3,7)]),endpoint)
        for k in (Q(0),Q(1,2**40),Q(1,4),Q(16),Q(128),Q(256)):
            u,w=(Q(2,3),Q(-1,5)),(Q(1,7),Q(4,9))
            got=target.evaluate(k,u,w)
            t=av(endpoint-ANCHOR);x=av(k);E=acb(0,2*x*t).exp()
            Phi=t*acb(0,x*t).exp()*(x*t).sinc()
            if k==0:
                N=t*t/2
            else:
                N=(Phi-t)/acb(0,2*x)
            contained(got['W'],E*cq(w)-av(Q(3,7))*Phi)
            contained(got['U'],cq(u)+Phi*cq(w)-av(Q(3,7))*N)
            passed.append(str(endpoint)+'/'+str(k))
        for k in (Q(0),Q(1,2**40),Q(1,4),Q(128),Q(256)):
            flow=EndpointFlow(polynomial_rows(k),endpoint)
            got=flow.evaluate(k,(Q(0),Q(0)),(Q(0),Q(0)))
            t=endpoint-ANCHOR
            contained(got['W'],cq((t,k*t*t)))
            contained(got['U'],cq((t*t/2,k*t*t*t/3)))
            passed.append('polynomial '+str(endpoint)+'/'+str(k))
        uncertain=rows([Q(0)])
        for r in uncertain:r['uniform_error']=Q(1,10**25)
        flow=EndpointFlow(uncertain,endpoint)
        got=flow.evaluate(Q(0),(Q(0),Q(0)),(Q(0),Q(0)))
        t=endpoint-ANCHOR
        contained(got['W'],cq((-Q(1,10**25)*t,Q(0))))
        contained(got['U'],cq((-Q(1,10**25)*t*t/2,Q(0))))
        passed.append('nonzero certified source uncertainty '+str(endpoint))
        odd_rows=rows([Q(0)])
        exact_W=Q(0);exact_U=Q(0)
        for j,r in enumerate(odd_rows):
            # Large, varying odd denominators model the exact-exp coefficient
            # error representation without constructing or sampling any source.
            exact=Q(1,10**1100+2*j+1);rounded=outward_error(exact)
            if not exact<=rounded<exact+Q(1,1<<512):raise RuntimeError('source error rounded inward')
            r['uniform_error']=rounded
            if r['right']<=endpoint:
                exact_W+=exact*Q(1,64);exact_U+=exact*Q(1,64)*(endpoint-r['center'])
        flow=EndpointFlow(odd_rows,endpoint)
        if flow.source_error_W<exact_W or flow.source_error_U<exact_U:raise RuntimeError('source uncertainty aggregation lost containment')
        json.dumps(flow.budget(),allow_nan=False)
        passed.append('large odd denominator source uncertainty serializes '+str(endpoint))
    for mutate in ('missing cell','wrong center','float coefficient','negative uncertainty','precision changed'):
        candidate=rows([Q(0)])
        try:
            if mutate=='missing cell':candidate.pop()
            elif mutate=='wrong center':candidate[0]['center']+=Q(1,128)
            elif mutate=='float coefficient':candidate[0]['coefficients']=(0.0,)+candidate[0]['coefficients'][1:]
            elif mutate=='negative uncertainty':candidate[0]['uniform_error']=Q(-1)
            if mutate=='precision changed':
                target=EndpointFlow(candidate,ENDPOINTS[0]);ctx.prec=BITS-1
                target.evaluate(Q(1),(Q(0),Q(0)),(Q(0),Q(0)))
            else:EndpointFlow(candidate,ENDPOINTS[0])
        except ValueError:passed.append('reject '+mutate)
        else:raise RuntimeError('mutation accepted: '+mutate)
    print(json.dumps({'status':'PASS_MANUFACTURED_ENDPOINT_CLOSED_FORM','checks':len(passed),'cases':passed,
                      'physical_source_evaluations':0,'retained_array_reads':0},sort_keys=True))

if __name__=='__main__':run()
