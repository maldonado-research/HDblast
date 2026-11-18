"""Check the time ordering and order of the propagator against a power series."""
import json
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('spectator',HERE/'spectator.py')
mod=module_from_spec(spec)
spec.loader.exec_module(mod)


def solution(a0,a1):
    # y''+(1+t)y=0, Taylor coefficients at the origin.
    coefficients=[Fraction(a0),Fraction(a1)]
    for j in range(48):
        val=coefficients[j]+(coefficients[j-1] if j else 0)
        coefficients.append(-val/((j+1)*(j+2)))
    return float(sum(coefficients)),float(sum(i*a for i,a in enumerate(coefficients)))


def matrix_at(steps,method):
    k=np.zeros(1)
    matrix=(np.ones(1),np.zeros(1),np.zeros(1),np.ones(1))
    for j in range(steps):
        def data(alpha):
            omega=np.sqrt(1+(j+alpha)/steps)
            return mod.State(1.,1.,1.,omega,0.,0.,omega)
        nodes=(data(.5),data(.5)) if method=='midpoint' else (data(mod.MINUS),data(mod.PLUS))
        matrix=mod.propagate(nodes,k,matrix,1/steps,method=='midpoint')
    return np.array([v[0] for v in matrix])


if __name__=='__main__':
    a,c=solution(1,0)
    b,d=solution(0,1)
    exact=np.array([a,b,c,d])
    modes={}
    for method in ('cf4','midpoint'):
        errors=[float(max(abs(matrix_at(n,method)-exact))) for n in (8,16,32,64)]
        ratios=[x/y for x,y in zip(errors,errors[1:])]
        if min(ratios)<(14 if method=='cf4' else 3.8):
            raise ValueError('Convergence benchmark failed: '+method)
        modes[method]={'errors':errors,'reduction':ratios}
    print(json.dumps(modes,indent=2))
