"""Mechanical NumPy-dispatch and no-mode counterterm preflight after run 001.

This does not evolve physical modes or change the registered numerical matrix.
It prevents recurrence of the ndarray-left operator dispatch failure.
"""
import argparse
import json
from pathlib import Path

import numpy as np

import smooth_frw_control as control


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    j=control.Jet.constant(np.ones((2,1)))
    array=np.arange(1,4,dtype=float)[None,:]
    operations={"array_plus_jet":array+j,"array_minus_jet":array-j,
                "array_times_jet":array*j,"array_divided_by_jet":array/j}
    dispatch={name:isinstance(result,control.Jet) and result.c.shape==(6,2,3)
              for name,result in operations.items()}
    assert all(dispatch.values()),dispatch
    k=np.array([0,.5,2,48,96])
    t=np.array([0,.125,.25,.375,.5,.625,.75,.875,1,1.25])
    cases=[]
    for amplitude in (0.,.2):
        ct=control.counterterms(k,t,amplitude)
        residual=ct['rho_prime']+3*ct['H'][:,None]*(ct['rho']+ct['p'])-ct['xp'][:,None]*ct['Q']/2
        maximum=float(np.max(abs(residual)))
        finite=all(np.all(np.isfinite(value)) for value in ct.values())
        shapes=all(ct[name].shape==(len(t),len(k)) for name in ('rho','p','Q','Q_prime','Q_second'))
        assert finite and shapes and maximum<2e-10,(amplitude,maximum)
        cases.append(dict(amplitude=amplitude,finite=finite,shapes_valid=shapes,
                          complete_subtraction_exchange_residual=maximum))
    static=control.counterterms(k,t,0.,static=True)
    w=np.sqrt(np.asarray(k,dtype=np.longdouble)**2+control.R)[None,:]
    expected={'rho':w/2,'p':k[None,:]**2/(6*w),'Q':1/(2*w)}
    static_errors={name:float(np.max(abs(static[name]-value))) for name,value in expected.items()}
    assert max(static_errors.values())<1e-15,static_errors
    a,x=control.background_jets(t,0.)
    h=x.c[0,:,None]-control.R
    hpp=2*x.c[2,:,None]
    expected_flat={
        'rho':w/2+h/(4*w)-h*h/(16*w**3),
        'p':k[None,:]**2/(6*w)-h*k[None,:]**2/(12*w**3)+h*h*k[None,:]**2/(16*w**5)
            -hpp*(2*k[None,:]**2/3+control.R)/(16*w**5),
        'Q':1/(2*w)-h/(4*w**3),
    }
    flat=control.counterterms(k,t,0.)
    flat_errors={name:float(np.max(abs(flat[name]-value))) for name,value in expected_flat.items()}
    assert max(flat_errors.values())<2e-14,flat_errors
    result=dict(status='PASS',physical_mode_evolution_executed=False,numpy_array_dispatch=dispatch,
                local_subtraction_cases=cases,static_errors=static_errors,flat_matching_errors=flat_errors,
                limitation='Selected numerical preflight cases; no new convergence or physical source result.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,allow_nan=False))


if __name__=='__main__':
    main()
