"""Independent quadrature of a coherence bound at fixed initial pressure and force."""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize

RELEASED_MODULES = {
    'squeeze.py': 'c15cdbccacf161493de964280054363d8b873e63bae4d5b549eb45e287dbfdc9',
    'transfer.py': '9bac2426b5633f7159859aa1844fc6889e2dae5a30f0fcf060fde50443eb1fa8',
}


def load_reference(repo):
    folder=repo/'research/HDBLAST_MOVING_PARENT_20261119/independent'
    for name, expected in RELEASED_MODULES.items():
        path=folder/name
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError(f"Mismatch in released moving-parent helper: {path}")
    spec=importlib.util.spec_from_file_location("positive_squares", folder/'squeeze.py')
    reference=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    return reference


def bounds(q, moments, initial_a, initial_mass, labels, order=12):
    # Composite Gauss quadrature of smooth response splines; record unit-radius support.
    x, w = np.polynomial.legendre.leggauss(order)
    left, width = q[:-1,None], np.diff(q)[:,None]
    k = (left + width*(x+1)/2).ravel()
    weight = (width*w/2).ravel()*3*k**2/8**3
    interp = CubicSpline(q,moments,axis=1)(k)
    d = weight@interp[labels.index('incoherent_excess')]
    absolute = d + weight@interp[labels.index('vacuum')]
    K = -(interp[labels.index('cross_real')]+1j*interp[labels.index('cross_imag')])/(np.sqrt(2)*d)
    energies = np.hypot(k/initial_a, initial_mass)
    h = np.stack((initial_mass/energies, (k/8)**2*initial_mass/energies))

    def fn(lam):
        residual = K - lam@h
        magnitude = abs(residual)
        slope = -(h@(weight*residual.real/magnitude))
        return weight@magnitude,slope

    # Any trial lambda is an upper certificate. A small dual slope checks tightness.
    opt = minimize(fn, [0.,0.], method='BFGS', jac=True, options={'gtol':1.e-12})
    upper, slope = fn(opt.x)
    # Project tiny quadrature residuals away, then restore the unit-disk radius.
    residual = K - opt.x@h
    witness = residual.conj()/abs(residual)
    gram = (h*weight)@h.T
    correction = np.linalg.solve(gram, h@(weight*witness.real))@h
    witness = (witness-correction)/(1+np.max(abs(correction)))
    attained = (weight@(K*witness)).real
    source_residual = np.max(abs(h@(weight*witness.real)))
    if np.max(abs(witness))>1+2e-14 or source_residual>5e-14 or upper-attained>3e-7:
        raise ArithmeticError('Primal and dual certificates disagree')
    # Split the absolute value at the spline zeros, integrating each signed
    # polynomial times q² exactly (up to floating-point roundoff).
    spline=CubicSpline(q, -moments[labels.index('cross_imag')]/(np.sqrt(2)*d))
    zeros=spline.roots(extrapolate=False)
    knots=np.unique(np.r_[q,zeros[np.isfinite(zeros)]])
    imag_width=np.diff(knots)[:,None]
    imag_points=knots[:-1,None]+imag_width*(x+1)/2
    imag_only=np.sum(3*imag_points**2/8**3*(imag_width*w/2)*abs(spline(imag_points)))
    return dict(diagonal_excess_contrast=d, diagonal_absolute_contrast=absolute,
                normalization='negative anomalous contrast divided by signed diagonal contrast',
                separable_bound=upper, quantum_bound=np.sqrt(2)*upper,
                separable_witness=attained, quantum_witness=np.sqrt(2)*attained,
                max_moment_residual=float(source_residual),
                imag_separable=imag_only, imag_quantum=np.sqrt(2)*imag_only,
                unconstrained_separable=weight@abs(K),
                separable_absolute_bound=upper*abs(d/absolute),
                quantum_absolute_bound=np.sqrt(2)*upper*abs(d/absolute),
                lambda_dimensionless=opt.x.tolist(),
                max_dual_gradient=float(max(abs(slope))),
                max_source_moment_discrepancy_aligned=float(max(abs(slope))))


def run(repo):
    sq=load_reference(repo)
    base=sq.base
    ref=json.loads(base.pinned(repo/base.REF,base.REF_SHA))['matter']['events'][0]['A_star']
    result = {}
    with zipfile.ZipFile(io.BytesIO(base.pinned(repo/base.WALL,base.WALL_SHA))) as archive:
        for steps, intervals in ((8192,32),(16384,64),(32768,128)):
            q=np.linspace(0,8,intervals+1)
            integrals={}
            for arm in sq.BRANCHES:
                wall=base.historical(archive,arm,ref)
                curves=base.modes(wall,q,steps)
                psi=sq.integrated(curves,q,base.WINDOWS['wide'],16)
                integrals[arm]=psi['full']['per_q'][0]
                initial_a=curves[2][0]
                initial_mass=abs(wall[2](curves[1][0]))
                del curves
            rows = {}
            for arm in sq.BRANCHES[1:]:
                diff=integrals[arm]-integrals['quiet']
                rows[arm]=bounds(q,diff,initial_a,initial_mass,sq.OUTPUTS)
                refined=bounds(q,diff,initial_a,initial_mass,sq.OUTPUTS,order=24)
                rows[arm]['max_factor_change_12_to_24']=max(
                    abs(rows[arm][name]-refined[name]) for name in
                    ('separable_bound','quantum_bound','separable_absolute_bound',
                     'quantum_absolute_bound','imag_separable','imag_quantum'))
            result[str(intervals)] = rows
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository',type=Path)
    print(json.dumps(run(parser.parse_args().repository),indent=2,allow_nan=False))
