"""Two evaluations of the electric Weyl curvature of a time-symmetric closed slice."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.special import eval_chebyu, eval_gegenbauer

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("closed_constraint", ROOT/"initial/constraint.py")
constraint = importlib.util.module_from_spec(spec)
spec.loader.exec_module(constraint)


def eigenvalue_difference(coefficients, amplitude, offset, branch, y):
    n = 2*np.arange(len(coefficients))
    y = np.asarray(y)
    sine2 = 1-y*y
    psi = eval_chebyu(n, y[:, None]) @ coefficients
    psi_y = (2*eval_gegenbauer(n[1:]-1, 2, y[:, None])) @ coefficients[1:]
    psi_yy = (8*eval_gegenbauer(n[1:]-2, 3, y[:, None])) @ coefficients[1:]
    p = psi_y/psi
    grad2 = 64*amplitude**2*y*y*sine2
    # Difference between the radial and either tangential mixed electric eigenvalue, in h^2 units.
    difference = (sine2*(-2*psi_yy/psi+6*p*p)-grad2/2)/psi**4
    q = amplitude*(4*y*y-1)+offset
    u, _ = branch.potential(q)
    # The 4D angular sectional curvature identity uses the areal radius and a single derivative.
    angular_weyl = (1+4*y*p-4*sine2*p*p+grad2/6)/psi**4-u/branch.u_top
    return difference, -1.5*angular_weyl, psi


def calculate(previous):
    chi, w = constraint.quadrature(256)
    y = np.cos(chi)
    branches = {}
    for label, branch in constraint.load_branches(previous).items():
        beta = branch.mu2/branch.u_top
        leading = (1-y*y)*(4*(beta-24)/35+16*(beta-10)*y*y/7)
        analytic_rms = float(np.sqrt(np.dot(w,leading**2)))
        series = {}
        for a in (0.001,0.0005):
            rows, profiles = [], []
            for sign in (1,-1):
                amp = sign*a
                coeff = constraint.solve_slice(branch, amp)
                diff, independent, psi = eigenvalue_difference(coeff, amp, 0, branch, y)
                residual = float(np.max(np.abs(diff-independent)))
                if residual>1e-12:
                    raise ArithmeticError("Electric and angular Weyl curvature disagree")
                w_proper = w*psi**6
                e2 = float(np.dot(w_proper,diff**2)/w_proper.sum())
                rows.append(dict(amplitude=amp,eigenvalue_difference_rms_over_h2=np.sqrt(e2),
                                 mean_weyl_squared_over_h4=(16./3)*e2,
                                 max_two_formula_difference=residual))
                profiles.append(diff)
            averaged = (profiles[0]+profiles[1])/(2*a*a)
            error = float(np.max(np.abs(averaged-leading))/np.max(np.abs(leading)))
            if error>5e-4:
                raise ArithmeticError("Signed profile inconsistent with quadratic coefficient")
            series[str(a)] = dict(signed=rows,even_profile_relative_error=error)
        z, independent, _ = eigenvalue_difference(constraint.solve_slice(branch,0),0,0,branch,y)
        if np.any(z!=0) or np.max(abs(independent))>2e-15:
            raise ArithmeticError("Homogeneous control has tidal curvature")
        branches[label] = dict(beta=beta,analytic_tidal_rms_over_h2_A2=analytic_rms,
                               analytic_mean_weyl_squared_over_h4_A4=16*analytic_rms**2/3,
                               amplitudes=series)
    return dict(description="Electric Weyl curvature on the time-symmetric initial slices",branches=branches)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True,help="Output JSON file")
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.previous),indent=2,allow_nan=False)+"\n")
