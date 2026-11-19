"""Time-symmetric radial Einstein–scalar constraint on the unit three-sphere."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import eval_chebyu


class Branch:
    def __init__(self, earlier, fitted, top, h, mu2):
        self.earlier, self.fitted = earlier, fitted
        self.top, self.h, self.mu2 = top, h, mu2
        self.u_top, self.du_top = earlier.potential_and_slope(fitted, top)

    def potential(self, q):
        """Corrected total u and Klein–Gordon force -3 du/dq, at offset q."""
        q = np.asarray(q)
        u, slope = self.earlier.potential_and_slope(self.fitted, self.top + q.reshape(-1))
        return u.reshape(q.shape) - self.du_top*q, -3*(slope.reshape(q.shape) - self.du_top)


def load_branches(previous):
    previous = Path(previous)
    spec = importlib.util.spec_from_file_location('previous_zonal', previous/'nonlinear/nonlinear.py')
    earlier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(earlier)
    hilltop, model, refs, tables = earlier.prior_model(previous/'prior')
    branches = {}
    for pm in (0.5, -0.5):
        key = f'phi_m={pm}'
        fitted = hilltop.BG.Potential(model, tables, pm, refs[key])
        top, h, mu2, _ = hilltop.field_top(fitted)
        branches[key] = Branch(earlier, fitted, top, h, mu2)
    return branches


def quadrature(count):
    chi = np.pi*np.arange(1, count+1)/(count+1)
    return chi, 2*np.sin(chi)**2/(count+1)


def basis(degrees, chi):
    return eval_chebyu(degrees, np.cos(np.asarray(chi))[..., None])


def evaluate(coefficients, chi):
    return basis(2*np.arange(len(coefficients)), chi) @ coefficients


def source(branch, amplitude, chi, offset=0.):
    q = offset + amplitude*(4*np.cos(chi)**2-1)
    g2 = (4*amplitude*np.sin(2*chi))**2
    total, _ = branch.potential(q)
    return g2, (total - branch.u_top)/branch.u_top


def residual(deviation, minus_lap, g2, du):
    # Evaluate the exactly rearranged residual about psi=1; naive subtraction
    # of several O(6) terms loses digits in small-amplitude convergence tests.
    s = deviation
    psi5_minus_psi = s*(4+s*(10+s*(10+s*(5+s))))
    return 8*minus_lap - 6*(psi5_minus_psi+(1+s)**5*du) - (1+s)*g2


def solve_slice(branch, amplitude, max_degree=24, nodes=128, offset=0.):
    """Galerkin–Newton solution on the near-round reflection-even branch."""
    if max_degree < 4 or max_degree % 2 or nodes <= max_degree:
        raise ValueError('Use an even degree >= 4 and more quadrature points than the maximum degree')
    degrees = np.arange(0, max_degree+1, 2)
    k2 = degrees*(degrees+2)
    chi, weights = quadrature(nodes)
    U = basis(degrees, chi)
    project = U.T*weights
    g2, du = source(branch, amplitude, chi, offset)
    c = np.zeros(len(degrees))
    for _ in range(15):
        s = U@c
        rhs = project @ residual(s, U@(k2*c), g2, du)
        if np.max(abs(rhs)) < 2e-16:
            c[0] += 1
            if c[0] - (degrees[1:]+1) @ abs(c[1:]) <= 0:
                raise ArithmeticError('Nonpositive conformal solution')
            return c
        derivative = 6 - 30*(1+s)**4*(1+du) - g2
        jacobian = np.diag(8*k2) + project @ (derivative[:, None]*U)
        c -= np.linalg.solve(jacobian, rhs)
    raise ArithmeticError('Constraint Newton iteration failed')


def second_order(branch):
    k2 = np.array([0, 8, 24])
    mass_ratio = branch.mu2/branch.u_top
    return (8 - k2/2 - mass_ratio)/(8*k2-24)


def centered_offset(branch, amplitude, max_degree=24, nodes=128):
    """Constant b near zero solving <b+A U2>_proper=0 with the metric re-solved."""
    if amplitude == 0:
        return 0.
    chi, weights = quadrature(383)
    U2 = 4*np.cos(chi)**2-1
    def mean(offset):
        c = solve_slice(branch, amplitude, max_degree, nodes, offset=offset)
        # The round integral of U2 vanishes exactly. Evaluate psi^6-1 using
        # log1p and expm1 to avoid subtracting two values close to unity.
        c[0] -= 1
        s = evaluate(c, chi)
        change = np.expm1(6*np.log1p(s))
        return offset+amplitude*(weights @ (U2*change))/(1+weights@change)
    return float(brentq(mean, -abs(amplitude), abs(amplitude), xtol=1e-23))


def check_constraint(branch, amplitude, c, offset=0.):
    # Lobatto grid includes the poles and differs from the projection nodes.
    chi = np.linspace(0., np.pi, 1537)
    n = 2*np.arange(len(c))
    d = c.copy()
    d[0] -= 1
    s = evaluate(d, chi)
    g2, du = source(branch, amplitude, chi, offset)
    scaled = residual(s, evaluate(d*n*(n+2), chi), g2, du)
    return dict(max_abs_equation_residual=float(max(abs(scaled))),
                max_abs_physical_Hamiltonian_over_u_top=float(max(abs(scaled)/(1+s)**5)),
                psi_min_sampled=float(1+min(s)), psi_max_sampled=float(1+max(s)),
                psi_positive_lower_bound=float(c[0]-(n[1:]+1)@abs(c[1:])))


def observables(branch, amplitude, c, offset=0.):
    chi, weights = quadrature(383)
    q = offset + amplitude*(4*np.cos(chi)**2-1)
    psi = evaluate(c, chi)
    volume_weights = weights*psi**6
    volume_ratio = np.sum(volume_weights)
    proper = volume_weights/volume_ratio
    u, force = branch.potential(q)
    # These calibrated profiles have a single stationary point on the field range.
    endpoints = sorted([offset-amplitude, offset+3*amplitude])
    top_in_range = np.clip(0., *endpoints)
    u_extrema, _ = branch.potential(np.r_[endpoints, top_in_range])
    if np.any(q*force < -1e-22):
        raise ArithmeticError('Fitted potential leaves monotone hilltop flanks')
    return dict(proper_volume_over_round=float(volume_ratio),
                proper_volume=2*np.pi**2*volume_ratio/branch.h**3,
                volume_radius_over_a_b=float(volume_ratio**(1/3)),
                proper_mean_q=float(proper@q),
                proper_mean_q_acceleration=float(proper@force+3*(proper@((u-branch.u_top)*(q-proper@q)))),
                min_local_H_dot_over_u_top=float(min(u_extrema)/branch.u_top),
                max_local_H_dot_over_u_top=float(max(u_extrema)/branch.u_top),
                volume_scale_acceleration_over_u_top=float(proper@u/branch.u_top))


def calculate(previous):
    degree, nodes = 24, 128
    result = dict(maximum_even_degree=degree, quadrature_points=nodes,
                  off_grid_points=1537, observable_quadrature_points=383, branches={})
    for key, branch in load_branches(previous).items():
        analytic = second_order(branch)
        records, pairs = {}, {}
        for amplitude in (0., .001, -.001, .0005, -.0005):
            c = solve_slice(branch, amplitude, degree, nodes)
            errors = check_constraint(branch, amplitude, c)
            fine = solve_slice(branch, amplitude, 32, 192)
            fine_delta = fine.copy()
            fine_delta[:len(c)] -= c
            norm_bound = np.arange(1, 2*len(fine), 2) @ abs(fine_delta)
            if errors['max_abs_equation_residual'] > 1e-11 or norm_bound > 1e-12:
                raise ArithmeticError('Constraint or angular refinement residual exceeded tolerance')
            records[str(amplitude)] = dict(psi_even_U_coefficients=c.tolist(), **errors,
                    **observables(branch, amplitude, c), fine_solution_sup_difference_bound=float(norm_bound))
            if amplitude == 0 and (c[0] != 1 or np.any(c[1:])):
                raise ArithmeticError('Zero-amplitude slice differs from round metric')
        centers = {}
        for amplitude in (.001, -.001, .0005, -.0005, 0.):
            offset = centered_offset(branch, amplitude)
            cc = solve_slice(branch, amplitude, offset=offset)
            center_values = observables(branch, amplitude, cc, offset)
            errors = check_constraint(branch, amplitude, cc, offset)
            if abs(center_values['proper_mean_q']) > 1e-17 or errors['max_abs_equation_residual'] > 1e-11:
                raise ArithmeticError('Centered data violate proper mean or Hamiltonian constraint')
            centers[str(amplitude)] = dict(offset=offset, psi_even_U_coefficients=cc.tolist(),
                    **errors, **center_values)
        for amplitude in (.001, .0005):
            pos, neg = (records[str(x)] for x in (amplitude, -amplitude))
            even = (np.array(pos['psi_even_U_coefficients'][:3]) + neg['psi_even_U_coefficients'][:3])/2
            even[0] -= 1
            error = even/amplitude**2 - analytic
            pairs[str(amplitude)] = dict(even_psi_change_over_A2=(even/amplitude**2).tolist(),
                    coefficient_errors=error.tolist(),
                    even_proper_volume_shift_over_A2=((pos['proper_volume_over_round']-1)+(neg['proper_volume_over_round']-1))/(2*amplitude**2),
                    even_volume_acceleration_shift_over_A2=((pos['volume_scale_acceleration_over_u_top']-1)+(neg['volume_scale_acceleration_over_u_top']-1))/(2*amplitude**2))
            if max(abs(error/analytic)) > 2e-4:
                raise ArithmeticError('Even data inconsistent with analytic leading geometry')
        refinements = []
        for nmax, nq in ((4, 32), (6, 48), (8, 64), (12, 96), (24, 128), (32, 192)):
            trial = solve_slice(branch, .001, nmax, nq)
            refinements.append(dict(maximum_degree=nmax, quadrature_points=nq,
                    off_grid_max_abs_equation_residual=check_constraint(branch,.001,trial)['max_abs_equation_residual']))
        result['branches'][key] = dict(top=branch.top, u_top=branch.u_top, a_b=1/branch.h,
                removed_u_slope=branch.du_top, mu2_over_u_top=branch.mu2/branch.u_top,
                analytic_psi_change_over_A2=analytic.tolist(),
                analytic_uncentered_proper_mean_q_over_A3=6*analytic[1],
                analytic_centering_offset_over_A3=-6*analytic[1],
                proper_mean_centered_runs=centers,
                analytic_proper_volume_shift_over_A2=6*analytic[0],
                analytic_volume_acceleration_shift_over_A2=-branch.mu2/(6*branch.u_top),
                runs=records, amplitude_pairs=pairs, plus_0p001_spatial_refinements=refinements)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, required=True,
                        help='Previous preparation archive study root with nonlinear/ and prior/')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.previous), indent=2, allow_nan=False)+'\n')
