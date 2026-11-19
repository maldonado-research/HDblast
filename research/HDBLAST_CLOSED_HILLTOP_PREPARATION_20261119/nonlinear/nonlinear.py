"""Classical zonal scalar evolution about a stationary closed de Sitter hilltop."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import eval_chebyu, eval_gegenbauer, hyp2f1

PROBE_X = 0.12
AMPLITUDES = (0.001, 0.0005, 0.00025, 0.0001, 0.00001)


def prior_model(prior):
    spec = importlib.util.spec_from_file_location("zonal_prior_hilltop", prior/"hilltop/hilltop.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    model, reference = module.BG.source(prior/"archive")
    return module, model, reference, model.build_tables()


def potential_and_slope(p, theta):
    phi, f = p.spline(theta).T
    dphi, df = p.spline(theta, 1).T
    z = 1 - phi/p.pm
    weight = 1 + p.b*phi
    denom = f*f*p.v0
    pot = z*z*weight/denom
    slope = z*(-2*weight/p.pm + z*p.b)*dphi/denom - 2*pot*df/f
    return pot, slope


def local_quadratic_force(p, top, direction):
    # Differentiate the quotient of the two local spline polynomials in an
    # offset about the top, taking the segment on each side of its nearby knot.
    idx = np.searchsorted(p.theta, top + direction*1e-6) - 1
    d = Polynomial([top-p.theta[idx], 1.])
    phi, f = (Polynomial(c[::-1])(d) for c in p.spline.c[:, idx, :].T)
    numerator = (1-phi/p.pm)**2 * (1 + p.b*phi)
    denominator = p.v0 * f**2
    series = []
    for power in range(4):
        coeff = numerator.coef[power] - sum(denominator.coef[j]*series[power-j]
                                             for j in range(1, power+1))
        series.append(coeff/denominator.coef[0])
    return -9*series[3]  # -3 u'''/2


def retarded_quadratic_integral(x, h, mu2):
    """Integrate the analytic homogeneous Green function times linear U_2 gain squared."""
    disc = math.sqrt(9 + 4*mu2/h**2)
    def basis(n, odd, t):
        sh, ch = math.sinh(h*t), math.cosh(h*t)
        a = n/2 + (3+disc)/4 + odd/2
        b = n/2 + (3-disc)/4 + odd/2
        return ch**n * (sh/h)**odd * hyp2f1(a,b,.5+odd,-sh*sh)
    displacement, velocity = basis(0,0,x), basis(0,1,x)
    def integrand(s):
        G = (velocity*basis(0,0,s) - displacement*basis(0,1,s)) * math.cosh(h*s)**3
        return G * basis(2,0,s)**2
    return quad(integrand,0.,x,epsabs=1e-14,epsrel=1e-13)[0]


def force_radii(p, top, mu2):
    radii = []
    grid = np.linspace(1e-5, .25, 1001)
    for sign in [-1, 1]:
        def residual(d):
            _, slope = potential_and_slope(p, top+sign*d)
            return np.abs(-3*slope/(mu2*sign*d)-1)-.1
        i = np.flatnonzero(residual(grid)>0)[0]
        radii.append(brentq(residual,grid[i-1],grid[i]))
    return radii


class Zonal:
    def __init__(self, p, top, h, mu2, maximum_n, quadrature_count, monitor_count, cap):
        self.p, self.top, self.h, self.mu2, self.cap = p, top, h, mu2, cap
        self.degrees = np.arange(0, maximum_n+1, 2)
        self.count = len(self.degrees)
        self.k2 = self.degrees*(self.degrees+2)
        chi = np.pi*np.arange(1, quadrature_count+1)/(quadrature_count+1)
        self.weights = 2*np.sin(chi)**2/(quadrature_count+1)
        self.basis = eval_chebyu(self.degrees, np.cos(chi)[:, None])
        self.project = self.basis.T * self.weights
        self.error_gram = np.max(abs(self.project @ self.basis - np.eye(self.count)))
        if self.error_gram > 1e-12:
            raise ArithmeticError("Radial basis normalization")
        monitor_angles = np.linspace(0., np.pi, monitor_count)[:, None]
        points = np.cos(monitor_angles)
        self.monitor = eval_chebyu(self.degrees, points)
        self.derivative = np.zeros_like(self.monitor)
        self.derivative[:, 1:] = (-2*np.sin(monitor_angles)
                                  * eval_gegenbauer(self.degrees[1:]-1,2,points))
        self.u_top, self.du_top = potential_and_slope(p, top)
        self.sided = [local_quadratic_force(p, top, sign) for sign in [-1, 1]]
        self.radii = force_radii(p,top,mu2)

    def potential(self, field):
        u, du = potential_and_slope(self.p, self.top + field)
        # The tiny residual fitted slope at the nominal stationary point is
        # removed consistently in the force and its antiderivative.
        return u - self.u_top - self.du_top*field, -3*(du-self.du_top)

    def energy(self, x, state):
        c, v = state[:2*self.count].reshape(2, self.count)
        inv_a2 = self.h**2/math.cosh(self.h*x)**2
        pot, _ = self.potential(self.basis @ c)
        return np.array([v @ v / 6, inv_a2*(self.k2*c @ c)/6, self.weights @ pot])/self.h**2

    def evolve(self, amplitude):
        n = self.count
        xh = self.h
        z0 = np.zeros(2*n+7)
        z0[1] = amplitude
        z0[2*n] = 1.  # linear U_2 displacement gain
        shape = self.basis[:, 1]
        coeff = self.weights @ (np.where(amplitude*shape < 0, *self.sided)*shape**2)

        def equation(x, y):
            c, v = y[:2*n].reshape(2, n)
            gain, gv, mean2, v2, mean_born, vb, work = y[2*n:]
            H = xh*math.tanh(xh*x)
            inv_a2 = xh**2/math.cosh(xh*x)**2
            _, force = self.potential(self.basis @ c)
            linear_field = amplitude*gain*shape
            _, force_linear_profile = self.potential(linear_field)
            mean_source = self.weights @ (force_linear_profile - self.mu2*linear_field)
            grad = inv_a2*(self.k2*c @ c)
            return np.r_[v, self.project @ force - inv_a2*self.k2*c - 3*H*v,
                         gv, (self.mu2-8*inv_a2)*gain - 3*H*gv,
                         v2, self.mu2*mean2 - 3*H*v2 + coeff*(amplitude*gain)**2,
                         vb, self.mu2*mean_born - 3*H*vb + mean_source,
                         H*(v @ v + grad/3)/xh**2]

        def field_stop(x, y):
            return self.cap - np.max(abs(self.monitor @ y[:n]))
        def energy_stop(x, y):
            return 0.1 - np.max(abs(self.energy(x,y)))
        def force_boundary(x, y):
            field = self.monitor @ y[:n]
            return min(self.radii[0]+min(field), self.radii[1]-max(field))
        force_boundary.direction = -1
        for event in (field_stop, energy_stop):
            event.direction = -1
            event.terminal = True
        # Stop the exactly stationary control at a specified comparison time.
        sol = solve_ivp(equation, (0., 0.8 if amplitude == 0 else 2.), z0, method="DOP853",
                        rtol=1e-11, atol=2e-15, max_step=.014, dense_output=True,
                        events=[field_stop, energy_stop, force_boundary])
        if not sol.success:
            raise ArithmeticError(sol.message)
        end = sol.t[-1]
        reason = next((name for name, ev in zip(("peak_field","averaged_energy"), sol.t_events)
                       if len(ev)), "control_endpoint")
        if amplitude and reason == "control_endpoint":
            raise ArithmeticError("Nonzero launch missed stop")
        energy0 = self.energy(0, z0)
        probe_points = np.linspace(0.,end,61)
        residual = []
        max_component = 0.
        for x in probe_points:
            y = sol.sol(x)
            E = self.energy(x,y)
            max_component = max(max_component, max(abs(E)))
            residual.append(sum(E-energy0) + y[-1])

        def record(x):
            y = sol.sol(x)
            c = y[:n]
            mean2, born = y[2*n+2], y[2*n+4]
            kinetic, grad, pot = self.energy(x,y)
            local_pot, _ = self.potential(self.monitor @ c)
            local_max = dict(kinetic=max((self.monitor @ y[n:2*n])**2)/(6*xh**2),
                             gradient=max((self.derivative @ c)**2)/(6*math.cosh(xh*x)**2),
                             abs_potential_shift=max(abs(local_pot))/xh**2)
            return dict(x=float(x), mean=float(c[0]), quadratic_mean=float(mean2),
                        full_force_on_linear_profile_mean=float(born),
                        mode_coefficients=c.tolist(), maximum_field=float(np.max(abs(self.monitor @ c))),
                        largest_local_energy_terms_over_Vstar=local_max,
                        energy_over_Vstar=dict(kinetic=kinetic, gradient=grad,
                                               potential_minus_top=pot, signed_sum=kinetic+grad+pot))
        return dict(amplitude=amplitude, stop=reason, initial=record(0.), at_probe=record(PROBE_X),
                    at_stop=record(end), weighted_sided_quadratic_source_coefficient=coeff,
                    x_first_local_10_percent_force_boundary=(float(sol.t_events[2][0])
                                                              if len(sol.t_events[2]) else None),
                    max_sampled_abs_energy_component_over_Vstar=max_component,
                    max_abs_energy_continuity_residual=max(abs(z) for z in residual),
                    energy_residual_over_largest_component=max(abs(z) for z in residual)/max(max_component,1e-30))


def calculate(prior):
    module, model, ref, tables = prior_model(prior)
    data = dict(probe_x=PROBE_X, primary_maximum_degree=16, primary_quadrature_nodes=96,
                primary_field_monitor_nodes=257, branches={})
    for key, cap in [("phi_m=0.5", .04), ("phi_m=-0.5", .08)]:
        p = module.BG.Potential(model,tables,float(key.split("=")[1]),ref[key])
        top, h, mu2, _ = module.field_top(p)
        primary = Zonal(p,top,h,mu2,16,96,257,cap)
        records = {}
        pairs = {}
        retarded = retarded_quadratic_integral(PROBE_X, h, mu2)
        retarded_errors = []
        control = primary.evolve(0.)
        if max(abs(z) for z in control["at_stop"]["mode_coefficients"]) != 0.:
            raise ArithmeticError("Stationary control moved")
        for A in AMPLITUDES:
            pair = []
            for amp in [A,-A]:
                run = primary.evolve(amp)
                independent = retarded * amp**2 * run["weighted_sided_quadratic_source_coefficient"]
                retarded_errors.append(abs(run["at_probe"]["quadratic_mean"]/independent-1))
                pair.append(run["at_probe"]["mean"])
                records[str(amp)] = run
                if np.sign(pair[-1]) != np.sign(primary.sided[0]):
                    raise ArithmeticError("Monopole sign disagrees with force")
                if run["energy_residual_over_largest_component"] > 1e-8:
                    raise ArithmeticError("Scalar continuity check failed")
            even = (pair[0] + pair[1])/2
            odd = (pair[0] - pair[1])/2
            pairs[str(A)] = dict(even_mean=even, odd_mean=odd, even_mean_over_A_squared=even/A**2,
                                 odd_to_even=odd/even)
        comparisons = []
        for A in [.001,.0005]:
            ratio = pairs[str(A)]["even_mean"]/pairs[str(A/2)]["even_mean"]
            if abs(ratio-4)>0.01:
                raise ArithmeticError("Incorrect leading amplitude scaling")
            comparisons.append(dict(larger_amplitude=A, even_mean_ratio=ratio))
        refinements = []
        for degree, nodes, monitors in [(8,48,129),(24,144,513)]:
            test = Zonal(p,top,h,mu2,degree,nodes,monitors,cap)
            for amp in [0.001,-0.001,1e-5,-1e-5]:
                other = test.evolve(amp)["at_stop"]
                base = records[str(amp)]["at_stop"]
                delta_x = other["x"]-base["x"]
                relative_mean = other["mean"]/base["mean"]-1
                if degree==24 and max(abs(delta_x),abs(relative_mean)) > 1e-6:
                    raise ArithmeticError("Angular refinement failed")
                refinements.append(dict(maximum_degree=degree, quadrature_nodes=nodes,
                                        field_monitor_nodes=monitors, amplitude=amp,
                                        stop_x_difference=delta_x, stop_relative_mean_difference=relative_mean))
        data["branches"][key] = dict(top=top,h=h,mu_squared=mu2,peak_amplitude_cap=cap,
                                     degrees=primary.degrees.tolist(), removed_force_at_top=float(-3*primary.du_top),
                                     sided_quadratic_force_coefficients=primary.sided,
                                     directional_10_percent_force_radii=primary.radii,
                                     basis_orthonormality_error=primary.error_gram,zero_control=control,
                                     quadratic_retarded_integral_at_probe_per_unit_source=retarded,
                                     max_relative_retarded_integral_error=max(retarded_errors),
                                     runs=records, amplitude_pairs_at_probe=pairs,
                                     halving_tests_at_probe=comparisons, spatial_refinements=refinements)
        if max(retarded_errors)>1e-9:
            raise ArithmeticError("Analytic retarded integral disagrees with forced mean")
    return data


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior),indent=2,allow_nan=False)+"\n")
