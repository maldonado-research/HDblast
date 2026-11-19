"""Audit Euclidean matching from a driven canonical parent plus a finite oscillator bath."""
import argparse
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import solve
from scipy.special import roots_legendre


def require(condition, message):
    if not condition:
        raise ValueError(message)


class FrozenBath:
    """Canonical proper-time Hamiltonian with a source term -a^3 J chi."""

    def __init__(self, a, p, mass, cutoff, count, coupling, gate, static, kinetic):
        nu, weight = roots_legendre(count)
        q = a * p
        upper = a * cutoff
        omega = q + (nu + 1) * (upper - q) / 2
        dw = weight * (upper - q) / 2
        g2 = omega * (omega**2 - q**2) * dw / (4 * math.pi**2)
        cf = coupling * gate
        c0 = np.sum(g2 / omega**2)
        c1 = np.sum(g2 / omega**4)
        inertia = 1 - kinetic * cf**2 * c1
        require(inertia > 0, "Nonpositive matched kinetic coefficient")
        stiffness = np.diag(np.r_[a**3 * (p*p + mass*mass) + static * a * cf**2 * c0,
                                  omega**2 / a])
        stiffness[0, 1:] = stiffness[1:, 0] = cf * np.sqrt(g2)
        kinetic_matrix = np.diag(np.r_[1 / (a**3 * inertia), np.full(count, 1/a)])
        zeros = np.zeros_like(stiffness)
        self.flow = np.block([[zeros, kinetic_matrix], [-stiffness, zeros]])
        self.force = np.zeros(2*(count+1))
        self.force[count+1] = a**3
        self.external_stiffness = p*p + mass*mass
        self.a, self.cf, self.omega, self.g2 = a, cf, omega, g2
        self.c0, self.c1, self.inertia = c0, c1, inertia

    def susceptibility(self, frequency):
        lhs = -1j * frequency * np.eye(len(self.force)) - self.flow
        return solve(lhs, self.force)[0]

    def reduction(self, frequency, static, kinetic):
        sigma = np.sum(self.g2 / (self.omega**2 - (self.a * frequency)**2))
        correction = (static * self.c0 - sigma)/self.a**2 + kinetic * self.c1 * frequency**2
        return self.external_stiffness - frequency**2 + self.cf**2 * correction

    def causal_drive(self, xi):
        # Three tests use a drive that rises from the distant past. With the
        # common final amplitude J=1, it samples the upper-half-plane resolvent.
        duration = 32 / xi
        def rhs(t, state):
            return self.flow @ state + self.force * math.exp(xi * (t - duration))
        out = solve_ivp(rhs, (0, duration), np.zeros_like(self.force), method="DOP853",
                        rtol=2e-11, atol=2e-14)
        require(out.success, out.message)
        return float(out.y[0, -1])


def frozen_test():
    p, m, cutoff = 1.5, 2.3, 11.
    y, gate = .4, .8
    xi = np.array([.25, .5, 1., 2.]) * p
    freq = [*list(1j*xi), 2.6 + .4j]
    calibrations = ((False, False, "bare"), (True, False, "static"), (True, True, "matched"))
    curves = {}
    result = {}
    for a in (1., 1.7):
        by_scheme = {}
        for static, kinetic, label in calibrations:
            model = FrozenBath(a, p, m, cutoff, 32, y, gate, static, kinetic)
            susceptibility = np.array([model.susceptibility(w) for w in freq])
            inverse = 1/susceptibility
            reduced = np.array([model.reduction(w, static, kinetic) for w in freq])
            error = max(abs(inverse-reduced))
            require(error < 3e-12, "Schur complement disagrees with full canonical response")
            dc = (1/model.susceptibility(0)).real - model.external_stiffness
            w_small = 1e-3
            # An uncalibrated bath adds +g^2/omega^4 to the inverse slope.
            slope = ((1/model.susceptibility(1j*w_small)).real -
                     (1/model.susceptibility(0)).real) / w_small**2
            slope_expected = 1 + (1-kinetic)*(model.cf**2 * model.c1)
            require(abs(slope-slope_expected) < 1.5e-8, "Independently estimated zero-frequency slope")
            require(abs(dc - model.cf**2*(int(static)-1)*model.c0/a**2) < 2e-12,
                    "DC shift disagrees with static completion")
            if label == "matched":
                response = (model.external_stiffness + xi**2 - inverse[:-1].real) / model.cf**2
                u = xi[:, None]**2
                omega2 = model.omega[None, :]**2
                spectral = (a**2 * u**2 * model.g2 / (omega2**2 * (omega2+a**2*u))).sum(axis=1)
                require(np.all(response > 0) and max(abs(response-spectral)) < 2e-12,
                        "Matched sign or frequency dependence")
            else:
                response = None
            curves[a, label] = susceptibility
            by_scheme[label] = dict(maximum_inverse_disagreement=float(error),
                                    static_inverse_shift=float(dc),
                                    finite_difference_inverse_slope=float(slope),
                                    expected_inverse_slope=float(slope_expected),
                                    kinetic_coefficient=float(model.inertia),
                                    positive_frequency_susceptibility_imag=float(susceptibility[-1].imag))
            if response is not None:
                by_scheme[label]["positive_remainder_from_inverse"] = response.tolist()
                by_scheme[label]["max_remainder_disagreement"] = float(max(abs(response-spectral)))
            require(susceptibility[-1].imag>0, "Retarded parent sign at the complex probe")
        result[str(a)] = by_scheme
    change = max(float(max(abs(curves[1., label]-curves[1.7, label]))) for _,_,label in calibrations)
    require(change < 1e-12, "Changing conformal normalization changed proper susceptibility")
    drive={}
    for st,kt,label in calibrations:
        model=FrozenBath(1.7,p,m,cutoff,10,y,gate,st,kt)
        value=model.causal_drive(1.15*p)
        reference=model.susceptibility(1j*1.15*p).real
        require(abs(value-reference)<2e-11, "Retarded drive disagrees with full first-order resolvent")
        drive[label]=dict(parent_endpoint=value, full_resolvent=float(reference),
                          abs_error=float(abs(value-reference)))
    return dict(proper_momentum=p, proper_mass=m, proper_cutoff=cutoff, gate=gate, coupling=y,
                euclidean_frequency_over_momentum=(xi/p).tolist(), complex_probe=[2.6,.4],
                modes=32, maximum_proper_response_change=float(change), scale_factors=result,
                retarded_drive=dict(start_from_zero=True, duration_in_e_folds=32,
                                    euclidean_frequency_over_momentum=1.15, modes=10, cases=drive))


def tail_check():
    # x = proper pair energy squared and t = cutoff^2 / x compactifies the tail.
    result = {}
    for p in (.019251120011914757, 1.5, 7.980748879988085):
        for r in (.25,.5,1.,2.):
            u = (r*p)**2
            q2 = p*p
            continuum = ((q2+u)*math.log1p(u/q2)-u) / (8*math.pi**2)
            independent = u**2/(8*math.pi**2*q2)*quad(
                lambda t:(1-t)/(1+u*t/q2), 0, 1, epsabs=3e-14)[0]
            require(math.isclose(continuum,independent,rel_tol=1e-13), "Continuum integral mismatch")
            for cutoff in (60.,120.,240.,480.):
                X = cutoff**2
                leading = u**2/(8*math.pi**2*X)
                ratio = quad(lambda t:(1-q2*t/X)/(1+u*t/X),0,1,epsabs=3e-14)[0]
                center = 1 - (q2+u)/(2*X)
                remainder_cap = u*(q2+u)/(3*X**2)
                require(center-4e-15 <= ratio <= center+remainder_cap+4e-15,
                        "Two-term tail enclosure failed")
                finite = u**2/(8*math.pi**2*q2)*quad(
                    lambda t:(1-t)/(1+u*t/q2), q2/X, 1, epsabs=3e-14)[0]
                require(math.isclose(finite+leading*ratio,independent,rel_tol=3e-14),
                        "Finite and omitted intervals fail to recover the whole spectrum")
                result[f"{p:g}/{r:g}/{cutoff:g}"] = dict(
                    positive_finite_remainder=finite, positive_continuum_remainder=continuum,
                    missing_fraction=leading*ratio/continuum, tail_over_leading=ratio,
                    rigorous_lower_over_leading=center, rigorous_upper_over_leading=center+remainder_cap)
    return result


def compact_quadrature(native):
    if native is None:
        # Native fork uses 24 Gauss-Legendre parent momenta on (0, 8).
        momenta=4*(np.polynomial.legendre.leggauss(24)[0]+1)
    else:
        momenta=np.asarray(json.loads((native/"FORK.json").read_text())["fine"]["q"])
        require(len(momenta)==24 and np.all((momenta>0) & (momenta<8)), "Native parent momenta")
    counts=(8, 12, 16, 24)
    quadratures={n:np.polynomial.legendre.leggauss(n) for n in counts}
    worst={n:0. for n in counts}
    with localcontext() as ctx:
        ctx.prec=52
        one=Decimal(1)
        for q in momenta:
            for cutoff in (60.,120.,240.,480.):
                beta=(q/cutoff)**2
                b=Decimal.from_float(float(beta))
                for r in (.25,.5,1.,2.):
                    k=r*r
                    kd=Decimal.from_float(k)
                    # Exact antiderivative of (1-t)/(1+kt), evaluated with guard digits.
                    ref=float(((one+kd)*((one+kd).ln()-(one+kd*b).ln())-kd*(one-b))/kd**2)
                    for n,(nodes,weights) in quadratures.items():
                        t=beta+(1-beta)*(1+nodes)/2
                        value=(1-beta)*(weights@((1-t)/(1+k*t)))/2
                        err=abs(value/ref-1)
                        worst[n]=max(worst[n],float(err))
    results={}
    for n in counts:
        # Worst case k=4, beta=0 gives distance of the pole A/B=3/2.
        # A Chebyshev relative residual bounds the integral at any allowed k,beta.
        delta=1/math.cosh((2*n-1)*math.acosh(1.5))
        exact_arithmetic_bound=2*delta/(1-delta)
        err=worst[n]
        require(err <= exact_arithmetic_bound+3e-15, "Mapped quadrature bound plus roundoff")
        results[str(n)] = dict(max_relative_discrepancy=err,
                               uniform_relative_bound_exact_arithmetic=exact_arithmetic_bound)
    require(results["16"]["max_relative_discrepancy"]<1e-12, "Compact frequency quadrature")
    return dict(native_momentum_count=len(momenta), cutoff_count=4, ratio_count=4, cases=384,
                integration_variable="t=(q/omega)^2", orders=results)


def compare_response_tables(path):
    """Validate archived finite-cutoff, continuum-limit and omitted-tail fields independently."""
    nodes,weights=np.polynomial.legendre.leggauss(24)
    with np.load(path,allow_pickle=False) as data:
        q,cutoffs,ratios = (data[key] for key in ("q","cutoffs","ratios"))
        require(data["orders"][-1]==2048 and data["analytic_matched"].shape==(4,24,4),
                "Unexpected archived response axes")
        relative={key:0. for key in ("analytic_matched","smooth_integral","finite_matched_2048",
                                     "infinite_matched","missing_tail")}
        fields={key:data[key] for key in ("analytic_matched","smooth_integral","finite_matched",
                                          "infinite_matched","missing_tail")}
        for ci,cutoff in enumerate(cutoffs):
            for qi,p in enumerate(q):
                beta=(p/cutoff)**2
                for ri,r in enumerate(ratios):
                    k=r*r
                    prefactor=(p*k)**2/(8*math.pi**2)
                    def integrate(lo,hi):
                        t=lo+(hi-lo)*(1+nodes)/2
                        return prefactor*(hi-lo)*(weights@((1-t)/(1+k*t)))/2
                    finite,tail,full=integrate(beta,1),integrate(0,beta),integrate(0,1)
                    index=(ci,qi,ri)
                    values=dict(analytic_matched=(fields["analytic_matched"][index],finite),
                                smooth_integral=(fields["smooth_integral"][index],finite),
                                finite_matched_2048=(fields["finite_matched"][(-1,*index)],finite),
                                infinite_matched=(fields["infinite_matched"][qi,ri],full),
                                missing_tail=(fields["missing_tail"][index],tail))
                    for key,(found,reference) in values.items():
                        relative[key]=max(relative[key],float(abs(found/reference-1)))
        for key,error in relative.items():
            limit=3e-9 if key=="finite_matched_2048" else 2e-12
            require(error<limit, f"Archived table field {key} differs from compact quadrature")
    return relative


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--native", type=Path, help="Optional native directory containing FORK.json")
    parser.add_argument("--tables", type=Path, default=Path(__file__).resolve().parents[1]/"response/TABLES.npz")
    args=parser.parse_args()
    payload=json.dumps(dict(frozen_hamiltonian=frozen_test(), continuum_tail=tail_check(),
                            compact_quadrature=compact_quadrature(args.native),
                            table_comparison=compare_response_tables(args.tables)),
                       indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload,end="")

if __name__=="__main__":
    main()
