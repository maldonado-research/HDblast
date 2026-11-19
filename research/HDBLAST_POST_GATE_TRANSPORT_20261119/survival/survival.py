"""Lower bounds on surviving spectral structure after an elastic scattering interval."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammainc, roots_legendre

ROOT = Path(__file__).resolve().parents[1]
Y = 0.1
SPIN = 2
END = 2.16


def library(directory, filename):
    path = ROOT / directory / filename
    spec = importlib.util.spec_from_file_location(directory + "_spectral_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def equilibrium(number, energy, edges):
    """Return binned unit-mass MB and FD probability, including the overflow bin."""
    theta_mb = energy / (3 * number)
    def cumulative(theta):
        return gammainc(3, edges / theta)
    classical = np.diff(np.r_[cumulative(theta_mb), 1.])
    # The series for z/(exp(x)+z) converges absolutely for the small z bracketed here.
    j = np.arange(1, 13, dtype=float)
    def series(z):
        return (-1.)**(j-1) * z**j
    def trial(z):
        terms = series(z)
        f3, f4 = np.sum(terms / j**3), np.sum(terms / j**4)
        theta = theta_mb * f3 / f4
        return SPIN / np.pi**2 * theta**3 * f3, theta
    z = brentq(lambda z: trial(z)[0] - number, 1.e-8, .05, xtol=1.e-16)
    _, theta = trial(z)
    coeff = series(z)/j**3
    fd_cdf = np.sum(coeff[:, None]*gammainc(3, j[:, None]*edges / theta), axis=0)/sum(coeff)
    fermi = np.diff(np.r_[fd_cdf, 1.])
    assert abs(trial(z)[0] / number - 1) < 1.e-10
    return classical, fermi, float(z), float(theta)


def excess_counts(survivors, target):
    return float(np.maximum(np.r_[survivors, 0.] - target, 0.).sum())


def sum_cells(x, factor):
    return x.reshape(-1, factor).sum(axis=1)


def target_energy_floor(q, omega, weights, duration, energy, inverse, y, order):
    """Necessary target column for ten percent of the original particles to collide."""
    nodes, w = roots_legendre(order)
    k = (omega[:, :, None] + q[:, None, None] * nodes) / 2
    probabilities = weights[:, :, None] * w / (2 * weights.sum())
    def fraction(column):
        return np.sum(probabilities * (-np.expm1(-8*np.pi*column/k)))
    column = brentq(lambda u: fraction(u)-.10, 0., 2., xtol=1.e-14)
    jensen = -np.log1p(-.1) * weights.sum() / (8*np.pi*inverse)
    mean_inv = column / duration
    floor = 2*np.pi**2/SPIN * mean_inv**2
    # At fixed inverse moment, the fully occupied ball minimizes energy.
    radius = np.sqrt(4*np.pi**2*mean_inv / SPIN)
    sphere_k = (nodes+1)*radius/2
    sphere = SPIN/(2*np.pi**2)*w*radius/2*sphere_k**2
    assert abs(np.sum(sphere/sphere_k)/mean_inv - 1.) < 2.e-12
    assert abs(np.sum(sphere*sphere_k)/floor - 1.) < 2.e-12
    assert 0 < jensen <= column and abs(fraction(column)-.10)<2.e-12
    return dict(target_tag_fraction=.10, necessary_integrated_inverse_moment=column,
                jensen_necessary_column=jensen, necessary_mean_inverse_moment=mean_inv,
                necessary_mean_target_energy=floor, emitted_per_sign_energy=y**2*energy,
                minimum_ratio_to_emitted_energy=floor/(y**2*energy))


def experiment(repo):
    daughter = library("daughters", "daughters.py")
    transport = library("transport", "transport.py")
    result = {"illustrative_y":Y, "spin_dof_per_sign":SPIN, "end":END, "histories":{}}
    for arm in daughter.ARMS:
        q, omega, weights = daughter.pair_nodes(repo, arm, "full")
        raw = daughter.moments(q, omega, weights)
        number, energy, inv = (raw[key] for key in ("0","1","-1"))
        mean = energy / number
        physical_count = Y**2 * number
        kf = (6*np.pi**2 * physical_count / SPIN)**(1/3)
        max_inv = SPIN*kf**2 / (4*np.pi**2)
        # Independently integrate the maximally occupied ball.
        p,w = np.polynomial.legendre.leggauss(8)
        k=(p+1)*kf/2
        measure=SPIN/(2*np.pi**2)*w*kf/2*k**2
        assert np.isclose(measure.sum(), physical_count, rtol=1.e-13)
        assert np.isclose((measure/k).sum(), max_inv, rtol=1.e-13)
        column = transport.column(repo, arm, END)
        target_floor = target_energy_floor(q, omega, weights, column["delta_eta"], energy, inv, Y, 40)
        check_nodes = daughter.pair_nodes(repo, arm, "full", 64, 3)
        lower_res = target_energy_floor(*check_nodes, column["delta_eta"], energy, inv, Y, 20)
        target_floor["relative_refinement"] = abs(target_floor["necessary_mean_target_energy"] /
                                                 lower_res["necessary_mean_target_energy"] - 1.)
        assert target_floor["relative_refinement"] < 3.e-6
        b = 8*np.pi*max_inv*column["delta_eta"]
        ncap = -np.expm1(-b * inv/number)
        ecap = -np.expm1(-b * number/energy)
        prefactor = 8/3*energy*column["J"]
        fine_h = 1/64
        fine = np.linspace(0,64, int(64 / fine_h) + 1)
        curves = daughter.pair_interpolants(repo, arm, 144, 5)
        counts = daughter.binned(q, "full", curves, fine)[0]/number
        clipped = abs(np.minimum(counts,0.).sum())
        assert clipped < 2.e-9 and abs(1-counts.sum())<1.e-8
        counts = np.maximum(0.,counts)
        def histogram_column(k_edge):
            def frac(u):
                active = k_edge > 0.
                hit = np.ones_like(k_edge)
                hit[active] = -np.expm1(-8.*np.pi*u/k_edge[active])
                return counts @ hit - .10
            return brentq(frac,1.e-8,2.)
        lower_col = histogram_column(fine[:-1])
        upper_col = histogram_column(fine[1:])
        assert lower_col <= target_floor["necessary_integrated_inverse_moment"] <= upper_col
        target_floor["independent_column_bracket"] = [lower_col,upper_col]
        target_floor["conservative_lower_energy_ratio"] = (
            2*np.pi**2/SPIN * (lower_col / column["delta_eta"])**2 / (Y**2*energy)
        )
        # Compare a separately evaluated coarse histogram; account for cumulative roundoff.
        curves2 = daughter.pair_interpolants(repo, arm, 80, 3)
        bins2 = daughter.binned(q, "full", curves2, fine[::16])[0]/number
        comparison = np.sum(abs(sum_cells(counts,16) - bins2))
        assert comparison < 1.e-7
        table={}
        for width in (1., .5, .25):
            m = int(round(width/fine_h))
            edges = fine[::m]
            mb, fd, z, theta = equilibrium(physical_count, Y**2*energy, edges)
            masses = sum_cells(counts,m)
            fd0 = excess_counts(masses, fd)
            mb0 = excess_counts(masses, mb)
            assert 0 < mb0 < 1 and 0 < fd0 < 1
            def contact(x, side):
                ends=fine[1:] if side=="lower" else fine[:-1]
                return excess_counts(sum_cells(counts*np.exp(-x*ends/mean),m),mb)
            thresholds={}
            for delta in (.10,.01):
                lo=brentq(lambda x:contact(x,"lower")-delta,0,32)
                hi=brentq(lambda x:contact(x,"upper")-delta,0,32)
                assert 0 < lo <= hi < 32
                thresholds[str(delta)]=dict(mean_depth_lower=lo,mean_depth_upper=hi,
                    coupling_product_necessary=lo / (mean*prefactor))
            # Before a tag's first collision, it redshifts freely. For the Pauli envelope
            # the smallest surviving fraction in a microbin is taken at its lower edge.
            survival=np.zeros(len(counts))
            survival[1:]=np.exp(-b/fine[1:-1])
            direct=excess_counts(sum_cells(counts*survival,m),fd)
            triangle=max(0.,fd0-ncap)
            table[str(width)]=dict(initial_TV_MB=mb0,initial_TV_FD=fd0,
                final_TV_FD_lower_from_fraction=triangle,final_TV_FD_lower_from_survivors=direct,
                final_TV_FD_lower=max(triangle,direct), contact=thresholds)
        result["histories"][arm]=dict(mean_k=mean,thermal_fugacity=z,thermal_comoving_scale=theta,
            fermi_ball_radius=kf,maximum_comoving_inverse_moment=max_inv,delta_eta=column["delta_eta"],
            upper_tag_exponent_numerator=b,upper_first_event_fraction_number=float(ncap),
            upper_first_event_fraction_initial_energy=float(ecap), bin_negative_roundoff=clipped,
            independent_bin_L1_discrepancy=float(comparison),bin_widths=table, target_energy_floor=target_floor)
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    args.output.write_text(json.dumps(experiment(args.repo.resolve()),indent=2,allow_nan=False)+"\n")
