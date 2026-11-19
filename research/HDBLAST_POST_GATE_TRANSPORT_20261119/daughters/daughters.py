"""Reconstruct the one-particle momentum distributions of massless Dirac daughters."""
import argparse
import csv
import importlib.util
import json
from math import comb
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.polynomial import Chebyshev
from scipy.fft import dct
from scipy.special import roots_legendre

PRIOR = Path("research/HDBLAST_QUANTUM_MEMORY_20261119")
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
BANDS = ("full", "low", "high")
ARMS = ("quiet", "plus", "minus", "far")
DEGREES = tuple(range(-1, 7))
EDGES = np.linspace(0, 64, 129)


@lru_cache(maxsize=2)
def previous(repo):
    path = Path(repo) / PRIOR / "quantum/spectrum.py"
    spec = importlib.util.spec_from_file_location("hdblast_ordered_spectrum", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=12)
def history(repo, arm):
    repo = Path(repo)
    native = repo / NATIVE
    seed = json.loads((native / "FORK.json").read_text())["fine"]
    model = previous(repo)
    quiet = model.wall(native, "quiet", seed)[1]
    loga, phi = model.wall(native, arm, seed)
    ends, solution, _ = model.fundamental(seed, loga, phi, quiet)
    return seed, ends, solution, phi, quiet, float(np.exp(loga(model.READOUT)))


@lru_cache(maxsize=16)
def drive(repo, arm, time_order):
    seed, ends, solution, phi, quiet, _ = history(repo, arm)
    return previous(repo).profile(seed, ends, solution, phi, quiet, time_order)


def band_bounds(q, band):
    if band not in BANDS:
        raise ValueError(f"Unknown band {band}")
    return np.full_like(q, 60.) if band == "high" else q, 60. if band == "low" else 120.


def number_density(repo, arm, omega, time_order):
    """Pair count dN/(y² dω), with the native radial q measure in each row."""
    seed = history(repo, arm)[0]
    q, wq = np.array(seed["q"]), np.array(seed["weight"])
    eta, measure, h, k = drive(repo, arm, time_order)
    # A common phase has no effect on the ordered quadratic form.
    centered_time = eta - .5*(eta[0]+eta[-1])
    X, C, P = (np.array(seed["states"]["retain"][name]) for name in ("X", "C", "P"))
    power = np.empty_like(omega)
    for j, frequencies in enumerate(omega):
        oscillation = np.exp(1j*frequencies[:, None]*centered_time)
        vh, vk = (oscillation @ (np.column_stack((h[j], k[j]))*measure[:, None]).astype(complex)).T
        # The matrix product has shape (Nω,2); its columns hold the amplitudes.
        power[j] = (X[j]*abs(vh)**2 + P[j]*abs(vk)**2
                    + 2*C[j]*np.real(vh.conj()*vk) - np.imag(vh.conj()*vk))
    if np.min(power) < -2e-16:
        raise ValueError("Negative ordered power")
    return wq[:, None]*(omega-q[:, None])*(omega+q[:, None])*power/(8*np.pi**2)


def pair_nodes(repo, arm, band="full", frequency_order=96, time_order=5):
    """Return native radial q, pair frequencies ω and one-sign comoving count weights.

    Each weight includes the q and ω quadrature and excludes y² H0³.
    """
    repo = Path(repo).resolve()
    q = np.array(history(repo, arm)[0]["q"])
    low, high = band_bounds(q, band)
    nodes, weights = roots_legendre(frequency_order)
    x = (nodes+1)/2
    exponent = 1 if band == "high" else 4
    omega = low[:, None] + (high-low)[:, None]*x**exponent
    widths = (high-low)[:, None]*exponent*x**(exponent-1)*weights/2
    return q, omega, widths * number_density(repo, arm, omega, time_order)


def moments(q, omega, pair_weights):
    """Raw one-sign moments, integrated analytically over the decay direction."""
    center = omega/2
    beta = q[:, None]/omega
    # atanh(beta) is evaluated from ω±q. For a node rounded to threshold
    # the pair weight vanishes and its log-weighted contribution is zero.
    d = omega-q[:, None]
    reciprocal = np.zeros_like(omega)
    np.log1p(np.divide(2*q[:, None], d, out=np.zeros_like(d), where=d>0), out=reciprocal)
    reciprocal /= q[:, None]
    values = {"-1": float(np.sum(pair_weights * reciprocal))}
    for degree in range(7):
        angular_mean = sum(comb(degree, j)*beta**j/(j+1) for j in range(0,degree+1,2))
        values[str(degree)] = float(np.sum(pair_weights * center**degree * angular_mean))
    return values


def angular_controls(q, omega, pairs):
    nodes, weights = roots_legendre(5)
    reference = moments(q, omega, pairs)
    errors = []
    for degree in range(7):
        value = sum(w/2*np.sum(pairs*((omega+q[:, None]*u)/2)**degree)
                    for u, w in zip(nodes, weights))
        errors.append(abs(value/reference[str(degree)] - 1))
    return max(errors)


def pair_interpolants(repo, arm, degree, time_order):
    """Interpolate dN/dω on each q→120 interval using Lobatto nodes."""
    q = np.array(history(repo, arm)[0]["q"])
    u = np.cos(np.arange(degree+1)*np.pi/degree)
    omega = q[:, None] + (120.-q)[:, None]*(u+1)/2
    density = number_density(repo, arm, omega, time_order)
    coefficients = dct(density, type=1, axis=1)/degree
    coefficients[:, [0,-1]] /= 2
    return [Chebyshev(row, domain=[qj,120.]) for row,qj in zip(coefficients,q)]


def binned(q, band, curves, edges):
    """Exact angular partial moments within each bin via cumulative pair integrals.

    Integrate over ω first. Along a given parent row the daughter k lies on the
    interval [(ω−q)/2,(ω+q)/2]. For each bin edge K the number and energy
    below K follow from an entire interval and one partially covered interval.
    """
    low, high = band_bounds(q, band)
    cumulative = np.zeros((2, len(edges)))
    for qj,lo,curve in zip(q,low,curves):
        upper_full = np.clip(2*edges-qj, lo, high)
        upper_partial = np.clip(2*edges+qj, lo, high)
        # z_j(x) = ∫ dN ω^j on [lo,x].
        omega_poly = Chebyshev([.5*(120+qj), .5*(120-qj)], domain=[qj,120])
        primitives = [(curve*omega_poly**j).integ() for j in range(3)]
        full = [p(upper_full)-p(lo) for p in primitives]
        part = [p(upper_partial)-p(upper_full) for p in primitives]
        cumulative[0] += full[0] + ((edges+qj/2)*part[0]-part[1]/2)/qj
        cumulative[1] += full[1]/2 + ((edges**2-qj**2/4)*part[0]+qj*part[1]/2-part[2]/4)/(2*qj)
    return np.diff(cumulative)


def log_spectrum_at_midpoints(q,band,curves,k):
    low, high = band_bounds(q, band)
    result = np.zeros_like(k)
    for qj, lo, curve in zip(q,low,curves):
        left = np.clip(2*k-qj,lo,high)
        right = np.clip(2*k+qj,lo,high)
        prim = curve.integ()
        result += (prim(right)-prim(left))/qj
    return k*result


def direct_bins(repo, arm, band, q, left, time_order):
    """Fold independent Fourier evaluations over the exact accepted angles in selected bins."""
    low, high = band_bounds(q, band)
    right = left+.5
    # The accepted angular interval changes form at the four interior candidates.
    shape = (len(q),len(left))
    boundaries = np.stack([np.broadcast_to(low[:,None],shape),
                           np.broadcast_to(high,shape),
                           2*left[None,:]-q[:,None],2*left[None,:]+q[:,None],
                           2*right[None,:]-q[:,None],2*right[None,:]+q[:,None]],axis=-1)
    boundaries = np.sort(np.clip(boundaries,low[:,None,None],high),axis=-1)
    nodes, weights = roots_legendre(24)
    widths = np.diff(boundaries,axis=-1) / 2
    omega = boundaries[:,:,:-1,None] + widths[:,:,:,None]*(nodes+1)
    lower = np.maximum(left[None,:,None,None],(omega-q[:,None,None,None])/2)
    upper = np.minimum(right[None,:,None,None],(omega+q[:,None,None,None])/2)
    fraction = np.maximum(upper-lower,0.)/q[:,None,None,None]
    density = number_density(repo,arm,omega.reshape(len(q),-1),time_order).reshape(omega.shape)
    numbers = widths[:,:,:,None]*weights*density*fraction
    return np.array([numbers.sum(axis=(0,2,3)),(numbers*(upper+lower)/2).sum(axis=(0,2,3))])


def occupancy_bound(repo,arm,time_order):
    """Frequency-independent upper bound from the norm of the ordered amplitude."""
    seed = history(repo,arm)[0]
    _, weights, h, k = drive(repo,arm,time_order)
    X, C, P = (np.asarray(seed["states"]["retain"][name])[:,None] for name in ("X","C","P"))
    instantaneous_variance = X*h**2+2*C*h*k+P*k**2
    amplitude_bound = np.sqrt(instantaneous_variance) @ abs(weights)
    return float(np.dot(seed["weight"],amplitude_bound**2))


def positive_bin_masses(value,tol):
    if np.min(value) < -tol:
        raise ValueError(f"Negative bin integrals: {np.min(value)}")
    return np.maximum(value,0.)


def experiment(repo):
    repo = Path(repo).resolve()
    archive = json.loads((repo/PRIOR/"quantum/SPECTRUM.json").read_text())["arms"]
    result = {"description":"Leading massless scalar Yukawa pair emission, one Dirac species",
              "readout_H0tau":2.10,"bands":{"full":["q",120],"low":["q",60],"high":[60,120]},
              "degeneracy_one_sign":2,"coupling_power_factored_out":2,
              "momenta":"Comoving k in H0 units; single-particle proper p=k/a at readout",
              "common_k_edges":EDGES.tolist(),"arms":{}}
    rows = []
    for arm in ARMS:
        a = history(repo, arm)[-1]
        curves = pair_interpolants(repo,arm,144,5)
        coarse_curves = pair_interpolants(repo,arm,80,3)
        bound = occupancy_bound(repo,arm,5)
        coarse_bound = occupancy_bound(repo,arm,3)
        bands = {}
        for band in BANDS:
            q, omega, weights = pair_nodes(repo,arm,band,96,5)
            raw = moments(q,omega,weights)
            raw_low = moments(*pair_nodes(repo,arm,band,64,3))
            refinement = max(abs(raw[j]/raw_low[j]-1) for j in raw)
            ang_error = angular_controls(q,omega,weights)
            integrated = binned(q,band,curves,EDGES)
            coarse = binned(q,band,coarse_curves,EDGES)
            selected = np.array([0,1,23,52,60,80,100,127])
            folded = direct_bins(repo,arm,band,q,EDGES[selected],5)
            bin_comparison = np.sum(abs(folded-integrated[:,selected]),axis=1)/[raw["0"],raw["1"]]
            if max(bin_comparison)>3e-9:
                raise ValueError("Direct angular folding")
            errors = np.sum(abs(integrated-coarse),axis=1)/[raw["0"],raw["1"]]
            total_errors = abs(integrated.sum(axis=1)/[raw["0"],raw["1"]]-1)
            quiet_ref = archive[arm]
            expected = (quiet_ref["retained_proper"]["emission"] if band=="high" else
                        quiet_ref["positive_frequency_cutoff_checks"]["60" if band=="low" else "120"]["archived"])
            proper_energy = 2*raw["1"]/a
            if abs(proper_energy/expected-1)>2e-7 or refinement>1e-6 or max(errors)>2e-7 or max(total_errors)>2e-8 or ang_error>5e-14:
                raise ValueError(f"Failed {arm}/{band}: {proper_energy} vs {expected}, {refinement}, {errors}, {total_errors}, {ang_error}")
            counts = positive_bin_masses(integrated[0],2e-12*raw["0"])
            energies = positive_bin_masses(integrated[1],2e-11*raw["1"])
            k = .5*(EDGES[:-1]+EDGES[1:])
            per_log = log_spectrum_at_midpoints(q,band,curves,k)
            occupied = counts / (np.diff(EDGES**3)/(3*np.pi**2)) # nψ with two spin states
            if np.max(occupied)>bound or abs(bound/coarse_bound-1)>1e-7:
                raise ValueError("Pointwise occupation envelope")
            for i,(left,right) in enumerate(zip(EDGES[:-1],EDGES[1:])):
                rows.append([arm,band,left,right,left/a,right/a,counts[i],energies[i],
                             per_log[i]/a**3, per_log[i]*k[i]/a**4,occupied[i]])
            bands[band] = {"one_sign_comoving_moments":raw,"pair_count_per_comoving_volume":raw["0"],
                           "both_sign_physical_number_density":2*raw["0"]/a**3,
                           "both_sign_proper_energy_per_comoving_volume":proper_energy,
                           "both_sign_physical_energy_density":2*raw["1"]/a**4,
                           "mean_daughter_proper_momentum":raw["1"]/raw["0"]/a,
                           "maximum_bin_averaged_occupation_coefficient":float(max(occupied)),
                           "max_occupation_bin_comoving_k":EDGES[np.argmax(occupied):np.argmax(occupied)+2].tolist(),
                           "max_raw_moment_relative_refinement":refinement,
                           "count_energy_binned_relative_L1_refinement":errors.tolist(),
                           "binned_count_energy_relative_closure":total_errors.tolist(),
                           "angular_moment_relative_error":ang_error,
                           "eight_direct_bins_absolute_disagreement_over_total":bin_comparison.tolist(),"relative_archived_energy_residual":proper_energy/expected-1}
        additivity = max(abs((bands["low"]["one_sign_comoving_moments"][str(j)] + bands["high"]["one_sign_comoving_moments"][str(j)])/bands["full"]["one_sign_comoving_moments"][str(j)]-1) for j in DEGREES)
        if additivity > 1e-7:
            raise ValueError("Band additivity failed")
        result["arms"][arm] = {"scale_factor":a,"bands":bands,"max_relative_moment_band_additivity":additivity,
                               "pointwise_per_spin_occupation_coefficient_upper_bound":bound,
                               "relative_bound_refinement":abs(bound/coarse_bound-1)}
    return result,rows


def save(study,rows,out):
    out.mkdir(parents=True,exist_ok=True)
    (out/"DISTRIBUTIONS.json").write_text(json.dumps(study,indent=2,allow_nan=False)+"\n")
    with (out/"SPECTRA.csv").open("w",newline="") as f:
        writer=csv.writer(f,lineterminator="\n")
        writer.writerow(["arm","band","k_low","k_high","p_low","p_high","one_sign_comoving_number",
                         "one_sign_conformal_energy","one_sign_dn_proper_dlnp_at_mid",
                         "one_sign_drho_proper_dlnp_at_mid","per_spin_bin_mean_f"])
        writer.writerows(rows)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    save(*experiment(args.repo),args.output)
