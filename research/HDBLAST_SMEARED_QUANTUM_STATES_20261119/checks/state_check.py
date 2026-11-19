"""Time-smeared minimum energy for the nonzero harmonics of unit-h closed de Sitter."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import roots_legendre


VOLUME = 2 * np.pi**2
MAX_L = 16384


def sample(T, order):
    """Positive half of the even measure w(z) dz, z=h*x, expressed in eta."""
    nodes, weights = roots_legendre(order)
    half_eta = np.arctan(np.sinh(T)) / 2
    eta = half_eta * (1 + nodes)
    weights = 2 * half_eta * weights
    c, s = np.cos(eta), np.sin(eta)
    z = np.arcsinh(np.tan(eta))
    one_minus_r2 = (1 - z/T) * (1 + z/T)
    bump = np.exp(-1 / one_minus_r2)
    normalization = np.dot(weights, bump / c)
    w = bump / normalization
    wz = -2*z/T**2 * w / one_minus_r2**2
    F = w*c
    Fprime = wz-w*s
    integrand = 2*F*np.cos(2*eta) + 0.5*Fprime*np.sin(2*eta)
    return eta, weights*w/c, weights*integrand


def gauss_forms(T, order, labels):
    eta, quadrature, by_parts = sample(T, order)
    c2 = np.cos(eta)**2
    mean_c4 = np.dot(quadrature, c2**2)
    mean_cos2c2 = np.dot(quadrature, c2*np.cos(2*eta))
    L = np.asarray(labels, dtype=float)
    diagonal = (2*L**2*mean_c4 - mean_cos2c2)/(4*L)
    anomalous = np.empty_like(L)
    for start in range(0, len(L), 64):
        block = L[start:start+64]
        anomalous[start:start+64] = (np.cos(2*block[:, None]*eta) @ by_parts)/(4*block)
    return diagonal, anomalous


def fourier_forms(T, count, L):
    """Periodic trapezoid on a conformal-time interval containing the extended bump."""
    eta = np.pi*np.fft.fftfreq(count)
    inside = abs(eta) < np.arctan(np.sinh(T))
    e = eta[inside]
    c, s = np.cos(e), np.sin(e)
    z = np.arcsinh(np.tan(e))
    one_minus_r2 = (1-z/T)*(1+z/T)
    bump = np.exp(-1/one_minus_r2)
    spacing = np.pi/count
    w = bump/(spacing*np.sum(bump/c))
    wz = -2*z/T**2*w/one_minus_r2**2
    envelope = np.zeros(count)
    envelope[inside] = 2*w*c*np.cos(2*e) + (wz-w*s)*np.sin(2*e)/2
    anomalous = spacing*np.fft.rfft(envelope).real[L.astype(int)]/(4*L)
    mean_c4 = spacing*np.sum(w*c**3)
    mean_cos2c2 = spacing*np.sum(w*c*np.cos(2*e))
    diagonal = (2*L**2*mean_c4-mean_cos2c2)/(4*L)
    return diagonal, anomalous


def state(diagonal, anomalous):
    minimum = np.sqrt((diagonal-anomalous)*(diagonal+anomalous))
    ratio = -anomalous/(diagonal+minimum)
    alpha = 1 / np.sqrt((1-ratio)*(1+ratio))
    beta = alpha*ratio
    # Difference avoids loss of the small gain for large L.
    gain = -anomalous**2/(diagonal+minimum)
    return alpha, beta, gain


def modes(L, z):
    eta = np.arctan(np.sinh(z))
    c, s = np.cos(eta), np.sin(eta)
    lam = L**2-1
    pre = np.exp(-1j*L*eta)/np.sqrt(2*L*lam)
    field = pre*(L*c+1j*s)
    velocity = -1j*lam*c**2*pre
    gradient = np.sqrt(lam)*c*field
    return field, velocity, gradient


def stress(L, alpha, beta, z):
    eta = np.arctan(np.sinh(z))
    c = np.cos(eta)
    _, velocity, gradient = modes(L, z)
    # Stable expression for the anomalous kinetic-plus-gradient combination.
    rho_cross = c**2/(4*L)*(np.cos(2*eta)+1j*L*np.sin(2*eta))*np.exp(-2j*L*eta)
    p_cross = (velocity**2-gradient**2/3)/2
    rho = beta**2*(abs(velocity)**2+abs(gradient)**2) + 2*alpha*beta*rho_cross.real
    p = beta**2*(abs(velocity)**2-abs(gradient)**2/3) + 2*alpha*beta*p_cross.real
    return L**2/VOLUME*rho, L**2/VOLUME*p


def direct_minimum(T, labels, alpha, beta, gain):
    """Check the complex field kinetic and gradient energies in a separate quadrature."""
    times, weights = roots_legendre(360)
    weights *= np.exp(-1/(1-times**2))
    weights /= weights.sum()
    z = T*times
    discrepancies = []
    for L, a, b, g in zip(labels, alpha, beta, gain):
        _, v, k = modes(L, z)
        def cost(r):
            squeezed_v = np.cosh(r)*v + np.sinh(r)*v.conj()
            squeezed_k = np.cosh(r)*k + np.sinh(r)*k.conj()
            return np.dot(weights, (abs(squeezed_v)**2+abs(squeezed_k)**2)/2)
        # Bracket the minimum with a derivative-free optimizer.
        found = minimize_scalar(cost, method="bounded", bounds=(-.2,.2),
                                options={"xatol":1e-11})
        c_cross = np.dot(weights, (v*v+k*k)/2)
        v_s = a*v+b*v.conj()
        k_s = a*k+b*k.conj()
        transformed_cross = np.dot(weights, (v_s*v_s+k_s*k_s)/2)
        discrepancies.append([abs(np.arcsinh(b)-found.x),
                              abs((cost(np.arcsinh(b))-cost(0))-g),
                              abs(c_cross.imag), abs(transformed_cross)])
    return dict(zip(("optimizer_squeeze_error", "objective_gain_error",
                     "reference_anomalous_imaginary", "stationary_anomalous_modulus"),
                    np.max(discrepancies, axis=0).tolist()))


def run_width(T):
    L = np.arange(2, MAX_L+1, dtype=float)
    alpha, beta, gain = state(*fourier_forms(T, 131072, L))
    _, beta2, _ = state(*fourier_forms(T, 65536, L))
    _, gauss_beta, _ = state(*gauss_forms(T, 1300, L[:4095]))
    checks = direct_minimum(T, L[:24], alpha[:24], beta[:24], gain[:24])
    checks["max_beta_quadrature_difference"] = float(max(abs(beta-beta2)))
    checks["max_gauss_fourier_beta_difference_through_n4095"] = float(max(abs(beta[:4095]-gauss_beta)))
    sampled = np.array([0., .1, .2, .3, .5, 1.])
    table = []
    symmetry, wronskian, continuity = 0., 0., 0.
    for z in sampled:
        rho,p = stress(L, alpha, beta, z)
        rho2,p2 = stress(L, np.hypot(1,beta2), beta2, z)
        cutoff = {}
        for max_label in (63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383):
            stop = max_label
            cutoff[str(max_label)] = [float(rho[:stop].sum()), float(p[:stop].sum())]
        table.append(dict(hx=z, sums_by_max_n=cutoff,
                          abs_density_pressure_quadrature_differences=[
                              abs(float((rho-rho2).sum())),abs(float((p-p2).sum()))]))
        rhominus,pminus = stress(L, alpha, beta, -z)
        symmetry = max(symmetry,max(abs(rho-rhominus)),max(abs(p-pminus)))
        q,qdot,_=modes(L,z)
        field = alpha*q+beta*q.conj()
        velocity = alpha*qdot+beta*qdot.conj()
        wronskian = max(wronskian,max(abs(1j*np.cosh(z)**3*(
            field.conj()*velocity-velocity.conj()*field)-1)))
        step = 1e-5
        totals = [stress(L, alpha, beta, z+i*step)[0].sum() for i in (-2,-1,1,2)]
        derivative = (totals[0]-8*totals[1]+8*totals[2]-totals[3])/(12*step)
        continuity=max(continuity,abs(derivative+3*np.tanh(z)*(rho+p).sum()))
    checks.update(max_time_reflection_shell_error=float(symmetry),
                  max_squeezed_wronskian_error=float(wronskian),
                  max_finite_difference_continuity_residual=float(continuity))
    thresholds = dict(optimizer_squeeze_error=1e-7, objective_gain_error=2e-12,
                      reference_anomalous_imaginary=1e-13, stationary_anomalous_modulus=2e-12,
                      max_beta_quadrature_difference=2e-12, max_gauss_fourier_beta_difference_through_n4095=2e-12,
                      max_time_reflection_shell_error=1e-13,
                      max_squeezed_wronskian_error=1e-12, max_finite_difference_continuity_residual=1e-8)
    for key, tol in thresholds.items():
        if checks[key] > tol:
            raise ArithmeticError(f"{T=}, {key} = {checks[key]:g} > {tol:g}")
    samples = {}
    for l in (2,3,4,5,8,16,32,64,128,256,512,1024,2048,4096):
        idx=l-2
        samples[str(l-1)] = dict(beta=float(beta[idx]), occupation=float(beta[idx]**2))
    return dict(hT=T, samples_by_n=samples, stress=table,
                sampled_density_difference=float(np.dot(L**2, gain)/VOLUME),
                reference_number_in_squeezed_state=float(np.dot(L**2,beta**2)), checks=checks)


def width_limit():
    """Check the pressure coefficient against an integral over the fixed profile."""
    nodes, weights = roots_legendre(256)
    g = np.exp(-1/(1-nodes**2))
    # Integral of k*g_hat(k) on the positive line = -integral of g'(y)/y.
    fourier_moment = 2*np.dot(weights, g/(1-nodes**2)**2)/np.dot(weights, g)
    pressure_coefficient = fourier_moment/(16*np.pi**2)
    rows = []
    for width in (.1, .05, .025, .0125, .00625, .003125):
        last_L = round(800/width)
        L = np.arange(2, last_L+1, dtype=float)
        fft_size = 2**(4*last_L).bit_length()
        alpha, beta, _ = state(*fourier_forms(width, fft_size, L))
        pressure_shell = L/(2*VOLUME)*(beta**2*(2*L**2/3-1)-alpha*beta*(4*L**2/3-1))
        density_shell = L/(2*VOLUME)*(beta**2*(2*L**2-1)+alpha*beta)
        density = float(density_shell.sum())
        row = dict(hT=width, max_n=last_L-1, bounce_density=density,
                   hT_squared_times_bounce_pressure=float(width**2*pressure_shell.sum()))
        if rows:
            row["density_log_slope_over_prediction"] = (
                (rows[-1]["bounce_density"]-density)*32*np.pi**2/np.log(2))
        rows.append(row)
    error = abs(rows[-1]["hT_squared_times_bounce_pressure"]/pressure_coefficient-1)
    if error > 1e-5 or abs(rows[-1]["density_log_slope_over_prediction"]-1) > 5e-4:
        raise ArithmeticError("Narrow-width pressure or density coefficient disagreed")
    return dict(pressure_coefficient=float(pressure_coefficient),
                density_logarithm_coefficient=float(-1/(32*np.pi**2)), rows=rows)


def zero_mode():
    # Unit example h=G=sigma_q=hbar=1 for a dimensionless identity check.
    sigma_q = 1.
    variance_p = 1/(4*sigma_q**2)
    z = np.linspace(-1.25,1.25,1601)
    shear = (1/np.cosh(z)*np.tanh(z)+np.arctan(np.sinh(z)))/2
    variance_q = sigma_q**2+shear**2*variance_p
    covariance = shear*variance_p
    determinant = variance_q*variance_p-covariance**2
    return dict(h_G_sigma_q_hbar_equal_one=True, bounce_density=variance_p/(2*VOLUME),
                bounce_fraction_of_de_sitter_density=variance_p/(2*VOLUME)*8*np.pi/3,
                max_covariance_determinant_error=float(max(abs(determinant-.25))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = dict(h=1, physical_density_pressure_units="h^4, hbar=1",
                  proper_time_variable="z=h*x", samples_weight="normalized exp[-1/(1-(z/hT)^2)]",
                  maximum_n=MAX_L-1, degeneracy="(n+1)^2",
                  stress_definition="squeezed state minus nonzero-mode exact reference, per physical volume",
                  runs=[run_width(width) for width in (.1,.2,.3)], narrow_window=width_limit(),
                  homogeneous_zero_mode=zero_mode())
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
