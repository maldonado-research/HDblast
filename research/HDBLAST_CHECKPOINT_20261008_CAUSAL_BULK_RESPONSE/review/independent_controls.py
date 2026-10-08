#!/usr/bin/env python3
"""Independent manufactured controls; no project data or producer code is read.

The nontrivial numerical tests are canonical spectral sum rules and convergence
of an independently terminated finite slab. These are diagnostic controls, not
interval proofs. All model values below are manufactured dimensionless units.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 65
OUT = Path(__file__).resolve().parent / "INDEPENDENT_CONTROLS_RESULT.json"
checks = []


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def exact(name, expression):
    residual = sp.simplify(expression)
    require(residual == 0, (name, residual))
    checks.append({"name": name, "kind": "exact_symbolic", "residual": "0"})


def numeric(name, got, expected, tolerance="1e-45"):
    error = abs(got - expected)
    require(error < mp.mpf(tolerance), (name, got, expected, error))
    checks.append({"name": name, "kind": "high_precision_diagnostic",
                   "got": str(got), "expected": str(expected),
                   "absolute_error": str(error), "tolerance": tolerance})


# A bound eigenfunction checks both variational signs and the residue by an
# independent Hilbert-space normalization, without differentiating D(z).
t, y = sp.symbols("t y", real=True)
omega = sp.sqrt(3)
q = sp.cos(omega*t)
phi = sp.exp(-y)*q/4
M, c, g, m2 = sp.Integer(2), sp.Integer(3), sp.Integer(1), sp.Rational(13, 4)
trace = phi.subs(y, 0)
exact("bound_mode_bulk_equation", sp.diff(phi, t, 2)-sp.diff(phi, y, 2)+M**2*phi)
exact("bound_mode_robin_equation", sp.diff(phi, y).subs(y, 0)-c*trace+g*q)
exact("bound_mode_brane_equation", sp.diff(q, t, 2)+m2*q-g*trace)
bulk_profile_norm = sp.integrate(sp.exp(-2*y)/16, (y, 0, sp.oo))
brane_weight = 1/(1+bulk_profile_norm)
exact("bound_residue_from_normalized_eigenvector", brane_weight-sp.Rational(32, 33))
Ebulk = sp.integrate((sp.diff(phi,t)**2+sp.diff(phi,y)**2+M**2*phi**2)/2, (y,0,sp.oo))
Eq = (sp.diff(q,t)**2+m2*q**2)/2
ER = c*trace**2/2
Eint = -g*q*trace
exact("bulk_energy_boundary_flux", sp.diff(Ebulk,t)+sp.diff(trace,t)*sp.diff(phi,y).subs(y,0))
exact("robin_reservoir_exchange", sp.diff(Ebulk+ER,t)-g*q*sp.diff(trace,t))
exact("brane_energy_exchange", sp.diff(Eq,t)-g*trace*sp.diff(q,t))
exact("total_energy_constant", Ebulk+Eq+ER+Eint-sp.Rational(99,64))
wrong_robin_residual = sp.simplify(sp.diff(phi,y).subs(y,0)+c*trace-g*q)
require(wrong_robin_residual != 0, "Wrong Robin sign unexpectedly accepted")
checks.append({"name":"wrong_robin_sign_detected", "kind":"negative_control",
               "rejected_residual":str(wrong_robin_residual)})
missing_interaction_drift = sp.simplify(sp.diff(Ebulk+ER+Eq,t))
require(missing_interaction_drift != 0, "Omitted interaction energy unexpectedly accepted")
checks.append({"name":"omitted_interaction_energy_detected", "kind":"negative_control",
               "rejected_energy_drift":str(missing_interaction_drift)})

# This static trial subspace exactly saturates the trace inequality. Its
# Schur complement establishes sharpness of the stability boundary; a finite
# Robin form term of the wrong sign fails this control.
a, q0, mm, cc, mass, gg = sp.symbols("a q0 mm cc mass gg", positive=True)
static_form = (mass+cc)*a**2 + mm*q0**2 - 2*gg*a*q0
exact("sharp_static_schur_complement",
      static_form.subs(a,gg*q0/(mass+cc)) - q0**2*(mm-gg**2/(mass+cc)))
exact("trace_saturation",
      sp.integrate(sp.diff(a*sp.exp(-mass*y),y)**2+mass**2*(a*sp.exp(-mass*y))**2,
                   (y,0,sp.oo))-mass*a**2)

# The full brane field is canonical. Independent spectral completeness
# requires integral rho_q = 1 and integral u*rho_q = m^2, including bound
# mass and residue. The positive numerical measure comes from outgoing
# scattering amplitudes directly in momentum p, not a copied rho(u) routine.
def continuum_moment(mass_squared, power):
    def integrand(p):
        u = 4+p*p
        outgoing_trace = 1/(3-1j*p)
        brane_response = 1/(mass_squared-u-outgoing_trace)
        scattering_weight = (2/mp.pi)*p*p/(9+p*p)*abs(brane_response)**2
        return u**power*scattering_weight
    # Resolve the manufactured above-threshold resonance in the no-bound
    # case explicitly; the tails are integrated, never silently truncated.
    return mp.quad(integrand, [0, mp.mpf('0.1'), 1, 2, mp.mpf('2.1'),
                               mp.mpf('2.2'), mp.mpf('2.3'), mp.mpf('2.4'), 3, 10, mp.inf])

bound_Z = mp.mpf(32)/33
bound_mass2 = mp.mpf(3)
for label, mm_value, z, weight in [
        ("bound", mp.mpf(13)/4, bound_mass2, bound_Z),
        ("no_bound", mp.mpf(9), mp.mpf(0), mp.mpf(0))]:
    for power, total_expected in [(0, mp.mpf(1)), (1, mm_value)]:
        numeric(f"canonical_{label}_spectral_moment_{power}",
                continuum_moment(mm_value,power)+weight*z**power, total_expected)

# A finite slab with Neumann condition at y=L has a boundary response
# 1/[c+s*tanh(s L)]. Derive it from its profile and boundary force, then
# check convergence to the outgoing half-space for an upper-half-plane
# frequency. This is independent of the continuum spectral integral.
s, length = sp.symbols("s length", positive=True)
profile = sp.cosh(s*(length-y))/sp.cosh(s*length)
exact("finite_slab_right_neumann", sp.diff(profile,y).subs(y,length))
exact("finite_slab_boundary_denominator",
      c*profile.subs(y,0)-sp.diff(profile,y).subs(y,0)-(c+s*sp.tanh(s*length)))
w = mp.mpc('1.1','0.7')
k = mp.mpf('0.6')
ss = mp.sqrt(k*k+4-w*w)
require(ss.real > 0, "Upper-half-plane root must decay into bulk")
halfspace = 1/(3+ss)
slab_errors = []
for slab_L in [1,2,4,8,16,32]:
    response = 1/(3+ss*mp.tanh(ss*slab_L))
    slab_errors.append({"L":slab_L,"absolute_error":str(abs(response-halfspace))})
require(all(mp.mpf(slab_errors[i+1]["absolute_error"]) < mp.mpf(slab_errors[i]["absolute_error"])
            for i in range(len(slab_errors)-1)), "Finite slab errors did not decrease")
numeric("finite_slab_outgoing_limit", 1/(3+ss*mp.tanh(ss*32)), halfspace)
checks.append({"name":"finite_slab_convergence_sequence", "kind":"high_precision_diagnostic",
               "frequency":str(w), "momentum":str(k), "values":slab_errors})

# Direct outgoing-wave work checks the physically decisive sign. All
# amplitudes are manufactured; the result is independent of the bound test.
ww, pp, QQ, coupling, robin = sp.symbols("ww pp QQ coupling robin", positive=True)
bulk_amplitude = coupling*QQ/(robin-sp.I*pp)
flux = sp.re((-sp.I*ww*bulk_amplitude)*sp.conjugate(sp.I*pp*bulk_amplitude)*(-1))/2
self_force_work = sp.re(coupling*bulk_amplitude*sp.conjugate(-sp.I*ww*QQ))/2
expected_flux = ww*coupling**2*pp*QQ**2/(2*(robin**2+pp**2))
exact("outgoing_noether_flux", flux-expected_flux)
exact("self_force_work_equals_negative_flux", self_force_work+expected_flux)

result = {"status":"PASS", "scope":"manufactured analytic controls only",
          "protected_project_data_read":False, "producer_code_imported":False,
          "precision_decimal_digits":mp.mp.dps,
          "versions":{"sympy":sp.__version__, "mpmath":mp.__version__},
          "numerical_results_are_interval_certificates":False,
          "checks":checks}
OUT.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"], "check_count":len(checks), "result_path":str(OUT)}))
