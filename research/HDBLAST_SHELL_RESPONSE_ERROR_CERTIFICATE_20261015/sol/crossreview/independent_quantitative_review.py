#!/usr/bin/env python3
"""Independent exact-rational cross-check; no production module is imported.

Re-derives adequate shell/collar boxes from the retained analytic endpoint
inequalities, then proves the published coarse interval with simple rational
constants. This checks consequences of analytic hypotheses, not those theorems.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
groups = []


def check(name, condition):
    if not condition:
        raise ArithmeticError(name)
    groups.append(name)


class R:
    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = self.lo if hi is None else Q(hi)
        if self.lo > self.hi:
            raise ArithmeticError("reversed interval")

    @staticmethod
    def r(x):
        return x if isinstance(x, R) else R(x)

    def __add__(self, x):
        x = R.r(x)
        return R(self.lo + x.lo, self.hi + x.hi)

    __radd__ = __add__

    def __neg__(self):
        return R(-self.hi, -self.lo)

    def __sub__(self, x):
        return self + (-R.r(x))

    def __rsub__(self, x):
        return R.r(x) + (-self)

    def __mul__(self, x):
        x = R.r(x)
        products = [u * v for u in (self.lo, self.hi) for v in (x.lo, x.hi)]
        return R(min(products), max(products))

    __rmul__ = __mul__

    def inv(self):
        if self.lo <= 0 <= self.hi:
            raise ArithmeticError("unsafe reciprocal")
        return R(1 / self.hi, 1 / self.lo)

    def __truediv__(self, x):
        return self * R.r(x).inv()

    def __rtruediv__(self, x):
        return R.r(x) * self.inv()

    def sq(self):
        return R(0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2),
                 max(self.lo**2, self.hi**2))

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def out(self):
        return [str(self.lo), str(self.hi)]


def background_polynomials(eta):
    # W is the superpotential here, not the master potential.
    w = Q(1, 3) + eta.sq() + eta.sq() * eta / 3
    w1 = 2 * eta + eta.sq()
    w2 = 2 + 2 * eta
    u1 = w1 * (w2 - Q(4, 3) * w)
    u2 = w2.sq() + 2 * w1 - Q(4, 3) * (w1.sq() + w * w2)
    return w, w1, u1, u2


c = Q(96714321425692, 161839258930731)
detuning = Q(1, 1000)
check("exact_model_decimal", c == 2 / Q("1.0357712571566784") - Q(4, 3))
check("retained_detuning_range", 0 < detuning <= Q(1, 500))
eta0 = -9 * c * detuning / 64
eta = R(eta0 - Q(799, 250) * detuning**2,
        eta0 + Q(799, 250) * detuning**2)
h20 = (1 + c) * detuning / 27 + ((1 + c)**2 / 36 - c**2 / 384) * detuning**2
h2 = R(h20 - Q(127, 500) * detuning**3,
       h20 + Q(127, 500) * detuning**3)
# A deliberately simple enclosure is enough; verify it by squaring.
rho = R(Q(12992, 100), Q(12993, 100))
check("rho_from_endpoint_curvature", 0 < h2.lo and rho.lo**2 < 1 / h2.hi
      and 1 / h2.lo < rho.hi**2)
w, w1, u1, u2 = background_polynomials(eta)
sigma = 2 * w + detuning * (1 + c + c * eta)
sigma1 = 2 * w1 + detuning * c
sigma2 = 4 + 4 * eta
a = sigma / 6
p = -sigma1 / 2
h = u1 / p
b = h - 4 * a + sigma2 / 2
mu = 1 / (rho * b)
shell_c = rho * (2 * h - 5 * a)
beta = mu * R(Q(-25, 4), Q(5, 2)) + shell_c / 2
check("regular_shell_signs", rho.lo > 0 and p.hi < 0 and b.lo > 0)
check("shell_mu_bounds", Q(1, 1000) < mu.lo and mu.hi < Q(1, 400))
check("shell_gamma_bound", (shell_c / 2).mag() < 250)
check("coherent_beta_bound", beta.mag() < 1000)
check("rho_squared_prefactor_bounds", 16641 < rho.sq().lo and rho.sq().hi < 17161)
check("p_squared_prefactor_bounds", Q(9, 625000000) < p.sq().lo
      and p.sq().hi < Q(49, 2500000000))

depth = Q(1, 1250000000)
outer = {
    "rho": R(129, 131), "eta": R(Q(-9, 100000), Q(-8, 100000)),
    "p": R(Q(-14, 100000), Q(-12, 100000)),
    "a": R(Q(111, 1000), Q(112, 1000)),
}
initial = {"rho": rho, "eta": eta, "p": p, "a": a}
rr, ee, pp, aa = (outer[k] for k in ("rho", "eta", "p", "a"))
_, _, uu1, uu2 = background_polynomials(ee)
rhs = {
    "rho": -rr.sq() * aa,
    "eta": -rr * pp,
    "p": -rr * (uu1 - 4 * aa * pp),
    "a": 1 / rr + rr * pp.sq() / 3,
}
first_exit = {}
for key in outer:
    flow = initial[key] + R(0, depth) * rhs[key]
    check("strict_first_exit_" + key,
          outer[key].lo < flow.lo <= flow.hi < outer[key].hi)
    first_exit[key] = flow.out()
hh = uu1 / pp
master_w = rr.sq() * (2 * hh.sq() - 10 * aa * hh + Q(35, 4) * aa.sq()
                       - uu2 - pp.sq() / 6) + R(Q(-35, 4), 0)
check("master_potential_collar_bound", master_w.mag() < 10**6)
check("alpha_positive_collar", pp.hi < 0 and p.hi < 0 and rr.lo > 0)

# Independently evaluate an intentionally coarser arithmetic chain than the
# production certificate. All upper bounds below are strict exact comparisons.
ell = Q(1, 10**10)
T = 3 * ell
L = 8 * ell
M = Q(10**6)
m = Q(1, 1000)
gamma = Q(250)
f_upper = Q(29)
g = 16 / ell + 1
energy_a = Q(1000001)
h_upper = Q(5060000000100)
check("energy_domain_inside_collar", 0 < T < L == depth)
check("bulk_residual_l2_bound_29", M**2 * (1 + T)**2 * L < f_upper**2)
check("energy_a_dominates_boundary_term", energy_a**2 * m > 1001**2)
check("forcing_root_upper", f_upper**2 + g**2 / m < h_upper**2)
x = energy_a * T / 2
exp_upper = 1 / (1 - x)
check("exponential_bound_budget", 0 < x < 1 and exp_upper < Q(1001, 1000))
energy_upper = Q(1520)
check("energy_norm_upper_1520", exp_upper * T * h_upper < energy_upper)
trace_factor_upper = Q(7, 1600)
check("trace_root_upper", T * (1 / m + gamma**2) < trace_factor_upper**2)
delta_upper = T * (f_upper + M * energy_upper) + trace_factor_upper * energy_upper
check("trace_error_below_eight", delta_upper < 8)

# 8-term arctangent intervals are enough for the coarse pi bounds needed here.
def atan_bracket(denominator):
    xx = Q(1, denominator)
    ss = sum(((-1)**j * xx**(2*j+1) / (2*j+1) for j in range(8)), Q(0))
    next_term = xx**17 / 17
    return R(ss, ss + next_term)

pi = 16 * atan_bracket(5) - 4 * atan_bracket(239)
check("coarse_machin_bracket", 3 < pi.lo <= pi.hi < Q(22, 7))
prefactor_lower = 3 * pi.lo**2 / (17161 * Q(49, 2500000000))
prefactor_upper = 3 * pi.hi**2 / (16641 * Q(9, 625000000))
check("simple_physical_prefactor", 80000 < prefactor_lower
      and prefactor_upper < 125000)
check("comparator_norm_lower", 141000**2 < 2 / ell and 141000 - 8 > 140000)
check("comparator_norm_upper", Q(1, 100000)**2 == ell
      and Q(1, 100000) * (16 / ell + 1) + 8 < 1600010)
c_lower = Q(80000) * 140000**2
c_upper = Q(125000) * 1600010**2
check("published_C_lower", c_lower > 10**15)
check("published_C_upper", c_upper < 4 * 10**17)
check("whole_correction_support_margin", Q(3) * ell + Q(5, 4) * Q(7, 2) * ell
      == Q(59, 8) * ell < L)
check("strict_evolution_margin", T < Q(7, 2) * ell)

paths = [
    "SHELL_RESPONSE_ERROR_CERTIFICATE.md", "certify_response.py", "check_response_certificate.py",
    "archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md",
    "archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json", "archive/derive_shell_box.py",
    "archive/sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/theory/JACOBI_ENDPOINT_MONOTONICITY.md",
    "archive/sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/bounds/RESCALED_BOUND_CERTIFICATE.json",
    "sol/RESIDUAL_TO_PHYSICAL_RESPONSE.md",
]
pins = {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}
check("primary_note_expected_pin", pins[paths[0]] ==
      "5f01be536c77123e9f667f4cd5514b52067641b2b389ba4e6eb655cb24654022")
check("retained_theorem_expected_pin", pins[paths[6]] ==
      "f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8")
check("retained_certificate_expected_pin", pins[paths[7]] ==
      "22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d")
check("original_sol_note_unchanged", pins[paths[8]] ==
      "b8c77c810c70aa286227c5dabe010a53942af84067eb03961f0186419a5ea96a")

print(json.dumps({
    "status": "PASS_CONDITIONAL_QUANTITATIVE_CROSS_REVIEW",
    "check_count": len(groups), "checks": groups,
    "production_module_imported": False,
    "physical_waveform_computed": False,
    "analytic_hypotheses_certified_by_program": False,
    "parameters_exact": {"ell": str(ell), "T": str(T), "L_energy": str(L)},
    "independent_shell_boxes": {"rho": rho.out(), "eta": eta.out(), "p": p.out(),
                                 "a": a.out(), "B": b.out(), "mu": mu.out()},
    "first_exit_flow_boxes": first_exit,
    "master_potential_collar_box": master_w.out(),
    "delta_upper_independent_exact": str(delta_upper),
    "C_lower_independent_strict": str(c_lower),
    "C_upper_independent_strict": str(c_upper),
    "pins_sha256": pins,
}, indent=2))
