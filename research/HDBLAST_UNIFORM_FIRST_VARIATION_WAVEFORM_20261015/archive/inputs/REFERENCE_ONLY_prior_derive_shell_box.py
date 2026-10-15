#!/usr/bin/env python3
"""Exact rational consequences of the retained analytic branch estimates.

This newly written program reads frozen references as data. It does not import
or execute any archived program, solve the background ODE, or certify the
analytic theorem on which the enclosures depend. Floats are display only.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
PROOF = "sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/theory/JACOBI_ENDPOINT_MONOTONICITY.md"
CERT = "sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/bounds/RESCALED_BOUND_CERTIFICATE.json"
PINS = {
    PROOF: "f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8",
    CERT: "22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d",
}


class I:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def box(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, other):
        o = self.box(other)
        return I(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.box(other)

    def __rsub__(self, other):
        return self.box(other) + -self

    def __mul__(self, other):
        o = self.box(other)
        candidates = [a * b for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(candidates), max(candidates))

    __rmul__ = __mul__

    def recip(self):
        if self.lo <= 0 <= self.hi:
            raise ValueError("division by interval containing zero")
        return I(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * self.box(other).recip()

    def __rtruediv__(self, other):
        return self.box(other) * self.recip()

    def sq(self):
        if self.lo <= 0 <= self.hi:
            return I(0, max(self.lo**2, self.hi**2))
        return I(min(self.lo**2, self.hi**2), max(self.lo**2, self.hi**2))

    def record(self):
        return {"exact": [str(self.lo), str(self.hi)],
                "display_only": [float(self.lo), float(self.hi)]}


def sqrt_bounds(q, bits=80):
    if q < 0:
        raise ValueError("negative square root")
    den = 1 << bits
    n = isqrt((q.numerator << (2 * bits)) // q.denominator)
    lo, hi = F(n, den), F(n + 1, den)
    if not lo * lo <= q <= hi * hi:
        raise ValueError("square-root enclosure failed")
    return lo, hi


def collar_bounds(shell):
    """Verified first-exit tube for the actual smooth branch, not a solver."""
    depth = F(8, 10**10)
    outer = {
        "rho": I(129, 131),
        "eta": I(F(-9, 10**5), F(-8, 10**5)),
        "p": I(F(-14, 10**5), F(-12, 10**5)),
        "a": I(F(111, 1000), F(112, 1000)),
    }
    initial = {"rho": shell["rho_b"], "eta": shell["eta_b"],
               "p": shell["p_b"], "a": shell["a_b"]}
    rho, e, p, a = (outer[k] for k in ("rho", "eta", "p", "a"))
    W = F(1, 3) + e.sq() + e.sq() * e / 3
    W1 = 2 * e + e.sq()
    U1 = W1 * (2 + 2 * e - F(4, 3) * W)
    U2 = 4 * (1 + e).sq() + 2 * W1 - F(4, 3) * (W1.sq() + 2 * (1 + e) * W)
    rhs = {"rho": -rho.sq() * a,
           "eta": -rho * p,
           "p": -rho * (U1 - 4 * a * p),
           "a": 1 / rho + rho * p.sq() / 3}
    first_exit = {}
    for name, box in outer.items():
        enclosure = initial[name] + I(0, depth) * rhs[name]
        if not (box.lo < enclosure.lo <= enclosure.hi < box.hi):
            raise ValueError("first-exit margin failed: " + name)
        first_exit[name] = {
            "outer_box": box.record(), "initial_box": initial[name].record(),
            "rhs_enclosure": rhs[name].record(),
            "flow_integral_enclosure": enclosure.record(),
            "strict_lower_margin": str(enclosure.lo - box.lo),
            "strict_upper_margin": str(box.hi - enclosure.hi),
        }
    h = U1 / p
    C = rho * (2 * h - 5 * a)
    Cz = rho.sq() * (2 * h.sq() - 10 * a * h + 5 * a.sq() - 2 * U2 - F(5, 3) * p.sq()) - 5
    spatial = rho.sq() * (2 * h.sq() - 10 * a * h + F(35, 4) * a.sq() - U2 - p.sq() / 6)
    potential = spatial + I(F(-35, 4), 0)
    mu = shell["mu"]
    boundary_beta = mu * I(F(-25, 4), F(5, 2)) + shell["C_b"] / 2
    bounds = {
        "mu_lower": F(1, 1000), "mu_upper": F(1, 400),
        "absolute_W_upper": F(10**6),
        "absolute_boundary_beta_upper": F(1000),
        "absolute_C_b_over_two_upper": F(250),
        "absolute_C_over_two_collar_upper": F(500),
        "rho_b_squared_lower": F(129**2),
        "rho_b_squared_upper": F(131**2),
        "p_b_squared_lower": F(12, 10**5)**2,
        "p_b_squared_upper": F(14, 10**5)**2,
    }
    def mag(i):
        return max(abs(i.lo), abs(i.hi))
    gates = {
        "positive_mu": bounds["mu_lower"] < mu.lo,
        "bounded_mu": mu.hi < bounds["mu_upper"],
        "bounded_W_for_all_real_t": mag(potential) < bounds["absolute_W_upper"],
        "bounded_boundary_beta_for_all_real_t": mag(boundary_beta) < bounds["absolute_boundary_beta_upper"],
        "bounded_C_b_over_two": mag(shell["C_b"] / 2) < bounds["absolute_C_b_over_two_upper"],
        "bounded_C_collar_over_two": mag(C / 2) < bounds["absolute_C_over_two_collar_upper"],
        "rho_b_squared_lower": bounds["rho_b_squared_lower"] < shell["rho_b"].sq().lo,
        "rho_b_squared_upper": shell["rho_b"].sq().hi < bounds["rho_b_squared_upper"],
        "p_b_squared_lower": bounds["p_b_squared_lower"] < shell["p_b"].sq().lo,
        "p_b_squared_upper": shell["p_b"].sq().hi < bounds["p_b_squared_upper"],
    }
    if not all(gates.values()):
        raise ValueError("coarse coefficient gate failed: " + repr(gates))
    return {
        "status": "PASS_EXACT_FIRST_EXIT_TUBE_AND_COEFFICIENT_BOUNDS",
        "inward_depth_exact": str(depth),
        "first_exit_inclusions": first_exit,
        "coefficient_boxes": {"C_collar": C.record(), "C_z_collar": Cz.record(),
                              "W_all_real_t_collar": potential.record(),
                              "boundary_beta_all_real_t": boundary_beta.record()},
        "simple_bounds_exact": {k: str(v) for k, v in bounds.items()},
        "gates": gates,
        "argument": "A smooth ODE solution starting in the initial box cannot first exit the outer box before this depth, because its time-integrated derivative is strictly interior on every such partial interval. The rational vector field is smooth near the compact outer box with rho>0. ODE continuation therefore covers the full depth. This is a conditional bound for the actual analytic branch, not a computed trajectory.",
    }


def main():
    for rel, expected in PINS.items():
        if sha256((ROOT / rel).read_bytes()).hexdigest() != expected:
            raise ValueError("source pin mismatch: " + rel)
    cert = json.loads((ROOT / CERT).read_text())
    primary = cert["primary_uniform_certificate"]
    const = primary["exact_constants"]
    if not primary["pass_all_sufficient_gates"]:
        raise ValueError("retained certificate did not pass")
    # The outward decimal constants in the analytic theorem, Eq. (29),
    # strictly exceed the exact retained bounds.
    if not (F(const["K_eta"]) < F(799, 250)
            and F(primary["K_H_best"]) < F(127, 500)):
        raise ValueError("outward theorem constants failed")
    c, delta = F(96714321425692, 161839258930731), F(1, 1000)
    if not 0 < delta <= F(primary["delta_star"]):
        raise ValueError("outside certified detuning range")
    e0 = -9 * c * delta / 64
    e = I(e0 - F(799, 250) * delta**2,
          e0 + F(799, 250) * delta**2)
    h0 = (1 + c) * delta / 27 + ((1 + c)**2 / 36 - c**2 / 384) * delta**2
    H2 = I(h0 - F(127, 500) * delta**3, h0 + F(127, 500) * delta**3)
    if H2.lo <= 0:
        raise ValueError("curvature interval is not positive")
    rho = I(sqrt_bounds(1 / H2.hi)[0], sqrt_bounds(1 / H2.lo)[1])
    W = F(1, 3) + e.sq() + e.sq() * e / 3
    W1 = 2 * e + e.sq()
    sigma = 2 * W + delta * (1 + c + c * e)
    sigma1 = 2 * W1 + delta * c
    sigma2 = 4 + 4 * e
    a, p = sigma / 6, -sigma1 / 2
    U1 = W1 * (2 + 2 * e - F(4, 3) * W)
    U2 = 4 * (1 + e).sq() + 2 * W1 - F(4, 3) * (W1.sq() + 2 * (1 + e) * W)
    h = U1 / p
    g = h - 4 * a
    B = g + sigma2 / 2
    C = rho * (2 * h - 5 * a)
    Cz = rho.sq() * (2 * h.sq() - 10 * a * h + 5 * a.sq() - 2 * U2 - F(5, 3) * p.sq()) - 5
    spatial_potential = rho.sq() * (2 * h.sq() - 10 * a * h + F(35, 4) * a.sq() - U2 - p.sq() / 6)
    if not (p.hi < 0 and B.lo > 0 and rho.lo > 0):
        raise ValueError("required regular-branch signs not established")
    boxes = dict(eta_b=e, H_squared=H2, rho_b=rho, sigma_b=sigma,
                 sigma_prime_b=sigma1, sigma_second_b=sigma2,
                 a_b=a, p_b=p, U_prime_b=U1, U_second_b=U2,
                 U_prime_over_p_b=h, p_y_over_p_b=g, B_b=B,
                 mu=1 / (rho * B), C_b=C, C_z_b=Cz,
                 W_at_t_zero_and_z_zero=spatial_potential)
    output = {
        "status": "PASS_EXACT_RATIONAL_CONSEQUENCES_OF_RETAINED_ANALYTIC_THEOREM",
        "scope": "Conditional analytic shell box; not a validated ODE trajectory or waveform.",
        "parameters_exact": {"delta": str(delta), "c": str(c)},
        "input_pins": PINS,
        "bounds_used": {"K_eta_upper": "799/250", "K_H_upper": "127/500"},
        "independent_square_root_enclosure_bits": 80,
        "required_signs": {"rho_b_positive": True, "p_b_negative": True,
                           "B_b_positive": True},
        "boxes": {k: v.record() for k, v in boxes.items()},
        "collar": collar_bounds(boxes),
    }
    path = ROOT / "ANALYTIC_SHELL_BOX_DELTA_1E3.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": output["status"],
                      "boxes_display_only": {k: v.record()["display_only"] for k, v in boxes.items()},
                      "collar_status": output["collar"]["status"],
                      "collar_depth": output["collar"]["inward_depth_exact"],
                      "simple_bounds_exact": output["collar"]["simple_bounds_exact"]}, indent=2))


if __name__ == "__main__":
    main()
