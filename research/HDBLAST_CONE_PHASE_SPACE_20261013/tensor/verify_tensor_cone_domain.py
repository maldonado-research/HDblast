#!/usr/bin/env python3
"""Exact algebra receipts for the tensor cone-domain note, not a PDE proof."""
import argparse
import json
from pathlib import Path

import sympy as s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--inject-failure", action="store_true")
    args = parser.parse_args()
    checks = {}
    controls = {}

    def zero(name, expression):
        residual = s.simplify(expression)
        if residual != 0:
            raise RuntimeError(f"{name}: residual {residual}")
        checks[name] = "passed"

    def reject(name, residual):
        witness = s.simplify(residual)
        if witness == 0:
            raise RuntimeError(f"wrong-formula control not rejected: {name}")
        controls[name] = str(witness)

    # Warped curvature in a tangential orthonormal frame. K is the slice
    # sectional curvature, here one. Write Ric_ij=B delta_ij, Ric_yy=A.
    r, r1, r2 = s.symbols("rho rho_prime rho_second", positive=True)
    A = -4 * r2 / r
    B = (3 - r * r2 - 3 * r1**2) / r**2
    scalar = A + 4 * B
    sectional = (1 - r1**2) / r**2
    zero("five_dimensional_tangential_weyl_trace_cancellation",
         sectional - 2 * B / 3 + scalar / 12)
    zero("five_dimensional_radial_weyl_cancellation",
         -r2 / r - (A + B) / 3 + scalar / 12)
    # An all-lowered tangential Weyl component has one rho^2 warp factor,
    # and passing to four orthonormal legs contributes rho^-4.
    C = s.symbols("slice_Weyl")
    zero("orthonormal_weyl_weight", r**2 * C / r**4 - C / r**2)
    reject("incorrect_bounded_weyl_weight", (C / r**2 - C).subs({r: 2, C: 1}))

    # Direct Bianchi I slice witness with trace-free diagonal perturbation.
    t, z = s.symbols("t perturbation")
    f = s.exp(-3 * t)
    polarizations = [s.Integer(1), s.Integer(-1), s.Integer(0)]
    hubbles = [1 + z * e * s.diff(f, t) / 2 for e in polarizations]
    hubble_sum = sum(hubbles)
    ric00 = -sum(s.diff(H, t) + H**2 for H in hubbles)
    zero("slice_linearized_Einstein_tt", s.diff(ric00 + 3, z).subs(z, 0))
    for i, H in enumerate(hubbles):
        ricii = s.diff(H, t) + H * hubble_sum
        zero(f"slice_linearized_Einstein_spatial_{i}",
             s.diff(ricii - 3, z).subs(z, 0))
    zero("slice_graviton_equation", s.diff(f, t, 2) + 3 * s.diff(f, t))
    electric0 = s.diff(-(s.diff(hubbles[0], t) + hubbles[0]**2) + 1,
                       z).subs(z, 0)
    zero("nonzero_slice_Weyl_component", electric0 + s.Rational(3, 2) * f)
    reject("zero_slice_Weyl_claim", electric0.subs(t, 0))

    # Exact conformal null map in the two-dimensional orbit space.
    u, v = s.symbols("u v", positive=True)
    T, R = (v - u) / 2, (v + u) / 2
    zero("conformal_Minkowski_square", R**2 - T**2 - u * v)
    zero("null_metric_uu", -s.diff(T, u)**2 + s.diff(R, u)**2)
    zero("null_metric_vv", -s.diff(T, v)**2 + s.diff(R, v)**2)
    zero("null_metric_uv_coefficient", 2 * (-s.diff(T, u) * s.diff(T, v)
                                            + s.diff(R, u) * s.diff(R, v)) - 1)

    # Quantitative prototype cutoff in the exact rho=y cone. The piecewise
    # linear cutoff is H1; smoothing its corners retains all power laws.
    y, eps = s.symbols("y epsilon", positive=True)
    norm = s.integrate(y**2, (y, 0, eps)) + s.integrate(
        y**2 * (y / eps - 2)**2, (y, eps, 2 * eps))
    radial_form = s.integrate(y**4 / eps**2, (y, eps, 2 * eps))
    zero("radial_cutoff_squared_norm", norm - s.Rational(13, 15) * eps**3)
    zero("radial_cutoff_form", radial_form - s.Rational(31, 5) * eps**3)
    inertial_cutoff = s.integrate(y**3 / eps**2, (y, eps, 2 * eps))
    zero("inertial_cutoff_energy", inertial_cutoff - s.Rational(15, 4) * eps**2)
    zero("inertial_angular_energy_tail", s.integrate(y, (y, 0, eps)) - eps**2 / 2)
    reject("nonvanishing_cutoff_form_limit", s.limit(radial_form, eps, 0) - 1)
    reject("wrong_radial_cutoff_power", s.diff(radial_form, eps, 3) - 0)

    # Smooth exceptional-sector translation vector coefficients.
    x, omega, nu = s.symbols("s omega nu")
    Omega = 1 + omega * x + nu * x**2
    inv = s.series(1 / Omega, x, 0, 3).removeO()
    Avec = 1 + omega * x + (nu - omega**2) * x**2 / 3
    Bvec = -2 * omega + s.Rational(4, 3) * (omega**2 - nu) * x
    zero("regular_translation_ODE_through_second_order",
         Avec - 2 * x * s.diff(Avec, x) - inv)
    zero("regular_translation_B_coefficient", (inv - Avec) / x - Bvec)
    # In a spatial plane with c=e1, ell0 vector's second component is
    # -x1*x2/r: derivative along x1 on the x1 axis is zero, while along a
    # diagonal ray xi2(lambda,lambda)/lambda=-1/sqrt(2), incompatible with
    # the linear derivative inferred on the two axes.
    reject("ell0_residual_is_C1_at_vertex", -1 / s.sqrt(2))
    # Full ambient continuity is stronger than fixed-finite-t boundedness.
    # A spacelike path tending toward a null direction has X1=eps,
    # X0=sqrt(eps^2-eps^6), hence r=eps^3 and xi1=r-X1^2/r.
    ambient_time = s.sqrt(eps**2 - eps**6)
    near_null_radius = s.sqrt(eps**2 - ambient_time**2)
    zero("near_null_path_radius", near_null_radius - eps**3)
    ell_component = near_null_radius - eps**2 / near_null_radius
    zero("near_null_residual_component", ell_component - (eps**3 - 1 / eps))
    scaled_limit = s.limit(eps * ell_component, eps, 0, dir="+")
    zero("near_null_scaled_divergence", scaled_limit + 1)
    reject("ell0_residual_is_ambient_bounded", scaled_limit)
    # Local p=k y: background cancellation integrates exactly to O(y^3).
    k = s.symbols("k")
    zero("exceptional_metric_cancellation", -s.integrate(k**2 * y**2 / 3, y)
         + k**2 * y**3 / 9)

    if args.inject_failure:
        zero("injected_failure", s.Integer(1))
    result = {
        "scope": "Exact algebra and prototype scaling checks; not a global PDE, gauge-completeness, or novelty certificate.",
        "positive_checks": checks,
        "positive_count": len(checks),
        "wrong_formula_controls": controls,
        "wrong_formula_count": len(controls),
        "cutoff_prototype": {
            "squared_radial_norm": str(norm),
            "radial_form": str(radial_form),
            "inertial_cutoff_energy": str(inertial_cutoff),
        },
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
