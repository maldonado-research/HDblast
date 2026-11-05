#!/usr/bin/env python3
"""Exact finite display audit; independent of the scientific producer."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json
import re

HERE = Path(__file__).resolve().parent
STEM = "PAST_SHELL_CONTACT"
REPORT_SHA = "4172579bd566e03be1e1ea09c6254dc3ef92e786d6b0aca6863bd528a9552f19"
PROOF_SHA = "d0fda29c840cc6aec840e2a832cc1f0d2c6b3666ae825b40e5e141519b463a79"


def main():
    ap = argparse.ArgumentParser()
    for name in ("report", "optimized", "proof", "output"):
        ap.add_argument("--" + name, required=True, type=Path)
    args = ap.parse_args()
    tex = HERE / (STEM + ".tex")
    source = tex.read_text()
    report = json.loads(args.report.read_text())
    checks = []

    def need(ok, name):
        if not ok:
            raise ValueError(name)
        checks.append(name)

    def val(path):
        obj = report
        for part in path.split("."):
            obj = obj[part]
        return F(obj)

    def interval(path):
        obj = report
        for part in path.split("."):
            obj = obj[part]
        return tuple(map(F, obj))

    def printed(token, name):
        need(token in source, name + " printed")

    need(sha256(args.report.read_bytes()).hexdigest() == REPORT_SHA,
         "frozen primary canonical report pin")
    need(args.report.read_bytes() == args.optimized.read_bytes(),
         "normal and optimized reports byte identical")
    need(sha256(args.proof.read_bytes()).hexdigest() == PROOF_SHA,
         "frozen mathematical proof pin")
    need(report["status"] == "PASS_CONDITIONAL_PAST_SHELL_CONTACT_SCALARS",
         "conditional actual-packet report passes")
    need(report["finite_checks"] == 27 == len(report["check_names"]),
         "primary finite-check roster")

    ell = F(1, 10**6)
    printed(r"\ell=10^{-6}", "registered width")
    printed(r"w_t(0,z)=e^{-z}q_\ell(z)", "actual nonzero initial velocity")
    printed(r"\frac{96714321425692}{161839258930731}", "fixed model parameter")
    printed(r"\sigma''=4\phi", "tension second derivative")
    printed(r"-w_z+\mu_b w_{tt}+\beta(t)w=0", "actual signed boundary law")
    printed(r"\mu_b=\frac1{\rho_b B_b}", "actual boundary mass definition")
    for key, lo, hi in [
            ("mu_b", F(1, 1000), F(1, 400)),
            ("C_b_over_two", F(189), F(300)-F(1, 160)),
            ("beta_all_real_t", F(180), F(300)),
            ("Lambda", F(1000), F(2500)),
            ("b_all_real_x", 180*ell, 300*ell)]:
        lower, upper = interval("raw_derived_boundary_intervals." + key)
        need(lo < lower <= upper < hi, key + " safe strict printed enclosure")
    printed(r"\frac1{1000}<\mu_b<\frac1{400}", "safe boundary mass interval")
    printed(r"189<C_b/2<300-\frac1{160}", "safe static boundary interval")
    printed(r"180<\beta(t)<300", "safe all-real-time spring interval")
    printed(r"1000<\Lambda<2500", "scaled mass interval")
    need(F(129)**2/3-F(35, 4) == val("energy.all_real_W_lower") > 5500,
         "all-real-time lower potential arithmetic")
    need(131**2*6 == val("energy.all_real_W_upper") < 103000,
         "all-real-time upper potential arithmetic")
    printed(r"\frac{22153}{4}>5500", "continuous potential lower bound")
    printed(r"102966<103000", "continuous potential upper bound")
    need(5500*ell**2 == val("energy.scaled_potential_lower"),
         "scaled continuous potential lower bound")
    printed(r"5500\ell^2<\mathcal V<103000\ell^2", "scaled potential enclosure")
    need(180*ell/2500 > 5500*ell**2, "boundary as well as bulk coercivity")
    need(F(5, 2) == val("prior_operator_credit.unscaled_local_constant_shift"),
         "retained operator constant shift")
    printed(r"b_0&=\ell[C_b/2+(5/2)\mu_b]", "spring receives same static shift")
    printed(r"b(x)&=b_0+\Lambda\omega(x)", "common scalar temporal shift")
    printed(r"\frac{-u_r(0)+b_0q}{\Lambda}", "operator endpoint sign and mass")
    need(not report["prior_operator_credit"]["auxiliary_extension_is_global_regular_cone_operator_claimed"],
         "no global auxiliary spectral identification")

    contraction = F(9, 2)*F(35, 4)*ell**2*(3*ell)**2
    need(contraction == F(2835, 8)*ell**4 ==
         val("energy.static_wave_Volterra_contraction_upper") < 1,
         "three-time wave-group contraction arithmetic")
    printed(r"\frac{2835}{8}\ell^4<1", "wave-group contraction display")
    epsilon = F(35, 4)*9*ell**4/(5500*ell**2)
    need(epsilon == F(63, 4400)*ell**2 == val("energy.relative_energy_change_upper"),
         "integrated temporal relative energy change")
    need(1-epsilon == val("energy.backward_lower_factor") and
         val("energy.backward_upper_factor") == 1, "two-sided energy comparison factors")
    printed(r"\epsilon=\frac{63\ell^2}{4400}", "energy loss display")
    printed(r"E'(x)=\frac12\int_0^\infty\mathcal V_xU^2", "complete bulk energy identity")
    printed(r"+\frac{b_x}2q^2", "complete boundary spring derivative")
    printed(r"\frac{\Lambda}{2}q_x^2+\frac{b(x)}2q^2", "both boundary energy summands")
    need((1+F(103001, 4)*ell**2)/2 == val("energy.initial_energy_upper_per_D_squared"),
         "basic initial energy coefficient")
    need(F(2, 1000) == val("energy.boundary_velocity_square_per_E0") and
         4*F(2, 1000) == val("energy.boundary_master_trace_square_per_E0_to_three"),
         "boundary velocity and integrated trace constants")
    printed(r"|q_p'(s)|^2\le E(0)/500", "boundary velocity bound")
    printed(r"\le E(0)/125", "integrated master trace bound")

    need(F(1) == val("graph_energy.unit_bump_L2_derivative_caps.Q") and
         F(16) == val("graph_energy.unit_bump_L2_derivative_caps.Q_r") and
         F(416) == val("graph_energy.unit_bump_L2_derivative_caps.Q_rr"),
         "analytic unit-bump derivative caps")
    need(4**4+2*4**3+2*4**2 == 416, "analytic second derivative triangle bound")
    printed(r"256+128+32=416", "second derivative norm arithmetic")
    twice_graph = (416+103000*ell**2)**2 + ell**2*((16+ell)**2+103000*ell**2)
    need(twice_graph == val("graph_energy.twice_initial_graph_energy_upper"),
         "all-channel initial graph energy exact majorant")
    need(1-epsilon == val("graph_energy.coercivity_factor") and
         twice_graph/(1-epsilon) < 417**2,
         "strict graph norm 417 from initial energy and coercivity")
    printed(r"\norm{A_0Z(x)}_{\HH}<417", "graph norm headline")
    need(val("graph_energy.A0_Z_H_norm_strict_upper") == 417,
         "graph norm headline agrees with report")
    need(not report["graph_energy"]["A0_half_Z_compact_support_claimed"] and
         report["graph_energy"]["all_endpoint_terms_retained_by_spectral_form_identity"],
         "nonlocal graph vector retains all endpoint terms")
    printed(r"Z^{[1]}=A_0^{1/2}Z", "actual square-root graph vector")
    printed(r"\norm{A_0Z_0}_{\HH}^2+a_0[f]", "spectral form identity retains endpoint terms")
    basic_upper = (16**2+(1+103000)*ell**2)/2
    need(basic_upper < val("graph_energy.initial_basic_energy_strict_upper") == 129,
         "strict basic energy 129")
    need(F(129, 125) < 2**2 and val("graph_energy.master_boundary_value_strict_upper") == 2,
         "strict actual master trace two")

    need(val("physical_scalar_trace.rho_b_lower_guard") == 129 and
         val("physical_scalar_trace.minus_p_b_lower_guard") == F(3, 25000),
         "physical denominator guards")
    need(F(1, 400)*ell == F(1, 20000)**2 ==
         val("physical_scalar_trace.sqrt_mu_b_ell_strict_upper")**2,
         "boundary graph mass square-root factor")
    printed(r"\sqrt{\mu_b\ell}<\frac1{20000}", "physical endpoint graph factor")
    normco = report["physical_scalar_trace"]["normalization_coefficients"]
    need(normco["graph_ell_power"] == -2 and F(normco["q_factor"]) == F(5, 2)
         and F(normco["scalar_factor"]) == 3, "physical trace normalization coefficients")
    printed(r"-\frac{(A_0Z)_b}{\ell^2}+\frac52q", "exact physical endpoint cancellation")
    printed(r"\frac{3\mathsf V_z(t,0)}{\rho_b^3p_b}", "actual physical scalar definition")
    physical = F(3)/(129**3*F(3, 25000))*(F(417, 20000)/ell**2+F(5, 2)*F(1, 400)*2)
    need(physical == F(1042500000000625, 4293378) ==
         val("physical_scalar_trace.derived_uniform_upper"),
         "exact graph-to-physical upper bound")
    need(physical < val("physical_scalar_trace.safe_strict_uniform_upper") == 250000000,
         "strict actual scalar bound 250000000")
    printed(r"\frac{417}{20000\ell^2}", "physical graph term scaling")
    printed(r"\frac52\,\frac1{400}\,2", "physical master term scaling")
    printed(r"&<250000000,\qquad -3\ell\le t\le0", "physical bound and exact time slab")
    need(report["physical_scalar_trace"]["actual_selected_first_linear_coefficient_bound"]
         and not report["physical_scalar_trace"]["finite_nonlinear_amplitude_window_supplied"],
         "actual linear coefficient scope")
    printed(r"\tau=\rho_b t", "background proper clock without extra amplitude factor")
    printed(r"H(\Omega)=n_1^2-n_2^2", "physical bounded quadrupole")

    eta = ell/4
    beta = F(103000)*(ell+eta)*eta/2
    need(eta == val("qualitative_activation.max_eta") and
         beta == val("qualitative_activation.max_beta"),
         "narrow-front auxiliary activation guards")
    need(1-eta-3*beta == val("qualitative_activation.positive_contradiction_numerator_lower") > 0,
         "positive contradiction numerator for backward velocity")
    need(not report["qualitative_activation"]["auxiliary_contradiction_lower_value_is_actual_response_bound"],
         "auxiliary contradiction supplies no actual amplitude lower bound")
    printed(r"\beta_*\le5M\ell^2/32", "auxiliary activation estimate")

    polygon = re.search(r"\\fill\[blue!10\]\s*(.*?)--cycle;", source)
    need(polygon is not None, "support figure polygon present")
    points = [tuple(map(F, xy)) for xy in re.findall(r"\(([0-9.-]+),([0-9.-]+)\)", polygon.group(1))]
    need(points == [(F(1), F(0)), (F(2), F(0)), (F(5), F(-3)),
                    (F(0), F(-3)), (F(0), F(-1))], "support figure exact front vertices")
    need(max(r for r, x in points) == val("energy.support_outer_depth") == 5 < 6 < 8,
         "actual support inside unmodified coefficient collar")
    for r, x in points:
        need(0 <= -x <= 3 and max(0, 1+x) <= r <= 2-x,
             "plotted vertex obeys backward finite-propagation enclosure " + str((r, x)))
    need(interval("energy.past_scaled_x_interval") == (F(-3), F(0)),
         "closed negative-time domain")
    need(report["selected_local_inputs_authenticated"] == 5,
         "five selected local input copies")
    need(all(report["new_activity"][key] == 0 for key in
             ["background_field_bump_or_source_queries", "legacy_checker_or_project_imports",
              "new_background_ODE_image_recomputations", "numerical_PDE_solves"]),
         "reported zero scientific queries, old imports, ODE images and PDE solves")
    printed("4 November 2026", "client date")
    printed("working manuscript", "working-paper status") if "working manuscript" in source else printed("Working manuscript", "working-paper status")
    printed(r"same $\ell=10^{-6}$ packet~\cite{waveforms}", "earlier actual future-waveform credit")
    need("/workspace/" not in source, "public-safe source paths")

    result = {
        "status": "PASS_DISPLAYED_PAST_SHELL_CONTACT_BOUNDS",
        "method": "stdlib Fraction; printed constants, frozen report and support diagram; no producer imports",
        "finite_checks": len(checks), "checks": checks,
        "audit_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "inputs": {name: {"filename": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}
                   for name, path in [("source", tex), ("primary_report", args.report),
                                      ("optimized_report", args.optimized), ("primary_proof", args.proof)]},
        "scope": "Finite displayed bounds only; retained background and reconstruction, linear wave construction and physical interpretation are reviewed separately."
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "finite_checks": len(checks)}, sort_keys=True))


if __name__ == "__main__":
    main()
