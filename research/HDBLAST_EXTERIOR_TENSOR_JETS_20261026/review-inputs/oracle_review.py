#!/usr/bin/env python3
"""Independent exact-formula audit of the new actual TT Fr oracle."""
import argparse
import hashlib
import importlib.util
import json
import sys
import tempfile
from fractions import Fraction as Q
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
FIXTURE_SHA256 = "a8beb5662b8f16f9f32ece9127b52531b985abdd090afdb750bafc3fa71f09fc"
PRODUCER_SHA256 = "d36f8aca0b257ed8bc9de6b2cbacbdcde24ebf3b447823e54d8af006638e08cb"
PRODUCER_REPORT_SHA256 = "079f8a186783245e7e3a26b77b85ddc1234062afbe5e22d68bc1f412dc6cbbca"
DEPENDENCY_REPORT_SHA256 = "112f7acec47be4b5c33a520c9862b5aaa40561aa83433f7d7c6804c4bb664464"


def need(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def contains(interval, exact):
    return interval.lo <= exact <= interval.hi


def output_guard_review(module, repo):
    rows = []

    def rejected(label, target, expected):
        try:
            module.guard_output(repo, target)
        except ValueError as error:
            need(expected in str(error), "WRONG_OUTPUT_GUARD_REASON:"+label)
            rows.append({"case":label, "status":"REJECTED", "reason":str(error)})
        else:
            raise ValueError("OUTPUT_GUARD_ACCEPTED:"+label)

    with tempfile.TemporaryDirectory(prefix="source-output-guards-") as directory:
        root = Path(directory)
        regular = root / "existing.json"
        regular.write_text("existing\n")
        rejected("existing_file", regular, "FRESH_OUTPUT_REQUIRED")
        symlink = root / "dangling.json"
        symlink.symlink_to(root / "missing.json")
        rejected("dangling_symlink", symlink, "OUTPUT_SYMLINK")
        need(len(module.PROTECTED_LABELS) == 9, "NINE_OLD_PROTECTED_OUTPUT_TREES")
        for label in module.PROTECTED_LABELS:
            rejected("old_frozen_tree:"+label, repo / "research" / label / "NEVER_CREATED_REVIEW_OUTPUT.json", "OUTPUT_INSIDE_FROZEN_CAPSULE")
        original_file = module.__file__
        try:
            # Only change the imported module's location marker; no protected
            # tree or frozen byte is created or altered for this guard test.
            capsule = root / "research" / "NEW_REVIEW_CAPSULE"
            module.__file__ = str(capsule / "inputs" / "actual_third_jets.py")
            rejected("own_packaged_capsule", capsule / "NEVER_CREATED.json", "OUTPUT_INSIDE_PACKAGED_CAPSULE")
        finally:
            module.__file__ = original_file
        module.guard_output(repo, root / "fresh.json")
        rows.append({"case":"fresh_external_output", "status":"ACCEPTED_WITHOUT_WRITE"})
        need(regular.read_text() == "existing\n" and not (root / "fresh.json").exists(), "GUARD_CHECKS_HAVE_NO_OUTPUT_SIDE_EFFECTS")
    return rows


def run(repo, injection, prior_path, output):
    path = HERE.parent / "inputs/actual_third_jets.py"
    need(sha(prior_path) == PRODUCER_REPORT_SHA256, "PRODUCER_REPORT_PIN")
    prior = json.loads(prior_path.read_text())
    before = sha(path)
    need(before == prior["source_module_sha256"] == PRODUCER_SHA256, "PRODUCER_REPORT_BINDS_REVIEWED_MODULE")
    module = load("independently_reviewed_actual_third_jets", path)
    module.guard_output(repo, output)
    guards = output_guard_review(module, repo)
    ctx = module.Context(repo)
    if injection == "wrong_negative_power_derivative":
        correct_power = module.Dual.__pow__

        def wrong_power(self, power):
            result = correct_power(self, power)
            if power < 0:
                result.d = -result.d
            return result
        module.Dual.__pow__ = wrong_power
    elif injection == "omit_yy_conversion":
        correct_physical = module.physical_expression

        def wrong_physical(amp, suffix, rho, a, ell):
            if suffix == "yy":
                return (1 / (ell * rho))**2 * amp.dr().dr()
            return correct_physical(amp, suffix, rho, a, ell)
        module.physical_expression = wrong_physical

    # Entire independent derivation from the collected source was performed in
    # dependency_review.py. Bind its normal and optimized exact results here.
    normal = json.loads((HERE / "DEPENDENCY_NORMAL.json").read_text())
    optimized = json.loads((HERE / "DEPENDENCY_OPTIMIZED.json").read_text())
    need(sha(HERE / "DEPENDENCY_NORMAL.json") == DEPENDENCY_REPORT_SHA256, "INDEPENDENT_DEPENDENCY_REPORT_PIN")
    need(normal == optimized and normal["status"] == "PASS_INDEPENDENT_TT_FR_INPUT_DEPENDENCY_REVIEW", "INDEPENDENT_NORMAL_OPTIMIZED_DEPENDENCY_IDENTITY")
    need(normal["source_sha256"] == sha(ctx.maps_path), "INDEPENDENT_REVIEW_BINDS_EXACT_SOURCE_MAP")
    need(ctx.dependency_audit["required_actual_master_total_order"] == 2 and ctx.dependency_audit["required_actual_gauge_total_order"] == 3, "PRODUCER_ORDER_AUDIT_MATCHES_INDEPENDENT")

    # Pinned expectations were generated from exact coordinate polynomials by
    # SymPy1.14.0. This final replay only imports Python's standard library.
    fixture_path = HERE / "EXACT_FIXTURE.json"
    need(sha(fixture_path) == FIXTURE_SHA256, "INDEPENDENT_EXACT_FIXTURE_PIN")
    fixture = json.loads(fixture_path.read_text())
    need(fixture["generator_sha256"] == sha(HERE / "generate_exact_fixture.py"), "EXACT_FIXTURE_GENERATOR_BINDING")
    need(fixture["source_map_sha256"] == sha(ctx.maps_path), "EXACT_FIXTURE_SOURCE_MAP_BINDING")
    need(fixture["sympy_version"] == "1.14.0", "EXACT_FIXTURE_SYMPY_VERSION")

    def taylor(name):
        return ctx.jets.T({tuple(map(int, key.split(','))): Q(value) for key,value in fixture["synthetic_normalized_Taylor_coefficients"][name].items()})

    amp_t, rho_t, a_t = map(taylor, ("amp", "rho", "a"))
    ell = Q(fixture["synthetic_ell"])
    synthetic_checks = []
    for test in fixture["synthetic_physical_expected"]:
        suffix, order, exact = test["suffix"], test["source_r_order"], Q(test["exact"])
        observed = module.physical_expression(amp_t, suffix, rho_t, a_t, ell)
        enclosure = observed.val() if order == 0 else observed.dr().val()
        need(contains(enclosure, exact), "INDEPENDENT_PHYSICAL_CHAIN:" + suffix + ":" + str(order))
        need(enclosure.hi - enclosure.lo < Q(1, 10**20), "SYNTHETIC_CHAIN_TEST_SHARPNESS")
        synthetic_checks.append({"suffix":suffix, "source_r_order":order, "exact":str(exact), "interval":enclosure.row()})

    # Compare production Dual AD against an independent symbolic gradient of
    # every one of the 78 rational monomials. Nonzero derivatives of c and rho
    # in these algebra tests exercise both negative and positive powers.
    terms = ctx.terms[5]
    names = sorted({name for term in terms for name in term["powers"]})
    algebra_rows = []
    for test in fixture["independent_full_source_symbolic_gradient_cases"]:
        case = test["case"]
        values = {name:Q(value) for name,value in test["values"].items()}
        rates = {name:Q(value) for name,value in test["rates"].items()}
        need(sorted(values) == names == sorted(rates), "ALL_SOURCE_SYMBOLS_EXERCISED")
        exact_value, exact_derivative = Q(0), Q(0)
        for term in terms:
            value = Q(term["coefficient"])
            for name,power in term["powers"].items():
                value *= values[name]**power
            exact_value += value
            exact_derivative += value * sum((power*rates[name]/values[name] for name,power in term["powers"].items()), Q(0))
        need(exact_value == Q(test["F"]) and exact_derivative == Q(test["Fr"]), "INDEPENDENT_LEIBNIZ_EQUALS_PINNED_SYMBOLIC_GRADIENT:"+str(case))
        observed = ctx.coeff.evaluate_encoded(terms, {name:module.Dual(values[name], rates[name]) for name in names})
        need(observed.v == exact_value, "INDEPENDENT_FULL_SOURCE_VALUE:"+str(case))
        need(observed.d == exact_derivative, "INDEPENDENT_FULL_SOURCE_LEIBNIZ:"+str(case))
        algebra_rows.append({"case":case, "F":str(exact_value), "Fr":str(exact_derivative)})

    # Actual probes replay the authenticated C3 enclosures and the separately
    # retained frozen value oracle. These are enclosure checks, not claims that
    # the unknown actual solution takes the interval midpoints.
    actual_rows = []
    for old in prior["probes"]:
        xx, rr = ctx.I(*old["x"]), ctx.I(*old["r"])
        value, radial = ctx.forcing_with_radial(xx, rr, audit=True)
        need(value.row() == old["F"] and radial.row() == old["Fr"], "ACTUAL_REPORT_REPLAY:"+old["label"])
        actual_rows.append({"label":old["label"], "F":value.row(), "Fr":radial.row()})
    ctx.finish()
    need(sha(path) == before, "REVIEWED_MODULE_CHANGED")
    return {
        "status":"PASS_INDEPENDENT_ACTUAL_TT_FR_ORACLE_REVIEW",
        "reviewer_source_sha256":sha(Path(__file__)),
        "producer_module_sha256":before,
        "producer_report_sha256":sha(prior_path),
        "dependency_report_sha256":sha(HERE / "DEPENDENCY_NORMAL.json"),
        "source_map_sha256":sha(ctx.maps_path),
        "exact_fixture_sha256":FIXTURE_SHA256,
        "fixture_generator_sympy_version":fixture["sympy_version"],
        "final_replay_dependencies":"Python standard library only",
        "independent_source_output_guard_checks":guards,
        "physical_conversion_exact_checks":synthetic_checks,
        "independent_full_source_symbolic_leibniz_cases":algebra_rows,
        "actual_probe_replay":actual_rows,
        "authenticated_C3_errors_retained":True,
        "C4_needed":False,
        "normal_and_optimized_dependency_reports_identical":True,
        "all_frozen_closures_reauthenticated_after_run":True,
        "scope":"Exact dependency, physical derivative and full-source AD formulas verified independently; unknown actual solution bounds retain the inherited authenticated C3 conditional premises.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--producer-report", type=Path, default=HERE.parent / "inputs/actual_source_radial_final_normal.json")
    parser.add_argument("--inject-error", choices=("none", "wrong_negative_power_derivative", "omit_yy_conversion"), default="none")
    args = parser.parse_args()
    need(not args.output.exists(), "FRESH_OUTPUT_REQUIRED")
    result = run(args.repo_root.resolve(), args.inject_error, args.producer_report.resolve(), args.output)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"], "physical_checks":len(result["physical_conversion_exact_checks"]), "full_source_cases":3, "actual_probes":len(result["actual_probe_replay"])}))


if __name__ == "__main__":
    main()
