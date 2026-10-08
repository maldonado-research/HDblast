#!/usr/bin/env python3
"""Present serialized certified outputs; never import numerical production code.

All display rounding is checked against exact Fractions. Floating-point
conversion is used only for drawing the display figure. No physical target or
source callback, original-array read, or likelihood computation occurs here.
"""
from __future__ import annotations

import csv
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import tempfile

_plot_cache = tempfile.mkdtemp(prefix="hdblast-presentation-cache-")
os.environ["MPLCONFIGDIR"] = _plot_cache
os.environ["XDG_CACHE_HOME"] = _plot_cache
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


BASE = Path(__file__).resolve().parent
WORK = BASE.parent
INPUTS = {
    "exact_data": WORK / "root/actual-normal005-001/DATA.json",
    "independent_display": WORK / "actual-independent-review/METRICS_DISPLAY_24.json",
    "science_summary": WORK / "root/actual-normal005-001/SCIENCE_SUMMARY.json",
    "target_budgets": WORK / "root/actual-normal005-001/TARGET_BUDGETS.json",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def pin(path):
    raw = path.read_bytes()
    return {"path": str(path), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


ROUND_CHECKS = []


def display(value, digits=10, upward=True):
    """Exact decimal-grid ceiling/floor, then Decimal scientific formatting."""
    x = Fraction(value)
    if x == 0:
        return "0E+0"
    ax = abs(x)
    exponent = len(str(ax.numerator)) - len(str(ax.denominator))
    if ax < Fraction(10) ** exponent:
        exponent -= 1
    require(Fraction(10) ** exponent <= ax < Fraction(10) ** (exponent + 1), "exact exponent")
    quantum_exp = exponent - digits + 1
    quantum = Fraction(10) ** quantum_exp
    scaled = x / quantum
    units = -(-scaled.numerator // scaled.denominator) if upward else scaled.numerator // scaled.denominator
    with localcontext() as ctx:
        ctx.prec = digits + 20
        value_dec = Decimal(units) * Decimal(10) ** quantum_exp
        result = format(value_dec, f".{digits - 1}E")
    rounded = Fraction(Decimal(result))
    require(rounded >= x if upward else rounded <= x, "outward decimal direction")
    require(abs(rounded - x) <= quantum, "outward decimal step")
    ROUND_CHECKS.append({"exact": str(x), "display": result, "direction": "CEILING" if upward else "FLOOR"})
    return result


def finite_interval(row, coord, digits=5):
    lo, hi = row["signed_intervals"][coord + "_error"]
    return f"[{display(lo, digits, False)}, {display(hi, digits, True)}]"


def main():
    initial_pins = {name: pin(path) for name, path in INPUTS.items()}
    docs = {name: json.loads(path.read_text()) for name, path in INPUTS.items()}
    data, review, summary = docs["exact_data"], docs["independent_display"], docs["science_summary"]
    rows = data["rows"]
    require(len(rows) == 24 and len({r["case"] for r in rows}) == 24, "24 unique cases")
    require(data["scope"] == "SAVED_ENDPOINT_MINUS_EXACT_FLOW_FROM_EXACT_REPRESENTED_INCOMING_STATE", "scope")
    require(data["nodes"] == 49152 and data["endpoint_comparisons"] == 98304, "coverage")
    require(data["manufactured"] is False, "actual serialized dataset")
    require(summary["status"] == "FINITE_RETAINED_ENDPOINT_ERRORS_ENCLOSED", "science status")
    require(summary["incoming_error_counted"] is False, "incoming component separate")
    by_case = {r["case"]: r for r in rows}
    require(set(by_case) == {r["case"] for r in review["rows"]}, "independent case identities")

    # Independent display values are approximate; exact DATA remains authority.
    crosschecks = []
    for r in review["rows"]:
        exact_row = by_case[r["case"]]
        comparisons = {
            "maximum_normalized_U_error": Fraction(exact_row["maxima"]["U_L1_upper"]),
            "maximum_normalized_W_error": Fraction(exact_row["maxima"]["W_L1_upper"]),
            "R_triangle_upper": Fraction(exact_row["sums"]["R_triangle_upper"]),
            "P_triangle_upper": Fraction(exact_row["sums"]["P_triangle_upper"]),
            "signed_R_error_midpoint": sum(map(Fraction, exact_row["signed_intervals"]["R_error"])) / 2,
            "signed_P_error_midpoint": sum(map(Fraction, exact_row["signed_intervals"]["P_error"])) / 2,
        }
        for name, exact in comparisons.items():
            approx = Decimal(r[name])
            quantum = Fraction(10) ** approx.as_tuple().exponent
            require(abs(Fraction(approx) - exact) <= quantum, "independent display exceeds final shown place")
            crosschecks.append({"case": r["case"], "field": name, "within_one_display_quantum": True})

    maxima = {coord: max(Fraction(r["maxima"][coord + "_L1_upper"]) for r in rows) for coord in ("U", "W")}
    radii = {coord: Fraction(summary["maximum_complete_endpoint_export_L1_radius"][coord]) for coord in ("U", "W")}
    sign_counts = {}
    for coord in ("R", "P"):
        intervals = [tuple(map(Fraction, r["signed_intervals"][coord + "_error"])) for r in rows]
        require(all(lo <= hi for lo, hi in intervals), "ordered signed intervals")
        sign_counts[coord] = {"positive": sum(lo > 0 for lo, hi in intervals), "negative": sum(hi < 0 for lo, hi in intervals), "contains_zero": sum(lo <= 0 <= hi for lo, hi in intervals)}
        require(sign_counts[coord]["contains_zero"] == 0, "all signed finite sums exclude zero")
    for r in rows:
        require(r["incoming_state_error"] == "ZERO_BY_EXACT_REPRESENTED_ANCHOR_DEFINITION", "anchor")
        require(r["contact_difference"] == "EXACT_ZERO_UNDER_IDENTICAL_COMPLETE_CONTACTS", "conditional contacts")

    # Exact rationals plus directed displays, covering all original row metrics.
    csv_rows = []
    for r in rows:
        source, grid = r["capsule"].split("/")
        out = {"case": r["case"], "source": source, "grid": grid, "K": r["K"], "eta": r["endpoint"], "node_count": r["node_count"], "last_k_exact": str(Fraction(r["last_k"]))}
        for group in ("maxima", "sums"):
            for name, val in r[group].items():
                out[name + "_exact"] = str(Fraction(val))
                out[name + "_decimal_ceiling"] = display(val)
        for coord in ("R", "P"):
            lo, hi = r["signed_intervals"][coord + "_error"]
            out[coord + "_signed_lower_exact"] = str(Fraction(lo))
            out[coord + "_signed_upper_exact"] = str(Fraction(hi))
            out[coord + "_signed_lower_decimal_floor"] = display(lo, upward=False)
            out[coord + "_signed_upper_decimal_ceiling"] = display(hi)
        out["scope"] = data["scope"]
        out["contact_difference"] = r["contact_difference"]
        out["incoming_state_error"] = r["incoming_state_error"]
        csv_rows.append(out)
    csv_path = BASE / "EXACT_ENDPOINT_METRICS_24.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    with csv_path.open(newline="") as f:
        serialized_rows = list(csv.DictReader(f))
    require(len(serialized_rows) == 24, "24 CSV data rows")
    require(all(a["case"] == b["case"] for a, b in zip(csv_rows, serialized_rows)), "CSV case order")

    # Figure: double precision is used only here, after exact decisions above.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 1, figsize=(13, 8.3), sharex=True, gridspec_kw={"hspace": 0.27})
    capsules = ["positive_B/coarse", "positive_B/fine", "signed_uB/coarse", "signed_uB/fine"]
    u_color, w_color = "#b65f00", "#1768a4"
    neg_color, pos_color = "#7546a4", "#12805c"
    for r in rows:
        group = capsules.index(r["capsule"]) * 3 + ["64", "128", "256"].index(r["K"])
        late = r["endpoint"] == "-7/2"
        x = group + (0.14 if late else -0.14)
        for coord, color, marker in [("U", u_color, "^"), ("W", w_color, "o")]:
            axes[0].plot(x, float(Fraction(r["maxima"][coord + "_L1_upper"])), marker=marker, color=color, markerfacecolor=color if late else "white", markersize=6, linestyle="none", markeredgewidth=1.3)
        lo, hi = map(Fraction, r["signed_intervals"]["P_error"])
        magnitude = (abs(lo) + abs(hi)) / 2
        lower_magnitude, upper_magnitude = sorted((abs(lo), abs(hi)))
        color = pos_color if lo > 0 else neg_color
        marker = "^" if lo > 0 else "v"
        axes[1].errorbar(x, float(magnitude), yerr=[[float(magnitude - lower_magnitude)], [float(upper_magnitude - magnitude)]], marker=marker, color=color, markersize=6, markerfacecolor=color if late else "white", linestyle="none", capsize=2, zorder=3)
        axes[1].plot(x, float(Fraction(r["sums"]["P_triangle_upper"])), marker="x", color="#343a40", markersize=7, linestyle="none", markeredgewidth=1.2, zorder=2)
    for ax in axes:
        ax.set_yscale("log")
        ax.grid(axis="y", which="major", color="#dddddd", linewidth=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        for boundary in [2.5, 5.5, 8.5]:
            ax.axvline(boundary, color="#c9c9c9", linewidth=0.8)
        for start in [0, 6]:
            ax.axvspan(start - 0.5, start + 2.5, color="#f4f6f8", zorder=-1)
    axes[0].set_title("A   Maximum node error: certified upper bounds", loc="left", fontweight="bold")
    axes[0].set_ylabel("Normalized Cartesian L1 upper")
    axes[0].legend(handles=[Line2D([], [], color=u_color, marker="^", linestyle="none", label="U = u₁/ε"), Line2D([], [], color=w_color, marker="o", linestyle="none", label="W = w₁/ε"), Line2D([], [], color="#555555", marker="o", markerfacecolor="white", linestyle="none", label="η = −4 (left)"), Line2D([], [], color="#555555", marker="o", linestyle="none", label="η = −7/2 (right)")], ncol=4, loc="upper left", fontsize=9, frameon=True)
    axes[0].set_ylim(4e-17, 2e-14)
    axes[1].set_title("B   Finite pressure sum: signed error magnitude and triangle bound", loc="left", fontweight="bold")
    axes[1].set_ylabel("Normalized finite pressure difference")
    axes[1].legend(handles=[Line2D([], [], color="#343a40", marker="x", linestyle="none", label="P triangle upper"), Line2D([], [], color=neg_color, marker="v", linestyle="none", label="|signed P error|, negative"), Line2D([], [], color=pos_color, marker="^", linestyle="none", label="|signed P error|, positive")], ncol=3, loc="upper left", fontsize=9, frameon=True)
    axes[1].set_ylim(6e-16, 4e-10)
    axes[1].set_xticks(range(12), ["64", "128", "256"] * 4)
    axes[1].set_xlim(-0.55, 11.55)
    for center, label in zip([1, 4, 7, 10], ["positive B / coarse", "positive B / fine", "signed uB / coarse", "signed uB / fine"]):
        axes[1].text(center, -0.145, label, ha="center", va="top", transform=axes[1].get_xaxis_transform(), fontsize=10)
    fig.suptitle("Later saved endpoints versus exact flow from the represented incoming state", fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.08, 0.04, "Ticks: retained momentum cutoff K. Open/left: η = −4; filled/right: η = −7/2. Pressure intervals are narrower than their symbols.\nDisplay only: exact rational CSV governs bounds and signs. Finite nodes/endpoints; no dense-path or continuum certificate.", fontsize=9, color="#444444")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.925, bottom=0.17)
    fig.savefig(BASE / "endpoint_comparison.png", dpi=260, facecolor="white")
    fig.savefig(BASE / "endpoint_comparison.pdf", facecolor="white", metadata={"Title": "Certified HDBLAST endpoint discrepancies", "Author": "HDBLAST research presentation", "Subject": "Serialized exact endpoint results; display only", "CreationDate": None, "ModDate": None})
    plt.close(fig)

    # Concise report: all numbers generated with exact outward rounding.
    worst = {coord: display(maxima[coord], 4) for coord in ("U", "W")}
    target_display = {coord: display(radii[coord], 4) for coord in ("U", "W")}
    for coord in ("U", "W"):
        require(maxima[coord] < Fraction(Decimal(worst[coord])), "strict reported endpoint upper")
        require(radii[coord] < Fraction(Decimal(target_display[coord])), "strict reported radius upper")
        require(radii[coord] < Fraction(1, 10**18), "reference radius threshold")
    table_lines = ["| Source / grid | η | U upper | W upper | Signed finite P interval | P triangle upper |", "| --- | ---: | ---: | ---: | --- | ---: |"]
    for r in rows:
        if r["K"] == "256":
            table_lines.append(f"| {r['capsule']} | {r['endpoint']} | {display(r['maxima']['U_L1_upper'], 4)} | {display(r['maxima']['W_L1_upper'], 4)} | {finite_interval(r, 'P', 5)} | {display(r['sums']['P_triangle_upper'], 4)} |")
    r_triangle = max(Fraction(r["sums"]["R_triangle_upper"]) for r in rows)
    p_triangle = max(Fraction(r["sums"]["P_triangle_upper"]) for r in rows)
    report = f"""# Certified later-endpoint errors and the physical-model limit

8 October 2026. The later retained states are now enclosed against the exact
scalar evolution from their **exact represented incoming state** at η=−9/2.
Coverage is 49,152 retained node occurrences and 98,304 endpoint comparisons at
η=−4 and −7/2. The two sources, two grids, three cutoffs K=64,128,256 and two
endpoints give 24 reported rows. Their largest normalized Cartesian L1 error
uppers satisfy **max ‖ΔU‖₁ < {worst['U']} and max ‖ΔW‖₁ < {worst['W']}**.

The target is the same prescribed scalar dynamics, U′=W and W′=2ikW−g, started
from the retained incoming coordinates. The discrepancy is **saved endpoint
minus exact flow**. It encloses the total subsequent integration/arithmetic
error at these snapshots, without separating unsaved historical step errors.
The earlier incoming-state discrepancy from the prescribed Bunch–Davies target
is a separate contribution: this comparison neither erases it nor counts it
again. No normalization projection or recalibration is applied.

## Precision and resolved finite differences

The maximum complete reference-enclosure L1 radii are exactly
`U: {radii['U']}` and `W: {radii['W']}`: respectively
less than {target_display['U']} and {target_display['W']}. The 1E−18 reference-radius
threshold checks target precision. It is not a tolerance for the saved endpoint
error, a physical accuracy claim, or the full pressure/contact gate.

All 24 signed finite density R intervals and all 24 signed finite pressure P
intervals exclude zero exactly. R has 13 positive and 11 negative cases; P has
six positive and 18 negative cases. These signed sums retain cancellation.
Their triangle uppers instead sum absolute contributions and need not be
attained. The largest triangle uppers are R < {display(r_triangle, 4)} and
P < {display(p_triangle, 4)}. Complete contact differences cancel only under
the declared identical-complete-contacts premise.

The table shows K=256; all **24 rows**, exact rational endpoints and outward
decimal displays are in [EXACT_ENDPOINT_METRICS_24.csv](EXACT_ENDPOINT_METRICS_24.csv).
U=u₁/ε and W=w₁/ε use the inherited normalized state convention. R/P use the
inherited a₀⁴/ε-scaled first-order finite stress convention, H=1, with the
retained positive momentum weights. No continuum quadrature accuracy is implied.

{chr(10).join(table_lines)}

The larger coarse-grid W and pressure bounds at K=256 locate a feature of these
finite retained outputs. They do not determine its cause or establish a
convergence order. Ratios of triangle bounds would not be ratios of realized
pressure errors.

![Certified endpoint comparison](endpoint_comparison.png)

The standalone [PDF figure](endpoint_comparison.pdf) and PNG show all 24 cases.
Panel A plots state-error uppers. Panel B distinguishes pressure triangle
uppers from signed-error magnitudes; triangle direction preserves the sign.
Within each cutoff the earlier endpoint is on the left and the later endpoint
on the right. Signed interval widths are narrower than the plotted symbols.
The drawing uses floating-point display coordinates; exact CSV fractions
govern every sign and bound.

## What the conditional bulk model adds

The independently accepted scalar half-space model derives a causal boundary
kernel, `G_R=1/[c+sqrt(k²+M²−(ω+i0)²)]`, from an explicit stable flat 4+1-dimensional
action with one spacelike extra direction. Its positive spectrum constrains
response and equilibrium noise jointly. A four-dimensional continuum of fields
reproduces the same reduced linear response and Gaussian measurement laws under
matched preparation and measurement protocol; this is no dimensional-origin test.

The accepted incoming-scattering extension adds an explicit energy supply as a
declared coherent incoming packet. It proves unit reflection, transient brane
excitation and eventual local decay for smooth continuum packets away from
threshold with zero initial bound projection in both position and velocity.
Merely setting q=q_dot=0 at a finite start is insufficient. Fixed-k statements
concern Fourier fibers; finite total brane energy requires smooth packets in k
as well. Their energy returns to the outgoing bulk. Bound components already
present persist. This free autonomous model derives no irreversible capture,
reheating or hot Big Bang.

The next physical priority is a specified gravitational/cosmological action and
an energy-transfer mechanism capable of depositing the incoming energy into
declared brane degrees of freedom, with source preparation, backreaction and
thermalization analyzed from that same model. Adding an arbitrary pulse does
not provide those ingredients. Separately, the numerical certificate still
needs control between saved times, continuous momentum integration, complete
contacts and the ultraviolet remainder. The unsaved original solver path is
not reconstructed by two snapshots.

Historical metric calibration remains **FAIL**; the full continuous
pressure/contact/time/momentum/UV certificate remains **UNRESOLVED**;
higher-dimensional Big Bang cause remains **NOT_ESTABLISHED**; external novelty
remains **NOT_ASSESSED**.

This is a presentation of independently certified serialized outputs, not a new
study or physical execution. `PRESENTATION_RECEIPT.json` records exact rounding
and independent-display checks; `MANIFEST.json` pins the four read-only inputs
and all presentation files. No original arrays or production numerical modules
were opened or imported by this presentation.
"""
    (BASE / "RESULTS.md").write_text(report)

    final_pins = {name: pin(path) for name, path in INPUTS.items()}
    require(initial_pins == final_pins, "serialized inputs unchanged")
    receipt = {
        "schema_version": 1,
        "status": "PASS_PRESENTATION_OF_CERTIFIED_SERIALIZED_OUTPUTS",
        "input_pins": initial_pins,
        "csv_rows": len(serialized_rows),
        "source_case_count": len(rows),
        "endpoint_comparisons": data["endpoint_comparisons"],
        "retained_node_occurrences": data["nodes"],
        "sign_counts": sign_counts,
        "exact_global_error_uppers": {c: str(v) for c, v in maxima.items()},
        "exact_maximum_target_L1_radii": {c: str(v) for c, v in radii.items()},
        "independent_display_crosschecks": crosschecks,
        "outward_rounding_checks": ROUND_CHECKS,
        "inputs_unchanged_after_render": True,
        "physical_source_callbacks": 0,
        "physical_target_evaluations": 0,
        "retained_array_reads": 0,
        "retained_array_decodes": 0,
        "observational_likelihood_evaluations": 0,
        "production_numerical_module_imports": 0,
        "plotting": {"library": "matplotlib", "version": matplotlib.__version__, "coordinates": "binary64 display only", "exact_authority": "DATA.json rationals, retained in CSV"},
    }
    (BASE / "PRESENTATION_RECEIPT.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    names = ["RESULTS.md", "EXACT_ENDPOINT_METRICS_24.csv", "endpoint_comparison.png", "endpoint_comparison.pdf", "render_results.py", "PRESENTATION_RECEIPT.json"]
    manifest = {"schema_version": 1, "status": "FROZEN_SCIENTIFIC_PRESENTATION", "input_pins": initial_pins, "files": [pin(BASE / n) for n in names], "scope": "Final presentation only; no new physical execution, retained-array decoding, novelty or origin claim."}
    (BASE / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": receipt["status"], "csv_rows": len(serialized_rows), "rounding_checks": len(ROUND_CHECKS), "independent_display_checks": len(crosschecks), "manifest_sha256": pin(BASE / "MANIFEST.json")["sha256"]}))


if __name__ == "__main__":
    main()
