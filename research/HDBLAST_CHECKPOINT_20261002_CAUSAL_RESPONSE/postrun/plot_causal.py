#!/usr/bin/env python3
"""Render registered saved values only; no scientific source functions imported."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/hdblast-causal-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/hdblast-causal-cache")
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter, NullLocator


COLORS = {"positive_B": "#176b91", "signed_uB": "#b44a2a"}
LABELS = {"positive_B": "Positive pulse B", "signed_uB": "Signed pulse uB"}
MARKERS = {"positive_B": "o", "signed_uB": "s"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(figure, output, stem):
    for extension in ("svg", "png", "pdf"):
        figure.savefig(output/(stem+"."+extension), dpi=200, bbox_inches="tight",
                       metadata={"Creator": "Post-run registered-data renderer"})
    plt.close(figure)


def timeline(axis, *, shade=True):
    if shade:
        axis.axvspan(-5, -3, color="#e9edf0", zorder=0)
    axis.axhline(0, color="#727b83", linewidth=.8, zorder=1)
    axis.set_xlim(-5.7, -1.3)
    axis.set_xticks([-5.5, -4.5, -4, -3.5, -2.5, -1.5])
    axis.set_xlabel(r"Conformal time $\eta$ ($H=1$)")
    axis.grid(axis="y", alpha=.2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError("Refusing to overwrite figures directory")
    data = json.loads(args.summary.read_text())
    if data["metrics"]["registered_source_observation_points"] != 12 or data["new_source_response_or_mode_evaluations"] != 0:
        raise RuntimeError("Unexpected presentation input scope")
    groups = {source: [r for r in data["registered_rows"] if r["source"] == source] for source in COLORS}
    for source, rows in groups.items():
        if [r["eta"] for r in rows] != [-5.5, -4.5, -4, -3.5, -2.5, -1.5]:
            raise RuntimeError("Registered observation grid changed")
    args.output.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11,
                         "figure.titlesize": 15, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.hashsalt": "causal-registration-a01d17e"})

    figure, axes = plt.subplots(1, 2, figsize=(11.4, 4.5))
    for source, rows in groups.items():
        times = [r["eta"] for r in rows]
        axes[0].scatter(times, [r["s_over_epsilon"] for r in rows], color=COLORS[source],
                        marker=MARKERS[source], s=53, label=LABELS[source], zorder=3)
        axes[1].plot(times, [r["continuum_y"] for r in rows], color=COLORS[source],
                     marker=MARKERS[source], markersize=6, linewidth=1.3, linestyle="--", label=LABELS[source], zorder=3)
    for axis in axes:
        timeline(axis)
    axes[0].set_title("Prescribed sources: recorded values")
    axes[0].set_ylabel(r"$f=s/\epsilon$")
    axes[0].text(-4, 1.10, "Exact source support", ha="center", color="#5b6770", fontsize=9)
    axes[0].set_ylim(-.5, 1.2)
    axes[0].legend(frameon=False, loc="lower right", fontsize=9)
    axes[1].set_title("Retarded variance: removed-cutoff result")
    axes[1].set_ylabel(r"$y=a^2\,\delta Q/\epsilon$")
    axes[1].text(-2, .0105, "Source is zero;\nresponse remains negative", ha="center", va="center", fontsize=9,
                 bbox={"boxstyle": "round,pad=.35", "facecolor": "white", "edgecolor": "#d9dfe2"})
    figure.suptitle("Causal memory at the twelve registered points", y=1.01)
    figure.text(.5, -.015, "Shading: −5 < η < −3. Markers: registered observations only. Dashed response lines guide the eye.\n"
                "Fixed geometry and incoming BD state; ε = 10⁻⁴. No stress, stability, or heating conclusion.", ha="center", fontsize=9, color="#4f5962")
    figure.tight_layout(w_pad=2.7)
    save(figure, args.output, "causal_memory")

    figure, axes = plt.subplots(1, 2, figsize=(11.4, 4.5), sharey=True)
    for axis, (source, rows) in zip(axes, groups.items()):
        selected = [r for r in rows if r["eta"] > -5]
        times = [r["eta"] for r in selected]
        points = [r["finite_k"][-1] for r in selected]
        axis.semilogy(times, [p["analytic_tail_bound"] for p in points], "D--", color="#505c64", markersize=5,
                       linewidth=1, label="Analytic combined UV tail bound")
        axis.semilogy(times, [p["primary_to_continuum_difference"] for p in points], "o", color=COLORS[source],
                       markersize=7, label=r"$|y_{K}-y_\infty|$: exact finite-K expression")
        axis.semilogy(times, [p["mode_to_continuum_difference"] for p in points], "^", markerfacecolor="none",
                       markeredgewidth=1.2, color="#1b2024", markersize=9, label=r"$|y_{K}-y_\infty|$: independent modes")
        axis.set_title(LABELS[source]+", K=256")
        axis.set_xlabel(r"Conformal time $\eta$")
        axis.set_xticks(times)
        axis.set_ylim(5e-15, 8e-6)
        axis.grid(axis="y", which="major", alpha=.25)
    axes[0].set_ylabel("Absolute difference or UV bound in y")
    axes[0].legend(frameon=False, loc="center left", fontsize=8.5)
    figure.suptitle("Observed regulator differences and analytic remainder bounds", y=1.01)
    figure.text(.5, -.015, "Exact pre-pulse zeros omitted on log scale. Mode and finite-K markers nearly coincide.\nOnly the UV remainder is analytically enclosed; numerical integration estimates remain separate.", ha="center", fontsize=9, color="#4f5962")
    figure.tight_layout(w_pad=2.3, rect=(0, .04, 1, 1))
    save(figure, args.output, "causal_agreement")

    figure, axes = plt.subplots(1, 2, figsize=(11.4, 4.5))
    cutoffs = [64, 128, 256]
    for source, rows in groups.items():
        maximum_delta = [max(r["finite_k"][i]["primary_to_continuum_difference"] for r in rows) for i in range(3)]
        maximum_tail = [max(r["finite_k"][i]["analytic_tail_bound"] for r in rows) for i in range(3)]
        maximum_cross = [max(r["finite_k"][i]["mode_to_primary_difference"] for r in rows) for i in range(3)]
        maximum_refinement = [max(r["finite_k"][i]["mode_refinement_difference"] for r in rows) for i in range(3)]
        color = COLORS[source]
        axes[0].loglog(cutoffs, maximum_delta, marker=MARKERS[source], color=color, label=LABELS[source]+": observed")
        axes[0].loglog(cutoffs, maximum_tail, marker=MARKERS[source], markerfacecolor="none", linestyle="--", color=color,
                       label=LABELS[source]+": analytic tail")
        axes[1].semilogy(cutoffs, maximum_cross, marker=MARKERS[source], color=color, label=LABELS[source]+": cross-route")
        axes[1].semilogy(cutoffs, maximum_refinement, marker=MARKERS[source], markerfacecolor="none", linestyle="--", color=color,
                         label=LABELS[source]+": refinement")
    for axis in axes:
        axis.set_xticks(cutoffs)
        axis.xaxis.set_major_formatter(ScalarFormatter())
        axis.xaxis.set_minor_locator(NullLocator())
        axis.xaxis.set_minor_formatter(NullFormatter())
        axis.set_xlabel("Registered comoving cutoff K")
        axis.grid(axis="y", alpha=.25)
        axis.legend(frameon=False, fontsize=8.5)
    axes[0].set_title("Regulator discrepancy and bound")
    axes[0].set_ylabel("Maximum over six registered times, in y")
    axes[1].set_title("Finite-K agreement and refinement")
    axes[1].set_ylabel("Maximum observed absolute difference, in y")
    figure.suptitle("Three frozen cutoffs, with no added evaluations", y=1.01)
    figure.text(.5, -.015, "Lines guide the eye between registered K values. Tiny observed cross-route/refinement differences are not certified accuracy bounds.",
                ha="center", fontsize=9, color="#4f5962")
    figure.tight_layout(w_pad=2.8)
    save(figure, args.output, "causal_cutoff")
    figures = {path.name: sha(path) for path in sorted(args.output.iterdir()) if path.suffix in (".svg", ".png", ".pdf")}
    provenance = {"purpose": "post-run rendering of registered saved observations only", "new_evaluations": 0,
                  "summary_sha256": sha(args.summary), "renderer_sha256": sha(__file__),
                  "matplotlib_version": matplotlib.__version__, "figures_sha256": figures}
    (args.output/"FIGURE_PROVENANCE.json").write_text(json.dumps(provenance, indent=2, sort_keys=True)+"\n")
    print("Rendered", len(figures), "figures from", len(data["registered_rows"]), "registered saved rows")


if __name__ == "__main__":
    main()
