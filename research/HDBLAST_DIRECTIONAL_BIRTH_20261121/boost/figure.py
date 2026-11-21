"""Plot the expanding branch of three matched wall evolutions through turnaround."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
CASES = (
    ("pole_cap", 1, "Flavor-specific vector resonance", "#00666a", "-"),
    ("prescribed_drag", 2, "Prescribed constant drag", "#3266a8", "--"),
    ("zero_drag", 0, "Zero drag", "#9d492b", "-"),
)


def replay(base, feedback, table, archive, jobs):
    previous = json.loads((base / "envelope/ENVELOPE.json").read_text())
    archived = json.loads(archive.read_text())["runs"]
    joined = runpy.run_path(str(base / "coupled/joined.py"))
    make_source = runpy.run_path(str(HERE / "evolve.py"))["make_source"]
    transform = runpy.run_path(str(base / "envelope/build.py"))["transform"]
    src = transform((base / "coupled/joined_forward.cpp").read_text(), centered=True)
    src = make_source(src, runpy.run_path(str(feedback / "feedback.py")),
                      (feedback / "angular.hpp").read_text())
    action = 5 * previous["screen"]["mass_4e12"]["elastic_drag_action_floor"]
    parameters = dict(previous["parameters"]["mass_4e12"], sample=1)
    with tempfile.TemporaryDirectory(prefix="hdblast-heating-figure-") as temp:
        dest = Path(temp)
        source, exe = dest / "figure.cpp", dest / "solver"
        source.write_text(src)
        subprocess.run(["g++", "-std=c++17", "-O3", "-I", str(HERE), str(source),
                        "-o", str(exe)], check=True)

        def simulate(case):
            name, mode, *_ = case
            wrapper = dest / name
            options = shlex.join(["2", repr(action), str(table.resolve()), str(mode)])
            wrapper.write_text('#!/bin/sh\nexec ' + shlex.quote(str(exe)) + ' "$@" ' + options + '\n')
            wrapper.chmod(0o755)
            result = joined["run"](wrapper, parameters)
            ref = archived[name]
            for key in ("maximum_radius", "maximum_tau", "peak_hot_fraction", "hot_efolds"):
                if not math.isclose(result.get(key, 0.), ref.get(key, 0.), rel_tol=2e-9):
                    raise ValueError(f"Trace disagrees with archived {name} {key}")
            columns = dict(zip(joined["COLUMNS"], np.asarray(result["trace"]).T))
            time = columns["t"]
            ratio = np.divide(columns["intrinsic"], columns["e"],
                              out=np.zeros_like(time), where=columns["e"] > 0)
            opening = result["bounce_tau"]
            end = ref["maximum_tau"]
            mask = (time >= opening) & (time < end)
            radius = np.r_[result["bounce_radius"], columns["a"][mask], ref["maximum_radius"]]
            fraction = np.r_[np.interp(opening, time, ratio), ratio[mask], np.interp(end, time, ratio)]
            return name, (radius, 100 * fraction, ref)
        with ThreadPoolExecutor(jobs) as pool:
            return dict(pool.map(simulate, CASES))


def draw(curves, output):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MultipleLocator

    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#899096", "axes.labelcolor": "#273139",
                         "xtick.color": "#454e52", "ytick.color": "#454e52"}):
        fig, ax = plt.subplots(figsize=(9, 5.9))
        fig.subplots_adjust(left=.11, right=.97, bottom=.24, top=.79)
        fig.text(.11, .94, "How much of the wall’s energy becomes random particle motion?",
                 ha="left", fontsize=13.8, weight="bold", color="#172b35")
        fig.text(.11, .887, r"Three expanding-wall simulations with the same source and daughter mass $mL=4\times10^{12}$",
                 ha="left", fontsize=10.5, color="#4b555a")
        ax.axhline(90, color="#59666c", linewidth=1, linestyle=(0, (2, 3)))
        ax.text(3.3, 91.8, "90% hot threshold", fontsize=10, color="#475357")
        for name, _, label, color, style in reversed(CASES):
            r, y, ref = curves[name]
            ax.plot(r, y, label=label, color=color, linestyle=style,
                    linewidth=1.9 if name == "pole_cap" else 1.35)
            ax.plot(r[-1], y[-1], marker="o", ms=4.5, mec=color, mfc="white", zorder=5)
        ref = curves["pole_cap"][2]
        start = ref["first_hot_radius"]
        end = start * ref["hot_radius_factor"]
        ax.annotate("", (end, 80), (start, 80),
                    arrowprops=dict(arrowstyle="|-|", color="#00666a", lw=1.1))
        ax.text((start+end)/2, 74.9, f'Radius grows by a factor of {ref["hot_radius_factor"]:.2f}',
                color="#00666a", ha="center", va="top", fontsize=10)
        ax.set(xlim=(1, 36), ylim=(-2, 103), xlabel="Wall radius in units of L",
               ylabel="Added wall energy in random particle motion (%)")
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.grid(axis="y", alpha=.17)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc="lower left", bbox_to_anchor=(-.01, 1.015),
                  ncol=3, frameon=False, fontsize=10, borderpad=0, columnspacing=1.5)
        fig.text(.11, .094, "Open circles mark the start of recollapse. The horizontal bracket marks the fixed pole’s",
                 fontsize=9.6, color="#4b555a")
        fig.text(.11, .055, "continuous post-exit interval above 90%. All three runs use the same cap on particle births.",
                 fontsize=9.6, color="#4b555a")
        fig.savefig(output, dpi=180, facecolor="white", metadata={"Software": "Matplotlib"})
        plt.close(fig)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--base", type=Path, default=HERE.parent / "prior")
    cli.add_argument("--feedback", type=Path, default=HERE.parent / "feedback")
    cli.add_argument("--table", type=Path, default=HERE / "BOOST_TABLE.dat")
    cli.add_argument("--archive", type=Path, default=HERE / "EVOLUTION.json")
    cli.add_argument("--jobs", type=int, default=2)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    curves = replay(args.base, args.feedback, args.table, args.archive, args.jobs)
    draw(curves, args.output)
