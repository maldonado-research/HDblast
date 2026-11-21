"""Draw the hot drag requirements and single-resonance ceilings from the drag and momentum receipts."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


HERE = Path(__file__).resolve().parent


def main():
    run = json.loads((HERE/"transport"/"DRAG.json").read_text())["runs"]["mass_4e12"]
    moment = json.loads((HERE/"replay"/"MOMENT.json").read_text())
    if run["mass_mL"] != 4e12 or moment["mass_mL"] != run["mass_mL"]:
        raise ValueError("Figure is labeled for daughter mass mL = 4 × 10¹²")
    momentum_bound = moment["grids"]["half_step"]["required_drag_fraction_at_action_upper"]
    ink = "#162c3a"
    muted = "#42545b"
    teal = "#14746e"
    orange = "#a44f1c"
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(9.3, 5.35))
    fig.subplots_adjust(left=.39, right=.96, top=.78, bottom=.23)
    fig.text(.037, .936, "The parent resonance falls short of the required drag",
             fontsize=17, weight="bold", color=ink)
    fig.text(.037, .873, "HDblast at " + r"$mL=4\times10^{12}$" +
             " with one spin-zero resonance and three equally coupled, unpolarized Dirac flavors.",
             fontsize=9.5, color=muted)

    items = [
        (3.73, "Same outgoing flavor", run["elastic_fraction_upper"], teal),
        (3.09, "All three outgoing flavors", run["inclusive_fraction_upper"], teal),
        (1.71, "After summing over directions", run["required_fraction_at_sweep_ceiling"], orange),
        (1.07, "At the published example action", run["required_fraction_at_published_action"], muted),
        (.43, "With zero mean birth momentum per stream", momentum_bound, orange),
    ]
    ax.set_xlim(0, 39.0)
    ax.set_ylim(-.03, 4.55)
    ax.set_xticks([0, 10, 20, 30])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}%"))
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#dbe2e4", linewidth=.7)
    ax.spines[:].set_visible(False)
    ax.tick_params(axis="both", length=0, labelsize=9, colors=muted)
    ax.tick_params(axis="y", pad=10, labelcolor=ink)
    ax.set_yticks([v[0] for v in items], [v[1] for v in items])
    ax.set_xlabel("Drag rate as a share of the reference elastic unitarity ceiling",
                  fontsize=9.3, labelpad=6, color=ink)
    for y, _, fraction, color in items:
        value = 100 * fraction
        ax.barh(y, value, height=.26, color=color)
        label = f"{round(value)}%" if abs(value-round(value)) < 1e-10 else f"≈{value:.2f}%"
        ax.annotate(label, (value,y), xytext=(6,0), textcoords="offset points",
                    va="center", color=color, weight="bold", fontsize=10.2)
    for y,label in [(4.32,"Maximum from the resonance"),(2.32,"Required by the hot fluid")]:
        ax.text(-.617,y,label,transform=ax.get_yaxis_transform(),color=ink,
                weight="bold",fontsize=10.5,ha="left",va="center",clip_on=False)
    fig.text(.037, .105, "The bounds of 16.46% and 33.60% use the available fermion states on a sharp birth shell. "
             "The latter also enforces", color=muted, fontsize=8.4)
    fig.text(.037, .071, "zero mean injected momentum within each charge stream. "
             "The 20% example uses the published action before this condition.", color=muted, fontsize=8.4)
    fig.text(.037, .030, "The reference ceiling uses hot, small-drift elastic scattering. The modeled drag rate is "
             + r"$\nu L=18$" + ".", color=muted, fontsize=8.4)
    fig.savefig(HERE/"hot-drag.png", dpi=190, facecolor="white")
    fig.savefig(HERE/"hot-drag.pdf", facecolor="white",
                metadata={"Title":"HDblast hot drag and single-resonance ceilings",
                          "CreationDate":None, "ModDate":None})
    plt.close(fig)


if __name__ == "__main__":
    main()
