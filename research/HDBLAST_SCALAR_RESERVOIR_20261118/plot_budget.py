"""Redraw the fine-grid matter budget and the late-time expansion comparison."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def plot():
    ledger = json.loads((ROOT / "ledger/RESERVOIR.json").read_text())["runs"]["fine"]["interval"]
    work = ledger["mass_work"] + ledger["filling_work"]
    keys = ["radiation_gain", "radiation_pressure_loss", "particle_gain", "particle_pressure_loss"]
    shares = [100 * ledger[key] / work for key in keys]
    labels = ["Radiation retained", "Radiation pressure work", "Particle energy retained", "Particle pressure work"]
    colors = ["#21676a", "#b46e2b", "#5377a3", "#857397"]
    with (ROOT / "clock/fine.csv").open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if float(r["s"]) >= 4]
    times = [float(r["s"]) for r in rows]
    reference = json.loads((ROOT/"clock/RESULTS.json").read_text())["fine"]
    normalization = float(reference["best_H_fit_to_a_minus_two"]["coefficient"])
    from math import exp, sqrt
    fig,(budget,history) = plt.subplots(1,2,figsize=(11,3.9),layout="constrained",
                                         gridspec_kw={"width_ratios":[1,1.13]})
    budget.barh(labels,shares,color=colors,height=0.55)
    budget.invert_yaxis()
    for i, number in enumerate(shares):
        budget.text(number+1.2,i,f"{number:.2f}%",va="center",fontsize=9)
    budget.set(xlim=(0,99),xlabel="Share of scalar work (%)",title="Where the work goes,  s = 1–4")
    history.plot(times,[float(r["H_over_H0"]) for r in rows],color="#172d42",lw=2.4,label="Simulated expansion")
    history.plot(times,[normalization*exp(-2*float(r["ln_a"])) for r in rows],
                 color="#c06a22",ls="--",lw=1.9,label="Fitted radiation-like expansion")
    history.plot(times,[sqrt(float(r["rad"])) for r in rows],color=colors[0],lw=2,
                 label="Material radiation term, square root")
    history.set(xlabel="Dimensionless proper time   $s=H_0\\tau$",ylabel="Rate in units of $H_0$",
                ylim=(0,.106),title="What drives expansion,  s = 4–6")
    history.legend(loc="center",fontsize=8,frameon=False)
    for ax in (budget,history):
        ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)
    history.grid(axis="y",alpha=.20)
    fig.savefig(ROOT/"scalar_budget.png",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    plot()
