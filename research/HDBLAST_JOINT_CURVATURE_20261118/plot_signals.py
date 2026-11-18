"""Draw the scalar and the two positive-weight null Ricci averages on the same clock."""
import csv
import json
from fractions import Fraction as F
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
certificate = json.loads((ROOT / "METRIC.json").read_text())
radius = float(F(certificate["absolute_mean_error"]) + F(1, 10**11))
ell = 1e-6
xs = [0.0, 1.0]
scalar_mid = [0.0, 0.0]
z = np.zeros(2)
speed = np.zeros(2)
curv_lo, curv_hi = [0.0, 0.0], [0.0, 0.0]

with (ROOT / "inputs" / "FUTURE_CELLS.csv").open() as stream:
    for row in csv.DictReader(stream):
        left, right = (int(row[key]) / 4096 for key in ("left_grid", "right_grid"))
        if right <= 1:
            continue
        if right > 25/16:
            break
        acc = np.array([int(row["H_lo"]) / 2**112 - 3e-6, int(row["H_hi"]) / 2**112 + 3e-6])
        dt = right-left
        z += speed*dt + acc*(dt**2/2)
        speed += acc*dt
        xs.append(right)
        scalar_mid.append(np.mean(acc) / np.cosh(ell*right)**1.5)
        curv_lo.append((16/27)*(z[0]-radius))
        curv_hi.append((16/27)*(z[1]+radius))

xs = np.array(xs)
curv_lo = np.array(curv_lo)
curv_hi = np.array(curv_hi)
mid = (curv_lo+curv_hi)/2
fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.15), layout="constrained")
teal, purple, warm = "#227973", "#57478a", "#b14325"
for ax in axes:
    ax.axvspan(5/4, 11/8, color="#e0e6e8", alpha=0.8)
    ax.axvspan(49/32, 25/16, color="#e0e6e8", alpha=0.8)
    ax.axhline(0, color="#8f9da2", lw=0.8)
    ax.set(xlabel="Elapsed model time  x", xlim=(0.9,25/16))
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", lw=0.5, alpha=0.4)
    ax.set_axisbelow(True)
axes[0].plot(xs,scalar_mid,color=purple,lw=2)
axes[0].set(title="Scalar quadrupole",ylabel="Normalized scalar  A",ylim=(-2.5,5.0))
axes[1].fill_between(xs,curv_lo,curv_hi,color=teal,alpha=0.30)
axes[1].fill_between(xs,-curv_hi,-curv_lo,color=warm,alpha=0.30)
axes[1].plot(xs,mid,color=teal,lw=2,label="Direction weighting  −")
axes[1].plot(xs,-mid,color=warm,lw=2,label="Direction weighting  +")
axes[1].set(title="Intrinsic null Ricci curvature",ylabel="Normalized first variation",ylim=(-0.17,0.17))
axes[1].legend(frameon=False,loc="upper left",fontsize=9)
fig.savefig(ROOT/"joint_signals.png",dpi=175,facecolor="white")
