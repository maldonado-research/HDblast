"""Reduce the passive B1 collar observer's full fields to native geometric and stress columns."""
import argparse
import csv
from pathlib import Path

import numpy as np


CUTS = (1.0, 4.0)
TARGETS = (-0.04, -0.02, 0.0)
LABELS = ("inner", "middle", "wall")


def potential(phi):
    # HDBLAST static_w.py, same field normalization as G_ab = T_ab.
    return ((((((-2/27)*phi)*phi + 17/18)*phi - 4/9)*phi - 5/3)*phi + 4/3)*phi - 1/6


def radial_weights(z, start, order):
    """Polynomial derivatives on the physical z grid, without ghost values."""
    weights = np.zeros((len(z)-start, len(z)))
    n = len(z)
    for i in range(start, n):
        if i == n - 1:
            # Five strictly interior samples for the main calculation. The wall value is excluded.
            stencil = np.arange(n - order - 1, n - 1)
        else:
            first = max(0, min(i - order//2, n - order))
            stencil = np.arange(first, first + order)
        dx = z[stencil] - z[i]
        scale = max(abs(dx))
        vandermonde = np.array([(dx / scale)**n for n in range(order)])
        derivative = np.zeros(order)
        derivative[1] = 1 / scale
        weights[i-start, stencil] = np.linalg.solve(vandermonde, derivative)
    return weights


def geometric_and_matter(fields, start, derivative):
    A, B, phi, at, phit = (fields[name] for name in ("A", "B", "phi", "A_T", "phi_T"))
    az = derivative @ A
    phiz = derivative @ phi
    a, b, p, at, phit = (col[start:] for col in (A, B, phi, at, phit))
    r4 = np.exp(4*a)
    inv_lapse = np.exp(-2*b)
    kinetic = inv_lapse*(phit**2 + phiz**2)/2
    mixed = inv_lapse*phit*phiz
    u = potential(p)
    mu = r4*inv_lapse*(at**2 - az**2)
    scale = (2/3)*r4
    return np.array([mu, scale*u*at, -scale*kinetic*at, scale*mixed*az,
                     scale*u*az, scale*kinetic*az, -scale*mixed*at])


def write_csv(path, header, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([format(float(x), ".17g") for x in row])


def interpolate(clock, records, target):
    if not clock[0] <= target <= clock[-1]:
        raise ValueError("Requested cut outside the sampled clock")
    i = max(1, int(np.searchsorted(clock, target)))
    f = (target-clock[i-1])/(clock[i]-clock[i-1])
    return records[i-1] + f*(records[i] - records[i-1])


def extract(collar, out, grid, order):
    files = sorted(collar.glob("sample_*.npz"))
    if not files:
        raise ValueError("Observer sample files are missing")
    times, clocks, records = [], [], []
    with np.load(files[0], allow_pickle=False) as fields:
        z = fields["z"]
    locations = [int(np.argmin(abs(z-q))) for q in TARGETS]
    if z[-1] != 0 or not np.all(np.diff(z)>0) or max(abs(z[i]-q) for i,q in zip(locations,TARGETS))>1e-8:
        raise ValueError("Expected the two nominal radii to be native nodes, within 1e-8")
    start = locations[0]
    deriv = radial_weights(z, start, order)
    for index, path in enumerate(files):
        if path.name != f"sample_{index:06d}.npz":
            raise ValueError("Observer cadence has a missing file")
        with np.load(path, allow_pickle=False) as fields:
            if not np.array_equal(z,fields["z"]):
                raise ValueError("Coordinate grid varies in time")
            clocks.append(float(fields["s"]))
            times.append(float(fields["T"]))
            records.append(geometric_and_matter(fields,start,deriv))
        if clocks[-1] >= CUTS[-1]:
            break
    if not all(np.all(np.diff(seq)>0) for seq in (times,clocks)) or not clocks[0] <= 1 < 4 <= clocks[-1]:
        raise ValueError("Incomplete window or nonmonotone observer clocks")
    clocks = np.array(clocks)
    times = np.array(times)
    records = np.array(records)
    if not np.isfinite(records).all():
        raise ValueError("Nonfinite geometric or stress term")
    chosen = (clocks > 1) & (clocks < 4)
    t = np.r_[interpolate(clocks,times,1),times[chosen],interpolate(clocks,times,4)]
    s = np.r_[1,clocks[chosen],4]
    lo, hi = (interpolate(clocks,records,cut) for cut in CUTS)
    data = np.concatenate(([lo], records[chosen], [hi]))
    indices = np.array(locations)-start
    tt = np.column_stack([s,t, *[data[:,j,k] for k in indices for j in (0,1,2,3)]])
    header = ["s","T"] + [col+"_"+label for label in LABELS for col in ("mu","potential_T","kinetic_T","mixed_T")]
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/(grid+"_time.csv"),header,tt)
    zz = np.column_stack([z[start:], *[cut[j] for cut in (lo,hi) for j in (0,4,5,6)]])
    labels = ["z"] + [col+"_"+when for when in ("initial","final") for col in ("mu","potential_z","kinetic_z","mixed_z")]
    write_csv(out/(grid+"_radius.csv"),labels,zz)
    print(f"{grid} native z surfaces " + ", ".join(format(z[i],".17g") for i in locations))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("collar", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("grid", choices=("coarse","fine"))
    p.add_argument("--order", type=int, choices=(3,5), default=5,
                   help="Number of polynomial stencil points. Public tables use five.")
    args=p.parse_args()
    extract(args.collar,args.output,args.grid,args.order)
