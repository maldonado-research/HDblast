"""Check temporal and radial five-dimensional planar mass balances from the public tables."""
import argparse
import csv
import json
import math
from pathlib import Path


LABELS = ("inner", "middle", "wall")
TERMS = ("potential", "kinetic", "mixed")
HERE = Path(__file__).resolve().parent


def read(path, coordinate, heading):
    with path.open(encoding="ascii", newline="") as stream:
        raw = csv.DictReader(stream, strict=True)
        if raw.fieldnames != heading:
            raise ValueError(f"{path.name}: unexpected header")
        rows = [{key: float(value) for key, value in row.items()} for row in raw]
    if len(rows)<3 or any(not math.isfinite(q) for row in rows for q in row.values()):
        raise ValueError(f"{path.name}: nonfinite or insufficient samples")
    if any(a[coordinate] >= b[coordinate] for a,b in zip(rows,rows[1:])):
        raise ValueError(f"{path.name}: nonincreasing coordinate")
    return rows


def trap(rows, coordinate, column):
    return math.fsum((b[coordinate]-a[coordinate])*(a[column]+b[column])/2
                     for a,b in zip(rows,rows[1:]))


def compare(rows, coordinate, dimension, label):
    keys = [term+"_"+dimension+"_"+label for term in TERMS]
    cols = dict(zip(TERMS,[trap(rows,coordinate,key) for key in keys]))
    integral = math.fsum(cols.values())
    mu = "mu_"+label
    delta = rows[-1][mu]-rows[0][mu]
    error = integral-delta
    thin = rows[::2]
    if thin[-1] is not rows[-1]:
        thin.append(rows[-1])
    thin_integral = math.fsum(trap(thin,coordinate,key) for key in keys)
    return {
        "points":len(rows),"mu_start":rows[0][mu],"mu_end":rows[-1][mu],
        "delta_mu":delta,"integrated_rhs":integral,"integrated_rhs_components":cols,
        "rhs_minus_delta":error,"fraction_of_absolute_delta":error/abs(delta) if delta else None,
        "every_other_node_integral":thin_integral,
    }


def summarize(folder, grid):
    time_cols = ["s","T"] + [c+"_"+w for w in LABELS for c in ("mu",*(v+"_T" for v in TERMS))]
    radial_cols = ["z"] + [c+"_"+w for w in ("initial","final") for c in ("mu",*(v+"_z" for v in TERMS))]
    times = read(folder/(grid+"_time.csv"),"T",time_cols)
    radial = read(folder/(grid+"_radius.csv"),"z",radial_cols)
    if not (times[0]["s"]==1 and times[-1]["s"]==4 and
            all(a["s"]<b["s"] for a,b in zip(times,times[1:])) and radial[-1]["z"]==0):
        raise ValueError("Expected shell clock interval and wall")
    tchecks = {name:compare(times,"T","T",name) for name in LABELS}
    collars = {}
    for target,label in zip((-0.04,-0.02),LABELS):
        i = min(range(len(radial)),key=lambda j:abs(radial[j]["z"]-target))
        if abs(radial[i]["z"]-target)>1e-8:
            raise ValueError("Wrong native surface")
        cut = radial[i:]
        init, final = [compare(cut,"z","z",when) for when in ("initial","final")]
        if not (init["mu_start"]==times[0]["mu_"+label] and final["mu_start"]==times[-1]["mu_"+label]
                and init["mu_end"]==times[0]["mu_wall"] and final["mu_end"]==times[-1]["mu_wall"]):
            raise ValueError("Geometric energies differ at shared corners")
        delta = final["delta_mu"]-init["delta_mu"]
        temporal_source = tchecks["wall"]["integrated_rhs"]-tchecks[label]["integrated_rhs"]
        radial_source = final["integrated_rhs"]-init["integrated_rhs"]
        collars[label] = {
            "inner_z":cut[0]["z"],"outer_z":cut[-1]["z"],"initial":init,"final":final,
            "change_in_geometric_mass_contrast":delta,
            "two_temporal_integrals_wall_minus_inner":temporal_source,
            "two_radial_integrals_final_minus_initial":radial_source,
            "temporal_minus_radial_integrals":temporal_source-radial_source,
        }
    return {"shell_clock_interval":[times[0]["s"],times[-1]["s"]],
            "coordinate_T_interval":[times[0]["T"],times[-1]["T"]],
            "fixed_z":tchecks,"radial_collars":collars}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tables",nargs="?",type=Path,default=HERE/"data")
    args = parser.parse_args()
    print(json.dumps({g:summarize(args.tables,g) for g in ("coarse","fine")},indent=2,allow_nan=False))
