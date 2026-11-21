"""Test whether any angular birth rate can preserve zero directed momentum on the heavy histories."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
import mmap
from pathlib import Path
import struct

HERE = Path(__file__).resolve().parent
STEP = struct.Struct("=12d")
R = Q(4, 3)


def primitives(u):
    """Integrals of 1+r*u and u*(1+r*u) with arbitrary lower limit."""
    return u + R*u*u/2, u*u/2 + R*u**3/3


def clipped(u):
    """Integrate |1+r*u| and u*|1+r*u| from -1 to u, with u >= 0."""
    start, cusp, end = primitives(Q(-1)), primitives(-1/R), primitives(u)
    return tuple(a - 2*b + c for a, b, c in zip(start,cusp,end))


def moment_bound():
    lo, hi = Q(0), Q(1)
    for _ in range(60):
        mid = (lo + hi)/2
        if clipped(mid)[1] < 0:
            lo = mid
        else:
            hi = mid
    lower, upper = clipped(lo)[0]/2, clipped(hi)[0]/2
    assert clipped(lo)[1] < 0 < clipped(hi)[1] and upper-lower < Q(1,10**17)
    # Separately integrate the original, unsimplified functions. A two-point
    # Gauss rule is exact through degree three on each side of the absolute value.
    cut = float((hi+lo)/2)
    mass = moment = 0.
    cusp = -1/float(R)
    for a,b in [(-1.,cusp),(cusp,cut)]:
        half,mid = (b-a)/2,(a+b)/2
        for u in [mid-half/math.sqrt(3),mid+half/math.sqrt(3)]:
            g = abs(1+float(R)*u)
            mass += half*g/2
            moment += half*u*g/2
    assert abs(moment)<1e-15 and abs(mass-float(upper))<1e-15
    return {
        "d_over_c": float(R), "angular_cut": cut,
        "cut_bracket_exact": [str(lo), str(hi)],
        "zero_momentum_capacity_over_c": float(upper),
        "capacity_over_c_upper_exact": str(upper),
        "quadrature_mean_moment_residual": moment,
    }


def find_witness(path, reference_peak):
    with path.open("rb") as stream, mmap.mmap(stream.fileno(),0,access=mmap.ACCESS_READ) as raw:
        crossings = []
        for row in STEP.iter_unpack(raw):
            t,dt,k0,k1,d0,d1,drest,f0,f1,c,k,h = row
            if d0*d1 < 0:
                if c<=0 or h<=0:
                    raise ValueError("Invalid shell speed or demand at an electric reversal")
                crossings.append((h/c,t,dt,d0/c,d1/c,f0,f1))
        digest = hashlib.sha256(raw).hexdigest()
    h_over_c,t,dt,d0,d1,f0,f1 = max(crossings)
    # Establish whether the highest h/c among reversing steps equals the worst
    # angle-integrated ratio anywhere on this published grid.
    if not math.isclose(h_over_c,reference_peak,rel_tol=5e-13):
        raise ValueError("Selected reversal differs from replayed summed peak")
    witnesses=[]
    for signed in [-float(R),float(R)]:
        portion = (signed-d0)/(d1-d0)
        if not 0<portion<1:
            raise ValueError("The reversing step does not reach the chosen signed field")
        witnesses.append({"signed_field_over_c": signed, "tau": t+portion*dt,
                          "cumulative_fraction": f0+portion*(f1-f0)})
    return {
        "trace_sha256": digest, "electric_reversal_count": len(crossings),
        "step_start_tau": t, "step_end_tau": t+dt, "step_dt": dt,
        "signed_field_over_c_endpoints": [d0,d1], "h_over_c": h_over_c,
        "witnesses": sorted(witnesses,key=lambda value:value["tau"]),
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--slices",type=Path,required=True,help="Directory retained by replay.py --slices")
    p.add_argument("--output",type=Path,default=HERE/"MOMENT.json")
    args=p.parse_args()
    base=json.loads((HERE/"results.json").read_text())["runs"]["mass_4e12"]
    pulse=json.loads((HERE.parent/"pulse/PULSE_MODES.json").read_text())
    drag=json.loads((HERE.parent/"transport/DRAG.json").read_text())["runs"]["mass_4e12"]
    trial=base["action"]
    if trial != pulse["actions"]["4e12"] or base["drag_floor"] != drag["reference_drag_floor"]:
        raise ValueError("Receipts use different normalizations")
    law=moment_bound()
    table={}
    for grid in ["native","half_step"]:
        source=args.slices / f"mass_4e12-{grid}.bin"
        peak=base["grids"][grid]["central99_998"]["witness"]["aggregate"]["peak"]
        w=find_witness(source,peak)
        ratio=w["h_over_c"]/law["zero_momentum_capacity_over_c"]
        summed_at_point=2*w["h_over_c"]/(float(R)+1/float(R))
        if not summed_at_point<1<ratio:
            raise ValueError("Summed and moment-constrained capacities fail witness inequalities")
        ceiling=trial/ratio
        min_drag=base["drag_floor"]/ceiling
        w.update({
            "summed_demand_over_capacity_at_witness": summed_at_point,
            "zero_momentum_demand_over_capacity": ratio,
            "necessary_action_upper_from_this_witness": ceiling,
            "required_drag_fraction_at_action_upper": min_drag,
            "elastic_pole_shortfall_factor": min_drag/drag["elastic_fraction_upper"],
            "inclusive_pole_shortfall_factor": min_drag/drag["inclusive_fraction_upper"],
            "pulse_band_fraction_of_parent_at_action_upper": pulse["band_energy_to_classical_parent"]["4e12"]*ratio,
        })
        table[grid]=w
    result={"mass_mL":4e12,"trial_action":trial,"angular_shape":law,"grids":table}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
