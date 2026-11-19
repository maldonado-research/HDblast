"""Check projected state changes and use the distributed response table independently."""
import csv
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def check_close(x, y, tolerance):
    if not math.isfinite(x) or not math.isfinite(y) or abs(x-y)>tolerance:
        raise ValueError((x,y,tolerance))


def integrate(y, dt):
    return dt/3*(y[0]+y[-1]+4*math.fsum(y[1:-1:2])+2*math.fsum(y[2:-1:2]))


def main():
    with (HERE/"RESPONSE.csv").open() as f:
        rows=[tuple(map(float,r)) for r in list(csv.reader(f))[1:]]
    report=json.loads((HERE/"KERNEL.json").read_bytes())
    if len(rows)!=2049:
        raise ValueError("Incomplete quiet kernel")
    dt=(2.1-.8)/2048
    for i,(s,k0,k2) in enumerate(rows):
        check_close(s,.8+i*dt,1e-13)
        if not math.isfinite(k0+k2):raise ValueError("Kernel samples")
    for mesh in report.values():
        for arm in mesh["arms"].values():
            linear=arm["components"]["linear"]
            exact=arm["components"]["exact"]
            check_close(linear["total_interference"],sum(linear[k] for k in
                        ("initial_basis","terminal_basis","mass_evolution","geometry_evolution")),1e-12)
            check_close(arm["linear_to_finite_ratio"],linear["total_interference"]/arm["delta_N"],1e-12)
            if abs(arm["linear_to_finite_ratio"]-1)>.0035: raise ValueError("Linear response accuracy")
            check_close(arm["delta_N"],exact["total_interference"]+exact["quadratic_amplitude"],3e-8)
            if "centered_derivative_at_step_0.01" in arm:
                check_close(arm["centered_derivative_at_step_0.01"],linear["total_interference"],2e-6)
        a,b=(mesh["arms"][n]["cumulative_linear_delta_N"] for n in ("near_plus","far_plus"))
        for term in ("mass","geometry"):
            check_close(a[term]["before_near_compensation"],b[term]["before_near_compensation"],4e-7)
    design=report["fine"]["designed_fixed_geometry_controls"]
    for c, target in zip(design["candidates"],design["trial_predictions"],strict=True):
        width=design["width"]
        weight=[math.sin(math.pi*((t-c)/width+.5))**4 if abs(t-c)<width/2 else 0. for t,_,_ in rows]
        check_close(integrate(weight,dt),design["area_per_unit_height"],1e-9)
        check_close(integrate([w*r[1] for w,r in zip(weight,rows)],dt),target,3e-8)
    for key,center in (("lowest",1.70),("highest",1.62)):
        sample=design["controls"][key]
        check_close(sample["center"],center,1e-12)
        if sample["predicted_delta_N_per_height"]*sample["direct_delta_N"]["1.0"]<=0:
            raise ValueError("Timing control sign")
    print("Verified the causal occupation response and 16 independent projections of the quiet kernel.")


if __name__=="__main__":
    main()
