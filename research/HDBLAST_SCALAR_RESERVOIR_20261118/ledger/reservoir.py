"""Replay both component budgets from derived CSVs using standard-library arithmetic."""
from bisect import bisect_left
import csv
import hashlib
import io
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PINS = {
    "coarse": "8c9fd6f7604fdd1e509c461f0a0c47e77805eb9f116f251f5ac44b4ffa22113a",
    "fine": "7dd4b6e3bc4b45840e05e5af1ba0ee6d7fe731053b522b183128f59e4d0949b0",
}
HEADERS = ["s","chi_energy","radiation_energy","mass_work_rate","filling_work_rate",
           "decay_rate","chi_pressure_rate","radiation_pressure_rate"]
CUTS = [1+3*i/8 for i in range(9)]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(label):
    raw = (ROOT / ("reservoir_"+label+".csv")).read_bytes()
    require(hashlib.sha256(raw).hexdigest()==PINS[label], "unrecognized "+label+" inputs")
    table = csv.reader(io.StringIO(raw.decode("ascii")))
    require(next(table)==HEADERS, "column definitions")
    rows = [list(map(float,row)) for row in table]
    require(len(rows)==155 and rows[0][0]==1 and rows[-1][0]==4, "sample set")
    for row in rows:
        require(len(row)==len(HEADERS) and all(map(math.isfinite,row)), "bad sample")
    require(all(a[0]<b[0] for a,b in zip(rows,rows[1:])), "proper times")
    return rows


def value(rows, clock, time):
    i = bisect_left(clock,time)
    if clock[i]==time:
        return rows[i]
    a,b=rows[i-1],rows[i]
    t=(time-a[0])/(b[0]-a[0])
    return [time]+[start+t*(end-start) for start,end in zip(a[1:],b[1:])]


def integrate(rows, low, high):
    clock=[r[0] for r in rows]
    selected=[value(rows,clock,low)]+[r for r in rows if low<r[0]<high]+[value(rows,clock,high)]
    terms=[math.fsum((b[0]-a[0])*(a[k]+b[k])/2 for a,b in zip(selected,selected[1:]))
           for k in range(3,len(HEADERS))]
    mass,fill,decay,chi_pressure,rad_pressure=terms
    chi=selected[-1][1]-selected[0][1]
    rad=selected[-1][2]-selected[0][2]
    gain=chi+rad
    chi_error=mass+fill-decay-chi_pressure-chi
    rad_error=decay-rad_pressure-rad
    total_error=mass+fill-chi_pressure-rad_pressure-gain
    return dict(interval=[low,high], particle_gain=chi, radiation_gain=rad, total_gain=gain,
                mass_work=mass, filling_work=fill, decay=decay, particle_pressure_loss=chi_pressure,
                radiation_pressure_loss=rad_pressure, particle_residual=chi_error,
                radiation_residual=rad_error, total_residual=total_error)


def analyze(label):
    rows=read(label)
    total=integrate(rows,1,4)
    bins=[integrate(rows,lo,hi) for lo,hi in zip(CUTS,CUTS[1:])]
    for k in total.keys()-{"interval"}:
        require(abs(math.fsum(cell[k] for cell in bins)-total[k])<1.e-12,"bin sum "+k)
    gross=total["mass_work"]+total["filling_work"]
    error=total["total_residual"]
    # All control residuals are signed integrated inputs minus measured gain.
    controls=dict(omit_mass_work=error-total["mass_work"], omit_filling_work=error-total["filling_work"],
                  omit_both_pressure_losses=error+total["particle_pressure_loss"]+total["radiation_pressure_loss"],
                  double_count_decay_as_external_input=error+total["decay"])
    thinned=[]
    for stride in (2,4):
        down=rows[::stride]
        if down[-1] != rows[-1]:
            down=down+[rows[-1]]
        result=integrate(down,1,4)
        thinned.append(dict(stride=stride, grid_points=len(down), signed_total_residual=result["total_residual"]))
    return dict(csv_sha256=PINS[label], interval=total, work_mass_fraction=total["mass_work"]/gross,
                work_filling_fraction=total["filling_work"]/gross, residual_over_gain=error/total["total_gain"],
                residual_over_work=error/gross, pressure_over_work=(total["particle_pressure_loss"]+
                    total["radiation_pressure_loss"])/gross, controls=controls, reduced_cadence=thinned,
                max_abs_subinterval_total_residual=max(abs(b["total_residual"]) for b in bins), subintervals=bins)


def main():
    runs={label:analyze(label) for label in PINS}
    big=runs["fine"]["interval"];small=runs["coarse"]["interval"]
    require(0<big["total_residual"]<0.00001 and 0.60<runs["fine"]["work_mass_fraction"]<0.61, "rounded results")
    output=dict(runs=runs, fine_vs_coarse_fractional_total_gain=big["total_gain"]/small["total_gain"]-1)
    print(json.dumps(output,indent=2))


if __name__=="__main__":
    main()
