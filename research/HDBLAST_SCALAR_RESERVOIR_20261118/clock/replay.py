"""Compare the expansion history with flat radiation expansion on the archived clock."""
import csv
import json
from decimal import Decimal as D, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
PI = D("3.1415926535897932384626433832795028841971693993751")


def best_scale(ratios):
    """Minimize max |A/q - 1| over positive q and A."""
    lo, hi = min(ratios), max(ratios)
    if lo <= 0:
        raise ValueError("The calibration requires expanding states with positive density")
    return {"coefficient": 2*hi*lo/(hi+lo), "max_fractional_error": (hi-lo)/(hi+lo)}


def evaluate(grid, spec):
    with (HERE / (grid+".csv")).open(newline="") as stream:
        data = [{key:D(value) for key,value in row.items()} for row in csv.DictReader(stream)]
    if data[0]["s"]!=1 or data[-1]["s"]!=6 or any(a["s"]>=b["s"] for a,b in zip(data,data[1:])):
        raise ValueError("Invalid clock grid")
    late = [r for r in data if r["s"] >= 4]
    if late[0]["s"]!=4 or len(late)<2:
        raise ValueError("Missing late interval")
    def h(r):return r["H_over_H0"]
    def fraction(r,key):return r[key]/h(r)**2
    start, end = late[0],late[-1]
    slope = (1/h(end)-1/h(start))/(end["s"]-start["s"])
    # Only the reference temperature uses g_star, with all radiation in one equilibrium bath.
    def reference_temperature(r):
        rho_over_H04 = D(str(spec["rho_b"]))*r["R"]/D(str(spec["b"]))
        return (30*rho_over_H04/(PI**2*D("10.75"))).sqrt().sqrt()
    early=[r for r in data if r["s"]<=4]
    peak=max(early,key=lambda r:r["R"])
    out={
        "sample_count_4_to_6":len(late),
        "best_H_fit_to_sqrt_R": best_scale([h(r)/r["R"].sqrt() for r in late]),
        "best_H_fit_to_a_minus_two": best_scale([h(r)*(2*r["ln_a"]).exp() for r in late]),
        "inverse_H_secant":slope,
        "effective_w_secant_if_flat_GR":2*slope/3-1,
        "radiative_fraction_of_H2_range":[min(fraction(r,"rad") for r in late),max(fraction(r,"rad") for r in late)],
        "weyl_fraction_of_H2_range":[min(fraction(r,"Wy") for r in late),max(fraction(r,"Wy") for r in late)],
        "vacuum_fraction_of_H2_range":[min(fraction(r,"vac") for r in late),max(fraction(r,"vac") for r in late)],
        "fractional_change_a4_R":(4*(end["ln_a"]-start["ln_a"])).exp()*end["R"]/start["R"]-1,
        "e_folds":end["ln_a"]-start["ln_a"],
        "peak_sample_clock_1_to_4":peak["s"],
        "peak_density_temperature_over_H0_if_gstar_10p75":reference_temperature(peak),
        "temperature_over_H0_at_4_if_gstar_10p75":reference_temperature(start),
    }
    return out


if __name__ == "__main__":
    specs=json.loads((HERE/"INPUTS.json").read_text())
    with localcontext() as context:
        context.prec=48
        result={k:evaluate(k,spec) for k,spec in specs.items()}
    print(json.dumps(result,indent=2,default=lambda value:format(value,".14g")))
