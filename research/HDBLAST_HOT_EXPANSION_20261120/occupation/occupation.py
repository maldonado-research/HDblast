"""Finite sum of excited photon-shell energies on the prescribed parity-odd wall trajectory."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

def from_gain(maximum):
    """Occupation reconstructed from reciprocal energy singular values."""
    return (maximum-1)**2/(4*maximum)

def even_budget(transfer):
    end = transfer["endpoint_reference"]
    denominator = 2*math.pi**2*end["comoving_sheet_energy"]
    shells = {}
    for label, stages in transfer["modes"].items():
        k = int(label)
        if 2 < k <= 1024:
            n = from_gain(stages["combined"]["max_energy_gain"])
            shells[label] = 2*(k*k-1)*k*n/denominator
    return dict(action_for_selected_shells_equal_coherent=math.fsum(shells.values()),
                max_wavenumber=1024, number_of_shells=len(shells), action_by_shell=shells)

def calculate(transfer, handoff, two_stage_transfer):
    endpoint = dict(zip(handoff["columns"], handoff["rows"][-1], strict=True))
    a, k_final = endpoint["a"], endpoint["K"]
    coherent = k_final*(endpoint["B"]**2+endpoint["X"]**2)/2
    shells = {}
    for label, helicities in transfer["modes"].items():
        k = int(label)
        if k == 2:
            continue  # The homogeneous background already uses the lowest multiplet.
        corrections = {}
        for h in ("-1","1","0"):
            # For a symplectic map between normalized free oscillator coordinates,
            # the vacuum occupation above the final zero point is (mean_gain-1)/2.
            n = from_gain(helicities[h]["max_energy_gain"])
            if n < -1e-12:
                raise ValueError("Passive rotation has negative vacuum occupation")
            corrections[h] = max(n,0.)*k*(k*k-1)/(2*math.pi**2*a**4)
        shells[label] = dict(odd_action_times_density=corrections["-1"]+corrections["1"],
                             even_control_action_times_density=2*corrections["0"])
    if any(int(k)<2 for k in shells):
        raise ValueError("S3 has no physical transverse mode at k=1")
    odd = math.fsum(s["odd_action_times_density"] for s in shells.values()) / coherent
    even = math.fsum(s["even_control_action_times_density"] for s in shells.values()) / coherent
    largest = max(shells,key=lambda k:shells[k]["odd_action_times_density"])
    return dict(two_stage=even_budget(two_stage_transfer),number_of_distinct_shells=len(shells), one_surface_gauge_field=True,
                end_scale_factor=a, coherent_sheet_density=coherent, shells=shells,
                action_for_odd_selected_shells_equal_coherent=odd,
                action_for_odd_selected_shells_at_one_percent=100*odd,
                action_for_even_control_selected_shells_equal_same_reference=even,
                largest_selected_odd_shell=int(largest),
                single_largest_odd_shell_action_for_equal=(
                    shells[largest]["odd_action_times_density"]/coherent))

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--helicity",type=Path,default=HERE.parent/"helicity/HELICITY.json")
    p.add_argument("--even",type=Path,default=HERE.parent/"helicity/EVEN.json")
    p.add_argument("--handoff",type=Path,default=HERE.parent/"axion/HANDOFF.json")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=calculate(json.loads(args.helicity.read_bytes()), json.loads(args.handoff.read_bytes()),
                     json.loads(args.even.read_bytes()))
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
