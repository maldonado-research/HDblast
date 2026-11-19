"""Cross-check the signed daughter response, pressure balance and analytic controls."""
import json
import math
from pathlib import Path
root = Path(__file__).resolve().parent
d = json.loads((root/"RESERVOIR.json").read_bytes())
s = json.loads((root/"SPLIT.json").read_bytes())
rows = []
for case, data in d["refine48/y1"].items():
    final = data["2.1"]
    ref = d["refine48/y1"]["quiet"]["2.1"]
    excess = 100*(final["daughter_physical_density"]/ref["daughter_physical_density"]-1)
    errors = []
    for n in ("1800","3600"):
        errors.append(max(abs(s[n][case][k]/final[k]-1) for k in s[n][case]))
    if not errors[1] < 9e-6 or errors[0]/errors[1] < 3:
        raise ValueError((case,errors))
    if not math.isclose(final["daughter_comoving_energy"]+final["chi_comoving_energy"],
                        final["source_work"]+final["chi_expansion_work"]+final["daughter_expansion_work"],
                        abs_tol=8e-6):
        raise ValueError("Energy balance")
    if any(x["a4_rho_daughter"] != 0 for x in d["fine/y0"][case].values()):
        raise ValueError("Stable field deposited radiation")
    rows.append(dict(case=case,physical_radiation=final["daughter_physical_density"],
                     change_percent=excess,split_relative_error=errors[1],
                     split_error_reduction=errors[0]/errors[1]))
print(json.dumps(rows, indent=2))
