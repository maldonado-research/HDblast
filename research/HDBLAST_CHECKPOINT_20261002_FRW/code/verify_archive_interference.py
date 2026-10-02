"""Independent SciPy Ci check of the leading archive jump integral."""
import json
from pathlib import Path
import numpy as np
from scipy.special import sici
out=Path("outputs")
jumps=json.loads((out/"primary_hermite_jumps.json").read_text())
reference=json.loads((out/"archive_interference.json").read_text())
eta=np.array([j["eta"] for j in jumps]);d=np.array([j["generic_first_jump"] for j in jumps])
gap=np.diff(eta);assert np.all(gap>0)
lo,hi=np.array([1e3,1e6])/np.min(gap)
i,j=np.triu_indices(len(eta),1);sep=eta[j]-eta[i]
cross=float(np.sum(2*d[i]*d[j]*(sici(2*sep*hi)[1]-sici(2*sep*lo)[1])))
diagonal=float(np.dot(d,d)*np.log(hi/lo))
ratio=(diagonal+cross)/diagonal
phase_bound=float(2*(1/lo+1/hi)*np.sum(np.abs(d[i]*d[j])/sep)/diagonal)
assert phase_bound<1e-3
difference=abs(ratio-reference["coefficient_ratio"])
assert np.isfinite(ratio) and abs(ratio-1)<=1e-3
assert difference<=1e-10
result=dict(status="PASS",method="SciPy sici independent cross-check",K=[float(lo),float(hi)],
    coefficient_ratio=ratio,phase_independent_relative_bound=phase_bound,JavaScript_ratio_absolute_difference=difference,
    scope="Leading archive amplitudes only; no exact high-k mode evolution.")
(out/"archive_interference_scipy.json").write_text(json.dumps(result,indent=2)+"\n")
print("HDBLAST_FRW_ARCHIVE_INTERFERENCE_JSON="+json.dumps(result,separators=(",",":")))
