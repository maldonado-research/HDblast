"""Compare a fresh independent RK4 replay with the fresh DOP853 calculation."""
import json
from pathlib import Path
import numpy as np

base=Path("outputs/replay")
read=lambda p: json.loads(Path(p).read_text())
r=read("outputs/actual_modes/actual_modes.json")
assert r["registered_controls_pass"] and not r["errors"]
for case in r["cases"]:
    for state in ("initial_state","out_state"):
        assert case[state]["ward_path_max_relative"] <= 1e-5
    assert case["paired_ward_path_relative"] <= 1e-5
p=r["cases"][0]
ind=read(base/"endpoint_basis_check.json")
four=read(base/"independent_rk4.json")
assert ind["status"] == four["status"] == read(base/"coherence_controls.json")["status"] == "PASS"
modes=ind["runs"][0]["rows"]
assert len(modes)==32
k=np.array([v["k"] for v in modes])
n=np.array([v["n"] for v in modes])
max_dn=float(np.max(np.abs(n-p["spectrum_n"])))
assert max_dn<1e-5
np.testing.assert_allclose(k,p["k"],rtol=1e-12,atol=1e-12)
# Derive the same endpoint independently from the archived knots and derivatives.
from scipy.interpolate import CubicHermiteSpline
g=np.array(read("data/independent_geometry.json"))
s=p["end"]
a=float(np.exp(CubicHermiteSpline(g[:,0],g[:,1],g[:,2])(s)))
phi=float(CubicHermiteSpline(g[:,0],g[:,3],g[:,4])(s))
F=np.array([sum(x*x for x in v["chi"]) for v in modes])
D=np.array([sum(x*x for x in v["chi_dot"]) for v in modes])
K=k*k/a**2
mass2=10000*(phi-.5)**2
om=np.sqrt(K+mass2)
F0=1/(2*a**3*om)
D0=om/(2*a**3)
w=np.array(p["weights"])
rho=float(w@(.5*(D-D0+(K+mass2)*(F-F0))))
pressure=float(w@(.5*(D-D0-(K/3+mass2)*(F-F0))))
current=float(w@(10000*(phi-.5)*(F-F0)))
actual=dict(rho=rho,p=pressure,J=current)
expected=dict(rho=p["paired_rho_final"],p=p["paired_p_final"],J=p["paired_current_final"])
errors={key:abs(actual[key]-expected[key]) for key in actual}
assert max(errors.values())<1e-5, errors
print(json.dumps(dict(status="PASS",max_occupation_absolute_difference=max_dn,
    endpoint_values=actual,endpoint_absolute_differences=errors,
    scope="independent 32-mode RK4 phase-sensitive endpoint stress/current"),indent=2))
