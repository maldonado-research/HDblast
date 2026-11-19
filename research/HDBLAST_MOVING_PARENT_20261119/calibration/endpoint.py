"""The readout Hamiltonian of a free tagged mode on the stored homogeneous shells."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline

WALLS = "research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip"
WALL_HASH = "7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15"
REF = "research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/{}/quiet/rolling_summary.json"
REF_HASH = {"fine": "abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810",
            "coarse": "91601211fce4a617838d3a709849a138ec29a23998cc7c3ba2416ec0b0b8aa8c"}
MOMENTA = np.array([0., 2., 4., 6., 8.])
CUT = 2.1


def source(path, expected):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Historical input has changed: {path}")
    return data


def csv(archive, path):
    with archive.open(path) as f:
        return np.genfromtxt(f, delimiter=",", names=True)


def evolve(data, reference, step):
    t = data["H0tau"]
    if np.any(np.diff(t) <= 0.) or not t[0] < CUT < t[-1]:
        raise ValueError("Shell proper-time coverage")
    loga = CubicHermiteSpline(t, data["A_b"] - reference, data["H_over_H0"])
    mu = CubicHermiteSpline(t, 100*(data["phi_b"]-.5), 100*data["v_over_H0"])

    def mass_factor(s):
        a = np.exp(loga(s))
        return a**3, MOMENTA**2/a**2 + mu(s)**2

    volume, e2 = mass_factor(t[0])
    inverse_width = volume*np.sqrt(e2)
    x0 = 1/np.sqrt(2*inverse_width)
    state0 = np.r_[x0, -1j*inverse_width*x0]

    def rhs(s, state):
        x, p = np.split(state, 2)
        volume, e2 = mass_factor(s)
        return np.r_[p/volume, -volume*e2*x]

    sol = solve_ivp(rhs, (t[0], CUT), state0, method="DOP853", max_step=step,
                    t_eval=[CUT], rtol=2e-10, atol=2e-12)
    if not sol.success:
        raise ValueError(sol.message)
    x, p = np.split(sol.y[:, -1], 2)
    volume, e2 = mass_factor(CUT)
    diagonal = abs(p)**2/volume + volume*e2*abs(x)**2
    anomalous = p**2/volume + volume*e2*x**2
    scaled_error = np.max(abs((diagonal**2-abs(anomalous)**2)/e2-1))
    if scaled_error > 5e-9:
        raise ValueError("Canonical Hamiltonian identity failed")
    return np.array([np.sqrt(e2), diagonal, abs(anomalous)]).T, float(scaled_error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository", type=Path)
    repo = parser.parse_args().repository
    output = {}
    max_refine = max_id = 0.
    with zipfile.ZipFile(io.BytesIO(source(repo/WALLS,WALL_HASH))) as z:
        for grid, digest in REF_HASH.items():
            reference = json.loads(source(repo/REF.format(grid),digest))["matter"]["events"][0]["A_star"]
            pre = csv(z, f"{grid}/y1_b1/pre.csv")
            output[grid] = {}
            for arm in ("quiet", "plus", "minus", "far"):
                data = np.concatenate([pre, csv(z, f"{grid}/y1_b1/{arm}.csv")])
                first, efirst = evolve(data,reference,.004)
                refined, efine = evolve(data,reference,.002)
                max_id = max(max_id,efirst,efine)
                max_refine = max(max_refine,float(np.max(abs(first-refined))))
                output[grid][arm] = [dict(q=float(q), physical_quantum=energy,
                                          normal=diag, anomalous_abs=abs_b,
                                          normal_over_quantum=diag/energy,
                                          anomalous_fraction=abs_b/diag)
                                     for q,(energy,diag,abs_b) in zip(MOMENTA,refined)]
    if max_refine > 2e-5:
        raise ValueError("Endpoint evolution refinement")
    print(json.dumps(dict(cut_s=CUT, units="H0", by_grid=output,
                          max_step_halving_abs=max_refine,
                          max_hamiltonian_identity_scaled_error=max_id),
                     allow_nan=False, indent=2))


if __name__ == "__main__":
    main()
