"""Read native RK stages during a quiet replay and test time-centering of the trace constraint."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parent


def replay(bulk, tables, grid, cache):
    spec = importlib.util.spec_from_file_location("live_stages", bulk / "live.py")
    live = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(live)
    native = live.evolve_b1.Run.rhs
    rows = []

    def rhs(run, t, y, matter):
        k, dm, f = native(run, t, y, matter)
        A, B, phi, AT, BT, phiT = y[:, -1]
        AZ, BZ, phiZ = (f[key][-1] for key in ("A_Z", "B_Z", "phi_Z"))
        AZZ = f["A_ZZ"][-1]
        na, nb, nf = (f[key] for key in ("gA", "gB", "gF"))
        e2b = math.exp(2*B)
        w = 1-phi+phi**3/3
        potential = .5*(phi**2-1)**2-(2/3)*w**2
        ham = (-2*e2b*potential + 6*AT*(AT+BT) - 12*AZ**2 + 6*AZ*BZ
               - 6*AZZ - phiT**2 - phiZ**2)
        fac = run.rb * math.exp(2*(A-run.matter.ref) - B)
        # Rates are d(a² h)/dT. The first uses the actual RK acceleration.
        actual = fac * (k[3, -1] + AT*(2*AT-BT))
        junction = fac*(na*(na+nb) + e2b*potential/3 - (nf**2+phiT**2)/6)
        gradient = fac*((AZ-na)*(AZ+na+BZ) + na*(BZ-nb) - (phiZ-nf)*(phiZ+nf)/6)
        constraint = -fac*ham/6
        rows.append([t, actual, junction, gradient, constraint])
        return k, dm, f

    with tempfile.TemporaryDirectory() as scratch, patch.object(live.evolve_b1.Run, "rhs", rhs):
        folder = Path(scratch)
        live.run_case(grid, "quiet", 0., 120, "static", tables, folder, cache)
        with np.load(folder / f"{grid}_quiet_zero.npz", allow_pickle=False) as d:
            end = float(d["s"][-1])
            size = len(d["T"])
            step = float(d["dt"])
            a = np.exp(d["A"]-float(d["Aref"]))
            net = float((a**2*d["h"])[-1]-(a**2*d["h"])[0])
            recorded_clock = d["T"]
    rows = np.array(rows)
    if len(rows) != 4*(size-1)+1:
        raise RuntimeError("Unexpected native RK stage count")
    stages = rows[:-1].reshape(-1, 4, 5)
    np.testing.assert_allclose(stages[:, :, 0], recorded_clock[:-1,None]+
                               step*np.array([0., .5, .5, 1.]), rtol=0, atol=2e-13)
    weights = np.array([1.,2.,2.,1.]) / 6
    four = step * np.einsum("nsk,s->k", stages[:, :, 1:], weights)
    keys = ("rate", "junction", "stencil", "hamiltonian")
    point = max(abs(rows[:,1]-rows[:,2]-rows[:,3]-rows[:,4]))
    stage_one = np.trapezoid(rows[::4,4], x=recorded_clock)
    if point > 1e-10 or abs(four[0]-net) > 1e-8:
        raise ValueError("Stage trace identity failed")
    return dict(endpoint_proper_time=end, observed_a2h_change=net, full_step_count=size-1,
                all_stage_integrals=dict(zip(keys, map(float, four), strict=True)),
                stage_one_only_hamiltonian_integral=float(stage_one),
                maximum_pointwise_algebra_error=float(point),
                four_stage_chain_rule_remainder=float(net-four[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bulk", type=Path, default=HERE.parent / "bulk")
    ap.add_argument("--output", type=Path, default=HERE / "STAGES.json")
    ap.add_argument("--cache", type=Path, help="Optional trusted native geometry cache used by the bulk replay")
    args = ap.parse_args()
    with tempfile.TemporaryDirectory() as scratch:
        cache = args.cache if args.cache is not None else Path(scratch)
        result = {g: replay(args.bulk, HERE.parent / "source/tables", g, cache)
                  for g in ("coarse", "fine")}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
