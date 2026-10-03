"""Separately registered K384 follow-up; preserves the original matrix failure."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

FROZEN_CORE_SHA256="32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a"
BASELINE_HASHES={
    "summary.json":"f6cc20865e61d8134c63ecc8a9fa6bf1695d9ffafd0eea0890ca9aeb471e0228",
    "A0.2_tight_K192/arrays.npz":"546c09231f504da70f3a28cb598fe40559d56d26af4b0e7b01a998926fc3ed0d",
    "A0.2_quadrature_K192/arrays.npz":"62cc23be671106b228f3bb2cb99870b4d15df41d1e00676ed7e35c9c5cce0994",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core(path):
    if sha256(path)!=FROZEN_CORE_SHA256:
        raise ValueError("Frozen repaired core hash mismatch")
    spec=importlib.util.spec_from_file_location("smooth_frw_frozen_core",path)
    core=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=core
    spec.loader.exec_module(core)
    return core


def old_curves(base,label,new_t):
    path=base/label/"arrays.npz"
    with np.load(path,allow_pickle=False) as archive:
        if not np.array_equal(archive["t"],new_t):
            raise ValueError("Baseline and new output times differ")
        result={name:archive[f"K192_{name}"].copy() for name in ("rho","p","Q")}
    if not all(np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError("Nonfinite preserved baseline")
    return result


def execute(args):
    if sha256(args.protocol)!=args.protocol_sha256:
        raise ValueError("Follow-up protocol hash mismatch")
    core=load_core(args.core)
    for relative,expected in BASELINE_HASHES.items():
        if sha256(args.baseline/relative)!=expected:
            raise ValueError(f"Frozen baseline artifact mismatch: {relative}")
    baseline=json.loads((args.baseline/"summary.json").read_text())
    if baseline["status"]!="FAIL" or baseline["amplitudes"]["0.2"]["status"]!="FAIL":
        raise ValueError("Expected preserved failed original curved matrix")
    if baseline["provenance"]["source_sha256"]!=FROZEN_CORE_SHA256:
        raise ValueError("Baseline core source mismatch")
    args.output.mkdir(parents=True,exist_ok=False)
    provenance=dict(created=datetime.now(timezone.utc).isoformat(),
        followup_protocol_sha256=args.protocol_sha256,wrapper_sha256=sha256(Path(__file__)),
        core_sha256=FROZEN_CORE_SHA256,baseline_artifact_sha256=BASELINE_HASHES,
        original_matrix_status=baseline["status"],python=sys.version,numpy=np.__version__,
        scope="Separate prospective K384 prescribed curved-background follow-up; original matrix remains FAIL.")
    core.finite_json(args.output/"provenance.json",provenance)
    summary=dict(status="IN_PROGRESS",provenance=provenance)
    core.finite_json(args.output/"summary.json",summary)
    # Only the bookkeeping set changes inside this process. The core file,
    # subtraction, physical model, integration settings and gates stay frozen.
    core.CUTOFFS=(24.0,48.0,96.0,192.0,384.0)
    tight=core.run_case(args.output,.2,core.TIGHT,cutoff=384.)
    quad=core.run_case(args.output,.2,core.QUADRATURE,cutoff=384.)
    half=core.run_case(args.output,.2,core.TAIL_STEP,cutoff=384.)
    static=core.run_case(args.output,0.,core.TAIL_STEP,cutoff=384.,static=True)
    cutoff=core.compare_observables(tight["integrated"]["192"],tight["integrated"]["384"])
    solver=core.compare_observables(tight["integrated"]["384"],half["integrated"]["384"])
    quadrature=core.compare_observables(tight["integrated"]["384"],quad["integrated"]["384"])
    continuity={
        "tight":core.compare_observables(old_curves(args.baseline,"A0.2_tight_K192",tight["t"]),tight["integrated"]["192"]),
        "quadrature":core.compare_observables(old_curves(args.baseline,"A0.2_quadrature_K192",quad["t"]),quad["integrated"]["192"]),
    }
    chosen=quad["diagnostics"]["cutoffs"]["384"]
    t=quad["t"]
    values=quad["integrated"]["384"]
    fine,fine_arrays=core.exchange_diagnostics(t,values,quad["a"],quad["H"],quad["xp"])
    coarse,coarse_arrays=core.exchange_diagnostics(t[::2],{key:value[::2] for key,value in values.items()},
                                                  quad["a"][::2],quad["H"][::2],quad["xp"][::2])
    time=core.difference_gate(fine_arrays["ledger"][::2],coarse_arrays["ledger"],atol=1e-5,rtol=.005)
    static_max=static["diagnostics"]["cutoffs"]["384"]["maxima"]
    static_pass=all(value<=2e-5 for value in static_max.values()) and static["diagnostics"]["wronskian_passed"]
    gates={"cutoff":core.all_passed(cutoff),"solver":core.all_passed(solver),
        "quadrature":core.all_passed(quadrature),"subset_continuity":all(core.all_passed(value) for value in continuity.values()),
        "time":time["passed"],"wronskian":all(run["diagnostics"]["wronskian_passed"] for run in (tight,quad,half,static)),
        "integrated_exchange":chosen["exchange_fine"]["ledger_passed"],
        "local_exchange":chosen["exchange_fine"]["local_passed"],
        "independent_ODE_work":chosen["ode_bare_work_passed"],
        "direct_trace":chosen["trace"]["direct_derivatives"]["passed"],
        "sampled_trace":chosen["trace"]["sampled_derivatives"][0]["passed"],
        "future_spectral_energy":chosen["future_spectral_energy"]["passed"],
        "static_negative_control":static_pass}
    summary.update(status="PASS" if all(gates.values()) else "FAIL",gates=gates,
        selected=quad["diagnostics"]["label"],amplitude=.2,cutoff=384,
        cutoff_K192_to_K384=cutoff,solver_phase_refinement=solver,quadrature_refinement=quadrature,
        subset_K192_continuity=continuity,time_refinement=time,
        exchange_fine=fine,exchange_coarse=coarse,selected_diagnostics=chosen,
        static_negative_control=dict(passed=static_pass,maxima=static_max,absolute_limit=2e-5),
        original_matrix_status="FAIL",original_curved_pressure_cutoff_gate=
            baseline["amplitudes"]["0.2"]["final_cutoff_gate"]["p"])
    core.finite_json(args.output/"summary.json",summary)
    print("HDBLAST_SMOOTH_FRW_FOLLOWUP_JSON="+json.dumps(summary,allow_nan=False),flush=True)
    return 0 if summary["status"]=="PASS" else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute",action="store_true",required=True)
    parser.add_argument("--core",type=Path,required=True)
    parser.add_argument("--baseline",type=Path,required=True)
    parser.add_argument("--protocol",type=Path,required=True)
    parser.add_argument("--protocol-sha256",required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    try:
        return execute(args)
    except Exception as exc:
        if args.output.exists():
            (args.output/"failure.json").write_text(json.dumps(dict(exception=type(exc).__name__,message=str(exc)),indent=2,allow_nan=False)+"\n")
        raise


if __name__=="__main__":
    raise SystemExit(main())
