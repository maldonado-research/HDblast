"""Authenticate the quantum-memory study or reproduce its sources, spectra and native walls offline."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
PRIOR = Path("research/HDBLAST_DISPERSION_MATCHING_20261119/gravity/GRAVITY.json")
BASE = "67b2e6ab85bb14d09f248415fdf31b7b66326c31"
PINS_SHA = "35217707f8b8a05d68ac61629ba581914784cf50bb5a9cb983c201e9dbe04b17"
CASES = ("fine_quiet", "fine_minus", "fine_far", "coarse_quiet", "fine_plus")
ORDERS = (0, 1, 2, 4, 8, 16, 32, 64)


def require(condition, description):
    if not condition:
        raise ValueError(description)


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def authenticate(repo):
    manifest = HERE/"SHA256SUMS"
    files = set()
    for line in manifest.read_text().splitlines():
        expected, name = line.split("  ", 1)
        path = (HERE/name).resolve()
        require(path.is_relative_to(HERE) and path != manifest and path not in files,
                f"Invalid or repeated capsule path {name}")
        require(digest(path) == expected, f"Checksum mismatch at {name}")
        files.add(path)
    inventory = {p.resolve() for p in HERE.rglob("*") if p.is_file()
                 and "__pycache__" not in p.parts and p != manifest}
    require(inventory == files, f"Capsule inventory mismatch {inventory ^ files}")
    pins = HERE/"SOURCES.json"
    require(digest(pins) == PINS_SHA, "Upstream manifest digest")
    sources = read(pins)
    require(sources["base"] == BASE and len(sources["sources"]) == 18, "Historical revision or inventory")
    paths = set()
    for row in sources["sources"]:
        path = (repo/row["path"]).resolve()
        require(path.is_relative_to(repo) and path not in paths, "Invalid historical path")
        require(path.stat().st_size == row["bytes"] and digest(path) == row["sha256"],
                f"Historical input disagrees at {row['path']}")
        paths.add(path)
    modes = read(HERE/"source/NATIVE_MODES.json")
    fork = read(repo/NATIVE/"FORK.json")
    require(digest(repo/NATIVE/"FORK.json") == modes["source_sha256"], "Archived modes reference")
    for grid in ("fine", "coarse"):
        require(all(modes[field] == fork[grid][field] for field in ("q","weight")), "Native momentum rule")
    print(f"Authenticated {len(files)} study files and {len(paths)} pinned historical files.", flush=True)


def targets(spectrum, packets, floor, gravity, previous):
    """Check physical identities and signed controls directly on the new outputs."""
    require(set(spectrum["arms"]) == set(packets["arms"]) == {"quiet", "minus", "plus", "far"},
            "Coverage of incoming histories")
    for arm, row in spectrum["arms"].items():
        quantum = row["retained_proper"]
        packet = packets["arms"][arm]
        require(abs(quantum["emission"]-packet["proper_comoving_energy"]) < 1e-8,
                f"Frequency-correlation trace and ordered emission for {arm}")
        minimum = floor["arms"][arm]
        fixed = minimum["preparations"]["retained"]["proper_band"]
        require(abs(quantum["emission"]-fixed) < 1e-8, f"Emission of the fixed parent state for {arm}")
        require(0 < minimum["analytic_proper_floor"] <= fixed and
                math.isclose(sum(minimum["per_mode_proper_floor"]), minimum["analytic_proper_floor"],
                             rel_tol=1e-10), f"Positive mode-summed energy floor for {arm}")
        budget = minimum["fair_budget"]
        require(minimum["analytic_proper_floor"] <= budget["primal_proper_band"] <= fixed and
                abs(budget["primal_proper_band"]-budget["dual_proper_bound"]) < 1e-10 and
                abs(budget["used_excess_over_vacuum"]-budget["allowed_excess_over_vacuum"]) < 1e-10,
                f"Optimal emission at the retained initial energy for {arm}")
        require(abs(quantum["emission"]-quantum["symmetrized"]-quantum["commutator"]) < 1e-10,
                f"Operator ordering for {arm}")
        require(abs(packet["dephased_energy"]-packet["proper_comoving_energy"]) < 1e-11,
                f"Energy preserved under dephasing for {arm}")
        for port in packet["ports"].values():
            require(0 <= port["witness"] <= port["trace_distance"] + 1e-12 <= .5+1e-11,
                    f"Distinguishability of equal-energy packets for {arm}")
    require(gravity["upstream_revision"] == BASE and gravity["native_runs"] == 24
            and gravity["exploratory_native_runs"] == 4, "Native gravity coverage")
    require(set(gravity["results"]) == {"fine_quiet", "coarse_quiet", "fine_minus"}, "Native wall histories")
    for case, steps in gravity["results"].items():
        for step, row in steps.items():
            endpoint = row["endpoint"]["h"]
            increment = endpoint["high_increment"]
            require(4.3e-7 < increment < 5.1e-7, f"Positive added-band expansion response on {case}")
            if case.endswith("quiet"):
                grid = case.split("_")[0]
                old = previous["grids"][grid]["steps"][step]["observed_60_minus_120"]
                require(abs(increment+old) < 4e-9, f"Independent static matching contrast {grid} {step}")
    for row in gravity["residual_control"]["runs"].values():
        require(row["h"] > 0 and row["h_from_trace"] > 0 and abs(row["h"]-row["h_from_trace"]) < 4e-11,
                "Sign and independent junction trace for the amplified small residual")


def tolerance(path):
    group, field = path[0], path[-1]
    if group == "gravity":
        if "residual_control" in path:
            if field == "relative_endpoint_to_high":
                return 0., 1e-4
            if field == "relative_history_to_high":
                return 0., 8e-6
            if field == "pair":
                return 2e-7, 1e-10
            return 0., 4e-11
        if field in ("relative_to_high_stride_two", "relative_to_full"):
            # These ratios normalize finite-difference noise. The numerator
            # and high-band norm are compared separately in physical units.
            if "h" in path:
                return 0., .0003
            if any(k in path for k in ("phi", "A", "v")):
                return 0., .00015
            return 2e-6, 5e-7
        if any(k in path for k in ("rho", "p", "E", "j")):
            return 1e-7, 2e-6
        if "pair" in path:
            return 1e-7, 5e-7
        if any(k in path for k in ("A", "phi")):
            return 1e-7, 2e-9
        return 2e-7, 2e-8
    if group == "quantum" and (field == "difference" or "refinement" in field
                               or "archive" in field or "difference" in field):
        return 2e-8, 2e-8
    if group in ("source", "sources", "csv"):
        return 3e-7, 1e-9
    if group == "theory" and "source_rule_comparison" in path:
        return 0., 4e-12
    if group == "theory" and field in ("rel_shortfall", "upper_bound", "analytic_upper_bound"):
        return 2e-8, 0.
    if group == "floor":
        if field == "quiet_optimum_relative_gap":
            return 2e-7, 1e-8
        if field in ("numeric_optimizer_max_relative_energy_error","max_relative_floor_change_on_joint_refinement"):
            return 2e-7, 1e-10
        if "per_mode_proper_floor" in path or field == "minimum_weighted_conformal_floor_per_mode":
            return 1e-7, 1e-22
    return 2e-8, 3e-12


def compare(new, old, path):
    if isinstance(old, dict):
        require(isinstance(new, dict) and new.keys() == old.keys(), f"Field inventory at {path}")
        for key in old:
            compare(new[key], old[key], (*path,key))
    elif isinstance(old, list):
        require(isinstance(new, list) and len(new) == len(old), f"Length at {path}")
        for index, (v,w) in enumerate(zip(new,old,strict=True)):
            compare(v,w,(*path,str(index)))
    elif isinstance(old, float):
        r, a = tolerance(path)
        require(math.isfinite(new) and math.isclose(new,old,rel_tol=r,abs_tol=a),
                f"Value at {'/'.join(path)}: computed {new}, expected {old}")
    else:
        require(type(new) is type(old) and new == old, f"Value at {path}")


def verify_rules(current):
    import numpy as np
    compare(read(current/"RULES.json"), read(HERE/"source/RULES.json"), ("rules",))
    with np.load(current/"RULES.npz",allow_pickle=False) as new, np.load(HERE/"source/RULES.npz",allow_pickle=False) as old:
        require(set(new.files) == set(old.files), "Frequency-rule array inventory")
        for field in old.files:
            if field in ("q","weight","orders"):
                np.testing.assert_array_equal(new[field], old[field])
            else:
                require(np.all(new[field] > 0), f"Positive quadrature array {field}")
                np.testing.assert_allclose(new[field],old[field],rtol=5e-10,atol=0,err_msg=field)


def verify_source_files(current, tables):
    checks = read(current/"CHECKS.json")
    require(len(checks["tables"]) == 80 and len(checks["comparisons"]) == 10, "Source coverage")
    compare(checks, read(HERE/"source/CHECKS.json"), ("source",))
    original = read(HERE/"source/SOURCES.json")
    rebuilt = read(current/"SOURCES.json")
    require(len(rebuilt["runs"]) == 40, "Source run count")
    for name, row in rebuilt["runs"].items():
        for scheme, entry in row["files"].items():
            published = original["runs"][name]["files"][scheme]
            new_hash, old_hash = entry.pop("sha256"), published.pop("sha256")
            require(digest(tables/entry["path"]) == new_hash and len(old_hash) == 64, f"Table digest {name} {scheme}")
    compare(rebuilt,original,("sources",))
    with (current/"COMPARISON.csv").open() as a, (HERE/"source/COMPARISON.csv").open() as b:
        new,old = csv.DictReader(a), csv.DictReader(b)
        require(new.fieldnames == old.fieldnames, "Source CSV fields")
        count = 0
        for left,right in zip(new,old,strict=True):
            require(all(left[k] == right[k] for k in ("case","scheme","bath","field")), "Source CSV identity")
            for k in ("high_normalized_linf","high_normalized_l1","readout_high_relative_error"):
                compare(float(left[k]),float(right[k]),("csv",k))
            count += 1
        require(count == 280, "Source CSV coverage")


def gravity_comparison(current, tables, work):
    original = read(HERE/"gravity/GRAVITY.json")
    archive = read(HERE/"source/SOURCES.json")["runs"]
    for name, value in current["source_sha256"].items():
        require(value == digest(tables/name), f"Exact source bytes for gravity {name}")
        base, scheme = Path(name).stem.rsplit("_",1)
        require(original["source_sha256"][name] == archive[base]["files"][scheme]["sha256"],
                f"Archived gravity input disagrees with archived sources at {name}")
    require(current["source_sha256"].keys() == original["source_sha256"].keys(), "Gravity source inventory")
    require(current["residual_control"]["auxiliary_source_sha256"] == digest(work/"sources/fine_quiet_120_difference.npz"),
            "Small residual source digest")
    for report in (current,original):
        report.pop("source_sha256")
        report["residual_control"].pop("auxiliary_source_sha256")
        for steps in report["results"].values():
            for row in steps.values():
                # These quotients are defined by separately published primitive
                # quantities. Reconstruct them exactly on each architecture.
                for d in row["endpoint"].values():
                    ratio = d.pop("absolute_error_over_high")
                    require(ratio == abs(d["four_minus_full"]/d["high_increment"]), "Endpoint relative error")
                for d in row["history_error"].values():
                    ratio = d.pop("relative_to_high")
                    require(ratio == d["error_norm"]/d["high_norm"], "Integrated relative error")
    compare(current, original, ("gravity",))


def replay(repo, jobs):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    def run(script,*args):
        print(f"Recomputing {script}",flush=True)
        subprocess.run([sys.executable,"-I","-B",str(HERE/script),*map(str,args)],env=env,check=True)
    with tempfile.TemporaryDirectory(prefix="hdblast-memory-") as tmp:
        root = Path(tmp)
        rules, tables, sources, wall = (root/n for n in ("rules","tables","sources","wall"))
        native = repo/NATIVE
        run("source/quadrature.py","--modes",HERE/"source/NATIVE_MODES.json","--output",rules)
        verify_rules(rules)
        reports = {}
        for script, name, arguments, group in (
                ("theory/check.py","theory/CHECK.json",["--rules",rules/"RULES.npz"],"theory"),
                ("quantum/spectrum.py","quantum/SPECTRUM.json",["--repo",repo],"quantum"),
                ("packets/packets.py","packets/PACKETS.json",["--repo",repo],"packets"),
                ("floor/floor.py","floor/FLOOR.json",["--repo",repo],"floor")):
            output = root/Path(name).name
            run(script,*arguments,"--output",output)
            reports[group] = read(output)
            compare(reports[group], read(HERE/name), (group,))
        run("source/source.py","--native",native,"--output",tables,"--cases",*CASES,
            "--high-orders",*ORDERS,"--jobs",jobs)
        run("source/verify.py","--native",native,"--tables",tables,"--output",sources)
        verify_source_files(sources,tables)
        output = root/"GRAVITY.json"
        run("gravity/gravity.py","--repo",repo,"--tables",tables,"--work",wall,"--jobs",jobs,"--output",output)
        gravity = read(output)
        targets(reports["quantum"],reports["packets"],reports["floor"],gravity,read(repo/PRIOR))
        gravity_comparison(gravity,tables,wall)
    print("Reproduced 80 source tables, 28 wall evolutions, ordered energies, frequency correlations and kernel bounds.",
          flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo",type=Path,default=HERE.parents[1],help="Root with the pinned historical inputs")
    ap.add_argument("--verify",action="store_true",help="Recalculate all numerical controls from the archived inputs")
    ap.add_argument("--jobs",type=int,default=2,help="Number of parallel source or native wall workers")
    options = ap.parse_args()
    require(options.jobs > 0, "Positive number of workers")
    repo = options.repo.resolve()
    authenticate(repo)
    targets(read(HERE/"quantum/SPECTRUM.json"),read(HERE/"packets/PACKETS.json"),
            read(HERE/"floor/FLOOR.json"),read(HERE/"gravity/GRAVITY.json"),read(repo/PRIOR))
    if options.verify:
        replay(repo,options.jobs)


if __name__ == "__main__":
    main()
