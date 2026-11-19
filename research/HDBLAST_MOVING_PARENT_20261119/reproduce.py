"""Verify the moving-parent capsule, or rebuild its numerical reports with --full."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
CASES = {
    "moving": "primary/MOVING.json",
    "checks": "primary/CHECK.json",
    "transfer": "independent/TRANSFER.json",
    "beam": "angular/BEAM.json",
    "flat": "calibration/FLAT.json",
    "endpoint": "calibration/ENDPOINT.json",
    "scale": "scales/POPULATION_SCALE.json",
    "coherence": "coherence/COHERENCE.json",
    "squeeze": "independent/SQUEEZE.json",
    "phase": "phase_feasibility/CERTIFICATES.json",
    "phase-refined": "phase_feasibility/REFINED.json",
    "coherence-refined": "coherence/REFINED.json",
}
SIGNS = {"plus": -1, "minus": 1, "far": 1}
WINDOWS = ("wide", "narrow")
BANDS = ("full", "screened")
TAGS = ("q0", "q2", "q4", "q6", "q8", "uniform", "tapered")
TRANSFER_KEYS = {
    "number": "number", "energy": "energy", "net_loss": "net_parent_loss",
    "annihilation_number": "ann_number", "stimulated_number": "stim_number",
    "annihilation_energy": "ann_energy", "stimulated_energy": "stim_energy",
}
HISTORICAL = {
    "research/HDBLAST_TAGGED_PARENT_20261119/primary/tagged.py",
    "research/HDBLAST_TAGGED_PARENT_20261119/primary/TAGGED.json",
    "research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip",
    *(f"research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/{grid}/quiet/rolling_summary.json"
      for grid in ("coarse", "fine")),
}
PCT = "difference_percent_of_incoherent_quiet"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(value, reference, label, rel=2e-8, absolute=1e-12):
    require(math.isfinite(value) and math.isfinite(reference)
            and math.isclose(value, reference, rel_tol=rel, abs_tol=absolute),
            f"{label}: found {value!r}, expected {reference!r}")


def authenticate(repository):
    manifest = json.loads((HERE / "SOURCES.json").read_bytes())
    covered = set()
    for row in manifest["sources"]:
        name = row["path"]
        path = (repository / name).resolve()
        require(not Path(name).is_absolute() and path.is_relative_to(repository)
                and name not in covered, f"Invalid historical source path {name}")
        payload = path.read_bytes()
        require(len(payload) == row["bytes"] and hashlib.sha256(payload).hexdigest() == row["sha256"],
                f"Historical input digest or size mismatch for {name}")
        covered.add(name)
    require(HISTORICAL <= covered, f"Historical inputs missing from manifest: {sorted(HISTORICAL-covered)}")
    return len(covered)


def primary_check(data, summary):
    close(data["cut_s"], 2.10, "proper readout")
    close(data["screen_pair_energy"], 60., "pair energy screen")
    close(data["upper_parent_momentum"], 8., "parent cutoff")
    counts = 0
    for grid in ("coarse", "fine"):
        arms = data["grids"][grid]
        for window in WINDOWS:
            quiet = arms["quiet"][window]
            for arm in ("quiet", *SIGNS):
                case = arms[arm][window]
                require(case["max_daughter_momentum_for_screen"] < 6.49e-6**(-1/3), "Daughter screen")
                require(case["max_wronskian_error"] < 1e-8 and case["max_parseval_error"] < 1e-8, "Mode identities")
                for tag in TAGS:
                    entry = case["spectrum"][tag]
                    close(entry["local_count_without_threshold"]+entry["positive_threshold_correction"],
                          entry["full"]["number"], "Timelike Parseval identity")
                    require(entry["positive_threshold_correction"] >= -1e-12, "Spacelike correction")
                    for band in BANDS:
                        v = entry[band]
                        close(v["ann_number"]-v["stim_number"], v["net_parent_loss"], "Signed parent loss")
                        for key in ("number", "energy"):
                            close(v["ann_"+key]+v["stim_"+key], v[key], "Inclusive "+key)
                            require(0 < entry["screened"][key] < entry["full"][key], "Screen ordering")
                        keys = ["number", "energy"]
                        if tag in ("uniform", "tapered"):
                            density = data["comoving_parent_density"][tag]/case["a_final"]**3
                            for key in ("number", "energy"):
                                close(v[key]*density, v[key+"_density"], key+" density")
                                keys.append(key+"_density")
                            close(3*v["isotropic_pressure"], v["energy_density"], "Isotropic daughter pressure")
                        if arm in SIGNS:
                            baseline = quiet["spectrum"][tag][band]
                            for key in keys:
                                require(SIGNS[arm]*(v[key]-baseline[key]) > 0,
                                        f"Incoherent sign {grid}/{arm}/{window}/{tag}/{band}/{key}")
                            counts += len(keys)
    require(counts == summary["pulse_sign_checks"] == 432, "Incoherent pulse-sign coverage")
    for name in ("parent_quadrature", "waveform_and_fourier_quadrature"):
        require(summary[name]["max_inclusive_contrast_change_percentage_points"] < 2e-5, name)
    require(summary["rest_relative_discrepancy"] < 5e-8, "Rest-mode limit")
    return counts


def rest_check(data, repository):
    rest_path = repository/"research/HDBLAST_TAGGED_PARENT_20261119/primary/TAGGED.json"
    rest = json.loads(rest_path.read_bytes())["grids"]
    worst = 0.
    for mesh, arms in data["grids"].items():
        for arm, windows in arms.items():
            for window, case in windows.items():
                for band in BANDS:
                    source = rest[mesh][arm][window][band]
                    computed = case["spectrum"]["q0"][band]
                    for key in ("number","energy","net_parent_loss"):
                        worst = max(worst,abs(computed[key]/source[key]-1))
    close(worst,data["max_relative_error_vs_published_rest"],"Published rest agreement",absolute=1e-12)
    require(worst<5e-8,"Published rest discrepancy")


def independent_check(reports):
    grids = reports["moving"]["grids"]
    transfer = reports["transfer"]["32768_128"]["arms"]
    for arm in ("quiet", "plus"):
        for window in WINDOWS:
            data = transfer[arm]["windows"][window]["data"]
            reference = grids["fine"][arm][window]["spectrum"]
            for band in BANDS:
                actual = data["screen" if band == "screened" else "full"]
                for tag in TAGS:
                    for key, label in TRANSFER_KEYS.items():
                        close(actual[tag][key], reference[tag][band][label],
                              f"Transfer {arm}/{window}/{band}/{tag}/{key}", rel=5e-7)
    for mesh in ("fine", "coarse"):
        for arm in ("quiet", *SIGNS):
            for window in WINDOWS:
                angular = reports["beam"]["results"][mesh][arm][window]
                case = grids[mesh][arm][window]
                for q in (0, 4, 8):
                    for band in BANDS:
                        v = angular[str(q)][band]
                        row = case["spectrum"][f"q{q}"][band]
                        close(v["energy"], row["energy"], "Beam energy", rel=3e-7)
                        close(v["longitudinal"]+2*v["transverse_each"], v["energy"], "Massless trace")
                        close(v["signed_flux"], q*row["net_parent_loss"]/case["a_final"], "Beam momentum", rel=3e-7)
                        if q == 0:
                            close(v["longitudinal_fraction"], 1/3, "Rest isotropy")
                        else:
                            require(v["longitudinal_fraction"] > 1/3, "Boosted pressure anisotropy")


def coherence_check(reports):
    primary = reports["coherence"]["results"]
    independent = reports["squeeze"]["refinements"]["128"]["contrasts_percent"]
    for mesh in ("fine", "coarse"):
        for arm, sign in SIGNS.items():
            for window in WINDOWS:
                rows = primary[mesh][arm][window]
                for label, row in rows.items():
                    field = label.removeprefix("screened_").removeprefix("physical_")
                    band = "screened" if label.startswith("screened_") else "full"
                    arm_rows = reports["moving"]["grids"][mesh]
                    for name, source in (("pulse",arm),("quiet","quiet")):
                        ref = arm_rows[source][window]["spectrum"]["uniform"][band][field]
                        close(row["incoherent"][name], ref, f"Coherence diagonal {label}", rel=2e-7)
                    low, high = row["minimizing_state"][PCT], row["maximizing_state"][PCT]
                    d = row["incoherent"][PCT]
                    require(low <= row["constant_phase_min_percent"] <= d <= row["constant_phase_max_percent"] <= high
                            and d*sign > 0, "Coherence extrema ordering")
                if mesh != "fine":
                    continue
                for band in BANDS:
                    direct = independent[arm][window]["screen" if band=="screened" else "full"]
                    for key in ("number", "energy"):
                        row = rows[("screened_" if band=="screened" else "")+key]
                        values = direct[key]
                        pairs = [(values["incoherent_percent"], row["incoherent"][PCT]),
                                 (values["squeezed_percent"], row["constant_phases"]["-0.5"][PCT]),
                                 (values["interval_percent"][0], row["constant_phase_min_percent"]),
                                 (values["interval_percent"][1], row["constant_phase_max_percent"])]
                        if band == "full":
                            density = rows["physical_"+key+"_density"]
                            pairs.extend([(values["local_squeezed_percent"],density["constant_phases"]["-0.5"][PCT]),
                                          (values["local_incoherent_percent"],density["incoherent"][PCT])])
                        for actual, reference in pairs:
                            close(actual, reference, f"Independent phase {arm}/{window}/{band}/{key}",
                                  rel=0., absolute=2e-5)
    # Count, energy and physical count differ in sign for the far control.
    for arm, signs in {"plus": (1,1,1), "minus": (-1,-1,-1), "far": (-1,1,1)}.items():
        for key, sign in zip(("number","energy","physical_number_density"), signs):
            value = primary["fine"][arm]["wide"][key]["constant_phases"]["-0.5"][PCT]
            require(sign*value > 0, f"Fixed-phase sign {arm}/{key}")


def quadratic_phase_check(report):
    coarse, refined = (report["refinements"][n] for n in ("64","128"))
    pairings = {
        "excess":"interval_all_phases_percent", "absolute":"absolute_interval_all_phases_percent",
        "local_excess":"local_interval_all_phases_percent",
        "local_absolute":"local_absolute_interval_all_phases_percent",
    }
    for arm in ("quiet", *SIGNS):
        for window in WINDOWS:
            for band in ("full","screen"):
                for key in ("number","energy"):
                    mode = refined["arms"][arm]["windows"][window][band][key]
                    require(mode["vacuum"]>=0 and mode["profile_squeezed"]>=0, "Absolute pair probability")
                    if arm=="quiet":
                        continue
                    row = refined["contrasts_percent"][arm][window][band][key]
                    other = coarse["contrasts_percent"][arm][window][band][key]
                    quiet = refined["arms"]["quiet"]["windows"][window][band][key]
                    ratio = (refined["arms"]["quiet"]["a_cut"]/refined["arms"][arm]["a_cut"])**3
                    norm = quiet["incoherent_excess"]
                    for profile, envelope in pairings.items():
                        field = "quadratic_phase_"+profile+"_percent"
                        sign = -SIGNS[arm]
                        require(sign*row[field]>0, f"Quadratic-phase sign {arm}/{window}/{band}/{key}/{profile}")
                        low, high = row[envelope]
                        require(low<=row[field]<=high, "Attainable phase interval "+field)
                        require(sign*(low+high)<0, "Number-diagonal midpoint sign "+field)
                        close(row[field],other[field],"Phase refinement "+field,rel=0.,absolute=2e-7)
                        def observable(state):
                            return state["profile_squeezed"] - (state["vacuum"] if "excess" in profile else 0.)
                        dilution = ratio if profile.startswith("local_") else 1.
                        value = 100*(dilution*observable(mode)-observable(quiet))/norm
                        close(value,row[field],"Direct amplitude "+field,rel=1e-8)
    quiet = refined["arms"]["quiet"]["windows"]["wide"]["full"]["number"]
    require(0<quiet["profile_squeezed"]<quiet["vacuum"], "Positive count and negative quiet excess")


def feasibility_check(report, coherence):
    labels = report["labels"]
    expected = {f"{mesh}/{win}/{band}/{arm}" for mesh in ("fine","coarse") for win in WINDOWS
                for band in ("number","screened_number") for arm in SIGNS}
    require(len(labels)==len(expected) and set(labels)==expected, "Phase constraints")
    require(len(report["incoherent"])==len(labels), "Phase normalizations")
    for label, value in zip(labels,report["incoherent"]):
        mesh, win, band, arm = label.split("/")
        reference = coherence["results"][mesh][arm][win][band]["incoherent"]["difference"]
        close(value, reference, "Diagonal normalization "+label, rel=2e-5)
    for group in report["groups"].values():
        rows = group["rows"]
        for kind in ("free", "common"):
            cert = group[kind]
            require(len(rows)==len(cert["dual"])==len(cert["margin_factors"]), "Certificate rows")
            require(all(x>=-1e-12 for x in cert["dual"]), "Negative dual weight")
            close(sum(cert["dual"]), 1., "Dual normalization")
            require(len(cert["z_real"])==len(cert["z_imag"])==(report["nodes"] if kind=="free" else 1),
                    "Witness dimensions")
            require(all(math.hypot(x,y)<=1+1e-9 for x,y in zip(cert["z_real"],cert["z_imag"])),
                    "Witness separability")
            close(min(cert["margin_factors"]), cert["lower"], "Primal bound")
            close(cert["lower"], cert["upper"], "Primal and dual closure", rel=0., absolute=2e-8)
    group=report["groups"]["all"]
    require(group["rows"]==list(range(len(labels))), "All 24 constraints")
    require(group["free"]["lower"]>1 and group["common"]["upper"]<1, "Common versus variable phase")
    smooth = report["smooth"]
    close(min(smooth["margin_factors"]),smooth["minimum"],"Explicit minimum margin")
    require(smooth["minimum"]>1 and len(smooth["margin_factors"])==len(labels)==len(smooth["contrasts"]),
            "Explicit smooth witness")
    for d,g,v in zip(report["incoherent"],smooth["margin_factors"],smooth["contrasts"]):
        require(d*v<0, "Smooth separable contrast sign")
        close(d*(1-g),v,"Smooth separable contrast")
    cap = report["common_angle_variable_modulus"]
    require(0<=cap["row"]<len(labels) and len(cap["row_upper_bounds"])==len(labels), "Shared-angle rows")
    close(min(cap["row_upper_bounds"]),cap["upper"],"Shared-angle ceiling")
    require(cap["upper"]<1, "Variable modulus phase bound")


def match_certificates(actual, expected):
    # Different primal or dual witnesses can attain the same finite-grid bound.
    for key in ("nodes","samples","padding","labels","momentum","incoherent"):
        match(actual[key],expected[key],"certificate/"+key)
    require(actual["groups"].keys()==expected["groups"].keys(),"Phase certificate schema")
    for name, group in actual["groups"].items():
        require(group["rows"]==expected["groups"][name]["rows"],"Certificate subset")
        for kind in ("free","common"):
            for edge in ("lower","upper"):
                close(group[kind][edge],expected["groups"][name][kind][edge],"Certificate "+edge,
                      rel=0.,absolute=5e-5)
    for section,key in (("smooth","minimum"),("common_angle_variable_modulus","upper")):
        close(actual[section][key],expected[section][key],section,rel=0.,absolute=5e-5)


def science(reports):
    count = primary_check(reports["moving"], reports["checks"])
    independent_check(reports)
    coherence_check(reports)
    coherence_check({**reports, "coherence":reports["coherence-refined"]})
    quadratic_phase_check(reports["squeeze"])
    feasibility_check(reports["phase"], reports["coherence"])
    feasibility_check(reports["phase-refined"], reports["coherence-refined"])
    flat = reports["flat"]["rows"]
    require(len(flat)==16, "Flat boost and duration coverage")
    for row in flat:
        work = row["net_switching_work"]
        require(work > 0, "Flat switching work")
        close(work, row["daughter_energy"]-row["signed_parent_depletion"], "Flat energy balance")
        close(row["pair_longitudinal_momentum_over_EL"], -row["equal_opposite_parent_transfer_over_EL"],
              "Flat momentum balance")
    require(reports["endpoint"]["max_hamiltonian_identity_scaled_error"] < 1e-8, "Endpoint Hamiltonian")
    scale = reports["scale"]
    close(scale["number"]*scale["scale_factor"]**3, scale["comoving_number"], "Comoving density")
    shifts = scale["fractional_shifts"]
    close(shifts["matter"]+shifts["force"], shifts["sum"], "Initial source corrections")
    require(shifts["force"] < shifts["sum"] < 0 < shifts["matter"], "Initial source signs")
    return count


def match(actual, expected, path):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), "Dictionary schema at "+path)
        for key in expected:
            match(actual[key], expected[key], path+"/"+key)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual)==len(expected), "List length at "+path)
        for index, value in enumerate(expected):
            match(actual[index], value, f"{path}/{index}")
    elif isinstance(expected, float):
        # Small differences of separately solved arms, residuals of numerical identities
        # and the phase of a tiny complex contrast lose precision before division.
        if path.startswith(("coherence/", "coherence-refined/", "squeeze/", "checks/")):
            close(actual, expected, path, rel=3e-5, absolute=3e-7)
        elif path.startswith(("endpoint/","transfer/")) or "error" in path:
            close(actual, expected, path, rel=3e-6, absolute=1e-8)
        else:
            close(actual, expected, path, rel=5e-7, absolute=3e-10)
    else:
        require(type(actual) is type(expected) and actual==expected, "Value at "+path)


def run(script, args, directory, name):
    output = directory/(name+".json")
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    print(f"Rebuilding {name}", file=sys.stderr, flush=True)
    with output.open("wb") as file:
        subprocess.run([sys.executable, "-I", "-B", str(HERE/script), *map(str,args)],
                       stdout=file, check=True, env=env)
    return output


def replay(repository, saved):
    # Temporary results feed the refinement analysis. The archived reports remain immutable.
    with tempfile.TemporaryDirectory() as temp:
        directory = Path(temp)
        files = {}
        for name, opts in (("moving", []), ("q-refined", ["--nodes",32]),
                           ("refined", ["--nodes",32,"--samples",2048,"--padding",64])):
            files[name] = run("primary/moving.py", [repository,*opts], directory, name)
        files["checks"] = run("primary/checks.py", [files[x] for x in ("moving","q-refined","refined")],
                              directory,"checks")
        runners = {"transfer":"independent/transfer.py", "beam":"angular/beam_stress.py",
                   "flat":"calibration/flat_moving.py", "endpoint":"calibration/endpoint.py",
                   "scale":"scales/population_scale.py", "coherence":"coherence/coherence.py",
                   "squeeze":"independent/squeeze.py"}
        for name, script in runners.items():
            files[name] = run(script, [] if name=="flat" else [repository], directory, name)
        files["coherence-refined"] = run("coherence/coherence.py",
            [repository,"--nodes",40,"--samples",2048,"--padding",32], directory,"coherence-refined")
        for name, opts in (("phase",[]),("phase-refined",["--nodes",48,"--samples",2048,"--padding",32])):
            files[name] = run("phase_feasibility/feasibility.py", [repository,*opts], directory, name)
        actual = {name:json.loads(files[name].read_bytes()) for name in saved}
        for name, reference in saved.items():
            if name in ("phase", "phase-refined"):
                match_certificates(actual[name], reference)
            else:
                match(actual[name], reference, name)
        rest_check(actual["moving"],repository)
        return science(actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=HERE.parent.parent,
                        help="HDblast checkout containing the pinned historical sources")
    parser.add_argument("--full", action="store_true", help="rerun the numerical drivers (NumPy and SciPy)")
    args = parser.parse_args()
    repository = args.repository.resolve()
    pins = authenticate(repository)
    saved = {name:json.loads((HERE/path).read_bytes()) for name, path in CASES.items()}
    rest_check(saved["moving"], repository)
    count = science(saved)
    if args.full:
        replay(repository, saved)
    print(f"Verified {pins} historical sources, {count} incoherent signs, independent spectra, "
          f"phase contrasts, smooth-state certificates, beam stress and energy checks. Full numerical replay: {args.full}.")


if __name__ == "__main__":
    main()
