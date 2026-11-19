"""Check the archive and reproduce the scalar states, photon dynamics and work balances."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
BASE = "579f71405079f9048704e2bcd51cde6a6865a24e"
PINS_SHA = "5a50fb005cc3016dc090b0e4792c32293f35d5257a0e2975fc03d82a25fb02aa"
PRODUCTS = {"round":"round/SLE.json", "patch":"patch/PATCH.json", "checks":"checks/STATE_CHECK.json",
            "radiation":"radiation/PHOTONS.json", "maxwell":"maxwell/MAXWELL.json"}
SCRIPTS = {"round":"round/sle.py", "patch":"patch/patch.py", "checks":"checks/state_check.py",
           "radiation":"radiation/photons.py", "maxwell":"maxwell/maxwell.py"}


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_files():
    manifest = HERE/'SHA256SUMS'
    seen = set()
    for line in manifest.read_text().splitlines():
        checksum, name = line.split('  ', 1)
        path = HERE/name
        require(not Path(name).is_absolute() and Path(name).as_posix() == name
                and path == path.resolve() and path.is_relative_to(HERE)
                and path != manifest and path not in seen, 'Invalid inventory path')
        require(digest(path) == checksum, f'Hash differs for {name}')
        seen.add(path)
    actual = {p for p in HERE.rglob('*') if p.is_file() and p != manifest
              and not (p.parent.name == '__pycache__' and p.suffix == '.pyc')}
    require(actual == seen, 'Archive coverage differs')
    pin_file = HERE/'SOURCES.json'
    require(digest(pin_file) == PINS_SHA, 'Unexpected source manifest')
    pins = read(pin_file)
    require(pins['base'] == BASE and pins['previous_doi'] == '10.5281/zenodo.24380589',
            'Historical revision differs')
    prior = set()
    for record in pins['sources']:
        path = HERE/record['bundled']
        require(path in seen and path not in prior and path.stat().st_size == record['bytes']
                and digest(path) == record['sha256'], 'Historical input mismatch')
        prior.add(path)
    require(prior == {p for p in seen if p.is_relative_to(HERE/'prior')}, 'Historical coverage differs')
    print(f'Verified {len(seen)} files and {len(prior)} historical objects.', flush=True)


def numerical_floor(path):
    """Absolute replay tolerances in each product's reported units."""
    name, field = path[0], path[-1]
    if name == "checks" and field == "optimizer_squeeze_error":
        # The residual of the independent bounded scalar search varies with BLAS.
        return 2e-8
    if name == "patch":
        return 1e-9 if field == "full_matrix_balance_error" else 2e-11
    if name == "round":
        if "density_difference_over_parent" in path:
            return 5e-14
        if field == "pressure_difference_over_H0_4":
            return 2e-6
        if field in {"rho_difference_over_H0_4","matter_measured_density_difference_over_H0_4"}:
            return 5e-9
        if "work" in path:
            return 5e-8
    if name == "maxwell" and ("low_occupations" in path or field in {"sector_occupation","excitation_energy"}):
        return 2e-14
    return 8e-11


def compare(actual, expected, path=()):
    label = '/'.join(path)
    if isinstance(expected, dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(), f'Different keys at {label}')
        for key in expected:
            compare(actual[key],expected[key],(*path,key))
    elif isinstance(expected, list):
        require(isinstance(actual,list) and len(actual)==len(expected), f'Different size at {label}')
        for i, (a,e) in enumerate(zip(actual,expected)):
            compare(a,e,(*path,str(i)))
    elif isinstance(expected, float):
        floor = numerical_floor(path)
        require(type(actual) in (int,float) and math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual,expected,rel_tol=2e-7,abs_tol=floor),
                f'{label}: expected {expected!r}, obtained {actual!r}')
    else:
        require(type(actual) is type(expected) and actual==expected, f'Scalar differs at {label}')


def verify_round_states(rolling, exact):
    widths = rolling["half_widths"]
    require(widths == [r["hT"] for r in exact["runs"]] == [.1,.2,.3], "Preparation windows differ")
    control = rolling["closed_de_sitter_check"]
    # The exact u-mode reference has q variance (n+1)/sqrt(n(n+2))
    # times instantaneous q variance. A real Bogoliubov beta adds 2 asinh(beta)
    # to its log variance. Compare two independently parameterized state solvers.
    for sampled in exact["runs"]:
        values = control["log_widths"][f"{sampled['hT']:.2f}"]["from_integrator"]
        require(len(values) == len(control["shells"]), "Missing exact reference mode")
        for n,u in zip(control["shells"],values):
            beta = sampled["samples_by_n"][str(n)]["beta"]
            expected = math.log((n+1)/math.sqrt(n*(n+2))) + 2*math.asinh(beta)
            require(abs(u-expected)<2e-11, "Closed de Sitter scalar preparations disagree")
        require(sampled["sampled_density_difference"] < 0, "Sampled state missed its energy minimum")
    for branch in rolling["branches"].values():
        reference = f"{rolling['reference_half_width']:.2f}"
        for name, coupling in branch["couplings"].items():
            for width, cutoffs in coupling["relative_to_width_0p20"].items():
                for row in cutoffs.values():
                    densities = row["sampled_density"]
                    require(densities[width] <= 5e-11 and densities[reference] >= -5e-11,
                            "A round trial failed its minimization inequalities")
                    for work in row["work"].values():
                        direct, metric = work["coupling_work_over_H0"], work["geometry_work_over_H0"]
                        error = work["comoving_energy_change_over_H0"] - direct - metric
                        require(abs(error)<1e-8, "Round scalar work failed to balance")
                        if name == "constant":
                            require(direct==0, "Constant coupling did direct work")


def verify_uneven_geometry(patch, maxwell):
    passage = read(HERE/"prior/prior/passage/PASSAGE.json")
    require(patch["angular_ell"] == maxwell["ell"] == 2
            and patch["widths"] == maxwell["windows"], "Uneven sector conventions differ")
    for seed, run in patch["runs"].items():
        em = maxwell["runs"][seed]
        original = passage["runs"][seed]
        require(abs(em["times"]["crossing"]-original["crossing_x"])<2e-8
                and abs(run["crossing"]-original["crossing_x"])<2e-8,
                "An uneven background missed the archived first passage")
        for record in run["radial_cutoffs"].values():
            for width, delta in record["differences_from_width_0p2"].items():
                energy = delta["energy"]
                work = delta["coupling_work"]+delta["metric_work"]
                require(abs(energy["after"]-energy["before"]-work)<5e-10,
                        "Uneven scalar energy difference failed to balance")
                require(min(delta["objective_penalties"].values())>=-1e-10,
                        "Uneven scalar trial failed the two sampling inequalities")
            require(record["full_matrix_balance_error"]<3e-9, "Uneven full-matrix balance failed")
        require(all(row["coupling_work"]==0 for row in
                    run["constant_Z_at_48"]["differences_from_width_0p2"].values()),
                "Uneven constant-coupling control failed")
        require(em["round_control_max_sector_occupation"]<1e-21, "Round Maxwell null failed")
        bound = em["bounce_to_first_crossing"]
        left, right = bound["classical_maxwell_energy_ratio"]
        e = em["primary"]["states"]["0.2"]
        ratio = e["crossing"]["added_first_launch_mode_energy"]/e["bounce"]["added_first_launch_mode_energy"]
        require(0 < left < ratio < right and
                abs(ratio-bound["window_0p2_added_mode_energy_ratio"])<1e-8,
                "A positive added Maxwell mode violates the redshift inequality")
        require(all(r["weyl_squared_spatial_average_over_h4"]>=0 for r in em["weyl"].values()),
                "Electric Weyl invariant was negative")
    for row in patch["round_closed_de_Sitter_control_n16"].values():
        require(max(row.values())<5e-10,"Matrix vacuum missed the round exact-mode comparison")


def verify_photons(thermal, rolling):
    historical = read(HERE/"prior/homogeneous/results.json")["branches"]
    for key, branch in thermal["branches"].items():
        old = historical[key]["runs"]["inverse_f_seed_0"]
        zero = branch["runs"]["0"]
        require(old["stop"]=="first_maximum","Historical zero-radiation maximum absent")
        for stage, prior_stage in [("first_crossing","first_minimum_crossing"),("first_maximum","at_stop")]:
            for field in ("x","a","q","h"):
                require(abs(zero[stage][field]-old[prior_stage][field])<2e-8,
                        "Zero-radiation historical event disagrees")
        require(abs(zero["first_crossing"]["x"]-rolling["branches"][key]["archived_crossing_x"])<2e-8,
                "Round state and thermal control used different first crossings")
        require(abs(4*branch["ui"]*branch["critical_d"]-1)<2e-12,"Bounce discriminant differs")
        require(abs(branch["thermal_states"][-1]["launch_fraction"]-.5)<2e-7,
                "Double-root radiation fraction differs")
        critical = branch["critical_check"]
        force = critical["force_at_release"]
        for endpoint in critical["critical"]:
            leading = -(force**2)*endpoint["x"]**4/24
            ratio = endpoint["log_radius"]/leading
            require(ratio>0 and abs(ratio-endpoint["quartic_coefficient_ratio"])<2e-10,
                    "Critical quartic contraction differs from its scaled endpoint")
        tightest = critical["subcritical"][-1]
        delta = math.sqrt(tightest["gap"])
        alpha = 2*branch["ui"]*delta/(1+delta)
        time_scale = math.sqrt(6*alpha)/abs(force)
        gain_scale = 1.5*(alpha/force)**2
        require(abs(tightest["time_of_maximum"]/time_scale-1)<1e-3 and
                abs(tightest["log_radius_gain"]/gain_scale-1)<1e-3,
                "Near-critical turnaround differs from its analytic scaling")
        for label, run in branch["runs"].items():
            d = float(label)*branch["critical_d"]
            ai = run["initial"]["a"]
            initial_fraction = d/ai**2
            require(0<=initial_fraction<.5,"Expanding launch exceeds the thermal boundary")
            for row in run["curve"]:
                ceiling = initial_fraction/row["a_over_ai"]**2
                require(abs(ceiling-row["ceiling"])<2e-9 and
                        -2e-9 <= row["photon_fraction"] <= ceiling+2e-9 <= initial_fraction+4e-9,
                        "Expanding free radiation exceeded its analytic ceiling")
            event = run["first_maximum"]
            require(abs(event["r_photon"]*event["a"]**4-d)<2e-9,"Photon energy failed to redshift")


def physics(products):
    verify_round_states(products["round"],products["checks"])
    verify_uneven_geometry(products["patch"],products["maxwell"])
    verify_photons(products["radiation"],products["round"])


def regenerate(destination, jobs):
    destination = destination.resolve()
    require(not destination.is_relative_to(HERE),'Use a separate output directory')
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(name):
        output = destination/PRODUCTS[name]
        output.parent.mkdir(parents=True,exist_ok=True)
        require(not output.exists(),f'Output already exists: {output}')
        command = [sys.executable,'-I','-B',str(HERE/SCRIPTS[name]),'--output',str(output)]
        print('Recomputing '+SCRIPTS[name],flush=True)
        proc = subprocess.run(command,cwd=HERE,env=env,stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,text=True)
        require(proc.returncode==0, f"{name} exited {proc.returncode}\n{proc.stdout[-4000:]}")
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        list(pool.map(run,PRODUCTS))
    data = {name:read(destination/path) for name,path in PRODUCTS.items()}
    for name,path in PRODUCTS.items():
        compare(data[name],read(HERE/path),(name,))
    physics(data)
    print('Verified all five numerical products and cross-product physical checks.',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true',help='Repeat the numerical experiments')
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--output-dir',type=Path,help='Keep regenerated files in this separate directory')
    args=parser.parse_args()
    require(args.jobs>=1 and (args.output_dir is None or args.verify),'Incompatible replay arguments')
    verify_files()
    if args.verify:
        if args.output_dir:
            regenerate(args.output_dir,args.jobs)
        else:
            with TemporaryDirectory(prefix='hdblast-states-') as tmp:
                regenerate(Path(tmp),args.jobs)
