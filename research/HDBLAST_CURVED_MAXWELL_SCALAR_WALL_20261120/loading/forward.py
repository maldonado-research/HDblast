"""Reproduce the forward Cauchy evolution with a pulse specified at the inner tube."""
import importlib.util
import math
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
COLUMNS = ("t", "a", "H", "psi", "pi", "rad", "sheet", "e", "V", "source_proper_time",
           "source_emitted", "source_returned", "source_radius", "K", "Bmag", "X", "wall_mass")


def waveform_constants():
    spec = importlib.util.spec_from_file_location("curved_wall_material", HERE / "loading.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stop, signal = module.waveform(40.)
    # At s=3 the cutoff windows are unity and sqrt(sigmoid(0))*sigmoid(0)=1/sqrt(8).
    return stop, (1/math.sqrt(8))/signal(3.)


def run(executable, constants, *, cells, depth=.002, dose=.123, late_cells=6, late_gauge=-1.,
        end=16., trace=False):
    stop, norm = constants
    args = [cells, depth, dose, end, norm, stop, .08, late_cells, .08, late_gauge, int(trace)]
    rows = subprocess.check_output([executable, *map(str, args)], text=True).splitlines()
    output, series = {}, []
    for row in rows:
        name, *values = row.split()
        if name == "@":
            if len(values) != len(COLUMNS):
                raise ValueError("Forward trace has the wrong number of columns")
            series.append(list(map(float, values)))
        else:
            output[name] = float(values[0])
    for key in ("steps", "postexit_em90_intervals"):
        output[key] = int(output[key])
    if series:
        output["time_series"] = series

    if min(output["min_lapse"], output["min_bulk_mass"], output["min_radial_gradient_norm"]) <= 0:
        raise ValueError("Forward integration lost the ordinary positive-mass branch")
    if output["tube_returned_proper_energy"] > output["tube_incident_proper_energy"] + 1e-12:
        raise ValueError("Forward source ledger returned more energy than was introduced")
    if max(output["friedmann_error"], output["surface_balance_error"], output["sheet_balance_error"]) > 1e-4:
        raise ValueError("Wall gravitational or energy identity failed")
    if late_gauge >= 0 and (output["end_radiation"] != 0 or output["fluid_comoving_radiation_produced"] != 0):
        raise ValueError("A coherent Maxwell test acquired a separate fluid component")
    return output


def calculate():
    constants = waveform_constants()
    scenarios = {
        "fluid_coarse": dict(cells=24), "fluid_medium": dict(cells=48),
        "fluid_primary": dict(cells=64, trace=True),
        "maxwell_coarse": dict(cells=24, late_gauge=13200),
        "maxwell_medium": dict(cells=48, late_gauge=13200),
        "maxwell_primary": dict(cells=64, late_gauge=13200, trace=True),
        "maxwell_no_return_coupling": dict(cells=24, late_gauge=0),
        "maxwell_weaker": dict(cells=24, late_gauge=9000),
        "maxwell_stronger": dict(cells=24, late_gauge=16000),
        "maxwell_trapped": dict(cells=24, late_gauge=22000),
        "change_depth": dict(cells=24, depth=.003),
        "change_depth_and_dose": dict(cells=24, depth=.003, dose=.124),
        "same_early_mesh": dict(cells=12), "no_coarsening": dict(cells=12, late_cells=12),
        "vacuum": dict(cells=6, dose=0, end=.5),
    }
    with TemporaryDirectory() as temporary:
        executable = Path(temporary)/"annulus"
        subprocess.run(["g++", "-O3", "-std=c++17", str(HERE/"forward.cpp"), "-o", str(executable)],
                       check=True)
        outputs = {name: run(executable, constants, **options) for name, options in scenarios.items()}
    fluid, wave = outputs["fluid_primary"], outputs["maxwell_primary"]
    if not (.60 < fluid["radiation_expanding_efolds"] < .62 and
            .09 < wave["longest_em90_expanding_efolds"] < .14 and wave["postexit_em90_intervals"] == 1):
        raise ValueError("The independent source benchmark failed its expanding-radiation test")
    if not (wave["exit_complete_tau"] < wave["maximum_tau"] and wave["end_H"] < 0 and
            wave["at_peak_comoving_late_scalar_to_sheet"] > 0 and
            outputs["maxwell_weaker"]["peak_exit_sheet_fraction"] < .9 and
            outputs["maxwell_stronger"]["peak_exit_sheet_fraction"] < .9):
        raise ValueError("The late scalar-to-Maxwell comparison failed")
    if not (abs(outputs["vacuum"]["end_radius"]-math.sqrt(4-.25)) < 1e-10 and
            abs(outputs["no_coarsening"]["radiation_expanding_efolds"] -
                outputs["same_early_mesh"]["radiation_expanding_efolds"]) < 1e-7):
        raise ValueError("Vacuum or late spatial coarsening control failed")
    return {
        "units": "AdS length=1; densities and pressures in units of tuned tension",
        "source": {"proper_depth": .002, "width": .001, "dose": .123, "calibration_dose": .12,
                   "cap": 40., "pulse_end_source_proper_time": constants[0]*.001,
                   "late_mesh_wall_time": .08},
        "trace_columns": list(COLUMNS),
        "parameters": scenarios,
        "runs": outputs,
    }


if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(), indent=2, allow_nan=False)+"\n")
