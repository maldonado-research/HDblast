"""Propagate two orthogonal coexact fields with their angular-mean Maxwell stress and scalar work."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = (
    "t",
    "a",
    "K",
    "X",
    "B",
    "inc",
    "D",
    "phi",
    "pi",
    "H",
    "u",
    "theta",
    "R",
    "c",
    "d",
    "target",
    "comoving_scalar",
    "scalar_em",
    "scalar_pressure",
    "parent_drive",
    "parent_mass",
    "emitted",
    "source_tau",
    "source_r",
    "emitted0",
    "emitted1",
    "emitted2",
    "probe_scalar_debit",
)
MODE = ("E", "B", "incoming", "net_radial", "gauge_work")


def build(parent, temp):
    prior = parent / "parent/parent/prior"
    gen = runpy.run_path(str(parent / "geometry/probe.py"))
    gen["build"](prior, temp, (5, 7))
    path = temp / "probe.cpp"
    src = path.read_text()

    def change(old, new):
        nonlocal src
        if src.count(old) != 1:
            raise ValueError("Unexpected inherited source at " + old)
        src = src.replace(old, new)

    change("constexpr double probe_lambda", "double probe_lambda")
    change(
        "constexpr int channels=",
        """bool backreact=true;
double alphas[3]={1.,0.,0.};
int source_profile=0;
double source_scale=1.;
constexpr int channels=""",
    )
    change(
        "PROBE_DRIVER=PROBE_NET+probe_count, COUNT=PROBE_DRIVER+probe_count",
        "PROBE_DRIVER=PROBE_NET+probe_count, EMIT0=PROBE_DRIVER+probe_count, EMIT1, EMIT2, PROBE_DEBIT, COUNT",
    )
    change("double i, o, a;", "double i, o, a; array<array<double,3>,unforced> more{};")
    # Pass every mode before forming the resulting wall and radial pressure.
    header = (HERE / "closure.hpp").read_text()
    extra, radial = header.split("// For each Maxwell mode")
    change("Values wall(", extra + "\nValues wall(")
    src = (
        src[: src.index("Geo radial(")]
        + "// For each Maxwell mode"
        + radial
        + "\n"
        + src[src.index("Geo shifted(") :]
    )
    change(
        "    return {sheet,",
        """    const auto wextra=extra_wall(y,a,k);
    const double extra_sheet=.5*k*(wextra.e2+wextra.b2);
    sheet+=extra_sheet; e+=extra_sheet; p+=extra_sheet/3.;
    press+=wextra.i2+wextra.o2-.5*wextra.b2;
    flux+=wextra.i2-wextra.o2;
    return {sheet,""",
    )
    change(
        "field[k] = {y[j + I], y[j + O], y[j + GAUGE]};",
        """field[k] = {y[j + I], y[j + O], y[j + GAUGE]};
                if (backreact) for(int m=0;m<unforced;m++) for(int c=0;c<3;c++)
                    field[k].more[m][c]=y[j+3+m+probe_count*c];""",
    )
    change(
        "m.a += mid[k][j] * field[j].a;",
        """m.a += mid[k][j] * field[j].a;
                    if(backreact) for(int m0=0;m0<unforced;m0++) for(int c=0;c<3;c++)
                        m.more[m0][c]+=mid[k][j]*field[j].more[m0][c];""",
    )
    # Source proper time is the time measured by the emitting tube, including the perturbed lapse.
    change(
        "const Geo (&geom)[P], double emitted, double gauge)",
        "const Geo (&geom)[P], const std::array<double,2>& extra_emitted, double gauge)",
    )
    change("(m==unforced ? 0. : emitted)", "(m==unforced ? 0. : extra_emitted[m])")
    change(
        "void probe_cell(", (HERE.parent / "fields/waveforms.hpp").read_text() + "\nvoid probe_cell("
    )
    change(
        "double emitted = incoming(y[SU], dose, norm, stop), outgoing = 0.;",
        """double primary=incoming(y[SU],dose,norm,stop);
        double emitted=alphas[0]*primary, outgoing=0.;
        std::array<double,2> extra_emitted={alphas[1]*primary,alphas[2]*source_scale*
            shaped(y[SU],source_profile,primary)};""",
    )
    change("probe_cell(y,dy,cell,geom,emitted,w.k);", "probe_cell(y,dy,cell,geom,extra_emitted,w.k);")
    change(
        "        probe_wall(y,dy,w,a,h);",
        """        double emit2=emitted*emitted, ret2=outgoing*outgoing;
        const double tube_weight=2.*g.n*r*r*r;
        dy[EMIT0]=tube_weight*emit2;
        dy[EMIT1]=tube_weight*extra_emitted[0]*extra_emitted[0];
        dy[EMIT2]=tube_weight*extra_emitted[1]*extra_emitted[1];
        if(backreact) for (int m=0;m<unforced;m++) {
            const int node=index(P*cells-1,3+m);
            const double inc=extra_emitted[m], out=y[node+probe_count];
            const double b=probe_lambda[m]*y[node+2*probe_count]/(r*r);
            fin+=inc*inc-out*out; pr+=inc*inc+out*out-.5*b*b;
            emit2+=inc*inc; ret2+=out*out;
        }
        const Extra v=extra_wall(y,a,w.k);
        const double electric_minus_magnetic=w.x*w.x-w.b*w.b+v.e2-v.b2;
        probe_wall(y,dy,w,a,h);""",
    )
    change(".5 * w.dk * (w.x * w.x - w.b * w.b)", ".5 * w.dk * electric_minus_magnetic")
    change(
        "w.dk*(w.x*w.x-w.b*w.b);",
        """w.dk*electric_minus_magnetic;
        dy[PROBE_DEBIT]=.5*a3*pi*w.dk*(v.b2-v.e2);""",
    )
    change("pi * (w.b * w.b - w.x * w.x)", "pi * (-electric_minus_magnetic)")
    change("r * r * r * emitted * emitted", "r * r * r * emit2")
    change("r * r * r * outgoing * outgoing", "r * r * r * ret2")
    change("a * a * a * iwall * iwall", "a * a * a * (iwall*iwall+v.i2)")
    change("a * a * a * ow * ow", "a * a * a * (ow*ow+v.o2)")
    change("if (argc != 36)", "if (argc != 43)")
    change(
        "double last_probe_print=-1.;",
        """double last_probe_print=-1.;
    backreact=std::stoi(argv[36]); alphas[1]=std::stod(argv[37]); alphas[2]=std::stod(argv[38]);
    alphas[0]=std::sqrt(1.-alphas[1]*alphas[1]-alphas[2]*alphas[2]);
    probe_lambda[0]=std::stod(argv[39]); probe_lambda[1]=std::stod(argv[40]);
    source_profile=std::stoi(argv[41]); source_scale=std::stod(argv[42]);""",
    )
    change(
        "if(t<.020 || t>6.24 || t-last_probe_print>=.0004)",
        "if(t+dt>=end-1.e-12 || t<.020 || t-last_probe_print>=.0003 || (y[PHI]>threshold_phi-field_width*.51 && gas.comoving_rest<.00020))",
    )
    change(
        "const double record[]={t+dt,a,wnew.k,wnew.x,wnew.b,y[s.index(0,I)],gas.comoving_rest,y[PHI]};",
        """const auto cs=oscillator(y);
            const double R=cs.mass_ratio, sp=std::sqrt(std::max(0.,R*R-1.));
            const double cp=sp>0 ? h*sp+R*(cs.mass_slope*y[PI]/(2.*daughter_mass))/sp : 0.;
            const double tp=sp>0 ? std::pow(std::acos(-1.),2)/6.*(feedback_action/std::pow(daughter_mass,4))*cs.power/(R*sp*sp) : 0.;
            const double driver=wnew.e-y[RAD]-wnew.sheet-cs.energy-gas.energy;
            const double record[]={t+dt,a,wnew.k,wnew.x,wnew.b,y[s.index(0,I)],gas.comoving_rest,y[PHI],
                y[PI],h,gas.u,gas.theta,R,cp,charge_mass*wnew.x,tp,a*a*a*driver,y[PSI_EM_WORK],
                y[PSI_PRESS_WORK],y[CHI_DRIVE_WORK],y[CHI_MASS_WORK],y[EMITTED],y[SU],y[VOLUME],
                y[EMIT0],y[EMIT1],y[EMIT2],y[PROBE_DEBIT]};""",
    )
    path.write_text(src)
    exe = temp / "backreaction"
    subprocess.run(
        ["g++", "-O3", "-std=c++17", str(path), "-I", str(prior / "boost"), "-o", str(exe)], check=True
    )
    return exe


def solve(exe, parent, outdir, temp, label, args, mesh, end, keep):
    prior = parent / "parent/parent/prior"
    metadata = json.loads((parent / "parent/history/markers.json").read_bytes())
    raw = temp / (label + ".bin")
    wrapper = temp / label
    ws = [2, metadata["action"], prior / "boost/BOOST_TABLE.dat", 1, 1.0, raw, *args]
    wrapper.write_text(
        "#!/bin/sh\nexec " + shlex.quote(str(exe)) + ' "$@" ' + shlex.join([str(x) for x in ws]) + "\n"
    )
    wrapper.chmod(0o755)
    cfg = json.loads((prior / "prior/envelope/ENVELOPE.json").read_bytes())["parameters"]["mass_4e12"]
    config = dict(
        cfg,
        sample=0,
        end=end,
        cfl=cfg["cfl"] * (0.5 if mesh == "half" else 1.0),
        early_cells=cfg["early_cells"] * (2 if mesh == "radial" else 1),
        late_cells=cfg["late_cells"] * (2 if mesh == "radial" else 1),
    )
    runner = runpy.run_path(str(prior / "prior/coupled/joined.py"))["run"]
    try:
        report = runner(wrapper, config)
    except subprocess.CalledProcessError as e:
        raise RuntimeError((label, e.stderr)) from e
    arr = np.fromfile(raw, dtype="=f8").reshape(-1, len(BASE) + 3 * len(MODE))
    raw.unlink()
    if not np.isfinite(arr).all():
        raise ValueError(("nonfinite state", label))
    b = {k: arr[:, i] for i, k in enumerate(BASE)}
    modes = arr[:, len(BASE) :].reshape(-1, 3, len(MODE)).transpose(1, 2, 0)
    t = b["t"]
    a = b["a"]
    W = 0.5 * a * a * a * a * b["K"] * (modes[:, 0] ** 2 + modes[:, 1] ** 2)
    budget_res = abs(W - modes[:, 3] - modes[:, 4])
    res = np.max(
        budget_res / (np.max(W + abs(modes[:, 3]) + abs(modes[:, 4]), axis=1)[:, None] + 1e-30)
    )
    if res > 3.0e-5:
        raise ValueError((label, res))
    if np.any(modes[-1]):
        raise ValueError("Unforced probe became nonzero")
    totals = {
        k: float(report[k])
        for k in (
            "friedmann_error",
            "driver_balance_error",
            "surface_balance_error",
            "sheet_balance_error",
            "tube_mass_error",
            "tube_radius_error",
            "gas_balance_error",
            "oscillator_balance_error",
            "min_lapse",
        )
    }
    message = dict(
        feedback=bool(args[0]),
        launch_coefficients=list(args[1:3]),
        curl_magnitudes=list(args[3:5]),
        second_profile=("early", "prepare", "release", "birth", "late_birth", "ramp", "bipolar")[
            int(args[5])
        ],
        second_scale=args[6],
        events={},
        mode_ledger_relative_max=res,
        checks=totals,
    )
    peak_phi = int(np.argmax(b["phi"]))
    message["maximum_phi"] = float(b["phi"][peak_phi])
    message["time_of_maximum_phi"] = float(t[peak_phi])
    message["backward_motion_after_final_release_start"] = bool(
        np.any((b["phi"] > 0.5265) & (b["pi"] < 0))
    )
    message["end_computed_time"] = float(report["end_tau"])

    message["maximum_K"] = float(np.max(b["K"]))
    # For a delayed channel, inspect the common past before its emitter is active.
    if args[5] == 2:
        before = np.searchsorted(b["source_tau"], 6.206) - 1
        if before >= 0:
            message["before_delayed_emission"] = {
                k: float(b[k][before]) for k in ("t", "source_tau", "phi", "pi", "K", "D", "emitted2")
            }

    message["mean_emitted_energy"] = float(sum(b[k][-1] for k in ("emitted0", "emitted1", "emitted2")))
    ranges = (
        abs(b["comoving_scalar"])
        + abs(b["scalar_em"])
        + abs(b["scalar_pressure"])
        + abs(b["parent_drive"])
        + abs(b["parent_mass"])
    )
    message["scalar_ledger_relative_peak"] = float(
        np.max(
            abs(
                b["comoving_scalar"]
                - b["scalar_em"]
                + b["scalar_pressure"]
                + b["parent_drive"]
                + b["parent_mass"]
            )
        )
        / max(ranges)
    )
    birth = int(np.searchsorted(b["D"], 0.0, side="right"))
    events = {"end": t[-1]}
    if birth < len(t):
        events.update(birth=t[birth - 1])
        # Compare at a fixed exposure from #197, irrespective of how the new total pulse work changes.
        absrate = 1.0e7 * np.abs(b["X"] * b["D"] * b["u"] / np.hypot(1.0, b["u"]))
        ac = np.r_[0.0, np.cumsum(0.5 * np.diff(t) * (absrate[1:] + absrate[:-1]))]
        target = 0.01 * metadata["runs"]["dedicated_half"]["absolute_work"]
        if ac[-1] > target:
            events["fixed_work"] = np.interp(target, ac, t)
    for name, ev in events.items():
        at = lambda v: float(np.interp(ev, t, v))
        message["events"][name] = {
            "time": float(ev),
            **{k: at(b[k]) for k in BASE if k not in ("t", "inc", "c", "d", "target")},
            "mode_E": [at(modes[m, 0]) for m in range(2)],
            "mode_B": [at(modes[m, 1]) for m in range(2)],
            "mode_energy_a4": [at(W[m]) for m in range(2)],
            "mode_gauge_work_a4": [at(modes[m, 4]) for m in range(2)],
            "mode_net_work_a4": [at(modes[m, 3]) for m in range(2)],
        }
    if keep and "fixed_work" in events:
        mask = (t >= events["birth"]) & (t <= events["fixed_work"] + 0.0001)
        saved = {
            k: b[k][mask]
            for k in ("t", "a", "H", "X", "B", "D", "u", "theta", "R", "c", "d", "target")
        }
        path = outdir / "source_traces"
        path.mkdir(exist_ok=True)
        np.savez_compressed(
            path / (label + ".npz"),
            **saved,
            kappa=np.array(args[3:5]),
            mode_E=modes[:2, 0][:, mask],
            mode_B=modes[:2, 1][:, mask],
            marker_time=events["fixed_work"],
        )
    return message


def validate(results):
    for label, row in results.items():
        if row["mode_ledger_relative_max"] > 3e-5 or row["scalar_ledger_relative_peak"] > 1e-7:
            raise ValueError(("Energy balance", label))
        if row["checks"]["friedmann_error"] > 1e-6 or row["checks"]["tube_mass_error"] > 4e-4:
            raise ValueError(("Mean geometry", label))
        if abs(row["end_computed_time"] - row["events"]["end"]["time"]) > 1e-10:
            raise ValueError(("Incomplete timing record", label))

    def birth(label):
        return results[label]["events"]["birth"]["time"]

    if {"baseline", "same_mode_rotation"} <= results.keys():
        if abs(birth("same_mode_rotation") - birth("baseline")) > 2e-5:
            raise ValueError("Equal-profile rotation of the lowest modes changed first birth")
    if {"baseline", "pair_0.04_active", "pair_0.08_active"} <= results.keys():
        if not birth("baseline") < birth("pair_0.04_active") < birth("pair_0.08_active"):
            raise ValueError("Unexpected order of the measured first births")
    for mesh in ("long", "half", "radial"):
        row = results.get("pair_0.20_" + mesh)
        if row and (
            "birth" in row["events"]
            or row["maximum_phi"] >= 0.53
            or not row["backward_motion_after_final_release_start"]
        ):
            raise ValueError(("Reflected scalar control", mesh))
    for mesh in ("native", "half", "radial"):
        row = results.get("delayed_matched_" + mesh)
        if row and (
            "birth" in row["events"] or abs(row["fixed_total_physical_emission_relative_error"]) > 2e-8
        ):
            raise ValueError(("Delayed equal-energy control", mesh))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parent", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--cases", nargs="+")
    p.add_argument(
        "--profiles",
        type=Path,
        required=True,
        help="FIELDS.json produced by the calibrated outer-profile solver",
    )
    opt = p.parse_args()
    base = opt.parent.resolve()
    out = opt.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    tables = json.loads(opt.profiles.read_bytes())["dedicated_half"]["channels"]
    profiles = ("early", "prepare", "release", "birth", "late_birth", "ramp", "bipolar")
    tasks = []

    def add(label, total, k1, k2, feedback, mesh="native", profile="early", end=6.36, withhold=False):
        # Match integrated emitting-tube energy on the reference geometry. The
        # perturbed run reports the actual energy after radius/lapse changes.
        multiplier = (
            0.0 if withhold else tables[str(k2) + "_" + profile]["flux_normalizing_multiplier"]
        )
        args = (
            int(feedback),
            float(total / np.sqrt(2.0)),
            float(total / np.sqrt(2.0)),
            k1,
            k2,
            profiles.index(profile),
            multiplier,
        )
        tasks.append((label, args, mesh, end))

    add("baseline", 0.0, 5, 7, True)
    add("same_mode_rotation", 0.08, 2, 2, True)
    for total in (0.04, 0.08):
        for feedback in (False, True):
            add(f"pair_{total:.2f}_" + ("active" if feedback else "passive"), total, 5, 7, feedback)
    for mesh in ("half", "radial"):
        add("pair_0.08_" + mesh, 0.08, 5, 7, True, mesh)
        add("pair_0.20_" + mesh, 0.20, 5, 7, True, mesh, end=8.0)

    for total in (0.10, 0.12, 0.14, 0.16, 0.20):
        add(f"pair_{total:.2f}_long", total, 5, 7, True, end=8.0)
    for total in (0.01, 0.025 * np.sqrt(2.0), 0.05, 0.08):
        add(f"delayed_{total:.4f}", total, 7, 7, True, profile="release", end=8.0)
    add("simultaneous_time_control", 0.025 * np.sqrt(2), 7, 7, True, end=8.0)
    add(
        "withheld_time_control",
        0.025 * np.sqrt(2),
        7,
        7,
        True,
        profile="release",
        end=8.0,
        withhold=True,
    )
    for mesh in ("native", "half", "radial"):
        add("delayed_matched_" + mesh, 0.025 * np.sqrt(2), 7, 7, True, mesh, "release", 8.0)
    for total in (0.005, 0.01, 0.02):
        add(f"prepare_{total:.3f}", total, 7, 40, True, profile="prepare", end=8.0)
    if opt.cases:
        wanted = set(opt.cases)
        tasks = [task for task in tasks if task[0] in wanted]
        if wanted != {x[0] for x in tasks}:
            raise ValueError("Unknown case")
    with tempfile.TemporaryDirectory(prefix="hd-mean-backreaction-") as wd:
        temp = Path(wd)
        exe = build(base, temp)

        def one(item):
            label, args, mesh, end = item
            print(label, flush=True)
            result = solve(exe, base, out, temp, label, args, mesh, end, True)
            if label.startswith("delayed_matched_"):
                # For a delayed waveform, a changed source radius alters physical energy
                # even at fixed launch coefficient. Adjust only that waveform to restore
                # the original integrated emitting-tube energy on this realization.
                target = json.loads(opt.profiles.read_bytes())[
                    "dedicated_" + ("spatial" if mesh == "radial" else mesh)
                ]["unit_emitted_energy"]
                energy = result["events"]["end"]
                for _ in range(3):
                    if abs(result["mean_emitted_energy"] / target - 1.0) < 2e-9:
                        break
                    fac = np.sqrt(
                        (target - energy["emitted0"] - energy["emitted1"]) / energy["emitted2"]
                    )
                    args = (*args[:-1], args[-1] * fac)
                    result = solve(exe, base, out, temp, label, args, mesh, end, True)
                    energy = result["events"]["end"]
                result["fixed_total_physical_emission_relative_error"] = (
                    result["mean_emitted_energy"] / target - 1.0
                )
                if abs(result["fixed_total_physical_emission_relative_error"]) > 2.0e-8:
                    raise RuntimeError("Physical source energy matching did not converge")
            return label, result

        with ThreadPoolExecutor(opt.jobs) as pool:
            results = dict(pool.map(one, tasks))
    validate(results)
    (out / "SOURCE.json").write_text(
        json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
