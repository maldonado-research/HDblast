"""Replay the tagged rest-parent calculation on the archived HDblast wall histories."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

CASES = {
    "primary": ("tagged.py", "TAGGED.json"),
    "independent": ("tagged_transfer.py", "TAGGED_TRANSFER.json"),
    "flat": ("gaussian.py", "GAUSSIAN.json"),
}
SIGNS = {"quiet": 0, "plus": -1, "minus": 1, "far": 1}
WINDOWS = {"wide": [1.15, 2.05], "narrow": [1.25, 1.95]}
BRANCHES = {"annihilation_number": "negative_frequency_number",
            "annihilation_energy": "negative_frequency_energy",
            "stimulated_number": "positive_frequency_number",
            "stimulated_energy": "positive_frequency_energy"}


def require(condition, description):
    if not condition:
        raise ValueError(description)


def close(x, y, label, rel=1e-8, absolute=0.0):
    require(math.isfinite(x) and math.isfinite(y) and
            math.isclose(x, y, rel_tol=rel, abs_tol=absolute), f"{label}: {x!r} vs {y!r}")


def match_saved(actual, expected, path):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), f"Schema mismatch at {path}")
        for key in expected:
            match_saved(actual[key], expected[key], f"{path}/{key}")
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), f"Length mismatch at {path}")
        for idx, item in enumerate(expected):
            match_saved(actual[idx], item, f"{path}/{idx}")
    elif isinstance(expected, float):
        if path.startswith("flat/"):
            close(actual, expected, path, rel=1e-8)
        elif path.endswith(("wronskian_error","wronskian","spectral_minus_time",
                            "final_eta_difference","relative_final_mode_error")):
            close(actual, expected, path, absolute=1e-8)
        elif "/contrasts" in path:
            close(actual, expected, path, rel=3e-5, absolute=1e-6)
        else:
            close(actual, expected, path, rel=3e-7, absolute=1e-12)
    else:
        require(type(actual) is type(expected) and actual == expected, f"Mismatch at {path}")


def check_primary(data):
    require(data["switches"] == WINDOWS and data["cut"] == 2.1 and data["omega_screen"] == 60,
            "Switch, screen or readout mismatch")
    require(set(data["grids"]) == {"fine", "coarse"}, "Missing mesh")
    close(data["flat"]["number"], data["flat"]["analytic_number"], "flat compact switch", rel=1e-10)
    tested = 0
    for mesh, arms in data["grids"].items():
        require(arms.keys() == SIGNS.keys(), f"Missing arm on {mesh}")
        for arm, sign in SIGNS.items():
            for window in WINDOWS:
                sample = arms[arm][window]
                require(sample["max_screened_single_daughter_physical_momentum"] < 6.49e-6**(-1/3),
                        "Screen exceeds the specified gravity scale")
                close(sample["time_identity"]["exact"], sample["full"]["number"],
                      f"{mesh}/{arm}/{window} Parseval", rel=2e-8)
                close(sum(sample["time_identity"]["terms"].values()), sample["time_identity"]["exact"],
                      f"{mesh}/{arm}/{window} time terms", rel=1e-10)
                for obs in ("number", "energy"):
                    fraction = 1-sample["screened"][obs]/sample["full"][obs]
                    require(0 < fraction < 1, "Unphysical screen fraction")
                    close(fraction,sample["fraction_above_screen"][obs],"Screen fraction",absolute=1e-13)
                for band in ("full", "screened"):
                    vals = sample[band]
                    close(vals["annihilation_number"]+vals["stimulated_number"],vals["number"],"Inclusive pairs")
                    close(vals["annihilation_energy"]+vals["stimulated_energy"],vals["energy"],"Daughter energy")
                    close(vals["annihilation_number"]-vals["stimulated_number"],vals["net_parent_loss"],
                          "Signed depletion")
                    for key in ("number", "energy", "net_parent_loss", *BRANCHES):
                        quiet = arms["quiet"][window][band][key]
                        percent = 100*(vals[key]/quiet-1)
                        close(arms[arm]["contrasts_percent"][window][band][key],percent,
                              f"{mesh}/{arm}/{window}/{band}/{key} stored contrast",absolute=1e-11)
                        require(percent*sign > 0 if sign else percent == 0,
                                f"Pulse sign for {mesh}/{arm}/{window}/{band}/{key}")
                        tested += bool(sign)
    return tested


def check_independent(primary, independent):
    require(set(independent) == {"4096","8192","16384","32768"}, "Independent mesh sequence changed")
    end = independent["32768"]
    worst = 0.0
    for arm in ("quiet","plus"):
        for window in WINDOWS:
            native = primary["grids"]["fine"][arm][window]
            other = end["arms"][arm]["windows"][window]
            spectrum = other["spectrum"]["32"]
            require(set(other["spectrum"]) == {"4","8","16","32"}, "Independent FFT sequence changed")
            pairs = [(native["full"][key],spectrum[frequency]) for key,frequency in BRANCHES.items()]
            pairs.extend((native["full"][key],spectrum[key]) for key in ("number","energy"))
            pairs.extend([
                (native["full"]["net_parent_loss"],
                 spectrum["negative_frequency_number"]-spectrum["positive_frequency_number"]),
                (native["time_identity"]["exact"],other["local_number"]),
                (native["time_identity"]["terms"]["field"],other["parent_number"]),
                (native["time_identity"]["terms"]["switch"],other["switch_number"]),
                (native["time_identity"]["terms"]["interference"],other["interference"]),
            ])
            for value, reference in pairs:
                close(value, reference, f"Independent {arm}/{window} coefficient", rel=5e-7)
                worst = max(worst, abs(value-reference)/abs(value))
            for kind in ("number","energy"):
                c = primary["grids"]["fine"]["plus"]["contrasts_percent"][window]["full"][kind]
                ci = end["contrasts"][window][kind+"_percent"]
                close(c, ci, f"Independent near+ {window}/{kind} contrast", absolute=1e-5)
    return worst


def check_flat(data):
    require([r["mT"] for r in data["data"]]==[1,2,4,8], "Gaussian durations changed")
    for row in data["data"]:
        u, a = row["mT"], row["analytic"]
        close(a["inclusive_pairs"],1+1/(2*u*u),"Gaussian inclusive count",rel=1e-13)
        close(a["annihilation"]+a["stimulated_creation"],a["inclusive_pairs"],"Gaussian count identity")
        close(a["annihilation"]-a["stimulated_creation"],a["signed_depletion"],"Gaussian depletion identity")
        close(a["daughter_energy"]-a["signed_depletion"],a["external_work"],"Gaussian work identity")
        require(0 < a["stimulated_creation"] < row["stimulated_creation_upper_bound"] and
                a["external_work"] > 0, "Gaussian tail or work")
        for key, value in a.items():
            close(row["quadrature"]["8192"][key],value,f"Gaussian u={u}, {key}",rel=2e-9)


def main():
    study = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository",type=Path,default=study.parent.parent,
                        help="HDblast checkout containing the published source inputs")
    parser.add_argument("--stored-only",action="store_true",help="verify the archived data without recomputing")
    args = parser.parse_args()
    result = {}
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for name, (script, output) in CASES.items():
        folder = study/name
        saved = json.loads((folder/output).read_bytes())
        if args.stored_only:
            result[name] = saved
        else:
            cmd = [sys.executable,"-I","-B",str(folder/script)]
            if name != "flat":
                cmd.append(str(args.repository.resolve()))
            run = subprocess.run(cmd,check=True,stdout=subprocess.PIPE,env=env)
            result[name] = json.loads(run.stdout)
            match_saved(result[name],saved,name)
    tested = check_primary(result["primary"])
    worst = check_independent(result["primary"],result["independent"])
    check_flat(result["flat"])
    mode = "archived numbers checked" if args.stored_only else "three calculations regenerated and checked"
    print(f"{mode}; {tested} signed pulse contrasts; maximum primary/independent coefficient "
          f"difference {worst:.3g} relative; four switched Gaussian controls passed.")


if __name__ == "__main__":
    main()
