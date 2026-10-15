#!/usr/bin/env python3
"""Hash-specific replay supporting the written independent Sol final review."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PINS = {
    "UNIFORM_FIRST_VARIATION_WAVEFORM.md": "6eada42c742aec6141925c244df21dc263e3555ae2a064e1e5d3fb7dcc9a82ad",
    "waveform_intervals.py": "f17f3a54d7aa34a0012b52eaccefb7e1e7c319c32ef51ed1286e0c53954f5593",
    "compute_waveform.py": "0eaa55feaf321c6ad57ef3334c9b351a18daba25f4ec1b737d7aaa19b61fa1b4",
    "check_waveform.py": "78d0789a84da840d5a7352d96df5c4c6c61787eca060bd33ea707ad8b531138e",
    "computed/WAVEFORM_ENCLOSURE.json": "2e9273951db60f0400d00dae036074bf886dad481e9d27cf807aff47e2645163",
    "computed/WAVEFORM_CELLS.csv": "bb4f4d659c920afd259c0538182623d3b40aae0f93b5f36a65d5d9405ea6e987",
}


def require(value, label):
    if not value:
        raise RuntimeError(label)


def call(script, args=(), optimized=False):
    cmd = [sys.executable, "-I", "-B"]
    if optimized:
        cmd += ["-O"]
    cmd += [str(ROOT / script), *map(str, args)]
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    return result


def main():
    pins = {}
    for path, expected in PINS.items():
        actual = sha256((ROOT/path).read_bytes()).hexdigest()
        require(actual == expected, "reviewed hash changed: " + path)
        pins[path] = actual
    report = json.loads((ROOT/"computed/WAVEFORM_ENCLOSURE.json").read_text())
    exact_sol = json.loads((HERE/"UNIFORM_WAVEFORM_CHECKS.json").read_text())["exact"]
    require(F(exact_sol["eta"]) == F(report["error_constants"]["eta_upper"]), "root and Sol error differ")
    require(F(exact_sol["eta_tail"]) == F(report["error_constants"]["eta_tail_upper"]), "root and Sol tail differ")
    require(F(exact_sol["support_aware_L2_error_upper"]) == F(report["error_constants"]["norm_error_upper"]), "root and Sol norm error differ")
    for path in ("sol/check_uniform_waveform.py", "sol/UNIFORM_WAVEFORM_CERTIFICATE.md",
                 "sol/UNIFORM_WAVEFORM_CHECKS.json", "sol/UNIFORM_WAVEFORM_CHECKS_OPTIMIZED.json",
                 "archive/PHYSICAL_GAIN_ENCLOSURE.json", "referee/independent_csv_audit.py",
                 "astra/peak_analysis/verify_peaks.py", "astra/peak_analysis/certificate.json"):
        pins[path] = sha256((ROOT/path).read_bytes()).hexdigest()
    replays = []
    for script in ("check_waveform.py", "referee/independent_csv_audit.py", "astra/peak_analysis/verify_peaks.py"):
        parsed = []
        for optimized in (False, True):
            result = call(script, optimized=optimized)
            require(result.returncode == 0, script + " replay failed: " + result.stderr)
            parsed.append(json.loads(result.stdout))
        require(parsed[0] == parsed[1], script + " optimization changed decisions")
        replays.append({"script": script, "normal_and_optimized": "PASS_IDENTICAL", "result": parsed[0]})
    generated = []
    with tempfile.TemporaryDirectory(prefix="hdblast-sol-final-review-") as temporary:
        for optimized in (False, True):
            dest = Path(temporary)/("optimized" if optimized else "normal")
            result = call("compute_waveform.py", ("--output-dir", dest), optimized)
            require(result.returncode == 0, "full generator replay failed: " + result.stderr)
            hashes = {}
            for name in ("WAVEFORM_ENCLOSURE.json", "WAVEFORM_CELLS.csv"):
                got = (dest/name).read_bytes()
                require(got == (ROOT/"computed"/name).read_bytes(), "full generator changed " + name)
                hashes[name] = sha256(got).hexdigest()
            generated.append({"optimized": optimized, "status": "PASS_BYTE_IDENTICAL", "outputs": hashes})
    failures = {
        "drop-determinant-coupling": "first majorant equation",
        "flip-boundary-damping": "boundary jet correct sign identity",
        "use-full-C": "wrong gamma equals full C rejected",
        "omit-velocity-input": "dropped velocity input rejected at peak",
        "delete-boundary-memory": "deleted boundary memory rejected",
    }
    injections = []
    forbidden = HERE/"INJECTED_SHOULD_NOT_WRITE.json"
    require(not forbidden.exists(), "failure injection output unexpectedly preexists")
    for name, label in failures.items():
        for optimized in (False, True):
            result = call("sol/check_uniform_waveform.py", ("--inject-error", name, "--output", forbidden.name), optimized)
            require(result.returncode != 0 and "failed check: " + label in result.stderr,
                    "injection did not fail at named gate: " + name)
            require(not forbidden.exists(), "failed injection wrote passing artifact")
            injections.append({"injection": name, "optimized": optimized, "status": "REJECTED", "failure_label": label})
    # The exterior peak gaps are rational lower bounds from the separately
    # replayed peak certificate.  Use the Sol published error here, avoiding
    # dependence on the slightly tighter Astra internal error assumption.
    eta = F(report["error_constants"]["published_eta"])
    ell = F(report["error_constants"]["ell"])
    gaps = (F("0.000000263559678293613961962519"), F("0.000000590303334992119176697884"))
    require(min(gaps) > 2*eta+ell, "actual peak exclusion fails with Sol error")
    heights = (F("4.235640422182980388580255"), F("4.235640422182980388580257"),
               F("-4.235640422086796048554158"), F("-4.235640422086796048554156"))
    require(F("4.2356403") < heights[0]-eta and heights[1]+eta < F("4.2356406"), "main maximum display not outward")
    require(F("-4.2356406") < heights[2]-eta and heights[3]+eta < F("-4.2356403"), "main minimum display not outward")
    acceptance = {
        "status": "PASS_HASH_SPECIFIC_SOL_FINAL_REVIEW",
        "pins": pins,
        "root_error_constants_match_independent_sol": True,
        "replays": replays,
        "full_4096_cell_generations": generated,
        "failure_injections": injections,
        "peak_exclusion_rechecked_using_sol_published_error": True,
        "review_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "Hash-specific analytic source review and exact-arithmetic replay. Conditional on retained background, reconstruction, exact-family and evolution hypotheses. Internal AI review, not external peer review or formal PDE certification.",
    }
    (HERE/"FINAL_WAVEFORM_REVIEW.json").write_text(json.dumps(acceptance, indent=2, sort_keys=True)+"\n")
    print("PASS: pinned main note, root evaluator/generator/checker, all 4098 cells, 10 explicit failure injections")


if __name__ == "__main__":
    main()
