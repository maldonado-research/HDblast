#!/usr/bin/env python3
"""Five focused input/convention rejection controls; no upstream proof replay."""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def producer():
    spec = importlib.util.spec_from_file_location("response_bound_producer", ROOT/"check_bounds.py")
    need(spec is not None and spec.loader is not None, "PRODUCER_LOADER")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def calculate():
    bound = producer()
    bound.authenticate()
    conventions = bound.load_json((ROOT/"CONVENTIONS.json").read_bytes())
    bound.validate_conventions(conventions)
    results = []
    with tempfile.TemporaryDirectory(prefix="response-bound-control-") as temp:
        target = Path(temp)/"source"
        shutil.copytree(ROOT, target)
        # copytree preserves frozen modes. Only this disposable fixture is
        # thawed: directory write/traverse permits cleanup, and only the
        # intentionally corrupted report needs file write permission.
        for directory in [target] + [p for p in target.rglob("*") if p.is_dir()]:
            directory.chmod(0o700)
        report = target/"inputs/ACTUAL_SOURCES_NORMAL.json"
        raw = report.read_bytes()
        report.chmod(0o600)
        report.write_bytes(b" " + raw[1:])
        try:
            bound.authenticate(target)
        except ValueError as error:
            need(str(error).startswith("INPUT_AUTHENTICATION:"), "SOURCE_CONTROL_WRONG_REJECTION")
            results.append({"name": "changed_source_report_byte", "status": "REJECTED_INPUT_AUTHENTICATION"})
        else:
            raise ValueError("SOURCE_CONTROL_ACCEPTED")
    changes = [
        ("dropped_retarded_half", "retarded_K_factor", "1"),
        ("wrong_derivative_weight", "derivative_weight", "3"),
        ("wrong_point_to_global_value_weight", "value_weight", "1"),
        ("wrong_lower_transfer", "lower_transfer", "1-q"),
    ]
    for name, key, value in changes:
        changed = copy.deepcopy(conventions)
        changed[key] = value
        try:
            bound.validate_conventions(changed)
        except ValueError as error:
            need(str(error) == "REGISTERED_CONVENTION_MISMATCH", "CONVENTION_WRONG_REJECTION:"+name)
            results.append({"name": name, "status": "REJECTED_REGISTERED_CONVENTION"})
        else:
            raise ValueError("CONVENTION_CONTROL_ACCEPTED:"+name)
    need(len(results) == 5, "CONTROL_COUNT")
    return {
        "status": "PASS_FIVE_TARGETED_REJECTION_CONTROLS",
        "scope": "Five in-process controls in this execution mode; convention registration and source-byte checks, not a proof of analytic premises or an exhaustive fault campaign.",
        "producer_sha256": hashlib.sha256((ROOT/"check_bounds.py").read_bytes()).hexdigest(),
        "controls_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "controls": results,
        "upstream_scientific_producers_replayed": False,
    }


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,
         "USE_PYTHON_I_S_B")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw = json.dumps(calculate(), indent=2, sort_keys=True)+"\n"
    if args.output:
        output = args.output.resolve()
        need(not output.is_relative_to(ROOT), "OUTPUT_INSIDE_SOURCE")
        with output.open("x") as handle:
            handle.write(raw)
    else:
        sys.stdout.write(raw)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        sys.exit("RESPONSE_BOUND_CONTROLS_FAILED:"+str(error))
