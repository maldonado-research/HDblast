#!/usr/bin/env python3
"""Data-only, bounded artifact retention on the read-only model job."""
import argparse
import importlib.util
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    for name in ("round-root", "context", "out"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    # Explicit trusted sibling; isolated Python omits cwd, environment and user-site imports.
    spec = importlib.util.spec_from_file_location("trusted_round_control", Path(__file__).with_name("round_control.py"))
    control = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(control)
    control.safe_files(args.round_root)
    context = control.read(args.context)
    if Path(args.round_root).name != Path(context["round_path"]).name:
        raise ValueError("Incorrect selected round path")
    destination = Path(args.out)
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(args.context, destination / "context.json")
    shutil.copytree(args.round_root, destination / "round")
    control.safe_files(destination / "round")


if __name__ == "__main__":
    main()
