"""Build the working note in a fresh external directory; no mathematics runs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

NAME = "CHARACTERISTIC_ODD_TRACE"
HERE = Path(__file__).resolve().parent

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir
    need(out.is_absolute() and out.resolve() == out and not out.exists(),
         "Fresh canonical absolute output directory required")
    need(all(not p.is_symlink() for p in [out, *out.parents]), "Output symlink rejected")
    need(not out.is_relative_to(HERE), "Build output must be outside source directory")
    source = HERE / (NAME + ".tex")
    need(source.is_file() and not source.is_symlink(), "Regular editable source required")
    out.mkdir(parents=True)
    shutil.copyfile(source, out / source.name)
    env = os.environ.copy()
    env["SOURCE_DATE_EPOCH"] = "1793750400"
    env["FORCE_SOURCE_DATE"] = "1"
    compiler = subprocess.check_output(["pdflatex", "--version"], text=True).splitlines()[0]
    passes = []
    for number in range(1, 4):
        proc = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", source.name],
                              cwd=out, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (out / ("pass-" + str(number) + ".txt")).write_bytes(proc.stdout)
        need(proc.returncode == 0, "pdflatex failed on pass " + str(number))
        passes.append({"pass": number, "exit_code": proc.returncode,
                       "stdout_sha256": digest(proc.stdout)})
    log = (out / (NAME + ".log")).read_text(errors="replace")
    need(not re.search(r"undefined references|Citation .*undefined|Reference .*undefined|Overfull \\[hv]box", log),
         "Unresolved reference or overflowing final layout")
    pdf = out / (NAME + ".pdf")
    info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    match = re.search(r"^Pages:\s+(\d+)\s*$", info, re.M)
    need(match is not None, "PDF page count required")
    receipt = {"schema": "HDBLAST_ODD_TRACE_MANUSCRIPT_BUILD_V1",
               "status": "PASS_THREE_PASS_WORKING_NOTE_BUILD",
               "compiler": compiler, "pages": int(match.group(1)), "passes": passes,
               "source": {"filename": source.name, "bytes": source.stat().st_size,
                          "sha256": digest(source.read_bytes())},
               "pdf": {"filename": pdf.name, "bytes": pdf.stat().st_size,
                       "sha256": digest(pdf.read_bytes())},
               "final_undefined_references_or_overflow": False,
               "mathematical_producers_executed": 0}
    (out / "COMPILATION.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "pages": receipt["pages"],
                      "pdf_sha256": receipt["pdf"]["sha256"]}, sort_keys=True))

if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, OSError, ValueError, subprocess.SubprocessError) as error:
        print(json.dumps({"status": "FAIL", "witness": str(error)}, sort_keys=True))
        sys.exit(2)
