#!/usr/bin/env python3
"""Two fresh deterministic three-pass PDF builds; reject final warnings."""
from hashlib import sha256
from pathlib import Path
import json
import os
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
STEM = "INCOMING_PREPARATION_LEDGER"
EPOCH = "1793750400"  # 4 November 2026, 00:00:00 UTC.


def main():
    source = HERE / (STEM + ".tex")
    (HERE / "review").mkdir(exist_ok=True)
    env = dict(os.environ, SOURCE_DATE_EPOCH=EPOCH, FORCE_SOURCE_DATE="1", TZ="UTC")
    builds = []
    for name in ("fresh-build-1", "fresh-build-2"):
        dest = HERE / name
        dest.mkdir(exist_ok=False)
        passes = []
        for number in range(1, 4):
            command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                       "-output-directory=" + name, source.name]
            run = subprocess.run(command, cwd=HERE, env=env, capture_output=True, text=True)
            (dest / ("pass-" + str(number) + ".txt")).write_text(run.stdout + run.stderr)
            if run.returncode:
                raise RuntimeError("pdflatex failed in " + name + " pass " + str(number))
            log = (dest / (STEM + ".log")).read_text()
            bad = re.findall(r"^.*(?:Warning|Overfull|Underfull|Missing character|^!).*$",
                             log, re.MULTILINE)
            passes.append({"pass": number, "returncode": run.returncode,
                           "warnings_or_box_messages": bad})
            if number == 3 and bad:
                raise RuntimeError("final pass has warnings: " + repr(bad))
        pdf = dest / (STEM + ".pdf")
        info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
        page_count = int(re.search(r"^Pages:\s+(\d+)", info, re.MULTILINE).group(1))
        if not 10 <= page_count <= 13:
            raise RuntimeError("page count outside manuscript scope")
        builds.append({"directory": name, "passes": passes, "pages": page_count,
                       "pdf_bytes": pdf.stat().st_size,
                       "pdf_sha256": sha256(pdf.read_bytes()).hexdigest(),
                       "final_log_sha256": sha256((dest / (STEM + ".log")).read_bytes()).hexdigest()})
    if builds[0]["pdf_sha256"] != builds[1]["pdf_sha256"]:
        raise RuntimeError("fresh deterministic PDFs differ")
    final_pdf = HERE / (STEM + ".pdf")
    shutil.copyfile(HERE / "fresh-build-1" / final_pdf.name, final_pdf)
    subprocess.run(["pdftotext", "-layout", str(final_pdf),
                    str(HERE / "review" / "PDF_READBACK.txt")], check=True)
    receipt = {
        "status": "PASS_TWO_FRESH_DETERMINISTIC_BUILDS",
        "source_date_epoch": EPOCH, "force_source_date": "1", "timezone": "UTC",
        "source_sha256": sha256(source.read_bytes()).hexdigest(),
        "build_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "pdf_sha256": builds[0]["pdf_sha256"], "pages": builds[0]["pages"],
        "pdfs_byte_identical": True, "builds": builds,
        "toolchain": subprocess.check_output(["pdflatex", "--version"], text=True).splitlines()[0],
        "scope": "Compilation checks only; mathematical checks and visual review are separate."
    }
    (HERE / "COMPILATION_RECEIPT.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: receipt[k] for k in
                      ["status", "pdf_sha256", "pages", "pdfs_byte_identical"]}, sort_keys=True))


if __name__ == "__main__":
    main()
