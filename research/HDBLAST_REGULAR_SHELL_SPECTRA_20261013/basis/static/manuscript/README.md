# Standalone technical note

`scalar-response-monotonicity.pdf` is the 10-page technical note by Ricardo
Maldonado, Independent Researcher, ORCID 0009-0009-3937-6527, dated 11 October
2026. The document explicitly identifies AI-assisted mathematical research and
internal computational review; it does not claim external peer review,
dynamical stability, cosmological evidence, or global novelty.

The complete mathematical body is derived from the byte-for-byte frozen source
`../theory/JACOBI_ENDPOINT_MONOTONICITY.md`, SHA256
`f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8`.
The renderer verifies exact preservation of all 63 display-math payloads,
240 inline-math payloads, and 31 equation tags. It changes heading levels and
allows line breaks in repository paths; it does not alter the mathematics.

To rebuild from the frozen Markdown, install Python 3, Pandoc, XeLaTeX, and the
DejaVu fonts, then run:

```sh
python3 manuscript/build_note.py
```

The command works from the capsule root; an absolute path to the script also
works from another directory. It authenticates the proof, creates the
self-contained TeX, and compiles twice in a temporary directory with shell
escape disabled. The document
date is fixed to 2026-10-11, and `SOURCE_DATE_EPOCH` fixes the PDF metadata date.

Alternatively, compile the already generated, self-contained TeX directly:

```sh
cd manuscript
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error scalar-response-monotonicity.tex
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error scalar-response-monotonicity.tex
```

`COMPILATION_RECEIPT.json` records source, TeX, PDF, and builder hashes, tool
versions, math-preservation checks, and compiler diagnostics.
`RENDERING_REVIEW.json` separately records text and page-layout checks of the
distributed PDF. These are typesetting checks, not external mathematical review.
To retain detailed engine diagnostics, add `--diagnostics-dir /path/to/output`
to the builder command. The default leaves no auxiliary files in the capsule.
PDF text extraction can reorder mathematical fractions and scripts; use the
PDF or TeX when checking formulas.
