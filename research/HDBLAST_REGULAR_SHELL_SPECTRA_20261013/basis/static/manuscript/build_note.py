#!/usr/bin/env python3
"""Render the frozen mathematical note with Pandoc and compile with XeLaTeX.

Run from any directory with Python 3, Pandoc, and XeLaTeX installed.
The proof Markdown is authenticated and never modified. The resulting TeX is
self-contained and can subsequently be compiled without Python or Pandoc.
"""
from __future__ import annotations

import hashlib
import json
import os
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'theory' / 'JACOBI_ENDPOINT_MONOTONICITY.md'
EXPECTED = 'f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8'
STEM = 'scalar-response-monotonicity'

PREAMBLE = r'''\documentclass[11pt]{article}
\usepackage[letterpaper,margin=0.85in]{geometry}
\usepackage{fontspec}
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont[Scale=MatchLowercase]{DejaVu Sans Mono}
\usepackage{amsmath,amssymb,mathrsfs}
\usepackage{longtable,booktabs,array,calc}
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue!45!black,urlcolor=blue!45!black,
 pdftitle={Regular-cone Jacobi identities and controlled scalar boundary response},
 pdfauthor={Ricardo Maldonado},
 pdfsubject={Static Einstein-scalar branch monotonicity and exact sufficient bounds},
 pdfkeywords={Einstein-scalar, de Sitter brane, Jacobi field, static boundary problem}}
\setlength{\emergencystretch}{3em}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.45em}
\setcounter{secnumdepth}{-1}
\allowdisplaybreaks[2]
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\title{Regular-cone Jacobi identities and\\controlled scalar boundary response}
\author{Ricardo Maldonado\\Independent Researcher\\
\small\href{https://orcid.org/0009-0009-3937-6527}{ORCID: 0009-0009-3937-6527}}
\date{11 October 2026}
\begin{document}
\maketitle
\begin{abstract}
We study a static Einstein--scalar boundary problem for two identical regular
interiors joined at a de Sitter brane. A regular-cone Jacobi identity gives an
exact curvature factor in the derivative of the geometric junction residual.
The existing weighted contraction conditions imply strict bulk convexity;
this in turn proves strict residual monotonicity and removes a previously
restrictive sufficient inequality for local branch uniqueness. A Riccati
comparison sharpens the endpoint estimate. For the specified exact-decimal
polynomial model, exact arithmetic verifies all remaining sufficient conditions
for every $0<\delta\leq0.002$, together with explicit scalar and curvature
remainders and lower curvature relative to a separate model with
scalar-independent detuning. Uniqueness is restricted to the stated tube and
$T\geq1/k$. The note includes all conditions and the mathematical proof.
\end{abstract}

\textbf{Research and review status.}
This is an AI-assisted mathematical research note with internal mathematical
and computational review; it has not received external peer review.
OpenAI Codex assistance was used for derivation, checking, programming, and
drafting. The author attribution does not by itself assert a personal rerun of
every computation. The results are conditional static mathematical theorems
under the explicit hypotheses below. They do not establish dynamical stability,
a cosmological origin mechanism, or evidence for a higher-dimensional origin
of the Big Bang. No global priority or novelty claim is made.

\textbf{Source preservation.}
The accompanying \href{../theory/JACOBI_ENDPOINT_MONOTONICITY.md}{frozen proof source}
in the capsule's \texttt{theory} directory is the
exact frozen proof source. This typeset note preserves its complete mathematical
body, including all eleven sufficient inequalities and the registered-model
corollary. The source hash and reproducible compilation receipt accompany the
TeX and PDF.
'''


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args, *, input_text=None):
    environment = os.environ.copy()
    environment['SOURCE_DATE_EPOCH'] = '1791676800'  # 2026-10-11 00:00:00 UTC
    return subprocess.run(args, cwd=HERE, input=input_text, text=True,
                          capture_output=True, check=True, env=environment)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--diagnostics-dir', type=Path,
                        help='Optional directory for compilation transcripts and the engine log.')
    args = parser.parse_args()
    if digest(SOURCE) != EXPECTED:
        raise ValueError('Frozen mathematical source SHA256 mismatch')
    source_text = SOURCE.read_text()
    heading, body = source_text.split('\n', 1)
    if heading != '# Regular-cone Jacobi identities and a stronger uniform branch estimate':
        raise ValueError('Unexpected proof heading')
    fragment = run(['pandoc', '--from=markdown+tex_math_single_backslash',
                    '--to=latex', '--wrap=preserve', '--shift-heading-level-by=-1'],
                   input_text=body).stdout
    # Give repository paths legal line-break points, without changing their text.
    def linked_path(match):
        path = match.group(1).replace(r'\_', '_')
        return r'\href{' + path + r'}{\nolinkurl{' + path + '}}'
    fragment = re.sub(r'\\texttt\{(\.\./[^{}]+)\}', linked_path, fragment)
    source_tags = re.findall(r'\\tag\{([^}]*)\}', source_text)
    if re.findall(r'\\tag\{([^}]*)\}', fragment) != source_tags:
        raise ValueError('Equation tags changed during conversion')
    math_counts = {}
    for label, pattern in [('display', r'\\\[(.*?)\\\]'),
                           ('inline', r'\\\((.*?)\\\)')]:
        original = re.findall(pattern, source_text, re.S)
        converted = re.findall(pattern, fragment, re.S)
        if original != converted:
            raise ValueError(f'{label} mathematical payload changed during conversion')
        math_counts[label] = len(original)
    tex = HERE / (STEM + '.tex')
    tex.write_text(PREAMBLE + '\n' + fragment + '\n\\end{document}\n')
    passes = []
    pdf = HERE / (STEM + '.pdf')
    if args.diagnostics_dir is not None:
        args.diagnostics_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='scalar-note-latex-') as scratch:
        for index in (1, 2):
            result = run(['xelatex', '-no-shell-escape', '-interaction=nonstopmode',
                          '-halt-on-error', '-output-directory', scratch, tex.name])
            output = result.stdout + result.stderr
            if args.diagnostics_dir is not None:
                (args.diagnostics_dir / f'compile-pass{index}.log').write_text(output)
            passes.append(dict(pass_number=index, exit_code=result.returncode))
        log = (Path(scratch) / (STEM + '.log')).read_text()
        if args.diagnostics_dir is not None:
            (args.diagnostics_dir / (STEM + '.log')).write_text(log)
        shutil.copyfile(Path(scratch) / pdf.name, pdf)
    receipt = {
        'status': 'COMPILED_PDF_PASS',
        'source_relative_path': str(SOURCE.relative_to(HERE.parent)),
        'source_sha256': EXPECTED,
        'tex_sha256': digest(tex),
        'pdf_sha256': digest(pdf),
        'build_script_sha256': digest(Path(__file__)),
        'equation_tag_count': len(source_tags),
        'equation_tags_preserved': True,
        'all_mathematical_payloads_preserved_exactly': True,
        'math_payload_counts': math_counts,
        'source_date_epoch': 1791676800,
        'typesetting_only_adjustments': [
            'The source H1 is represented by the document title; H2 headings become sections.',
            'Repository paths use breakable monospace formatting and relative PDF links without changing their text.',
        ],
        'passes': passes,
        'pandoc_version': run(['pandoc', '--version']).stdout.splitlines()[0],
        'xelatex_version': run(['xelatex', '--version']).stdout.splitlines()[0],
        'overfull_boxes': re.findall(r'Overfull \\[hv]box[^\n]*', log),
        'undefined_references': bool(re.search(r'undefined references|Reference .* undefined', log)),
        'missing_characters': re.findall(r'Missing character:[^\n]*', log),
        'scope': 'Typesetting validation, not an independent mathematical proof or external peer review.',
    }
    (HERE / 'COMPILATION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: receipt[key] for key in ('status', 'equation_tag_count',
          'overfull_boxes', 'undefined_references', 'missing_characters')}, indent=2))


if __name__ == '__main__':
    main()
