# Build the technical draft

From the capsule root, with Pandoc and TeX Live available:

```sh
pandoc manuscript/SHELL_NOETHER.md --standalone --from markdown+tex_math_single_backslash --to latex --include-in-header=manuscript/preamble.tex -V geometry:margin=0.9in -V fontsize=10pt -o /tmp/SHELL_NOETHER.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp /tmp/SHELL_NOETHER.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp /tmp/SHELL_NOETHER.tex
```

The packaged PDF has seven pages and31 numbered equations. Its exact bytes are recorded in BUILD_AND_REVIEW.json and the capsule manifest. A new build need not be byte-identical because PDF creation timestamps and installed typesetting tools can differ. Algebra replay does not require a TeX installation.
