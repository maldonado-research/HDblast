# Integrated technical continuation

The reviewed technical draft is
`COVARIANT_FLUX_AND_PHYSICAL_DOMAINS.md`, with its synchronized standalone
LaTeX and compiled 10-page PDF. It is frozen within this capsule for GitHub distribution; it has not been submitted to APS.

The candidate includes the exact current and mixed-spectral shell balance,
conditional registered cone limit, exact flat smooth-flux comparator,
radial-versus-inertial energy comparator, tensor smoothness/completion
distinction, and local smooth-gauge vertex-charge estimate. The companion
source reports and prior-art register provide detailed derivations and checks.

From this directory, regenerate the LaTeX with:

```sh
pandoc COVARIANT_FLUX_AND_PHYSICAL_DOMAINS.md \
  --from=markdown+tex_math_single_backslash --standalone --to=latex \
  --variable fontsize=11pt --variable geometry=margin=0.9in \
  --variable colorlinks=true --include-in-header=preamble.tex \
  --output=COVARIANT_FLUX_AND_PHYSICAL_DOMAINS.tex
pdflatex -interaction=nonstopmode -halt-on-error COVARIANT_FLUX_AND_PHYSICAL_DOMAINS.tex
pdflatex -interaction=nonstopmode -halt-on-error COVARIANT_FLUX_AND_PHYSICAL_DOMAINS.tex
```

The generated `.tex` contains the preamble directly and is independently
compilable with `pdflatex`. It includes `../figures/cone-limits.pdf` if
present, and otherwise displays a clearly labeled placeholder.

The final candidate build has no LaTeX warnings, errors, overfull boxes or
underfull boxes. Extracted text confirms the author, affiliation, ORCID,
email and technical-draft status. All ten rendered pages were subsequently inspected by the root integration agent; see `../reviews/VISUAL_REVIEW_RECEIPT.json`.

Candidate hashes and build properties are in `BUILD_STAGE_RECEIPT.json` (the historical build-stage status precedes the final capsule freeze). Any
subsequent edit requires new hashes and affected review/build checks.
