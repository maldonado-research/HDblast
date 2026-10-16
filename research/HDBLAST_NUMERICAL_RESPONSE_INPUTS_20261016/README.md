# Certified numerical inputs for the selected formal second-order shell response

Checkpoint 20261016 UTC. This continues the conditional formal construction in
[HDBLAST_COMMON_GAUGE_FORCING_20261016](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_COMMON_GAUGE_FORCING_20261016).
The same retained smooth regular background, actual first solution, aligned
quadrupole `Y=n1²−n2²`, shell law, gauge and second initial prescription are used.
Earlier manuscripts and frozen capsules are unchanged.

The new calculation supplies:

* **Actual continuous approximation errors through scaled order three** for
  explicit first-master and gauge comparisons over `0<=x<=3,r>=0,x+r<=8`.
  These come from differentiated true integral/PDE identities, retaining the
  master dynamic boundary, reflection terms and nonzero stabilized initial jets.
* **A shared affine background approximation with scaled radial C4 errors**,
  retaining the inherited regular-branch shell uncertainty, plus whole-cell
  enclosures for the actual second initial profiles and their full one-jet.
* **A complete exact 23-amplitude homogeneous bulk/boundary operator**, with
  achieved harmonic-gauge, scalar, matter, shell and initial source projections.
  Both TT features and the scalar `DC[k]` coupling are retained.
* **The sharp selected-space pointwise constant `sqrt(321/19)`**, approximately
  `4.11032`, proved by exact sphere moments and global polynomial bounds.
* **A refined global majorant `||U||X<1.03e33`**, about 990 times smaller than
  the former rounded upper bound. This is not an evaluated response amplitude.

The full numerical second-order waveform remains unfinished. The geometric
quadratic Ricci/reduced metric-source adapter, source approximation-error
propagation, evaluated integrals and continuous second-response residual still
require computation. No exact nonlinear family, finite-amplitude remainder,
invariant force or evidence for a higher-dimensional Big Bang origin is claimed.

## Readable derivations

The five-page [technical working note](manuscript/NUMERICAL_INPUTS_NOTE.pdf) has
[editable LaTeX source](manuscript/NUMERICAL_INPUTS_NOTE.tex). It is a separate
computational note, not a replacement for the frozen shell-Noether manuscript
or a journal submission.

* `numeric/APPROXIMATION_ERROR_PROOF.md`: actual first-field error hierarchy.
* `numeric/INITIAL_PROFILE_PROOF.md`: continuous interval evaluation and the
  shared background approximation.
* `angular/EXACT_ANGULAR_OPERATOR.md`: exact homogeneous operator matrices.
* `angular/source-audit/SOURCE_INVENTORY.md`: achieved source identities and
  the explicit remaining metric-source adapter.
* `angular/SELECTED_SPACE_CLOSURE.md`: natural quadratic origin and Schur closure.
* `evaluation-bound/EVALUATION_BOUND.md`: sharp pointwise evaluation proof.
* `response/REFINED_RESPONSE_BOUND.md`: refined norm join and exact tail counts.

Independent internal reviews are in `numeric-review/` and `angular-review/`.
They are AI-assisted computational reviews, not external specialist peer review.
The literature comparison records specifically inspected primary results and
their limits; it does not establish external novelty.

## Quantities and conventions

Use `ell=1e−6`, `x=t/ell`, `r=z/ell`, `dy=−rho dz`; expansions use `epsilon²/2`.
The second field is `U=(k,w)` in the accepted parallel frame and

    ||U||X = max(sup|U|H/3, sup|U_x|H, sup|U_r|H).

H uses probability-normalized sphere L2 and a positive Frobenius metric norm.
The accepted contraction is unchanged: `q=0.002109`. The complete angular space
still has 287 sections. For this selected quadratic state, equivariance places
the formal response in E_v, which has 23 amplitudes (5 Dirichlet/18 Robin).
The sharp bound applies to E_v, not to arbitrary fields or to a candidate/error
that leaves this subspace. The universal 287-space factor 17 is unchanged.

The old `1.02e36` and new `1.03e33` X bounds concern the same selected formal
coefficient. Refining a bound does not change the solution. Using `4.111` as a
strict rational upper bound on the selected evaluation constant gives

    sup|U|F < 1.270299e34,
    sup|U_x|F, sup|U_r|F < 4.23433e33.

These large conservative numbers reflect the narrow packet and selected local
initial prescription. No nonlinear remainder or physical instability follows.
The initial prescription is not unique. A graph of its metric/scalar angular
means is in `manuscript/figures/INITIAL_PROFILE_ENCLOSURES.pdf` (and PNG).
Bands enclose entire radial cells; floating-point plotting is display-only.
The scalar mean has signed actual witness intervals with magnitude exceeding
`7e30`; six alternating signed cells prove at least five distinct initial zeros.
These are initial coefficient zeros, distinct from the earlier coordinate shell
source zeros and from an evolved physical observation.

## Reproduction and scope of verification

The source archive includes all executable scientific inputs used here, with
byte/SHA-256 pins. Inherited certificates are authenticated declared premises;
their PDE existence proofs are not silently re-proved by these new programs.
Every acceptance comparison uses exact Fraction arithmetic or outward rational
intervals. Display-only decimals do not decide acceptance.

`MANIFEST.json` closes the full source inventory. `REPLAY_ROSTER.json` explicitly
lists each positive normal/optimized process and direct rejection process, its
arguments, complete golden report and specific rejection reason. The runner
checks the source before/after, compares complete typed JSON and records raw
byte agreement. Wrong-formula/input processes must fail for their declared
reason and leave no output. Nested controls are recorded separately; process
counts are not theorem counts. The roster has 25 top-level positive
processes (24 scientific normal/optimized runs and one corruption harness),
104 direct rejection processes, and 32 nested rejections inside positive jobs:
16 initial-algebra controls, 14 source-adapter controls, and two isolated
changed-master-input checks. Historical aggregate control reports are
authenticated provenance and are not added again to these fresh counts.

From a checkout of the repository, use the separately pinned verifier:

    python -I -S -B hdblast/continuation/2026.10.16-numerical-response-inputs/verify_and_reproduce.py --repo-root . --pins hdblast/continuation/2026.10.16-numerical-response-inputs/PUBLICATION_PINS.json --output-dir /tmp/hdblast-numerical-inputs-fresh-check

Choose a new external output path. The verifier authenticates source and ZIP,
checks all member bytes, extracts a fresh copy and executes the registered
scientific roster. All acceptance programs need only Python's standard library.
Rebuilding the optional plot uses Matplotlib; rebuilding the optional PDF uses
pdfLaTeX. Those artifacts are included and byte-pinned, but typesetting/plot
reproducibility is separate from the standard-library mathematical replay.

## Next numerical target

Transform the certified first comparisons and shared background approximation
into the common-gauge metric/scalar two-jets; propagate their correlated errors
through every quadratic source, complete the metric-source amplitude adapter,
and evaluate the causal iterations in E_v. Certify the complete continuous
integral residual, including initial, shell, reflection and interpolation errors.
The retained estimate is

    ||U−V||X <= ||V−b−TV||X/(1−q).

Exact Neumann tail counts alone do not certify computed integrals. A declared
fixed-gauge scalar mean can be a numerical diagnostic; an invariant shell
observation needs its own geometric definition and error estimate.

No new Zenodo version or journal submission accompanies this checkpoint.
Private editorial history and raw web retrievals remain outside the supplement.
