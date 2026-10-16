# Independent internal manuscript review

Accepted with no remaining blocking finding. This review binds the final seven-page manuscript and its numerical table to the accepted geometric, projection, actual-source and shell records. It is an internal mathematical/computational review, not external specialist peer review or an originality assessment.

The final TeX has SHA256 `b859ff182e93689fbe8efc1fdc8fd7ec089a944660edaf68278ae1afc5b441cb`; the table has SHA256 `28f049a3ab1aa2baa26e82b2967a233a9c1f18ac476ab35e662faf9537eaaf6a`; the seven-page PDF has SHA256 `afa359b73629a075f37252d301f2ac1212965ba9b2d9ca31b19d9832ec7c50dd` and 298590 bytes. Every PDF page was visually inspected and the extracted text was read. Equations, table, theorem statements and references are legible and complete.

## Mathematical scope

The Ricci Hessian, matter normalization, second-parameter-derivative convention and scalar equation agree with the independent derivation in `review/INDEPENDENT_DERIVATION.md`. The actual metric-lowered gauge covector is distinguished from background lowering. The off-shell reduced Hessian includes the `2 P C1` term, and the scalar off-shell terms are displayed. Their harmonic specialization is applied to the actual constrained field. The frame amplitude formulas and the forcing pair `(-2 ell^2 rho^2 QEred, +ell^2 rho^2 Qphi)` agree with the retained conventions.

The angular assertion uses exact polynomial reconstruction and symbolic independence, with 20 metric and 3 scalar amplitudes. The complete 81-cell causal-domain cover is distinguished from the 64-cell radial profile and 6-cell shell cover. The three displayed table intervals are outward enclosures of their exact recorded endpoints. The shell table entry is explicitly the one-copy outward-normal scalar junction Hessian. The full forcing norm bound is 17387513156573.

For the intrinsic shell scalar, the material pullback and joint field/embedding transformation are explicit. The vanishing background tangential gradient makes its zeroth and first coefficients vanish, so its quadratic coefficient is invariant under the stated gauge and relabelling transformations. Its angular mean uses the fixed background time foliation. The monotonicity argument establishes exactly two simple mean zeros in each named bracket; it does not establish a global count of four. The signs are signs of a Lorentzian squared gradient at quadratic order, not of an invariant force.

The evolved scalar-mean corollary follows from the retained conditional derivative bound and authenticated initial mean intervals. The mean functional has norm one. On the complete stated rectangles, the drift is less than 5.15e30 and the strict sign margin is 1.85e30. The intermediate-value consequence supplies at least one radial mean zero separately at every stated time. It asserts neither uniqueness nor a continuous zero branch, and is not a new waveform calculation or a gauge-invariant observable.

The manuscript accurately retains the conditional smooth background/master/gauge premises. It does not claim a numerical second-response residual, a finite-amplitude nonlinear family, a physical force, an origin theorem, or external novelty. It cites the established curvature and gauge framework and states the remaining residual-validation work.

## Corrections resolved before acceptance

The final version describes background radial jet enclosures through order four; identifies the table's shell Hessian normalization; defines the base triangle separately from its product with the sphere; and defines the full forcing pair and its norm. These were wording and notation corrections without changes to the accepted source maps or interval calculations.

## Replay and evidence

`check_manuscript_arithmetic.py` independently checks the exact table endpoints, cover counts, all 64 negative profile means, accepted forcing norm, response margin arithmetic, time conversion and radial-zero scope. Isolated standard-library normal and optimized runs produced byte-identical reports with SHA256 `0a8fb6a1cc50e4ee0be364c74d1e5b8839932084ddb569b34391c8f54413cef7`. The checker binds the exact TeX, table, PDF and referenced records. Its arithmetic replay does not replace the formula/proof reading or seven-page visual review recorded here.

`MANUSCRIPT_REVIEW_RECEIPT.json` records the exact accepted bindings. `PUBLIC_FILES.json` is the complete public allowlist for this directory, and `REPLAY_ROSTER.json` supplies the two portable positive replays. There are no direct or nested negative controls in this manuscript arithmetic checker. Existing geometric and numerical control receipts remain separately applicable. The sealed `review/` files and all producer files were unchanged by this manuscript review.
