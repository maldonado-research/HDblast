# Internal audit and execution record

This review is internal and AI-assisted; it is not external peer review.

The first numerical protocol was committed at 99dd818755c493127ae4b4e0cad7bbf21fd141de before any quantum-mode integration in this round. Classical inputs, earlier Gaussian estimates and crossing locations were already known and disclosed. Code first ran at c71054d10c4ecd6fd993251946a0e9d02119903b.

[GitHub Actions run 36938696104](https://github.com/maldonado-research/HDblast-archive/actions/runs/36938696104), job 110624900993, completed successfully. Source compilation, pinned input SHA-256 checks, ten registered mode variants, curve/figure export and artifact upload completed. The fetched complete JSON has zero case errors and all registered control flags true. Its source SHA-256 is recorded in outputs/actual_modes.json.

Independent reviewers checked the action variations, current sign, clock mappings, backward-state work rebasing, Bogoliubov signs, coherent factors, endpoint energy relation and conformal coupling. No blocking algebra/sign error was identified. The work-integral independence and path-versus-endpoint distinctions were explicitly reviewed; all reported path residuals also fall below tolerance.

Root independently executed 162 local stress/current/coherence cases in V8, with maximal stress/current absolute discrepancy 4.55e-13, plus four conformal moving-cutoff controls. Those are algebraic controls, not the actual mode experiment.

A separate RK4 implementation then integrated four selected actual-history modes at 100k/200k steps, and all 32 modes at 200k steps. Source and outputs are supplied. Endpoint energy, pressure and scalar-current reconstruction agree with DOP853 within 6.1e-8 in the declared units; occupation agrees within 1.9e-11. This checks phase-sensitive observables as well as occupation.

The endpoint WKB projection, full 32-node RK4 comparison, heat-kernel matching discussion and illustrative flat-space loop numbers are post-registration exploratory additions. They do not retroactively change the primary protocol or its ten cases.

The endpoint work normalizer is |E(t)|+|E(initial)|+|scalar work|+|pressure work|+1e-30, with quantities appropriately differenced/rebased for the pair. It is not the integral of absolute power. These are auxiliary quadratures sharing the mode integrator, not separate numerical quadrature of the final saved curves.

The source uses fixed comoving weights; no moving-cutoff energy is silently injected. The state pair and metric are the same within each case. Coarse-background and alternative-spline cases are comparisons between different reconstructed geometries, not state subtractions across them.

An independent EFT audit confirmed the missing quartic potential and intrinsic-curvature operators. Flat MS-bar illustrations are scheme-dependent vacuum matching pieces, not radiation predictions, and they use the older A1 endpoint explicitly identified in their output.

The final public runner additionally replays the independent JavaScript controls and validates the curated download. Its actual status and link are recorded in the publication pull request, without replacing the pinned first-run evidence above.
