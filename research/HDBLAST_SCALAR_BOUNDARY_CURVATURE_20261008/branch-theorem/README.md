# Regular branch existence and uniform curvature error

This continuation addresses the main mathematical assumption in the uploaded APS research note. The full argument is in `UNIFORM_BRANCH_THEOREM.md`. It constructs a regular doubled-interior branch for every sufficiently small positive detuning, proves uniqueness within a stated weighted neighborhood, and bounds both scalar and curvature remainders uniformly as the shell radius diverges.

The proof uses a bounded weighted inverse on every finite interval, a Volterra equation for the metric, and an ordinary implicit-function argument only at positive detuning. It does not apply an implicit-function theorem at an infinite-radius solution.

The proof snapshot's final section records its status when written. Subsequently `../independent-math/UNIFORM_PROOF_REVIEW.json` independently checked 22 implications and accepted the exact source SHA256 `544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51`. This is an internal mathematical review, not external peer review or a novelty determination.

`../numerical-audit/uniform-constants/EXACT_UNIFORM_CONSTANTS.json` evaluates all 12 sufficient inequalities with exact rational arithmetic. Its polynomial and transcendental outer bounds make the positive radius explicit. For the registered model with the stated decimal parameter interpreted as an exact rational, the certified interval reaches approximately `1.0915542913337286e-7`. This is roughly 4,581 times smaller than the smallest detuning in the uploaded numerical sweep. The existence proof therefore does not turn those 48 floating-point configurations into certified finite-detuning solutions.

For the simple exact model

    k=1, W(eta)=3+(5/2)eta^2, f(eta)=1+eta,

the sufficient inequalities hold through `delta=1/1000`. It is a separate stated model and was not one of the uploaded benchmark families.

`derive_curvature_suppression.py` reads the exact constant certificate and derives the sign consequence recorded in `CURVATURE_SUPPRESSION_COROLLARY.json`. Throughout each model's entire certified interval,

    H^2 < k f0 delta/3 + f0^2 delta^2/36.

The reference on the right is the exact constant-scalar result with the scalar boundary coupling removed while its value at the stationary point is retained. For nonzero coupling it generally fails the scalar junction of the original model. The corollary compares static curvatures and supplies neither stability nor a time-dependent relaxation mechanism.

These results can strengthen a focused mathematical paper. They do not establish external originality, APS's editorial threshold, a physical energy scale, reheating, an observational prediction, or a higher-dimensional origin of the Big Bang.
