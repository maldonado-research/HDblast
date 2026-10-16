# A certified first relative-harmonic gauge wave

16 October 2026. Conditional research checkpoint for the retained Einstein–scalar background and one quadrupole state. This checkpoint constructs the **first-order coordinate gauge vector** on a positive-time causal domain. It does not solve the second-order Einstein–scalar boundary problem.

The starting perturbation solves the retained linear equations in longitudinal coordinates. To use it in the relative-harmonic formulation, its gauge defect must be removed by an actual compatible vector wave. Specifying its initial acceleration alone does not construct that vector at later times.

For \(Y=n_1^2-n_2^2\), write the raw contravariant components as

\[
\xi^y=fY,\qquad \xi^t=bY,\qquad \xi^i=dD_\kappa^iY.
\]

The vector solves \((\Box_0+\mathrm{Ric}_0)\xi=-C_1^{\mathrm{long}}\), has zero initial value and velocity, and retains the shell conditions \(f=0\), \(b_y=d_y=0\). The tangential values remain free and reparametrize the shell. In scaled coordinates \(x=t/\ell\), \(r=z/\ell\), \(\ell=10^{-6}\), the certified domain is

\[
0\le x\le3,\qquad r\ge0,\qquad x+r\le8.
\]

The coupled physical half-line retarded operator is a contraction in a declared norm of values and first derivatives, with factor below **0.011412**. Its convergent series uniquely constructs the actual vector. Compatible differentiated iterates converge through third derivatives; boundary jump and corner terms are retained. Coefficients are bounded on the whole collar by a fresh first-exit calculation, rather than importing bounds that apply only to the initial annulus.

An explicit auxiliary free-wave formula has a certified continuous residual and error bound relative to this actual forced vector. Its radial-amplitude error is below **\(5.47\times10^{-9}\)** throughout the triangle. It is used as an approximation, with the retained master’s initial-velocity, boundary and potential contributions included in the error. The actual master remains the solution of its registered integral equation. Five actual signed point enclosures include a point on the reflected characteristic; these signs concern the coordinate gauge amplitude.

The continuous raw-amplitude bounds are

\[
|f|<7.47\times10^{-8},\qquad |b|<1.14\times10^{-12},\qquad
|d|<7.49\times10^{-20}.
\]

Each first scaled derivative is bounded by \(1.14\times10^{-6}\), \(7.59\times10^{-13}\), and \(5\times10^{-20}\), respectively. The separate higher-jet certificate bounds all mixed second and third derivatives. A physical \(t,z\) derivative of total order \(j\) divides the corresponding scaled bound by \(\ell^j\). A higher proper-\(y\) derivative also differentiates its background conversion factor; it is not obtained by simply dividing by a constant radius.

The transformed first fields are \(h_1^{\mathrm{wave}}=h_1^{\mathrm{long}}+\mathcal L_\xi g_0\) and \(\varphi_1^{\mathrm{wave}}=\varphi_1^{\mathrm{long}}+\mathcal L_\xi\varphi_0\). Their linear relative-harmonic defect vanishes by the vector equation. This coordinate transformation preserves the geometric first-order shell laws by pullback. Coordinate amplitudes and their bounds do not describe an invariant force or a physical second response.

## Proof and verification map

- `operator/GAUGE_OPERATOR_AND_CONSTRUCTION.md`: geometric operator, angular closure, scaled matrices and shell rows.
- `coefficients/README.md`: global coefficient/source intervals and their derivatives, with exact rational checks.
- `regularity/MIXED_WAVE_REGULARITY.md`: half-line kernels, boundary traces, contraction and actual C3 convergence.
- `solver/GAUGE_WAVE_SOLUTION_AND_DEFECT.md`: explicit approximation and continuous residual-to-error estimate.
- `higher-norms/`: quantitative bounds for all mixed gauge-vector derivatives through order three.
- `fresh-review/`: separate internal mathematical review and independent finite checks.
- `literature/LINEAR_GAUGE_WAVE_PRIMARY_REVIEW.md`: bounded primary-reference comparison and applicability limits.

Finite algebra checks and deliberately wrong controls supplement the all-point analytic arguments. They do not replace their hypotheses. The separate review is internal AI-assisted review, not external specialist peer review.

## Reproduce this checkpoint

The source uses Python's standard library. `REPLAY_ROSTER.json` declares the exact programs, execution modes, arguments and complete expected reports. `reproduce.py` authenticates the full source inventory before and after running that roster and writes results outside the source directory.

Obtain the manifest SHA256 from the separately published `hdblast/continuation/2026.10.16-first-gauge-wave/PUBLICATION_PINS.json` or the dedicated CI workflow. Then run `python -I -S -B reproduce.py --manifest-sha256 THE_EXTERNALLY_AUTHENTICATED_HASH --output-dir /tmp/hdblast-gauge-replay`. Choose a new output directory. The separate publication verifier also authenticates the complete ZIP, its member bytes and the runner before extracting and replaying it. These checks replay the new checkpoint while treating the authenticated earlier background/master results as declared premises; their full proofs remain in the linked earlier public checkpoints.

## Inherited premises and remaining work

The result assumes the retained regular background branch and the smooth actual first-master/contact solution. Their authenticated reports and the earlier [interior master proof](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_INTERIOR_MASTER_JETS_20261015) and [initial-data/gauge checkpoint](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015) identify those premises explicitly. Enclosing interval cells does not define a differentiable solution by arbitrary choices inside the cells.

Earlier checkpoints certified the first master’s interior derivatives and constructed compatible second-order physical initial data. The remaining task is to transform the complete first metric–scalar pair and its quadratic bulk and shell sources consistently, then solve and certify the compatible second-order evolution. A first-order gauge vector alone does not complete that task. No finite-amplitude nonlinear family or nonlinear remainder bound is established here.

This is conditional mathematics. It does not establish a higher-dimensional origin of the Big Bang, observational evidence, mathematical priority, or journal readiness. Earlier frozen public capsules and the manuscript remain preserved. No journal submission or new Zenodo version accompanies this checkpoint.
