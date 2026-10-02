# E1: curved-source readiness and ultraviolet regularity

Registered 1 October 2026 in America/Los_Angeles (2 October UTC, used in the directory name). This protocol is committed before executing this round's archive audit or oscillator controls. Earlier finite-band and flat-pulse results are known. The analytic possibility of interpolation/state ultraviolet artifacts has already been derived; this is not a blind discovery claim.

## Inputs and scope

Use the exact fine and coarse source-free tuned A1 archives already published in HDBLAST_CHECKPOINT_20261001_QFT: delta=0.1,Y=0,G=100,phi_star=0.5, s=H0*tau, with scale factor a(0)=1. Their byte hashes are preserved in INPUTS.json. The original archive commit is 8f67197b730d4e1c43554b86f224c29cc72629eb; the inherited public checkpoint is pinned at 4d81e28a8f1ffb7fa6fff37316bfc6837a62c7c1. Use the prior loader's reliable-data and duplicate-time rules, without replacing bad values or extrapolating the solution.

This is a readiness audit of the reconstructed background and state, not a renormalized curved-shell source, new five-dimensional evolution, or physical rejection/confirmation of HDBLAST. The original delta=0.001 model and previous finite-band verdicts retain their scope.

## Analytic quantities to audit

For the minimally coupled scalar, u=a chi obeys u''+[k²+U]u=0 with U=a²[x-Hdot-2H²], x=G²(phi-phi_star)². Primes denote conformal time and dots s derivatives. The C1 cubic-Hermite reconstruction of ln a has potential jumps DeltaU=-a² DeltaHdot. For an admissible incoming ultraviolet state and otherwise smooth evolution, a knot contributes beta=DeltaU/(4k²) times a phase. Distinct interior knots yield the logarithmic particle-energy coefficient sum(DeltaU²)/(32pi² a_f⁴). Report this coefficient, all jumps and dominant knots on the reliable interval ending at s=6.9, with earlier-window information when useful.

Compute one-sided polynomial limits; do not estimate jumps with a macroscopic finite difference. Check their equivalent closed-form Hermite expressions independently. State/basis boundary mismatches are separate: the prior amplitude-WKB initial data omit a''/a and generically give a k^-2 beta tail; the physical-Hamiltonian reference generically gives k^-1 beta. Report their leading coefficients separately from knot effects. The latter is a basis/state completion issue, not a measured particle bath.

Audit the existing C2 cubic spline and two smoother interpolating reconstructions: quintic C4 and septic C6 B-splines. Report interpolation fidelity to stored values and derivative fidelity to stored H and phi_dot, plus sensitivity of U, curvature and curvature-squared stress ingredients. Do not classify a smoothing choice as a physical solution or a Hadamard completion. Any additional fit/grid choice must be recorded in a pre-execution addendum or explicitly labeled exploratory.

## Independent oscillator controls

Use the separate smooth frequency step U(eta)=1+1.5[1+tanh(eta/epsilon)] with epsilon in {0.1,0.3,1} and k in {0.5,2,8}. Exact asymptotic occupation is
n=sinh²[pi epsilon(omega_plus-omega_minus)/2]/[sinh(pi epsilon omega_plus)sinh(pi epsilon omega_minus)],
omega_minus=sqrt(k²+1), omega_plus=sqrt(k²+4).
Integrate from -18epsilon to +18epsilon with the asymptotic incoming mode. Compare a complex DOP853 implementation against an independently written four-real-component Radau implementation with analytic Jacobian. Primary tolerances: rtol=2e-11, atol=2e-13, max_step=min(epsilon/20,0.2/omega_plus). Required maximum occupation absolute error 1e-9; relative error <=1e-4 when n_exact>1e-9; Wronskian error <=1e-8. Report all nine cases, including tiny occupations where relative accuracy is meaningless. Record finite-endpoint preparation limits.

The abrupt control uses exact frequency-step matching, n=DeltaU²/[4omega_plus omega_minus(omega_plus+omega_minus)²], and checks k^4 n tending to DeltaU²/16. This is an algebraic asymptotic control; it is not a full archived high-k mode integration. Smooth-step UV suppression is assessed against its exact formula, with no claim to have evolved the entire archive at unbounded k.

## Decision and reporting

A resolved nonzero Hermite curvature jump blocks promoting that surrogate unchanged into a continuum absolute-source calculation. A nonzero leading initial-state ultraviolet mismatch is a separate blocker. Smoother interpolation alone cannot clear both requirements. Readiness remains unestablished without common-action finite curvature matching, controlled background derivatives and admissible state preparation.

Separate algebraic/control success from the physical readiness verdict. Keep numerical failures, source/input hashes, all cases and measured reconstruction sensitivities. Do not relax thresholds after observing results. Internal AI-assisted reviews are not external peer review; no novelty or Nobel-level claim follows from the audit.
