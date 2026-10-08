# Independent review of the conditional half-space scalar model

8 October 2026. Verdict: **ACCEPT_CONDITIONAL_ANALYTIC_TESTBED** for the exact
source bytes pinned in `REVIEWED_SOURCE_HASHES.json`. There is no remaining
mathematical blocker in the declared scalar domain. This verdict is not a
physical-origin result, empirical validation, novelty judgment, or authorization
for any protected physical-source, trajectory, or likelihood execution.

The accepted source manifest SHA-256 is
`d82d09d1786b1317dbd5fb8852e618a4ecd63b26bb70b729e9be4f041a5bc15c`.
The theorem SHA-256 is
`f2d695e131252452b82fd060258ae808005b5abed78ff3b1ccda96a6eec0b384`.
The final contract SHA-256 is
`f3542522a33980cc0c8b9ca29ccc5dd9cad90c157b3f717c5382edafb14dbf8c`.
Six source/receipt files were checked against the manifest; their exact bytes
and the manifest are preserved read-only under `reviewed-source/`. Any later
source change requires a new comparison and receipt. The producer's inherited
scope-document list is retained as provenance, not independently re-audited
by this review.

## Accepted mathematical content

**Variation and stability.** With the bulk occupying y>=0 and signature +----,
integration of -Phi_y delta Phi_y produces +Phi_y(0) delta Phi_0. The stated
boundary equation, Phi_y(0)=c Phi_0-gq, and brane force +g Phi_0 follow from the
same action. The Robin contribution to the energy is +c Phi_0^2/2, and the
interaction contribution is -gq Phi_0. There is no second half-space factor.

The trace inequality is sharp: the transverse profile exp(-My) saturates it.
On that trial subspace the zero-momentum spatial form is
`(c+M) A^2 - 2g A Q + m^2 Q^2`. Its Schur complement is
`m^2-g^2/(c+M)`. This independently establishes both the proposed stable side
and sharpness of its boundary. The displayed relative-form bound with r<1
controls the trace coupling on the stated H^1 direct-sum domain. Equivalence
to the positive uncoupled form gives a closed positive form, its self-adjoint
operator, and energy-space wave evolution. Very slowly varying brane profiles
turn a negative transverse trial direction into a negative full spatial sector.
The equality case is correctly described as a generalized k=0 zero mode,
rather than a normalizable constant function on the infinite brane.

The accepted contract deliberately takes c>0, M>0, m^2>0 and g nonzero. It is
narrower than an optional c>=0 extension; no c=0 bound-mode classification is
accepted by silently taking expressions containing 1/c at face value.

**Retarded branch, causality and modes.** The decaying upper-half-plane solution
has Phi proportional to exp(-s y), Re s>0. At positive propagating frequency its
boundary value is s=-ip, producing exp(+ipy) with the declared exp(-i omega t)
convention. Thus Im G_R>0, and the self-energy subtracts g^2 G_R from the brane
inverse response. The positive spatial operator excludes upper-half-plane
instability poles. The Klein-Gordon spectral representation supplies causal
support in the brane light cone; analyticity alone would not be enough to
establish that spatial support, but the source supplies the stronger argument.

The bare c>0 Robin bath has no normalizable bound state. The coupled system may
have exactly one below threshold. D(z) is strictly decreasing there, D(0)>0,
and its threshold limit has the sign given in the source. The threshold-tuned
case is not an L^2 bound state. A coupled bound state carries a separate positive
delta spectral weight and cannot be erased by incoming continuum damping.

**Spectral measure and FDT.** The Robin eigenfunctions have boundary weight
`(2/pi) p^2/(c^2+p^2)` in dp. The transformation u=M^2+p^2 divides this by 2p,
yielding exactly the proposed rho_b(u), with no missing factor of two or pi.
For one canonical oscillator of energy E, the imaginary part of its retarded
response is `pi/(2E)[delta(omega-E)-delta(omega+E)]`; its symmetrized covariance
has `pi/(2E)coth(beta E/2)` multiplying the sum of those deltas. This verifies
the signs and factors in both noise formulas under the stated one-half
anticommutator convention.

The full coupled ground/KMS state supplies the bound oscillator's occupation.
The incoming eta noise remains a continuum object; any homogeneous bound
sector is added to q and its covariance. The continuum identity
`S_q=|chi_R|^2 N_eta` cannot reconstruct that missing delta sector. The final
contract now states this explicitly. The quantum state is understood as a
Gaussian state on the smeared observable algebra in infinite volume; the
result does not require a trace-class global Gibbs density matrix.

The |omega|^-1 force-noise tail gives a logarithmically divergent raw equal-time
variance even at fixed k. Smooth time smearing, as included in the final
spacetime-smearing premise, controls it. Spatial smearing alone does not. The
retarded Stieltjes integral still converges, and this model makes no finite
unsmeared stress or vacuum-energy assertion.

**Energy.** The conserved quantity includes bulk, brane, Robin and interaction
energies. Omitting the interaction term fails even an exact stable bound-mode
solution. The positive-y Noether flux is -Phi_dot Phi_y; its outgoing harmonic
average is positive and equals minus the brane self-force work. This response
describes driven brane energy escaping into the bulk. An independently prepared
incoming state or dynamical source would be needed to supply bulk-to-brane
energy; the equilibrium state premise does not derive it.

**Rival statements.** At fixed k, elimination from a finite polynomial matrix
in frequency gives a rational response. Finite derivative order, translation
invariance, a well-defined linear probe and rational readout are material
hypotheses. The square root has simple zeros at distinct nonzero frequencies
when M>0, so it cannot equal a rational function on a regular open interval.
Local continuation of the selected boundary branch handles intervals lying
on the continuum. This argument establishes exact functional inequality only.
It excludes neither interacting multiparticle continua nor a finite-bath fit
to a finite noisy filtered experiment. Real-axis delta spectral measures do
not converge pointwise to the continuous absorptive density; the final source
restricts approximation claims to appropriate smoothing or upper-half-plane
domains.

The Robin spectral transform is unitary and transfers the Robin energy along
with the bulk energy. Although its boundary coupling vector is not unweighted
L^2 in p, it is a valid form coupling because
`int dp u_p(0)^2/(M^2+p^2)=1/(c+M)`. The transformed positive quadratic action
is therefore an exact four-dimensional infinite-species replica. With matched
full state, source, contacts and measurement protocol it gives the same q
quantum correlations and operational probabilities. A geometric extra
dimension is not identified. Classical Gaussian versions may be compared as
smeared random distributions. Noncommuting q(t) operators do not acquire an
ordinary classical joint path probability.

## Independent controls and their limits

`independent_controls.py` imports no producer code and reads no project data.
It passed under `python -O`, using explicit failure checks that remain enabled.
The 22 recorded checks comprise 14 exact symbolic checks, two negative
controls, and six high-precision numerical diagnostics, including a convergence
sequence. The numerical diagnostics use 65 decimal working precision and are
not interval certificates.

- For the manufactured M=2, c=3, g=1, m^2=13/4 example, q=cos(sqrt(3)t) and
  Phi=exp(-y)q/4 solve both boundary/brane equations and the bulk equation.
  The independent Hilbert-space norm of the bulk profile is 1/32. Normalizing
  the full eigenvector gives Z_b=32/33 without differentiating the inverse
  response. All four energy exchanges agree, and total energy per unit brane
  volume is exactly 99/64. These plane-wave statements are Fourier controls,
  not claims of finite total energy on R^3.
- Two negative controls detect a reversed Robin sign and the nonzero energy
  drift caused by dropping interaction energy from this exact solution.
- Canonical spectral completeness requires `int dmu_q(u)=1` and
  `int u dmu_q(u)=m^2`. These are independent constraints from the norm and
  quadratic form of the pure q vector. Continuum scattering-amplitude
  integrals plus the normalized bound weight satisfy both moments for the
  bound example. They also satisfy both for M=2, c=3, g=1, m^2=9 with no bound
  state. In the bound example the continuum moments are 1/33 and 15/44;
  omitting the bound sector would fail substantially. Reported numerical zero
  discrepancy means equality at working precision, not an exact quadrature
  certificate. The second mass-squared moment diverges; no stronger moment
  assumption is used.
- An independently terminated slab with Neumann condition at y=L has the
  boundary response `1/[c+s*tanh(sL)]`. Its differential boundary conditions
  were checked exactly. At the manufactured k=0.6, omega=1.1+0.7i, the slab
  response converges to the outgoing half-space response; discrepancy falls
  from about 0.00327 at L=1 to 1.13e-55 at L=32. This tests the outgoing branch
  through a different boundary-value problem, not through the producer's
  spectral-integral routine.
- A direct complex-amplitude Noether calculation independently checks the
  outgoing flux and its equality to minus the brane self-force work.

The frozen producer control script was separately inspected and replayed from
the read-only source snapshot under `python -O`. Its output must agree byte
for byte with the frozen producer receipt; final replay verification is recorded
in `INDEPENDENT_REVIEW_RECEIPT.json`. Producer algebra controls complement the
functional-analysis review but cannot replace it.

The delegated spectral/noise/rival audit in `noise-rival-audit/REVIEW.md` was
performed on a provisional theorem hash. This main review reread the final
frozen theorem and contract, and verified that the quantum-probability,
spacetime-smearing and classical-distribution qualifications are resolved.
The delegated provisional report is supporting analysis, not the final hash
acceptance receipt.

## Corrections resolved and remaining work

The final source resolves the requested operational quantum-probability
qualification, finite-bath approximation qualification, time-smearing
precision, and separation of the bound q covariance from incoming force noise.
No unresolved blocker remains within the accepted conditional theorem.

A useful next analytic calculation is to fix smooth probe/readout filters and
bound finite-species quadrature and UV-tail errors in the same filtered
response/covariance norm. Without these choices, exact branch-cut separation
does not establish experimentally useful discrimination. A physical test still
needs a calibrated probe, physical parameter units and support, a preparation,
measurement errors, named rivals and a likelihood/holdout contract. A Big Bang
origin claim additionally needs gravity, source dynamics, energy supply,
backreaction and thermalization. These are missing physical ingredients, not
conclusions established by this review.

No retained numerical array was decoded or read, and no physical source,
physical target or likelihood evaluation occurred. No repository or remote
mutation occurred. The current HDBLAST de Sitter scalar track stays separate:
metric calibration **FAIL**; full continuous certificate **UNRESOLVED**;
higher-dimensional cause **NOT_ESTABLISHED**; external novelty **NOT_ASSESSED**.
