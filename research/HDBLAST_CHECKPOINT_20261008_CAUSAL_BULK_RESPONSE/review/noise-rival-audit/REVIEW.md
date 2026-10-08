# Spectral, noise and continuum-rival audit

8 October 2026. Independent, read-only analytic review of
`physical-model/HALF_SPACE_RESPONSE_THEOREM.md`. No physical target, likelihood,
source, or protected trajectory calculation was run. No source edit was made.

The provisional source was changing during review. The latest inspected content
had SHA-256 `275277b3e4296ed02524e82e0cd088923008e1649a1d84374fa8fc0e193ca682`.
This report covers its spectral/FDT and rival assertions; the parent review
handles the classical energy proof and numerical controls independently.

## Verdict

Bounded acceptance: no remaining mathematical blocker was found in the reviewed
claims under the declared quadratic, positive, canonical-normalization,
spacetime-smearing, fixed-state and fixed-measurement premises. This is not a
physical validation, novelty judgment, or authorization to evaluate protected
targets. It does not show that a finite experiment can distinguish the model
from finite-dimensional approximations.

An initial substantive wording objection has been resolved in the latest text:
the continuum rival now gives identical smeared quantum states/correlations and
measurement probabilities for a specified common protocol. It explicitly avoids
assigning a classical joint path probability to noncommuting q operators at
different times. That qualification is necessary.

## Checks and reasoning

1. **Robin spectral normalization.** The stated modes are delta-normalized with
   respect to dp. Their boundary weight is `(2/pi)*p^2/(p^2+c^2)`. Dividing this
   weight by `du/dp=2p` gives exactly
   `rho_b(u)=sqrt(u-M^2)/(pi*(c^2+u-M^2))`. Thus there is no missing factor of two
   or pi in (4). The high-u density is of order `u^(-1/2)`; its resolvent
   integral converges. No UV subtraction is needed for this particular response.

2. **FDT convention.** One oscillator with `E^2=k^2+u` has
   `Im [E^2-(omega+i0)^2]^-1 = pi/(2E) [delta(omega-E)-delta(omega+E)]`.
   Its half-anticommutator spectrum is
   `pi/(2E) coth(beta E/2) [delta(omega-E)+delta(omega+E)]`.
   This verifies the factor and sign in (14), with the declared Fourier and
   canonical source conventions. Applying the same calculation to the bath
   boundary weights verifies (15). Positivity concerns the positive spectral
   measure or the even symmetrized noise; `Im G_R` itself is odd in frequency.

3. **State and bound-sector completeness.** A full coupled KMS state supplies
   the thermal covariance of the bound oscillator when it exists. Stationary
   incoming bath noise determines only the scattering continuum. The text
   correctly restricts `S_q=|chi_R|^2 N_eta` to that continuum and retains the
   bound delta contribution separately. The threshold-tuned case may give an
   integrable threshold singularity, which does not supply an extra normalizable
   bound delta. Gaussian Wick reconstruction must use the appropriate ordered
   quantum two-point functions; the symmetrized function together with the
   response/commutator determines them in this model.

4. **UV and distributions.** The stated `N_eta ~ g^2/|omega|` vacuum-frequency
   tail at fixed k gives a logarithmically divergent equal-time force variance.
   This is consistent with a positive Gaussian covariance on smooth spacetime
   test functions and with a convergent retarded kernel. It does not define
   ordinary pointwise noise or a finite stress tensor. A quantum ground/KMS
   state here is an algebraic Gaussian state in infinite volume; no trace-class
   global Gibbs density matrix is required or established.

5. **Finite-free-sector theorem.** At fixed k a finite local quadratic system
   with finite derivative order gives a finite polynomial matrix in frequency;
   any well-defined linear response obtained by eliminating fields is rational.
   A rational readout preserves this property. If its response equalled the
   declared response on a regular open interval, nonzero g would make the
   square root rational. The simple zeros of `k^2+M^2-omega^2` prohibit that.
   Local analytic continuation justifies the same argument on a regular segment
   of the continuum. The theorem appropriately excludes interacting continua,
   nonlinear/composite readouts, nonlocal filters and arbitrary numbers of
   effective hidden degrees of freedom. It proves exact functional inequality,
   not a nonzero separation in a finite noisy experiment.

6. **Exact four-dimensional continuum rival.** The Robin transform is unitary
   on the bulk one-particle spatial Hilbert space and diagonalizes both the bulk
   gradient energy and its Robin boundary energy. The q coupling transforms to
   the displayed integral of `u_p(0) chi_p`. Although the coupling vector is not
   square-integrable without a weight, it is well-defined as a form because
   `int dp u_p(0)^2/(M^2+p^2)=1/(c+M)`. The declared strict stability inequality
   therefore controls it. Matching the full state and the smeared q measurement
   algebra indeed transfers all q correlations and operational probabilities.
   This is an exact infinite-species replica; it supplies no finite-species,
   gravitational, finite vacuum-energy, or independent UV-completion theorem.

## Optional precision improvements

- Say **smooth spacetime smearing** explicitly for force noise. Spatial smearing
  alone cannot remove the fixed-k equal-time logarithmic divergence.
- Replace the remaining phrase "Classical Gaussian versions have identical q
  process laws" with "Classical Gaussian versions have identical laws for the
  smeared q random distributions" if pointwise classical processes are not
  separately constructed.
- An optional normalization cross-check is the q sum rule
  `int rho_q(u) du + Z_b = 1`, with `Z_b=0` when there is no bound state. It
  follows from canonical normalization and the large-imaginary-frequency
  behavior `chi_R(iR,k) ~ 1/R^2`; the sum rule is not an additional assumption.

These refinements are not objections to the corrected conditional result.
