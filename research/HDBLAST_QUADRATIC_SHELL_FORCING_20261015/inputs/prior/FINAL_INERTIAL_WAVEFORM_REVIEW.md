# Independent review of the wider inertial shell waveform

**Accepted as conditional linear first-variation mathematics.** The accepted
source bytes are pinned in `FINAL_INERTIAL_WAVEFORM_REVIEW.json`. No substantive
mathematical or computational blocker remains in the reviewed checkpoint.

The comparison retains the positive shell inertia and the prescribed nonzero
initial velocity. Its integration-by-parts kernel has the correct signs,
including the `-1` in the leading postpulse coefficient. Odd reflection
cancels the constant-potential displacement contribution after the pulse;
the nonzero velocity moment, coefficient variation, boundary lift, bulk
scattering and preceding pulse memory remain in the residual estimate.
The resulting derivative trace bound is obtained directly from incoming
characteristics, without inferring a pointwise normal trace from an energy norm.

The independent exact rational derivation proves global comparison error
below `1.022e-7` and postpulse error below `1.645e-10`. The separately reviewed
continuous comparator interval, enlarged by the published `1.7e-10`, proves
`3.2e-8 < H(x) < 5.3e-8` at every real `2 <= x <= 3`. Its physical quadrupole
coefficient has negative sign because the retained gain is negative. This
does not assign one sign to the angle-dependent quadrupole harmonic.

Every one of the 4353 exported cells was independently reconstructed in both
Python modes. The audit imports the separately reviewed free comparator;
it does not import the parent generator or its error evaluator. It uses the
pinned independent error derivation, exact rational gain/clock products,
whole-cell square sums and the `1/ell` proper-time norm Jacobian. It reproduces
the fresh wider-state interval
`6.06e11 < C_(3ell) < 7.44e11`, including the positive tail contribution.
No narrow-state response coefficient is transferred.

Independent signed endpoint calculations establish at least one zero in
each of `(1.49998,1.49999)` and `(1.96,1.97)`. No zero uniqueness, total zero
count or uniqueness of extrema follows. The normalized trace and the
physical coefficient have the same zeros and opposite signs.

The weighted matrix Taylor recurrence, its degree-eight remainder and the
global short-time positivity argument justify the comparison kernel ranges.
The earlier local-interval positivity concern was corrected by a global
derivative bound before acceptance. The manuscript's displayed norm interval
now rounds outward. Reproduction dependencies and the externally supplied
manifest digest are described accurately in the README.

The public residual-energy methods note has the correct boundary work sign,
positive boundary kinetic term, both coefficient derivative terms, coercivity
and trace constants. Its finite growth bound is correct. It controls a
boundary value under its premises and explicitly leaves the derivative
observable and nonlinear defect map to additional arguments. No theorem or
private research material is transferred by that note.

The figure was visually inspected: the physical lobe signs, negative tail
inset, model-time axis, continuous bands and midpoint-guide qualification
agree with the certificate. Floating rendering is not used for certification.

Executed checks and fault controls are reported by component in the JSON
receipt, without treating repeated or overlapping groups as new theorems.
The review concerns the local scientific source and completed calculations;
it does not certify a future publication, deployment or remote readback.
The retained background, reconstruction, inverse, reflection, gauge, support
and evolution hypotheses remain explicit inputs. There is no effective
nonlinear epsilon/remainder bound, finite-amplitude nonlinear waveform,
observational calibration or cosmological-origin conclusion. Internal
independent AI review and finite arithmetic are distinct from external peer
review and formal verification of the analytic PDE theorem.

The earlier conditional complete second Taylor coefficient through shell
contact is preserved. This checkpoint's quantitative next target is a
validated second-order observable for the wider state and an effective
nonlinear defect-to-error estimate; it does not revoke that earlier
conditional construction.
