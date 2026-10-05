# Independent review of the primary finite-band/error transport result

The independent action/bilinear derivation and its first exact symbolic
receipt were completed before reading the new primary implementation or
derivation. The full approximate-state residual and closed band kernels
were also sent to the root agent before reading the primary coefficients.

Reviewed:

- `../verify_finite_band.py`, standard-library exact sparse polynomials.
- `../EXACT_CHECKS_INITIAL.json`, 33 passed checks.
- `../theory/FINITE_BAND_ERROR_TRANSPORT.md`.

The independently derived full raw residual, signed f-g term, canonical
drift, closed kernels, Pi convention and K=64/128/256 component bounds agree.
The primary's exact joint-complex-residual coefficient improves the simple
componentwise bound and is valid by Cauchy--Schwarz followed by
`integral k sqrt(k^2+L^2) dk`. The fixed-background L is positive; a general
background would replace L^3 by |L|^3 in the lower-endpoint term.

The envelope-reflection first-moment bound alpha/2 is valid: the inherited
Cauchy envelope is the same even function about every equal-width panel
center, so its global reflection about -4 preserves its value. Its weighted
full-domain first moment is therefore one half of its total mass. For every
prefix, positivity and t-s<=b-s give the same upper bound. This is a calculus
argument, not a consequence of the finite 33 polynomial checks alone. The
primary document correctly makes that distinction.

The bare kernel formulas have apparent poles only in their closed quotient
expressions. Their finite-interval sine/cosine integral definitions prove
entireness, independently of checking ten Taylor coefficients. The finite
Taylor checks verify the implemented low-order coefficients and limits;
they are not claimed to formalize the full analytic theorem.

No blocking algebra or scope defect was found. The conditional contacts,
exact source-polynomial coefficients and finite-band integrability premises
are explicit. No source callback or archive array is evaluated. The new
pressure bound is kept separate from the nine-probe source-operator gate;
the twelve-case pressure/contact certificate remains UNRESOLVED.

One strengthening was sent to the root: the elementary constant real
`U=A,W=0` counterexample to Ward-only state certification is unprojected and
does not preserve the first-order Wronskian. A stronger example preserves it:

\[
U=A(k)e^{2ik(t-a)},\qquad W=2ikU,
\]

with A real, nonnegative and smooth, supported strictly between momentum
probes. Then `c_R=Re U-Im W/(2k)=0`, the Ward defect vanishes, all probe
values vanish, and `R_m(a)=3L(a)^2 A(k)/(2k)` has a nonzero band integral.
This example is a first-order Bogoliubov perturbation; an exact normalized
state can include the usual second-order correction to its positive-frequency
coefficient. Thus even a checked first-order Wronskian and all registered
probe values do not determine the continuum stress or incoming state.

Assessment: **PASS conditional analytic theorem and truncation component**.
Assessment of a full physical pressure/contact certificate: **UNRESOLVED**.
External mathematical novelty: **NOT_ASSESSED**.
