# New finite-collar check: scalar stress, both boundaries, and geometry work

The instrumented reruns support a numerical scalar-stress balance on a fixed finite collar. They also confirm the sign found from the archived junction: the shell supplies net energy through the collar's wall boundary during the selected interval. The collar nevertheless loses energy because its inner-boundary outflow is larger. This is a local, observer-specified stress balance, not depletion of a positive global five-dimensional energy reservoir.

This experiment is separate from the archived trace analysis in `BOUNDARY_WORK_AND_RESERVOIR.md`. It uses new, densely recorded six-field profiles. The archived NPZs alone did not contain those data.

## Fixed diagnostic and independent derivatives

The two original B1 cells with phi_star=0.5, G/H0=100, y=1 and b=6.49e-6 were rerun with profile output at computational spacing Delta T approximately 0.01. The geometry lane selected Z in [-0.04,0] and shell time s=H0 tau in [1,4] before reading the new snapshot values; `COLLAR_PROTOCOL.md` records the choice. The saved profile extension reaches at least Z=-0.08. Both grids supply 155 time rows on the clipped interval.

The diagnostic reads A,B,phi,A_T,B_T,phi_T. It computes phi_Z and B_Z with its own local polynomial stencils, without importing the producer or using its ghost/Neumann derivative. At Z=0 the derivative stencil uses only strictly interior samples; the wall field value is excluded from that derivative fit. Degree-two and degree-four stencils are reported separately. The fields themselves are the evolving state of a solver that imposes the junction; this is an independent diagnostic of that numerical solution, not an independent physical realization.

Spatial and temporal integration both use the declared composite trapezoid. Polynomial degree four does not make this a fourth-order integral method. Endpoint times are interpolated from the saved monotone (s,T) clock, and already weighted/integrated columns are interpolated separately at those times. The energy normalization is common a(0)=1 comoving volume, energies multiplied by kappa5-squared, per one retained interior.

## Exact identity used

For the explicit observer field xi=e^(-B) partial_T, define

\[
E=e^{3(A-A_0)-B}\left[\frac{\phi_T^2+\phi_Z^2}{2}+e^{2B}U\right],
\quad F=-e^{3(A-A_0)-B}\phi_T\phi_Z,
\]
\[
G=e^{3(A-A_0)-B}\left[
\frac{3A_T-B_T}{2}\phi_Z^2-\frac{3A_T+B_T}{2}\phi_T^2
+(3A_T+B_T)e^{2B}U+B_Z\phi_T\phi_Z\right].
\]

The scalar equation implies partial_T E+partial_Z F=G. `FINITE_COLLAR_IDENTITY.md` derives it, including the observer field and all potential/geometry terms. It needs the scalar equation on the evolving metric; it does not by itself test the Einstein equations. For the fixed collar,

\[
\Delta E_{\rm coll}=\int F(-0.04)\,dT-\int F(0)\,dT+
\int\!\int_{-0.04}^{0}G\,dZ\,dT.
\tag{1}
\]

Positive F(0) means energy leaving the interior toward the shell. At that boundary F(0)dT=-a^3 wv d tau. Twice the wall integral therefore matches the orientation and normalization of the total shell influx; it is not multiplied by another lapse.

## Results, with the cancellations exposed

The degree-four derivative results are:

| Term, per retained copy | Coarse grid | Fine grid |
|---|---:|---:|
| Initial collar energy | 1.09178485297 | 1.09133527175 |
| Final collar energy | -14.5568261772 | -14.5614762015 |
| Measured Delta E_coll | -15.6486110302 | -15.6528114733 |
| Signed inner influx, integral F(-0.04) | -57.8092605798 | -57.8288108335 |
| Outward wall flux, integral F(0) | -37.7829107823 | -37.7958661738 |
| Metric work, double integral G | +4.37684961098 | +4.37976806431 |
| RHS of (1) | -15.6495001865 | -15.6531765955 |
| Endpoint minus RHS residual | +0.000889156324 | +0.000365122169 |
| Residual / max(abs(Delta E),abs(RHS)) | 5.68169e-5 | 2.33258e-5 |
| Residual / sum of absolute integrated terms | 7.69049e-6 | 3.15693e-6 |

The smaller last-row numbers arise from cancellation between larger terms and must not replace the preceding error scale. These finite discrepancies are not continuum error bounds. The energy becomes negative because the actual bulk potential is included; replacing it by kinetic energy alone would change the observable and invalidate (1).

With degree-two derivatives, the residuals are -0.00301412064 and -0.00050868033, respectively 1.92614e-4 and 3.24978e-5 of the energy-change scale. Retaining every other time row changes the degree-four predicted energy change by approximately 2.35e-4 and 2.34e-4 relative, larger than the full-cadence residuals. Thus small balance residuals do not certify a converged integration error or the precision of each term.

The terms have a clear signed interpretation. On the fine grid, the wall supplies about 37.796 units to this collar, but about 57.829 units leave through its inner face. The declared observer's geometry-work contribution is about +4.380. Their sum gives the negative collar energy change. This does not demonstrate that the matter source was supplied by net depletion of the entire bulk.

## Independent wall-gradient comparison and failures of omitted terms

The degree-four interior-gradient measurement gives two-copy wall integrals -75.5658215646 and -75.5917323476. Evaluating the scalar junction on the same reruns gives -75.5658292477 and -75.5917341773. This comparison uses measured interior gradients on one side and the imposed junction expression on the other; it is stronger than substituting the junction on both sides. It still concerns the same numerical solution, not separate observational evidence.

Dropping metric work from the balance leaves residuals about +4.37774 and +4.38013. Reversing the wall-flux sign leaves about +75.56671 and +75.59210. These are meaningful large failures compared with the complete finite budget. The potential and metric-work terms are part of the specified stress current, not adjustable corrections selected to fit an endpoint.

## Replaying the limited result

The four `collar_data/*.csv` files contain only time, collar energy, the two fluxes, metric work, and the separate junction prediction. The public-size replay needs only Python's standard library:

```sh
python replay_collar.py --output COLLAR_BALANCE.json
```

The private extractor `collar_extract.py` implements the independently specified spatial derivative stencils and records SHA256 pins for every source snapshot. It never uses a producer Neumann derivative. The extraction receipts retain provenance for the rerun arrays. The optional SymPy script `check_collar_identity.py` verifies the finite displayed identity and nonzero residuals for four incorrect formulas. None of these checks turns the numerical experiment into a theorem of global energy conservation, a positive-reservoir bound, or a successful reheating scenario.
