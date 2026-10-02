# Proposed diagnosis of the failed metric Ward ledger

The completed third round retains scientific **FAIL**. Root reports all four
mode runs and twelve primary rows completed within budget, with 59 internal
failures in Ward endpoint/refinement tests and no physical-quantity refinement
failures. That pattern supports a targeted test of the ledger integrator. It
does not establish cross-route accuracy, operator consistency or a successful
metric calibration. This document proposes a new, explicitly registered
saved-data diagnostic. It has not been executed. Its companion verifier
checks 22 exact symbolic identities and rejects six mutations normally and
with optimization; no physical values are loaded or calculated.

## The discrete residual and the sampling problem

Write R=a0^4 delta_rho/epsilon and P=a0^4 delta_p/epsilon. At fixed physical
mass and fixed comoving cutoff, the continuous identity is

```
R' = F,
F = L*(R-3P)-3*h'*a0^4*(rho0,K+p0,K),  L=-1/eta.
```

The finite-band baseline term is essential during the source. For a uniform
history with spacing Delta, define the signed two-cell defect

```
e_m = R_(2m+2)-R_(2m)
      -Delta/3*(F_(2m)+4*F_(2m+1)+F_(2m+2)).
```

Its sum is exactly the endpoint difference minus the composite Simpson
ledger. Preserve these signs: sums of absolute cell defects lose the
cancellations that determine the endpoint failure. No finite difference or
continuous conservation assumption is required for this telescoping identity.

For a single harmonic R(t)=exp(i*Omega*t), theta=Omega*Delta, the exact
two-cell defect about t=0 is

```
2i*[sin(theta)-theta*(2+cos(theta))/3].
```

The Simpson/exact integral transfer, where sin(theta) is nonzero, is
`A(theta)=theta*(2+cos(theta))/(3*sin(theta))`. Its expansion is
`1+theta^4/180+O(theta^6)`. At theta=pi the exact endpoint increment is zero,
whereas the defect is `-2i*pi/3`: the ledger can accumulate a spurious
contribution from an alternating sampled derivative. The closed defect
formula is valid at this point even though the transfer ratio is singular.

Canonical stress oscillations contain Omega=2k. At K=256, theta ranges up
to 4 for Delta=1/128 and 2 for Delta=1/256. The coarse history therefore samples
k>64*pi above the ordinary Nyquist frequency; the fine history's canonical
phase lies below Nyquist but is not in a uniform small-theta regime. The
source and background envelopes are not strictly band-limited, so the latter
fact is not a sampling-accuracy theorem. The local GL8 forcing integrator is
different: mapping one time cell to[-1,1] gives phase parameter k*Delta, at
most 2 or 1. Its polynomial exactness through degree 15 does not supply a
physical error bound without bounds on derivatives of the full forcing
integrand. A failure of Simpson does not imply a failure of that mode solver.

The earlier maxima, 0.02107360235728546 and 1.9552444892255006e-05, occur at
different epochs. Their ratio is not a convergence order. Even a ratio at
the same epoch/cutoff would confound changed momentum and mode-time grids
and oscillatory cancellations; a factor 16 is only an asymptotic expectation
under suitable smoothness and fixed-integrand assumptions. A late failure
does not locate its origin: the global ledger retains errors accumulated
inside the source. Small Wronskian/energy diagnostics do not bound stress or
ledger accuracy.

## An exact source-free ledger from retained modes

Both eta_a=-2.5 and eta_b=-1.5 have retained complex modes and lie outside
the compact source. On this interval h and all its derivatives vanish.
All directional metric subtraction contacts and the Ward baseline term
therefore vanish; the response stresses are the bare mode variations. Keep
the actual unprojected u_a,w_a and the actual momentum quadrature. For each k,

```
Omega=2k, E(t)=exp(i*Omega*(t-eta_a)),
d=w_a/(i*Omega), c=u_a-d,
u(t)=c+d*E(t), w(t)=w_a*E(t).
```

Do not discard Re(c). It contains any retained linear normalization defect;
setting it to zero would be a projection. The direct per-mode response
densities, before multiplying by k^2/(2*pi^2) and the momentum weights, are

```
R_k(t) = [(2k^2+3L^2)*Re(c)
           +Re((3L^2*d-L*w_a)*E)]/(2k*epsilon),
P_k(t) = [(2k^2/3-L^2)*Re(c)
           +Re((-(4k^2/3+L^2)*d-L*w_a)*E)]/(2k*epsilon),
F_k(t) = 3L^3*Re(c)/(k*epsilon)
           +Re(((2k^2*L+3L^3)*d+L^2*w_a)*E)/(k*epsilon).
```

These formulas follow directly from the unrestricted complex bilinears and
the free ODE. They do not assume a numerically vanishing Wronskian. The pure
verifier also checks the full bare metric Ward identity with arbitrary u,w
and all metric contacts during the source; the identity holds for the ODE
vector field at any phase-space point. Consequently, an instantaneous Ward
algebra check alone cannot establish accuracy of a discrete trajectory.

An independently integrated free ledger uses the moments
`J_n=integral L(t)^n E(t) dt`, with

```
J2=[L*E]_a^b-i*Omega*J1,
J3=([L^2*E]_a^b-i*Omega*J2)/2.
```

In `d*(2k^2*J1+3*J3)+w_a*J2`, w_a=i*Omega*d causes J1 to cancel exactly.
No special-function integral is needed. The resulting exact ledger is

```
I_k = 3*Re(c)*(L_b^2-L_a^2)/(2k*epsilon)
      +Re(d*(3*[L^2*E]_a^b-i*Omega*[L*E]_a^b))/(2k*epsilon).
```

This is a phase-aware antiderivative of the independently expanded F_k.
Its inputs are the initial raw complex modes and background; it does not
use the measured final R or define a stress through conservation. As it
must, the antiderivative equals the correctly propagated free stress
increment. It is a free-flow/operator/ledger consistency diagnostic,
not a third independent physical theory implementation.

## Fixed prospective diagnostic and interpretation

Before running this diagnostic, freeze its source, the immutable third-round
archive/result hashes, the corrected-reader evidence, these formulas and
the following criteria in a separate registration. It is an openly post hoc
diagnosis of an already failed experiment, with outcomes still uncomputed.
Preserve all prior routes, ledgers, tolerances and FAIL receipts.

Use both sources, both resolutions and K=64,128,256: twelve fixed tests on
the single preselected interval[-2.5,-1.5]. Evaluate the analytic free
profiles and ledger from retained initial modes using the same saved momentum
nodes/weights. Preserve long-double values through exact binary-rational
conversion with roundtrip guards; use fixed 80/100-digit calculations and a
maximum 1e-12 precision-gap/serialization diagnostic. This is an empirical
arithmetic check, not an interval certificate. Keep the 900-second and 262144-KiB
resource limits and process momentum/time values in bounded blocks.

Reconstruct the two endpoint modes through the exact free flow and report
`delta_u=u_b-u_a-(E_b-1)*w_a/(i*Omega)` and `delta_w=w_b-E_b*w_a`.
Their signed contribution to final density is exactly

```
E_flow,k = [(2k^2+3L_b^2)*Re(delta_u)/k
             -Im(delta_w)-L_b*Re(delta_w)/k]/(2*epsilon).
```

A triangle bound replaces the real/imaginary parts by magnitudes with
coefficients `(2k^2+3L_b^2)/k` and `1+L_b/k`; report this possibly loose bound
separately. It bounds the mathematical projection of the measured free-flow
defects, not the inherited initial-state error or total numerical error.

From the unchanged stored history calculate the signed reset Simpson ledger
S_ab=S_b-S_a and the endpoint increment DeltaR=R_b-R_a. Independently calculate
I_ab by the free phase formula and form

```
D_S=DeltaR-S_ab,
D_cont=DeltaR-I_ab,
E_Q=I_ab-S_ab,
D_S=D_cont+E_Q.
```

The final identity is exact; numerical closure must meet 1e-12. Require every
sampled post-support R and P to agree with the analytic free profile within
2e-7, the original cross-stress scale; require both |D_cont| and the weighted
signed |E_flow|<=2e-7 in all twelve tests. These new diagnostic criteria
cannot grant the old calibration PASS. Compare fine-history Simpson at its
native step and at double step using `F_fine[::2]`: this contrast changes only
ledger weights on one fixed trajectory, eliminating the earlier resolution
confounding. Do not use a factor 15 Richardson estimate as a certified error.

Call a *post-support ledger error demonstrated at the original gate scale*
only if all the consistency criteria pass and at least one reset |D_S|>2e-6
has |D_cont|<=0.1*|D_S| and |E_Q|>=0.9*|D_S|. Otherwise report which check
fails, or that no post-support gate-scale attribution was demonstrated.
No fallback tolerance, alternative interval, state, cutoff or algorithm is
selected after seeing these diagnostic results.

A successful attribution supports replacing Simpson for the source-free
segment in a later separately registered numerical implementation. It does
not authorize changing this failed run or establish the active-source ledger.
If the reset defect is small, signed cell defects can localize the earlier
accumulation, but scalar sampled histories alone cannot certify the continuous
integral. A general active-source replacement needs separately specified
phase-aware quadrature/dense mode information and operator checks. The current
data justify the bounded exact free diagnostic; they do not yet select that
general replacement.

## The independent reader error is separate

For C=a0''/a0=2L^2, `C^(n)=2*(n+1)!*L^(n+2)` because L'=L^2. The old
reader's `(n+2)!*L^(n+2)` agrees only at n=0 and multiplies the next derivatives
by (n+2)/2. The producer's retained C0,C1,C2 are 2L^2,4L^3,12L^4, which are
correct. Preserve the reader's original failed evidence and qualify corrected
saved-data checks separately. Reader metadata and literal-subtraction ULP
failures are not evidence of a new physical discovery. In particular this
post-support proposal evaluates no directional local contacts, so it does
not silently repair or reuse that wrong metadata formula.

Actual shifted-root propagators, lapse/gauge controls, bulk/boundary/state
completion, coupled dynamics and stability remain independent prerequisites.
The proposed identities and diagnostic are established numerical analysis
and field-theory consistency tools, not a fundamental mathematical novelty.
