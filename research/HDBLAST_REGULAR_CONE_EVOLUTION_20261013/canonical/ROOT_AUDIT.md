# Independent background and canonical-variable audit

13 October 2026. Internal calculation, not external peer review.

The root calculation independently verifies the regular-master conjugation
and exact non-gauge Cauchy jets. It uses arbitrary symbolic coefficient
functions, not a fitted background or numerical evolution.

For `f=Phi_s`, `L=Omega_s/Omega+2f_s/f`, the positive density
`A=1/(Omega f²)` has `A_s/A=-L`. Writing `psi=sqrt(A)F` gives

    (Box_eta - Vcan) psi = 0,
    Vcan = 3L - 2s L_s + sL² + (8/3)s f².

The cancellation of nonradial first derivatives is checked separately from
the radial identity, so the calculation applies to arbitrary Cartesian
fields. This auxiliary scalar equation gives finite propagation and local
energy estimates; its action is not identified with the gravitational
symplectic action.

Independently expanding the background equations gives

    rho(y) = y + rho3 y³ + ... ,  rho3 = -U_h/36,
    phi(y) = phi_h + kappa y²/2 + phi4 y⁴ + ...,
    kappa = U'_h/5,
    phi4 = kappa (U''_h/56 + U_h/126).

The exact conformal-radius correction matters: `y=r+(rho3/2)r³+...`.
Consequently

    Omega_s(0) = -U_h/24,
    L(0) = U''_h/7 - 5U_h/56,
    Vcan(0) = 3U''_h/7 - 15U_h/56.

These values also agree with the separately derived curved minimally
coupled master `P/(rho³p)`: its potential at the cone is
`-3U''_h/7-5U_h/14`, and conformal conjugation produces the same `Vcan(0)`.
The finite coefficients do not remove the hypothesis `f(0)!=0` from the
master reconstruction. In particular, the stationary detuning limit is
not covered uniformly.

For the exact annular `l=2` data in the manuscript, the cone Wronskian is
`14C²ab H_2²` and the doubled angular-integrated current density at
`v=2t0` is `28t0³C²ab`. The two independent gauge-invariant relational
scalars have nonzero coefficients whenever `|t0²L(0)|<1`.
The checker also retains `lambda=l(l+2)` symbolically before specializing.

The wrong-formula controls reject a canonical mass-sign reversal, omission
of the conformal-radius correction, the wrong null Euler-derivative
factor, omission of angular dependence and a copy-count error. An explicit
injected failure is supported and remains effective under optimized Python.

These symbolic identities are components of the analytic proof. They do
not certify global shell evolution, a gravitational energy completion,
external novelty or a cosmological mechanism.
