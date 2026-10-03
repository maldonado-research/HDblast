# Independent active-source action check

This is pure symbolic preparation. No saved physical array, physical source
value, response, or source integral was loaded or evaluated. The final verifier
passes 26 exact identities and rejects 12 symbolic mutations under normal
Python and `python -O`. Use `NORMAL_FINAL.json`, `OPTIMIZED_FINAL.json`, and
`FORMULAS_FINAL.json`; the earlier receipts describe the earlier development
version and remain historical evidence.

## Action and unrestricted complex amplitudes

For the minimally coupled scalar, the canonical mode action gives

```
v'' + (k^2 + M - C)*v = 0,
M = a^2*m^2,  C = a''/a,  ell = a'/a,
R_bare = (|v'-ell*v|^2 + (k^2+M)*|v|^2)/2,
P_bare = (|v'-ell*v|^2 - (k^2/3+M)*|v|^2)/2.
```

The first check differentiates these two separate quadratic operators at an
arbitrary complex phase-space point. With `M'=2*ell*M` and
`ell'=C-ell^2`, it proves `R_bare'=ell*(R_bare-3*P_bare)`. No state
normalization or Wronskian value is imposed. Pressure comes from the action
operator, not from this Ward relation.

Fix `m^2=2`, `xi=0`, `H=1`, `a0=L=-1/eta`, and fixed physical phi. For
`a=a0*(1+epsilon*h)` to first order,

```
delta(M-C)/epsilon = g = 4*L^2*h - 2*L*h' - h'',
v0=exp(-ik*eta)/sqrt(2k), delta_v=v0*u, w=u',
u'=w,  w'=2ik*w-epsilon*g.
```

Here `u=X+iY` and `w=Z+iT` are unrestricted. The first-order density and
pressure, scaled by `a0^4/epsilon`, follow directly by Taylor expanding the
quadratic operators and `a^-4`:

```
R0 = (2*k^2+3*L^2)/(4*k),
P0 = (2*k^2/3-L^2)/(4*k),
Rmode = ((2*k^2+3*L^2)*X/k - T - L*Z/k)/(2*epsilon),
Pmode = ((2*k^2/3-L^2)*X/k - T - L*Z/k)/(2*epsilon).

Rcontact = L*h'/(2*k) - 2*k*h - 2*L^2*h/k,
Pcontact = L*h'/(2*k) - 2*k*h/3,
R=Rmode+Rcontact, P=Pmode+Pcontact.
```

The contacts combine three separately checked terms: the operator variation
`delta ell=h'`, the fixed-mass variation `delta M=4*L^2*h`, and the
`a^-4` prefactor variation `-4*h*(R0,P0)`. There is no scalar mass-law
contact when physical phi is fixed.

## During the source

The mode contribution alone obeys

```
Rmode' = L*(Rmode-3*Pmode) + L*g/(2*k).
```

The direct contact operators obey

```
Rcontact' = L*(Rcontact-3*Pcontact)
            -3*h'*(R0+P0) - L*g/(2*k).
```

Thus the full bare identity is

```
R'=F, F=L*(R-3*P)-3*h'*(R0+P0).
Fcontact=-2*k*h'-(5/2)*L^2*h'/k-2*L^3*h/k.
```

Multiply each term by the same fixed `k^2/(2*pi^2)` measure and momentum
weights. After subtracting separately derived stress subtraction operators,
the baseline becomes the same finite-band bare-minus-subtraction baseline,
not a substituted continuum value. The verifier proves this subtraction
transfer as an exact conditional algebra identity. It does not independently
recompute all WKB coefficients or their actual momentum integrals.

The ODE also gives phase-aware coordinates without using a stress endpoint:

```
d=w/(2ik),       d'=2ik*d+i*epsilon*g/(2k),
c=u-w/(2ik),     c'=-i*epsilon*g/(2k),
Re(c)'=0,
[w*exp(-2ik*(t-a))]'=-epsilon*g*exp(-2ik*(t-a)).
```

These identities specify the forced-mode information that an independent
active-source ledger needs. They do not make a sampled scalar history
sufficient to integrate its continuous oscillatory stress accurately.

## What the identity can diagnose

Let continuous trajectory defects be

```
e_u=u'-w,
e_w=w'-2ik*w+epsilon*g.
```

The exact off-shell bare Ward defect is

```
R'-F = ((2*k^2+3*L^2)*Re(e_u)/k
        -Im(e_w)-L*Re(e_w)/k)/(2*epsilon).
```

This is a signed projection. A small projection does not bound every
trajectory error. The on-shell symbolic Ward identity holds at any
phase-space point and therefore cannot establish physical trajectory
accuracy. A numerical active-source test still needs a prospectively fixed
trajectory/quadrature strategy, independent source integration, endpoint
operator checks, unchanged prior failure evidence, and explicit precision
and resource criteria. This proof supplies established action/ODE consistency
mathematics; it is not a fundamental novelty claim or evidence for the
higher-dimensional blast hypothesis.
