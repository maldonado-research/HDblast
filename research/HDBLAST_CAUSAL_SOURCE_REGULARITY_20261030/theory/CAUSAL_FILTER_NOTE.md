# A past-time filter for a distinct numerical candidate

This note audits an auxiliary reconstruction, not the actual differentiated
physical source. Let `FH=17387513156573` be the authenticated uniform
full23 forcing bound in the positive H norm. This is a pointwise amplitude
bound for the complete physical forcing, including its metric/scalar
normalizations. The different response quantity `X(KF+LG)` is not used as FH.
All source blocks and their Gram normalization are rechecked in the executable.

## Definition and domain

Extend the physical half-line amplitude vector F by zero for negative scaled
time. For a width `epsilon>0`, define

```
k_epsilon(s)=30 s^2(epsilon-s)^2/epsilon^5, 0<=s<=epsilon,
F_epsilon(x,r)=integral_0^epsilon k_epsilon(s) F(x-s,r) ds.
```

The fixed default is `epsilon=1/1024`. This is component smoothing in the
registered coordinate/frame representation, not a covariant physical average.
Couple the complete physical source before applying the output parity
reflection. Time convolution then commutes with that fixed output reflection;
no parity/matrix commutation is assumed.

The domain is `0<=x<=3/2`, `max(0,x-5/4)<=r<=7/4-x`. If `x'=x-s>=0`,

```
7/4-x'-r=(7/4-x-r)+s >=0,
r-x'+5/4=(r-x+5/4)+s >=0.
```

Thus every nonnegative past sample at fixed r remains in the domain; negative
samples use the specified zero extension. This domain lies in the inherited
actual-source triangle `0<=x<=3,r>=0,x+r<=8`. If `x+r<=1`, every past sample
also has `x-s+r<=1`, so the exact zero wedge is preserved. A future-time
filter would fail both these arguments; the checker includes counterexamples.

## Kernel and time derivatives

For the unit kernel `k(u)=30u^2(1-u)^2`, direct rational integration gives

```
integral k=1, integral u k=1/2, integral u^2 k=2/7,
k(0)=k(1)=k'(0)=k'(1)=0,
||k'||_L1=2 k(1/2)=15/4.
```

The kernel is nonnegative, and `k'(u)=60u(1-u)(1-2u)` has the stated sign
on the two half intervals. Since the compact zero extension of k is C1
and has integrable second derivative, convolution against bounded F is C2
in time. Translation continuity in L1 of the kernel derivatives justifies
the continuous derivatives without differentiating the unknown F. Therefore

```
||F_epsilon||_H <= FH,
||F_epsilon-F||_H <= 2 FH,
||(F_epsilon)_x||_H <= 15 FH/(4 epsilon),
||(F_epsilon)_xx||_H < 24 FH/epsilon^2.
```

For the last optional bound, `k''=60(1-6u+6u^2)` has its roots at
`(3+-sqrt(3))/6`. The extrema of k' have modulus `10 sqrt(3)/3`, giving
`||k''||_L1=40 sqrt(3)/3<24`. The first derivative bound at epsilon=1/1024
is `3840 FH=66768050521240320`.

Time-only convolution supplies no new radial derivative estimate by itself.
Candidate second jets require the separate principal-wave/initial/boundary
argument. The actual source's qualitative continuity remains an inherited
premise where that separate argument needs it.

## Explicit candidate change and C1 charge

Keep the actual target U, actual F and G, constant actual Robin matrix M,
homogeneous coefficient operator La, first field, and the quadratic-continued
smoothed initial data fixed. Define the distinct auxiliary candidate

```
Z_epsilon=P[Bs_quadratic+K(La W0_quadratic+F_epsilon)+LG],
P=(I-L M PN)^(-1).
```

The actual target still uses F. The exact difference from the published
quadratic candidate is `Delta=Z_epsilon-W1_quadratic=P K(F_epsilon-F)`.
Using the registered `X(K f)<= (3/2) sup||f||_H` and `||P||<=1/(1-b)` gives

```
X(Delta) <= 3 FH/(1-b) < 5.217e13.
```

The checker adds this rational quantity to the exact previous target-error
majorant, rather than adding it to the rounded readable bound. The result
is still strictly below `3.541e25`. This retains the same readable C1 upper
certificate for a different numerical candidate. It is not a claim that
the unknown actual errors are equal, or that a large upper certificate is
a small relative error.

Unlike the changed outer initial continuation in the preceding checkpoint,
this filter can affect earlier causal regions. The published selected-point
tables are not inherited unchanged; any use requires a perturbation charge
or new computation. The generic C1 comparison above supplies such a charge
for fields and first jets, but does not itself produce an updated table.

## Initial face and physical shell corner

For bounded F, the filter and its first two time derivatives tend to zero
on the initial face: the kernel, its first derivative and its second
derivative yield integration estimates of orders `x^3,x^2,x`, respectively.
In particular, `F_epsilon(0,r)=0` and `(F_epsilon)_x(0,r)=0`. The two
candidates have the same displacement and velocity. Their initial mixed
and radial second-jet differences also vanish, while their principal
equations give

```
Delta_xx(0,r)=-F(0,r).
```

Consequently the error forcing for `e=U-Z_epsilon` is

```
r_epsilon=La(Z_epsilon-W0_quadratic)+F-F_epsilon,
r_epsilon(0,r)=F(0,r).
```

The previous identity `r1(0,r)=0` cannot be reused. The actual F(0,r) is not
assumed zero in the interior active packet. Its magnitude is bounded by FH.
At the physical corner `(0,0)`, however, the inherited open zero collar of
initial data and sources is preserved by this past filter. Hence the
new initial source discrepancy does not break the shell corner matching.
The actual affine shell law and its boundary derivative remain unchanged
because M,G and the initial data are unchanged.

## Limits and verification

The rough discrepancy bound `2 FH` has no decay with epsilon. Even if actual
F is continuous for nonnegative time, negative-time zero extension can create
an initial jump, so uniform approximation to F on the closed initial face
does not follow as epsilon tends to zero. Away from that face an actual
modulus of continuity would be needed for a quantitative convergence rate.
The derivative cost grows as `1/epsilon`; smaller width is not an automatic
improvement of a later residual estimate.

This route does not bound actual F_x, C2 error to U, or the full natural
Einstein/scalar residual by itself. Replacing F by its auxiliary approximation
inside a numerical candidate also prevents inheriting exact bulk Bianchi or
tangential-constraint preservation. A separate residual calculation must
retain `F-F_epsilon`. No exact nonlinear spacetime, finite-amplitude
realization, cosmological prediction or higher-dimensional Big Bang evidence
is established.

`check_causal_filter.py` runs using only the standard library, authenticates
the complete source and inherited reports, verifies exact kernel identities,
reconstructs the complete forcing norm, charges the candidate change, and
rejects eight wrong formulas. Normal and optimized Python produce identical
reports. This audit has no network or publication side effects.
