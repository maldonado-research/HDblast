# Independent enclosure of the boundary reconstruction

This checker bounds the exact integral-defined boundary reconstruction `W`
by its difference from the previously certified first-feedback candidate
`V=B0+B1`. It imports neither the new trace-moment evaluator nor its formulas
or output. It shares the authenticated old candidate point and old field
certificates. Its independent part is the residual-resolvent enclosure of the
entire boundary correction; it is not an independent replay of every earlier
field calculation.

The reconstruction uses the actual constant shell coefficient `M`, distinct
from the frozen rational candidate matrix `M0`, and the actual affine source
`G`. The source is enclosed by inherited uniform magnitude bounds. A source
point value or sign is not resolved here.

## Residual correction and the actual boundary equation

Use inward scaled radial derivative and `PN` for the eighteen Robin rows.
The old candidate's residual against the actual law is

```
R_V = PN V_r-M PN V-G
    = -M0 PN B1-(M-M0)PN V-G.
```

Let `Jq(t)=integral_0^t q(s)ds`. Write the correction as the characteristic
lift `C=Lq`, whose shell trace is `-Jq`, inward derivative is `q`, and scaled
time derivative is `-q`. Then

```
R_(V+C) = R_V+q+M Jq.
```

The Volterra equation `q+M Jq=-R_V` therefore makes `W=V+C` obey
`PN W_r=M PN W+G` exactly. Its odd components vanish, so the five homogeneous
Dirichlet traces and their scaled-time derivatives are preserved. The lift
satisfies `C_xx-C_rr=0` on the physical side under the inherited regularity
and zero-corner assumptions. Hence `W_xx-W_rr=f0`, where `f0` still uses the
free wave `B0` and frozen bulk matrices. It does not follow that
`f0=Aactual W_x+Bactual W_r+Dactual W+Factual`: the full bulk problem is
not made self-consistent by this correction.

This is the same exact boundary reconstruction obtained by solving its trace
equation for `Z`, because subtracting the old candidate gives the displayed
correction equation with the same zero corner. No matrix inverse or
commutation with the bulk matrices or parity projector is used.

## Uniform shell-history bounds

The endpoint is `x=9/8`, and the active boundary history has length `d=1/8`.
The old nominal data vanish on the zero collar after smoothing. The actual
source's vanishing through shell time `1` is an **inherited conditional
property**: the previous formal actual `C1` response is zero on `x+r<1`,
so its shell value and inward derivative vanish there; its exact actual Robin
law then gives `G=0`. Continuity extends this to `x=1`. A magnitude bound on
`G` alone would not establish this support property.

The checker authenticates the earlier profiles, continuous initial
interpolation certificate and formal-response report. It recomputes local
initial caps from P1 nodal convexity plus certified interpolation errors over
every cell intersecting `[0,9/8+max(delta_u,delta_v)]`. Positivity of the
smoothing kernel transfers these caps to the smoothed profiles and their
first derivatives. It also verifies frozen matrix block norms and all actual
minus frozen Robin coefficient error bounds. This reproduces the prior local
cone caps in a separate arithmetic implementation.

For every shell time `0<=t<=9/8`, the free trace is bounded by

```
T0 = U0+d Vcap,
Fcap = (nA+nB)(U1+Vcap)+nD(U0+(9/8)Vcap),
B1cap = ((9/8)^2-1) Fcap/2 + d nM0 T0.
```

The active physical cone area is `(t^2-1)/2`, bounded by the displayed area at
`9/8`. It is not `(t-1)^2/2`. All factors are nonnegative, making these
endpoint cap formulas uniform over the entire preceding shell history.
The old endpoint's scalar residual interval alone could not justify a
uniform-history correction bound.

With `deltaM>=||M-M0||`, `nM=nM0+deltaM` and uniform `Gcap>=||G||`,

```
Rcap = nM B1cap+deltaM T0+Gcap,
sup||q|| <= Rcap+nM d sup||q||,
Qcap = Rcap/(1-nM d).
```

The checker verifies `0<nM d<1/10000`; the denominator is positive. The exact
rational consequences include

```
Qcap < 2.616505e24,
d Qcap < 3.270631e23.
```

Thus `||W-V||<=d Qcap` in value at this point and its coordinate first
derivative differences are bounded by `Qcap` in the registered Gram norm.
All odd-row corrections, including the radial derivative correction, are
exactly zero. Adding these component norm errors to the authenticated old
point table encloses all 23 reconstructed candidate jets without evaluating
any of the new higher trace moments.

Both TT rows admit a sharper elementary enclosure: their actual and frozen
Robin matrix rows and old candidate coefficients vanish exactly. Hence
`R_V,TT=-G_TT`, `q_TT=G_TT`, and
`W_TT=-J G_TT`, `W_TT,x=-G_TT`, `W_TT,r=G_TT`.
The report uses each inherited TT source coefficient cap directly for its
derivatives, and `d` times that cap for its value. These symmetric intervals
do not assert that a TT signal is nonzero.

The resulting positive scalar-mean candidate intervals imply the outward
rounded bounds

```
8.2903759e30 < W_w:0 < 8.2904081e30,
1.9746216e32 < partial_x W_w:0 < 1.9746223e32,
2.1515e27 < partial_r W_w:0 < 2.1568e27.
```

These are formal coefficients in the registered gauge and scaled coordinates.
Their signs are not physical forces, invariant observables or evidence of
a higher-dimensional Big Bang origin.

## Conditional transfer to the actual formal response

The old whole-domain `X(U-V)<2.184e28` estimate and the correction bound give
`X(U-W)<2.184e28+Qcap<2.184262e28` on the part of the retained cone whose
boundary history satisfies `x-r<=9/8`. For `x<r` the boundary correction
vanishes. This is a **local-history** transfer, not the whole-cone bound
through shell time `5/4`; that requires the separate global calculation.

At the shell, the shared exact actual boundary law gives

```
PN (U-W)_r = M PN(U-W).
```

The affine source cancels. The radial error is at most `nM d` times this
local `X` error, instead of the generic derivative error. Applying this to the
reconstructed candidate gives a positive conditional actual scalar-mean
radial interval approximately

```
2.15077139088e27 < partial_r U_w:0 < 2.15750034910e27.
```

The value and scaled-time derivative are also positive under the independent
coarse transfer. Every exact rational endpoint is in the JSON report. These
bounds inherit the earlier first-field, background, coefficient, source,
zero-wedge and formal-response assumptions; they are not a new numerical C2
error certificate or a full nonlinear spacetime construction.
For the zero TT Robin rows, this shared-law radial error is identically zero:
both the actual formal response and the reconstruction have inward derivative
equal to the same unknown actual `G_TT`, enclosed by the inherited source cap.

## Reproduction

```
python check_robin_coarse.py --repo-root /path/to/HDblast
python -O check_robin_coarse.py --repo-root /path/to/HDblast
```

By default the authenticated prior capsule is
`research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030` beneath the root. The optional
`--contact-capsule` argument selects a separate local copy. All checks are
explicit, survive optimized Python and produce byte-identical output. The
program imports only the Python standard library and does not execute any
scientific dependency code. Hash authentication preserves the identity of the
previous inputs; it does not independently re-prove their conditional field
certificates. No external specialist or proof-assistant review is claimed.
