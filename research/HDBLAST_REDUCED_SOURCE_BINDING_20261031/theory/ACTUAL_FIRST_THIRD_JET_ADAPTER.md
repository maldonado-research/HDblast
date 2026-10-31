# Actual first-field ordinary third-jet adapter

31 October 2026. This additive implementation returns derivative enclosures
of the same retained first fields. The raw actual master/gauge higher-jet
certificates and their analytic solution/regularity premises already exist.
The old numerical adapter and every earlier source map and candidate remain
unchanged. No smaller approximation error or new regularity theorem is claimed.

## Result and interface

`Adapter(repo_root).actual_jet(x,r,clip_causal=False)` returns the current
metric profiles `A,B,C,D,E,L,M` and scalar profile `s` with ordinary physical
suffixes

```
'', y, t, yy, yt, tt, yyy, yyt, ytt, ttt.
```

There are 70 metric and ten scalar slots. Lower two-jets are returned directly
from the authenticated original `actual_first_jets.actual_jet`, with their
complete interval endpoint bytes unchanged. The 32 third slots are coarse
uniform magnitude enclosures of derivatives of those same actual fields on

```
0 <= x <= 3, r >= 0, x+r <= 8,
x=t/ell, r=z/ell, ell=10^-6, dy=-ell*rho*dr.
```

The new finite part is the eight pure-radial `yyy` entries and a reusable
complete third-jet interface. All 24 time-containing third entries match the
previous original-bulk-Fx certificate exactly. Scalar third derivatives are
also supplied because the existing actual master-fourth certificate is
available; only scalar second derivatives are needed by the current natural
and reduced Bianchi source identity.

These are component enclosures for the actual retained functions, conditional
on their stated construction. Choosing unrelated values in the returned
intervals need not produce coherent jets, solve a PDE or define a continuous
field. The implementation does not differentiate independently selected
lower-jet interval endpoints.

## Exact physical chain rule

For any smooth composite profile F(x,r), put `Lrad=1/(ell*rho)`. The same
background satisfies `rho_y=a*rho` and `a_y=-rho^-2-p^2/3`. Thus

```
F_y   = -Lrad F_r,
F_yy  = Lrad^2 F_rr + a Lrad F_r,
F_yyy = -Lrad^3 F_rrr - 3a Lrad^2 F_rr
           +(a_y-a^2)Lrad F_r,
F_yyt = (Lrad^2 F_xrr + a Lrad F_xr)/ell,
F_ytt = -Lrad F_xxr/ell^2,
F_ttt = F_xxx/ell^3.
```

The coefficient `-3a` includes both differentiation of `Lrad^2` and the
derivative of the existing `a Lrad F_r` term. The final coefficient contains
both `a_y` and `-a^2`. Ordinary y and t derivatives commute because the
background coefficients are static. Here t is the same de Sitter coordinate
time used in the original adapter; conversion to an orthonormal time
derivative still includes the existing rho factor and frame connection in
the geometric operator. No covariant angular derivatives are returned here.

The new checker derives all ten formulas independently as operators
`D_y=-(ell*rho)^-1 D_r` and `D_t=ell^-1 D_x`, using
`rho_r=-ell*rho^2*a` and `a_r=-ell*rho*a_y`, then compares the collected
Laurent coefficients with the actual adapter conversion function. It also
checks every relevant permutation of the mixed y/t derivative words.

## Actual raw inputs and derivative inventory

The source expressions are the unchanged registered first frame amplitudes
from `HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py`:

```
u=R*vartheta*W, v=B*vartheta*(W_r/ell+kW),
A=-4u-2f_r/(ell*rho), Bmix=f_x/(ell*rho)+b_r/ell,
C=f/(rho*c)-c*d_r/ell, D=-2(u+a*f+b_x/ell),
E=c*d_x/ell-b/c, L=2(u+a*f+H*b), M=2d,
s=v+p*f.
```

The implementation calls that authenticated function directly with a new
fourth-order raw normalized-jet interface. It imports no old bulk-Fx
producer. Its raw interval slots are enclosures of actual derivatives from
the retained certificates:

| Output | Sufficient raw actual derivatives |
| --- | --- |
| Metric profiles through order three | Master through three; raw gauge through four |
| Scalar profile through order three | Master through four; raw gauge through three |
| Reconstruction/background coefficients | Ordinary spatial and temporal jets through three |

The actual raw master-fourth report is
`HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json`, SHA256
`d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e`.
The actual raw gauge-fourth report is
`HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/GAUGE4.json`,
SHA256 `94579087d39d9c849e9264d40bd1d71d64892a668fe3e2e08aade9ca4296b873`.
Their statuses, identical scaled domain/ell and common input hashes are
explicitly joined. Master orders zero to three use the same retained
`HIGHER_JETS_NORMAL.json` exact caps. The gauge report supplies all 45 raw
componentwise mixed derivative caps through order four.

Background coefficient derivatives come from the existing differentiated
actual ODE recurrences in `interval_core.background`. They are scaled by
`ell^n` when lifted from z to r jets. Temporal coefficient derivatives come
from the actual analytic cosh/tanh/vartheta recurrences already implemented
in `actual_first_jets.temporal`. These are differentiated exact-function
identities, not derivatives of a numerical remainder. The scalar's k-third
coefficient can use the existing exact polynomial potential through its
fourth derivative; no new potential trajectory is assumed.

Only coefficient jets through order three enter the returned amplitude jets.
Temporary coefficient-order-four placeholders are not used: the checker
replaces them by a nonzero sentinel and verifies all 32 third intervals are
unchanged. Raw gauge/master order-four inputs are used where required, and
a corruption replacing their caps by zero fails an authenticated shared-slot
comparison. Normalized jet coefficients are divided by `i!j!`; derivative
extraction restores those factors. No fourth amplitude output or missing
fifth raw derivative is used.

Combining uniform raw actual derivative caps over the causal triangle with
coefficient enclosures over the full x/r rectangle merely widens the output
enclosure. It does not assert the actual raw solutions exist outside the
triangle. With `clip_causal=True`, the original lower-jet assertion and the
new third-jet assertion both concern only the nonempty intersection with the
same triangle.

## Authentication, checks and limits

Before any historical numerical module is imported, the adapter authenticates
all 168 members of the closed geometric-source capsule against its fixed
manifest, rejecting symlinks, extra/missing files and changed bytes. The
older fourth-jet source capsules have no complete manifests; their explicitly
listed higher-jet reports and source scripts are authenticated by fixed hashes.
The geometric capsule contains the complete executed historical import/input
closure. The older master/gauge/Fx producers are not imported or replayed by
this adapter. Their higher-jet certificates remain stated inputs.

The checker authenticates the new adapter before import, verifies exact
chain-rule coefficients, normalized Taylor factors, all 24 shared actual
third intervals, all 32 coefficient-order-four sentinel invariances, and
lower-two-jet equality at interior, one-sided physical-shell and clipped
causal-edge interfaces. The latter cases test the interface; the all-point
enclosure follows from the uniform actual caps and exact finite-jet
composition, rather than sampling those cases. Invalid time/domain cells
are rejected. No smooth reflection of physical coefficients is assumed.

Ten deliberately wrong subprocesses per execution mode must fail the named
mathematical identity and leave no success report: pure-radial sign, omitted
middle/background radial terms, wrong `a^2` sign, omitted mixed rho term,
wrong mixed sign, missing time ell, missing Taylor factorial, omitted new
radial slot, and falsely discarded actual fourth-derivative inputs.
Normal and optimized complete report bytes must agree. The checks use only
standard-library exact rational arithmetic and explicit failures, not assert.

This adds derivative-availability plumbing for future source differentiation.
It neither evaluates a new bulk-source norm nor improves an approximation
error, proves a second-order evolution, constructs a finite-amplitude
nonlinear family, establishes a physical force, or supplies evidence for a
higher-dimensional Big Bang origin.
