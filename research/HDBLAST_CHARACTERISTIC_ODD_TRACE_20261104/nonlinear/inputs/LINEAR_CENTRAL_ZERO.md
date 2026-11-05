# A continuous-cell refinement of the conditional past zero bracket

4 November 2026, America/Los_Angeles. This component is separate from the
frozen past-shell-contact paper. It uses the companion past–future theorem
and a retained continuous-cell certificate; it does not execute that old
certificate's producer or independently reprove its analytic premises.

## Statement

Let `ell=10^-6` and let `S(t)` be the embedding-included first physical
scalar amplitude of the same registered quadrupole packet. Retain the
same-state, clock, branch, regularity and linear uniqueness premises of
the [companion proof](../primary/PROOF.md), including

```
|S(-ell*x)-S(+ell*x)| < B,
B = 264000000001/13738809600 < 20,       0<=x<=3.
```

Retain the actual continuous-cell enclosures of the wider inertial
waveform certificate. Then

\[
 S(-6145\ell/4096)>63,\qquad S(-6143\ell/4096)<-53.
\]

Consequently, continuity gives at least one zero of the actual first
physical scalar **amplitude** in

\[
 \boxed{(-6145\ell/4096,\,-6143\ell/4096)}.                 \tag{1}
\]

The open bracket has width `ell/2048`, exactly 1024 times narrower than the
previous two-lobe bracket of width `ell/2`. Its scaled endpoints are
`-1.500244140625` and `-1.499755859375`. This comparison is to that coarse
past bracket, not to the older, narrower future zero bracket. No
uniqueness, optimality or count of all zeros is asserted.

## 1. The retained cells enclose the actual physical amplitude

The selected old primary proof's Section 4 says its CSV enclosures apply
on complete closed time cells, rather than only sampled times. Its
actual-minus-comparator error has already been included before the
physical gain and `cosh(t)^(-3/2)` factors are multiplied outward.

The old JSON report binds the selected CSV by exact byte count and SHA256.
It declares 4353 contiguous cells covering `0<=x<=3`, with time denominator
`D_t=4096` and value denominator `D_v=2^112`. For each CSV row, the stored
integer pair `(ell_S_lo,ell_S_hi)` therefore means

\[
 \ell S(\ell x)\in
 [\mathrm{ell\_S\_lo}/D_v,\mathrm{ell\_S\_hi}/D_v]
 \quad\hbox{for every real }x\hbox{ in that closed row cell}.       \tag{2}
\]

These are actual coupled first-variation enclosures, not the `H0`
comparator values. The wider certificate has precisely the registered
`ell=10^-6`, initial displacement `Q(z/ell)` and initial velocity
`exp(-z)Q(z/ell)`. The proper clock is `tau=rho_b*t`, anchored at the same
section. The companion proof supplies the complete same-state join.

In particular, the raw `ell_S` integers must first be divided by both
`D_v` and `ell`. The uniform difference bound is in model units of `S`,
not in units of `ell*S`.

## 2. Uniform interval transport preserves two nearby cell signs

If a retained row gives `S(ell*x)` in `[L,U]` throughout its closed cell,
the companion difference theorem gives

\[
 L-B<S(-\ell x)<U+B                                      \tag{3}
\]

throughout the reflected closed cell. The following two selected rows
have strict sign margins after this widening:

| CSV row, zero based | Future scaled-time cell | Safe strict enclosure for `S(-ell*x)` |
|---|---|---|
| 2047 | `[3071/2048,6143/4096]` | `(-187,-53)` |
| 2050 | `[6145/4096,3073/2048]` | `(63,197)` |

All endpoints in this table follow from exact integer/Fraction arithmetic
on the stored actual `ell_S` intervals. The canonical report records the
full rational enclosures. Their approximate endpoints are
`[-186.278268,-53.884096]` and `[63.097981,196.439531]`, respectively;
these decimals are explanatory only, not decision inputs.

The first row is negative also at its closed right endpoint
`x=6143/4096`. The second is positive also at its closed left endpoint
`x=6145/4096`. Time reflection reverses their order. Hence the two strict
amplitude signs in the statement hold, and the intermediate value theorem
proves (1). No derivative or monotonicity estimate is used.

The script additionally scans all retained cells to identify these as the
nearest strict-sign cells on the two sides of `x=3/2` under this particular
whole-cell widening rule. The two intervening whole-cell intervals contain
zero after transport. This is a statement about those retained row bounds,
not an optimal localization theorem or an exclusion of sharper methods.

## 3. What this refinement does and does not establish

The first amplitude is continuous under the retained smooth dynamic-shell
development. The angular perturbation is `Phi_1=S(t)(n_1^2-n_2^2)`; its
angular factor has zero mean and nodal directions at other times. The
zero in (1) is a temporal zero of the amplitude after possible past shell
contact. It is not a claim that a background scalar, spacetime or physical
state disappears.

This bracket remains wider than the old future central zero bracket.
The old narrow future witness brackets and small negative late tail do
not follow from the uniform `B` transport; their margins are too small.
No globally signed angular field, sign of a second coefficient, physical
energy flux, incoming source, SI scale, finite nonlinear parameter window,
or cosmological-origin evidence is supplied.

The positive-operator/reflection/interval arguments are established
methods. This refinement is selected-certificate postprocessing for one
fixed packet. Broader originality and significance are not settled by
these checks. Earlier conditional formal future second Taylor-coefficient
existence results are not contradicted or replaced; quantitative particular
past data and nonlinear development remain separate questions.

## 4. Reproduction and retained premises

The three selected public sources are included in `inputs/` and listed in
`SOURCE_PINS.json`, with exact original repository paths and Git blob IDs.
`verify.py` authenticates them, checks the report's grid/value scales and
complete CSV roster, widens actual physical intervals exactly and verifies
the strict signs and bracket width. Its normal and optimized reports must
agree byte for byte. It imports no project scientific module and executes
no legacy producer, checker campaign, numerical PDE or field/source/bump
oracle.

Byte authentication does not prove the old continuous-cell theorem, its
background estimates, coupled reconstruction or the companion reflection
bound. Those remain explicit conditional mathematical premises. Replaying
this postprocessor is not a fresh certification of the older full
calculation campaign.
