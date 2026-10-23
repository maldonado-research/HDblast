# One-point boundary feedback: independent analytic audit

**ACCEPTED AS A CONDITIONAL ANALYTIC FOLLOW-ON FORMULA.**
`PROPOSED_UNVERIFIED_FOLLOWON`: no boundary source integral or response value
has been evaluated, and no sign, physical constraint, or novelty claim follows.

Audited file: `ONE_POINT_FEEDBACK_ENCLOSURE_DRAFT.md`, SHA256
`86caf742fc1c8b766f5c3b14ad63e146c22b09d6600a0c7d65dfb9f73be9f96e`.
Audit date: 23 October 2026. The proof below independently checks the displayed
constants against the retained operator definitions; a pure rational-arithmetic
check verified every displayed draft equality. No producer was run, no source
values were evaluated, and no frozen input was edited.

Use the exact roots and input pins recorded in `NEXT_BOUNDARY_CALCULATION.md`.
In particular C denotes the retained `FULL_SECOND_RESPONSE_CONSTRUCTION.md`
and B denotes `FORCED_INCREMENT_BOUND.md`.

## Domain, operator, and trace audit

At (a,b)=(3/2,1/4), let d=5/4 and
`D+={0<=s<=a,u>=0,|u-b|<=a-s}`. If `(t,v)` belongs to the backward
physical cone of `(s,u)` in D+, then
`|v-b|<=|v-u|+|u-b|<=a-t`; thus it lies in D+. A boundary history used
at `(s,u)` has time tau<=s-u<=a-b=d. The backward cones of all those
boundary points also lie in D+. In particular the reflected-source integrals
and recursive Robin traces are fully covered. All points satisfy s+u<=7/4<8.

Restrict the SAME actual solution and fixed operator to this domain. With
`||V||_X=max(sup|V|_H/3,sup|V_x|_H,sup|V_r|_H)`, the reflected K has
derivative bound `a sup|f|_H` and value/3 bound `(a^2/6)sup|f|_H`.
Since a^2/6<=a, its C0-to-X norm is <=a. The reflection J=P_N-P_D is
an isometry; no extra reflection multiplicity is needed in that bound.
The negative lift L has derivative bound sup|g| and value/3 bound
`(d/3)sup|g|`, hence norm <=1. The boundary value trace of V has norm
<=3||V||_X (C:208–255).

The full bulk combination is bounded by
`(10^-4+10^-4+3*10^-6)||V||_X=sigma||V||_X`, where sigma=203/10^6.
The Robin map has norm <=m=1/2000. Consequently

    q_D=a sigma+3m=3609/2000000,
    ||S||_X,D <= a F_H+G_H=52162547172519/2,
    B_D=(a F_H+G_H)/(1-q_D)
       =52162547172519000000/1996391.

Here F_H=17387513156573 and G_H=3851400 are the retained full-state source
caps (B:11–47); restricting their domain preserves them. This argument does
not assume modal decoupling, a new source error estimate, or a smaller global
norm for the original solution. It applies the retained norm estimates to a
closed dependency region.

## Point allowance and exact support audit

For `f=A W_x+B W_r+D0 W` and `g=M P_N W|boundary`, we have
`|f|_H<=sigma B_D` and `|g|_H<=3m B_D`. The sum of the unfused
direct and reflected cone areas is a^2. Therefore

    |W(a,b)-S(a,b)|_H
      <=[(a^2/2)sigma+3dm]B_D
       =877739181271977213/15971128.

The exact fixed zero wedge x+r<1 is a standing retained premise, not a
numerical observation (C:219–223; B:25–27). W and its first derivatives
vanish on its interior, and continuity supplies the boundary values needed
for the estimates. Thus f=0 when s+u<1 and g=0 for s<1.

Because d>1, the triangle `s>=0,u>=0,s+u<1` is contained in D- and hence
in D+. Its area is 1/2. Removing it from each unfused cone subtracts total
area 1, so K's remaining norm weight is `(a^2-1)/2=5/8`. The boundary
lift's remaining length is d-1=1/4. Thus the draft's sharpened bound is valid:

    E_D=[(5/8)sigma+(3/4)m]B_D
       =209432626897663785/15971128
       =13113201954.030033... .

The minus sign of L and the sign of M matter when integrating S and the
feedback itself. Triangle-inequality bounds depend only on their magnitudes;
this audit does not change either sign or place an exponential inside L.
All Robin feedback remains in T.

The covariant yy probability mean has H-dual norm one. Accordingly an actual
source interval `[S_minus,S_plus]` implies
`lambda(W(a,b)) in [S_minus-E_D,S_plus+E_D]`. No further factor three or
angular point-evaluation factor belongs here: E_D is already a value estimate.

## Optional improvement using the same premises

The actual W is C1 even though it is not numerically reconstructed. On the
boundary segment, W(1,0)=0 and |W_x|_H<=B_D, so the fundamental theorem
of calculus gives

    |W(s,0)|_H <= (s-1)B_D,  1<=s<=d.

Integrating this bound yields the smaller Robin value allowance
`m(d-1)^2 B_D/2`, instead of `3m(d-1)B_D`. Keeping every other estimate
unchanged gives the also-valid optional bound

    E_trace=[(5/8)sigma + m(d-1)^2/2]B_D
           =14866325944167915/3992782
           =3723300181.2189884... .

The segment [1,d] at u=0 lies in D+, so no outside derivative estimate is
used. This improvement requires only the already-assumed C1 actual solution,
its exact zero wedge, and the same B_D. It needs no source derivatives or
numerical C1 reconstruction. It is a further ordinary norm consequence,
not a newly evaluated response.

## Proof delta relative to the first follow-on specification

The earlier four-term plan correctly needs rigorous C1 reconstruction or
equivalent continuous error control to evaluate feedback iterates. Such a
reconstruction is NOT required for the cheaper first point test
`actual S(a,b) interval + analytic bound on all feedback`. The point test
can use exact clipped source-cell areas and shell interval lengths, with
all inherited first-field and source uncertainties retained. Its result may
contain zero; no sign is promised. Existing conditional analytic completion
of the physical total U remains distinct from this auxiliary W enclosure
and from an evaluated full numerical waveform.
