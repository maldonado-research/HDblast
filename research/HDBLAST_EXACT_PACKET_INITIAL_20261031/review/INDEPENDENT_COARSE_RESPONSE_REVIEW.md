# Independent response-premise audit

The displacement-only change preserves the mathematical premises used by the
retained physical-side C2 construction. Its completed displacement is C2,
its unchanged smoothed P1 velocity is at least C4, and both have the flat
initial shell collar. The proof below concerns the specified integral-defined
finite candidate and the retained formal coefficient target. It does not
assert that auxiliary numerical initial data solve the physical constraints
exactly.

The independent review read the original derivation and regularity review:

- `HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/theory/ORIGINAL_TARGET_C2_NOTE.md`
- `HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/review/INDEPENDENT_C2_REGULARITY_REVIEW.md`
- `HDBLAST_ROBIN_RECONSTRUCTION_20261030/theory/DIRECT_ERROR_NOTE.md`

The same physical-side linear maps, actual coefficient operator La, original
F/Fx, and affine G/M are retained. Put L0 for the frozen constant bulk
operator, P=(I-L M PN)^(-1), Dfield=W0-Bs and C=W1-W0. Linear subtraction,
without commuting any matrix through reflection, gives

```
r0=La Dfield+(La-L0)Bs+F,
C=P K r0,
r1=La C.
```

The construction uses common initial displacement and velocity in Bs, W0
and W1. Thus Dfield and C have zero initial displacement, time derivative
and tangential radial derivative. Since La is first order, r1(0)=0 follows
exactly. It requires neither F(0)=0 nor zero initial acceleration of a
correction. In fact

```
Dfield_xx(0)=L0 Bs(0),
C_xx(0)=r0(0).
```

The scalar test C(x,r)=x^2/2 independently exhibits the distinction: C and
its first jets vanish at x=0, whereas C_xx(0)=1. Differentiation of Kf
therefore retains its initial-velocity free-wave term. The independent
arithmetic checker retains f00 and r00 and rejects actual copied scripts
that drop either term by exact failure of the old complete differentiated
ledger reconstruction.

Bs satisfies the homogeneous free reflection law. Subtract it before
differentiating Dfield's shell law: the source is M PN Bs+G, and its
derivative is M PN Bs_x+G_x because M_x=0. C has zero affine shell source
because W1 and W0 use the same G/M. The target error also has zero affine
shell source. Thus G_x belongs in Dfield_x and vanishes from C_x and the
final error boundary law. The checker retains M*Bs and G_x in the former,
reproduces the retained exact old F1 ledger, and rejects a copied script
that drops the differentiated boundary source.

The exact exponential q is flat at r=1, and the nominal denominators are
uniformly separated from zero, so every required completed displacement
jet vanishes there. The old P1 velocity has zero first two node vectors and
is convolved with radius 1/2048, leaving it exactly zero through 2079/2048,
which exceeds one. Hence both initial functions vanish through r=1. Causal
K and L, local coefficient multiplication, and the Volterra inverse P
preserve the wedge x+r<=1. On the physical-side domain the same 3/4 value
factor, 1/4 boundary history factor and full contraction constant remain
available. The outer displacement extension is not convolved, so it does
not spread support below one.

The response construction needs displacement through rr and velocity
through r, together with the original continuous F_x/G_x and coefficient
commutator. It does not differentiate the new completed displacement a
third time globally. The independent regularity calculation proves C2 at
the outer join and explicitly declines global C3. The prescribed-forcing
time derivative constructions produce physical-side C2 candidate jets;
the same differentiated Picard argument gives summable n*q^n tails for the
target error. The principal equation then recovers rr continuously.
One-sided shell/domain traces suffice. No C2 extension through the
artificial reflected shell is assumed.

The new absolute candidate jets and the paired initial errors enter
different equations. `check_coarse_response.py` uses only the independent
48-leaf displacement certificate, independently reconstructs absolute
unchanged velocity caps from the pinned complete P1 table and final slope,
and retains complete unchanged velocity errors from the old C2 ledger.
It reads no new producer report. It first reproduces the retained exact
F1 and each retained C2 cap, then computes fresh global Dfield and correction
bounds and the differentiated ledger. Its conservative independent caps
still pass all response thresholds: approximately E=1.42466e24,
Bxx=Bxr=2.37685e26, Brr=2.38528e26. These are certified error upper bounds,
not lower bounds on the unknown true error or pointwise physical ratios.

Canonical current response artifacts are
`portable-response-normal/RESPONSE.json` and
`portable-response-optimized/RESPONSE.json`; they are byte-identical.
Current regularity artifacts are
`portable-regularity-normal/REGULARITY.json` and
`portable-regularity-optimized/REGULARITY.json`; they are byte-identical.
The source map itself is frozen by a hardcoded SHA256 in the checker.
Actual copied-script formula controls are in `response-controls-isolated-normal/`
and `response-controls-isolated-optimized/`; all six fail explicitly without a
successful report, and their canonical control reports are byte-identical.

Fresh finite-response point/quadrature/constraint tables remain a separate
required downstream calculation. The old shell point (9/8,0) cannot be
copied because removal of displacement convolution changes its initial
footprint. These analytic and arithmetic checks do not independently
evaluate such tables or establish invariant observables, a nonlinear
shell-contact family, or a cosmological interpretation. The previously
retained conditional exact initial/precontact theorem is not declared
absent.
