# An explicit monopole preparation with a conditional exact-family completion

1 November 2026. Portable explicit preparation; retained source files are unchanged. The public source commit is `b51908d7e62e7ec80cef82438edd00923ff94723`.
The public source pins and the published canonical monopole derivation pin are
in `../gauge/SOURCE_PINS.json`. The analytic weighted inverse, background and shell
evolution results remain retained conditional inputs. This is a new choice
of nonlinear completion, not recovery of the previous uninstantiated weights.

**Assessment:** the construction works, with one momentum normalization
correction. If `tau2=pi2/rho0=diag(C,D,D,D)`, choose
`C=0, D=-QM/(6a)`, rather than `D=QM/(3a)`. It avoids a fourth-order weighted
boundary-value solve for the second monopole, while still using the retained
weighted IFT to construct the full exact constraint family. Other angular
second components and higher orders can continue to depend on the weights.

## 1. Specify the bump and the explicit second datum

Use exactly the prior first seed, background, outward proper-radial coordinate
`y`, probability-sphere average and derivative convention
`x_e=x0+e v_seed+e^2 B/2+...`. The normalized canonical mean Hessian sources
are the prior `QH=<D^2 H[v_seed,v_seed]>` and
`QM=<D^2(M/rho0)_y[v_seed,v_seed]>`. Thus the required mean equations are
`Pg h=-QH`, `Pm tau=-QM`; neither source is a half-coefficient.

Fix the bump explicitly, with no new shape parameter:

```
b(y)=Q(z(y)),
Q(z)=exp(4-ell^2/[(z-ell)(2ell-z)]) for ell<z<2ell,
Q=0 otherwise, ell=10^-6, dy=-rho dz.
```

Let `yL=y(2ell)<yU=y(ell)`. Define the finite radial integrals

```
I=integral_(yL)^(yU) rho^4 QH dy,
J=integral_(yL)^(yU) rho^2 b dy >0,
f=(2I/(9J))b,
C[f]=-3(f''+3af')/2-Vf,
xi(y)=-(rho^4 p^2)^-1 integral_(yL)^y rho^4(QH+C[f])ds.
```

Extend the integral/fields by the displayed compact-support construction.
The retained background gives `rho>=129`, `|p|>=3/25000` and
`111/1000<=a<=14/125`, so both divisions by `p^2` and `a` are legitimate.
The bump is strictly positive inside its support interval and flat at the
endpoints. The selected first packet is nonzero, so the previous strict
moment result also gives `I<0`, although only finiteness of I is needed here.

Set the radial canonical second datum `B=(h,0,rho0 tau,0)` by

```
h=f gamma0+Lie_(xi partial_y)gamma0,
h_rr=f+2xi', h_angular=f+2a xi,
tau_rr=0, tau_angular=-QM/(6a).
```

The scalar and scalar momentum-density sectors are zero. Every quantity is
a fixed smooth functional of the actual initial seed/background and the
explicit bump; a numerical quadrature may later evaluate it, but a singular
weighted radial inverse is not needed to define this mean datum.

## 2. Exact constraints, support and regularity

The previously derived radial identities are

```
Pg(f gamma0)=C[f],
Pg(Lie_(xi partial_y)gamma0)=rho^-4(rho^4 p^2 xi)',
integral rho^4 C[f]dy=-(9/2)integral rho^2 fdy.
```

Our f gives `integral rho^4(QH+C[f])dy=I-(9/2)(2I/(9J))J=0`.
Hence the primitive defining xi vanishes on both exterior sides of the
packet interval. The integrand is smooth and flat at both endpoints; xi
therefore extends smoothly by zero. Both h and tau have support in
`[ell,2ell]` in z, strictly inside the correction annulus
`A'={ell/2<z<5ell/2}`. All finite Sobolev norms required by the theorem are
finite, despite possibly large inverse powers of p in the numerical size.

Differentiating the primitive gives `Pg h=-QH`. For the normalized momentum,

```
Pm(C,D)=-2[C'+3a(C-D)],
Pm(0,-QM/(6a))=-QM.
```

The suggested `D=QM/(3a)` would instead produce `+2QM` and fail. Because a
radial scalar constraint is a scalar and a radial momentum constraint is a
radial one-form, these formulas establish `P B=-P0 Q` as full tensors, where
`Q=D^2 Chat[v_seed,v_seed]` and P0 is the rotational monopole projector.
No angular mean momentum component is left unaccounted for: rotationally
invariant tangential one-forms on S3 vanish.

## 3. The full weighted IFT leaves this second monopole unchanged

Keep any fixed admissible radial defining function and weights of the
retained theorem; do not assert their unspecified values have been recovered.
For its same fixed canonical correction operator A and isomorphism
`L=w^-2 P A`, replace the affine trial by

```
F(e,U)=w^-2 Chat(x0+e v_seed+e^2 B/2+A U),
U(0)=U'(0)=0,
L U''(0)=-w^-2(Q+P B).
```

Adding a fixed smooth compactly supported polynomial trial term preserves
the theorem's nonlinear mapping estimates: its finitely many scaled-cover
seed constants simply enlarge the local bound. The derivative in U at zero
is still L. The fixed-k IFT therefore applies, after possibly reducing the
parameter interval. The same support and zero-extension conclusions hold.

The needed averaging argument does **not** claim that each rotation preserves
the chosen coordinate-reflection invariant subspace. Instead:

1. On the full weighted tensor spaces, define P0 by Haar averaging the
   natural SO(4) pullback. The radial background and radial weights make P,
   its formal adjoints, A and L rotationally equivariant. Thus they commute
   with P0 on compactly supported smooth fields.
2. Haar averaging is bounded in the weighted Sobolev norms. One can use the
   background covariant norms, where the action is isometric; equivalence
   with the theorem's scaled-cover norms gives the same bounded map. The
   concentric annulus is rotation invariant, so averaging preserves compact
   support and extends to the dotted closures by density.
3. P0 maps the reflection-invariant spaces into themselves because its
   output is fully rotationally invariant, hence invariant under the selected
   reflections as well. For radial scalar functions and radial one-forms,
   reflection invariance is immediate. The corresponding tensor/density
   statement uses their natural pullback and the fixed background density.
4. Consequently `L P0 U''=P0 L U''=0`, since
   `P0(Q+P B)=P0 Q+P B=0` and multiplication by `w^-2` commutes with P0.
   The retained L is injective on the reflection-invariant dotted domain,
   hence `P0 U''=0`. No adjoint-kernel term can be added: the source theorem
   already removes the noninvariant boost/rotation kernels and its invariant
   inverse is an isomorphism. Equivalently, the radial adjoint kernel is zero.

It follows that the exact curve satisfies

```
P0 x''(0)=B,  x'(0)=v_seed,
phi2=varpi2=0 in every angular sector.
```

The entire second monopole is independent of the remaining radial weight
choice. This assertion is only about this deliberately modified trial and
its second derivative. It does not make the full finite-parameter family,
its higher derivatives or nonmonopole second data weight independent.

## 4. Physical initial data and the same TT transfer

Let `s=ut/rho`, `j=vt/rho`, as in the prior mean derivation. Its exact
canonical conversion now gives the explicit mean physical velocity data

```
<Pi2>=-uj/3,
<K2>rr=QM/(3a)+6us, <K2>angular=-us,
<partial_t gamma2>rr=2rho[QM/(3a)+16us/3],
<partial_t gamma2>angular=-(4/3)rho us.
```

The selected initial lapse `N_e=rho(1+e uY)` and zero shift are retained.
Thus `phi2=partial_t phi2=0`, `<g2_tt>=-rho^2u^2/3`, and `g2_ti=0` initially,
as before. The remaining gauge velocities are supplied by the same exact
relative-harmonic completion, rather than chosen independently.

The added metric h is radial with an isotropic angular block, and the added
canonical momentum is likewise angular trace. They have zero angular TT
projection. Every weighted correction has zero angular TT projection by
the inherited adjoint argument. The quadratic canonical conversion depends
only on the unchanged first tangent at `pi0=0`; its angular term is pure
trace. Consequently both `P_TT gamma2` and `P_TT partial_t gamma2` remain
zero. This proves the hypotheses needed for the retained TT forced-response
identification, not equality of all second initial sectors with the old family.

The initial exact family still equals the background on the same open shell
collar, since all trial and correction supports remain inside A'. The same
harmonic construction supplies initial gauge and corner compatibility. The
same strict domain arithmetic, state neighborhood and prescribed-slab
theorem apply after shrinking the parameter interval. The full first field
is unchanged by linear uniqueness. Both copies can use the identical curve.
Thus, conditional on the same retained existence/regularity assumptions,
the registered TT transfer and first-only observables survive this new
preparation. No effective parameter radius is supplied.

## 5. Limits and next use

This explicitly selects mathematical second monopole data; it is not a
physical generator of the packet. It neither evaluates the monopole after
contact nor proves a shell-volume sign, heating or radiation history. A
subsequent calculation can propagate these stated initial mean profiles
using the retained coupled scalar response equations and evaluate the
anchored shell observables. Such a calculation must identify this new
completion by its explicit bump choice, rather than calling it the old
weighted-only initial family.

The gauge exact checkers in `../gauge/checks/` verify the new algebra, factor correction,
compact-support integral cancellation and physical mean conversion. The
SO(4)/dotted-space/IFT argument is analytic and is not certified by finite
symbolic checks. This bounded feasibility review is internal AI-assisted
work, not external specialist review or a current-publication claim.

The canonical source derivation is the repository-relative public file
`research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/weighted-mean/WEIGHTED_MEAN_INITIAL_DATA.md`.
It is used for the first tangent, canonical Hessian and radial identities,
not for selecting the new second mean. Source keys and byte hashes appear
in the adjacent gauge SOURCE_PINS.json.
