# Independent screen magnetic Weyl audit

The feasibility gate passes for the explicitly registered coordinate
prescription below. The full symmetric screen tensor has zero first
variation for every first scalar harmonic metric perturbation. The
registered standing quadrupole does allow a quadratic magnetic TT term.
The complete quadratic coefficient agrees with an independent direct
five-dimensional curvature calculation, including the normalized time
observer and screen volume corrections. This audit supplies no geometric
clock transport, invariant observable, full nonlinear family, detector
interpretation, or cosmology claim.

## Prescription and exact Codazzi reduction

Use signature `(-++++)`, the curvature convention
`R^a_bcd=partial_c Gamma^a_db-partial_d Gamma^a_cb+...`, and positive
extrinsic curvature `K=(dot gamma-L_beta gamma)/(2N)`. The background is

```
g0=dy^2+rho^2[-dt^2+cosh(t)^2 kappa],
R=rho cosh(t), theta=tanh(t), H=theta/rho, a=rho_y/rho.
```

On the perturbed metric take `U` to be the future unit normal of registered
`t`; within that time slice take `n` to be the positive unit spatial
normal of registered `y`. Explicitly,

```
U=N^-1(partial_t-beta^i partial_i),
n=gamma^{iy} partial_i/sqrt(gamma^{yy}).
```

Thus `U.n=0` and the screen is tangent to the coordinate spheres
`t,y=constant`. In general `n` is neither normalized `partial_y` nor the
normalized full spacetime gradient of `y`. The latter can have nonzero
inner product with `U`. These distinctions matter when shift is present.
This is a registered coordinate diagnostic; replacing it with a scalar
gradient and shell clock prescription requires a separate construction.

Let `s` be the screen metric, `chi_ab=gamma(D_a n,e_b)`, and
`v_a=K(n,e_a)`. In screen coordinate slots the exact identity is

```
H_ab = epsilon_(a^{cd} D^(4)_c K_{b)d}
     = curl_s(K^S)_ab+epsilon_(a^{cd} chi_{b)c} v_d.
```

Indeed `R_{bUcd}=D_c K_{bd}-D_d K_{bc}` with these conventions. Every
Schouten term in `C_{bUcd}` contains either `g_bU=0` or a screen metric
times `P_Ua`; its epsilon dual is antisymmetric in `a,b`. Symmetrization
therefore cancels the full Ricci and scalar subtraction exactly, without
using any field equation. The `chi v` term follows by decomposing the
four-dimensional spatial connection into its screen connection and
normal second fundamental form. The sign is positive for the stated
definition of `chi`.

An especially useful exact rewrite is

```
H=curl_s(K^S-Hs)+epsilon_(a^{cd}(chi-a s)_{b)c}v_d,
```

because `curl_s(f s)=0` for every scalar `f`, and
`epsilon_(a^{cd}s_{b)c}v_d=0` identically.

## Entire first-order tensor and second-order block isolation

The conformal product `dz^2+ds_dS4^2` is conformally flat, so `C0=0` for
every positive smooth `rho`. At first order,
`K^S-Hs=f(y,t)Y kappa+g(y,t) Hess_kappa Y` for an arbitrary scalar
harmonic first metric, including all radial/time scalar mixed amplitudes
and gauge-generated scalar amplitudes. The background `chi-a s` and `v`
both vanish. The symmetrized curl of a trace term vanishes; the curl of
the Hessian reduces, by the round-sphere commutator, to a metric times
`DY`, whose symmetric epsilon dual also vanishes. Consequently

```
H1=0 as the complete symmetric screen tensor.
```

This is stronger than zero TT projection. It needs no scalar or matter
equation and is not a claim that the full three-index observer magnetic
Weyl, or full `C1`, vanishes. Since `C0=0`, arbitrary first corrections of
the observer and projector cannot affect this first-order conclusion.
At second order their action on generally nonzero `C1` must be retained.

The linear second-metric map is an `SO(4)` intertwiner. Scalar harmonics
carry `(j,j)`, transverse-vector harmonics have `|jL-jR|=1`, and TT
harmonics have `|jL-jR|=2`. There is no intertwiner from the scalar or
vector blocks into TT. Thus every scalar and vector second metric
amplitude is removed by TT projection, including lapse, mixed radial
components and shift. Epsilon changes reflection parity without changing
these `SO(4)` labels. For the TT metric coefficient `u T`, the linear map is

```
L_H[uT] = u_t/(2 rho^2 cosh(t)) curl_kappa(T).
```

Only the first time derivative of the TT response is needed. No second
radial/time derivative of the second response is introduced.

## Independently derived complete scalar quadratic coefficient

Use the first metric amplitudes in the background orthonormal frame,

```
q_yy=A Y, q_y0=B Y, q_ya=C D_aY, q_00=D Y,
q_0a=E D_aY, q_ab=L Y kappa_ab+M Hess_abY.
```

For the registered `Y=n0^2-n1^2`, let
`V_a=D_aY`, `J_ab=Hess_abY`, and
`Aax_ab=epsilon_(a^{cd}V_c J_{b)d}`. Along the affine first metric path,
with second parameter derivative convention, the full result is

```
Q_H[q,q]_ab=(2F/R) Aax_ab,

F=[(D+4M-L) M_t-2M L_t]/(4rho)-E(D+L)/(2R)
  +(C/R-M_y/2)[R(C_t-theta C)/(2rho)-R E_y/2+B/2].
```

All amplitudes and their derivatives in this expression are physical
`y,t` derivatives. In particular `partial_t=ell^-1 partial_x`.

The apparent positive `B/2` in the final bracket is not the normal
extrinsic-curvature component. The independently derived latter is

```
v1_A=eta D_AY,
eta=R(C_t-theta C)/(2rho)-R E_y/2-B/2.
```

The bracket in `F` is `eta+B`. The extra `B` comes from the angular
`K2` connection contribution, which combines with the radial screen term.
This resolves the most plausible sign ambiguity in the reduction.

For further independent detail, write

```
alpha=(L_t+2theta L)/(2rho)+HD/2-aB,
beta =(M_t+2theta M)/(2rho)-E/R,
l=alpha-HL, m=beta-HM,
zH=C/R-aM-M_y/2, s0=(L-2M)/2.
```

With `K2` meaning the coefficient of `epsilon^2`, the four contributions
to `F` before cancellation are

```
curl(K2):              D beta/2+B zH-(E/R)(L-2M),
screen connection:    -s0 m,
epsilon variation:     M(-l+m),
radial screen term:    (C/R-M_y/2) eta.
```

Their sum reduces to the compact expression above. First amplitude `A`
and background radial expansion `a` cancel exactly. No assumption that
the quadratic term is zero is made. The epsilon term is essential even
though the complete symmetric `H1` vanishes: the unsymmetrized first
Codazzi contraction need not vanish.

## Exact angular conclusion

The actual ambient seed has eigenvalues `(1,-1,0,0)` and is not the
`SO(3)`-axisymmetric quadrupole with spectrum proportional to
`(3,-1,-1,-1)`. Its symmetric square contains the TT representations:

```
Sym^2(1,1)=(0,0)+(1,1)+(2,2)+(2,0)+(0,2).
```

An orientation reversing reflection fixes this seed and exchanges the
two frozen tensor modes. It requires opposite magnetic coefficients;
it does not require the magnetic response to vanish. Independent
all-point algebra in `symmetry/` proves

```
Pi_TT(DY tensor DY)=-(2/3)(Tplus+Tminus),
curl(DY tensor DY)=-2Tplus+2Tminus,
Aax=2(Tplus-Tminus).
```

The latter two are complete tensor identities, with zero trace and
divergence. The independent angular audit in `../angular/` finds
`curl(Tplus)=3Tplus`, `curl(Tminus)=-3Tminus` for outward orientation
`epsilon_ijk=n^A EPS_Aijk`, `EPS_0123=+1`; both norms are `1/6` under the
normalized sphere mean. Thus, for equal registered TT response `u`,

```
Pi_TT H2 = h(Tplus-Tminus),
h=3u_x/(2ell rho^2 cosh(t))+4F/R.
```

Reversing the fixed orientation reverses every magnetic sign. It leaves
nonvanishing and norm statements unchanged. No helicity interpretation
is inferred from the template names.

## Direct curvature verification and limits

`check_direct_curvature.py` has no project imports and uses
`fractions.Fraction` throughout. It builds a genuine coordinate metric
two-jet at a rational stereographic point with curved nonconstant `rho`,
nonzero background `theta`, and every one of the seven scalar metric
amplitudes and all six independent ordinary `y,t` two-jet entries
populated. It computes the five-dimensional inverse metric, Christoffel
symbols and their derivatives, Riemann, Ricci, scalar, Schouten and Weyl
directly through second parameter order. It separately normalizes `U`
and the induced screen epsilon through that order.

The affine-first fixture checks all 625 background Weyl slots, all
components of `H1`, exact symmetric Schouten cancellation, and exact
agreement of the full direct `H2` with `(2F/R)Aax`. The combined fixture
adds `epsilon^2 k/2`, with `k=R^2 u(Tplus+Tminus)` and nonzero
`u,u_t,u_y,u_yt,u_tt`. It reconstructs the frozen templates from their
ambient Weyl constants and pulls them back to stereographic coordinates.
Its full direct curvature agrees exactly with the quadratic expression
plus the linear coefficient `3u_t/(2rho^2 cosh(t))` multiplying
`Tplus-Tminus`. This verifies the orientation, time factor and second
parameter derivative convention independently of the Codazzi reduction.

Reasoned omission controls reject dropping `U1`, dropping the first
screen epsilon correction, and doubling the linear second metric term
by losing the `epsilon^2/2` normalization. Normal and Python optimized
reports agree byte for byte for both fixtures:

* `DIRECT_NORMAL.json`, `DIRECT_OPTIMIZED.json`;
* `DIRECT_COMBINED_NORMAL.json`, `DIRECT_COMBINED_OPTIMIZED.json`.

These direct checks are finite generic witnesses. The universal claims
rest on the exact Codazzi and scalar reduction above and the independent
all-point angular identities. No finite witness is presented as a proof
for arbitrary radial/time data. No actual source/response campaign or
frozen input was edited or rerun. A bounded actual coefficient still
requires the authenticated actual first jets and TT velocity enclosure.

The first-order identity would support a second-order gauge argument
only after a genuinely transforming geometric observer, radial
prescription and anchoring have been constructed. The present explicit
coordinate prescription is sufficient for a conditional local curvature
diagnostic and does not discharge that additional requirement.

The packaged symmetry checker accepts only a fresh, nonsymlink external
output with an existing real parent and creates it exclusively. Its
preflight rejects the source directory and any containing or targeted
`research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026` capsule. The
standalone control runner reproduces the two scientific reports, all six
scientific rejection controls, and twelve early output-guard controls
without writing into source. These receipts were regenerated from a
separate scratch copy of the source. The direct-curvature checker emits
only stdout. All review code and selected proof receipts are pinned in
`INDEPENDENT_AUDIT.json`.
