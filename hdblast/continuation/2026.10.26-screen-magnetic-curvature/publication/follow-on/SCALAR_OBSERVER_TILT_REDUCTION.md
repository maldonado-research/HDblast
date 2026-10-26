# Supplement: the scalar observer-tilt correction reduces to two local coefficients

This is a separate verification follow-on. It does not modify the frozen
60-file scientific capsule or the completed observer-anchoring note. It uses
no new actual-field numerical campaign and makes no new sign, full-solution,
or gauge-invariance claim.

For the same metric family at the same event, compare two oriented adapted
frames with the same background observer/radial pair. Restrict their relative
first variations to the standing-quadrupole scalar form

```
Delta U1=b0 Y n0+V v^a ea,
Delta n1=b0 Y U0+W v^a ea,
Delta(ea)1=V v_a U0-W v_a n0,
v_a=D_aY,  Y=n0_ambient^2-n1_ambient^2.
```

Here `b0,V,W` are arbitrary smooth functions of physical `y,t`; relative screen
rotation is matched. With the accepted orientation and normalized
symmetrization, the exact all-point reduction is

```
Delta H'' = 2(V P-W Q) A = 4(V P-W Q)(Tplus-Tminus),           (1)
A_ab=epsilon_(a^{cd}v_c D_{b)}D_dY,
P=(M_tt+theta M_t)/(2rho^2)-E_t/(rho R)+(D+L)/(2R^2),          (2)
Q=(M_yt+theta M_y)/(2rho)-E_y/(2R)
  -(C_t+theta C)/(2rho R)+B/(2R^2),                           (3)
R=rho cosh(t),    theta=tanh(t).
```

The radial/time boost coefficient `b0` drops out identically. These are physical
`y,t` derivatives of the same registered first metric amplitudes `A,...,M` as
in the frozen coefficient derivation. A generic scalar first metric has
nonzero `P,Q`: the correction has not been discarded on the grounds that
`H1=0`.

There is a further **exact correlated-amplitude identity** for the particular
registered first sector: substitution of its full scalar master plus gauge
representation gives `P=Q=0`, without using field equations or evaluating
actual first-field values. Consequently (1) is zero for that sector and this
specified scalar observer-tilt class. This conclusion is narrower than a
certificate of gauge invariance or of a completed geometric clock; its precise
scope is listed below.

## Reduction of the complete frame variation

Let `C1` be the full linear Weyl tensor in the common background orthonormal
frame, with `U=U0`, `n=n0`. Because `C0=0`, second frame variations do not enter.
The observer-change identity from the anchoring note expands to

```
Delta H''_ab = epsilon_(a^{cd} {
    b0 Y C1_{b)ncd}
   +V v^e C1_{b)ecd}+2V v_c C1_{b)UUd}
   -W v_{b)} C1_{nUcd}-2W v_c C1_{b)Und} } .                 (4)
```

The factor in (4) is already the full second parameter derivative: the cross
factor two cancels the one-half in the definition of `H`. All screen/epsilon
variations are included by using oriented adapted orthonormal frames.

For a linear scalar perturbation on this background, angular Codazzi and
constant-curvature Hessian commutation give

```
C1_{nUcd}=0,
C1_{bncd}=z(y,t)(kappa_bc v_d-kappa_bd v_c).
```

For example, the first identity follows from the antisymmetric angular curl
of `K_nd`: its scalar first variation is a gradient, while the remaining
first extrinsic-curvature products are symmetric in `c,d`. The radial Codazzi
identity gives the second form; its Weyl subtractions have the same
metric-times-gradient structure. Their symmetric epsilon contractions vanish
exactly, proving the `b0` cancellation and removing the term with `v_b`.

Every algebraic Riemann tensor restricted to a three-dimensional screen has
a unique screen Schouten form. Write

```
C1_{becd}=(kappa wedge S)_{becd},
S_bd=(trace coefficient)Y kappa_bd+S_H D_bD_dY,
C1_{bUUd}=(trace coefficient)Y kappa_bd+E_H D_bD_dY,
C1_{bUnd}=(trace coefficient)Y kappa_bd+Q D_bD_dY.
```

Direct epsilon contraction gives

```
epsilon_(a^{cd}v^e(kappa wedge S)_{b)ecd}=-2S_H A_ab,
epsilon_(a^{cd}v_c C1_{b)UUd}=E_H A_ab,
epsilon_(a^{cd}v_c C1_{b)Und}=Q A_ab.
```

Thus (4) is `2[V(E_H-S_H)-WQ]A`; `P=E_H-S_H`. The exact inherited angular
identity `A=2(Tplus-Tminus)` fixes the normalization. This is an all-point
contraction identity, not a fit to the direct fixtures below.

## Local Gauss–Codazzi calculation of P and Q

For a check of the metric-to-curvature coefficients, use positive time
extrinsic curvature and put

```
H=theta/rho,       a=rho_y/rho,
m=M_t/(2rho)-E/R,
j=M_y/2-C/R,
eta=R(C_t-theta C)/(2rho)-R E_y/2-B/2.
```

In the physical adapted frame, `m` is the Hessian coefficient of the first
screen time extrinsic curvature, `j` is the Hessian coefficient of the first
radial screen shape, and `K(n,ea)_1=(eta/R)v_a`.

The all-angular **Riemann** screen Schouten Hessian coefficient is

```
S_R,H=-L/(2R^2)-a j+H m.
```

The term `-L/(2R^2)` is the intrinsic screen-curvature variation; the `M HessY`
part of the angular metric is an intrinsic infinitesimal sphere
reparametrization and has no physical-frame intrinsic-curvature anisotropy.
The two remaining terms are the radial/time Gauss products.

The Riemann component with ordering `R_{bUUd}` has Hessian coefficient

```
E_R,H=m_t/rho+2H m+D/(2R^2)-a j.
```

The lapse Hessian supplies `D/(2R^2)-a j` with this index ordering.
Weyl Schouten subtractions cancel in `E_H-S_H`, including their metric
variation terms: the background tangential Schouten tensor is proportional
to the de Sitter metric. Consequently

```
P=E_R,H-S_R,H=m_t/rho+H m+(D+L)/(2R^2),
```

which expands to (2). The radial/time mixed Codazzi component gives

```
Q=m_y-eta/R^2+a m+H j,
```

which expands to (3); its relevant Weyl subtraction is a pure angular trace.
These calculations keep the complete seven scalar metric slots until the
cancellations. They are kinematic and use no first-field equation.

## Exact substitution of the registered first metric

To avoid confusing the gauge time amplitude with the relative boost, call
the registered gauge functions `f,b_g,d`. In ordinary physical derivatives,
the accepted first metric is

```
A=-4u+2 f_y,
B=f_t/rho-rho (b_g)_y,
C=f/R+R d_y,
D=-2(u+a f+(b_g)_t),
E=cosh(t) d_t-b_g/cosh(t),
L=2(u+a f+theta b_g),
M=2d.                                                       (5)
```

The `u` terms are the registered scalar master block. The other terms are its
registered scalar Lie-gauge block. Their common functional dependence matters;
independently widened numerical boxes should not be treated as independent
functions that violate these exact identities.

For `P`, use

```
D+L=2(theta b_g-(b_g)_t),
M_tt+theta M_t=2(d_tt+theta d_t),
E_t=cosh(t)(d_tt+theta d_t)
    +(theta b_g-(b_g)_t)/cosh(t).
```

Substitution into (2) cancels both pairs exactly. For `Q`, use

```
M_yt+theta M_y=2(d_yt+theta d_y),
E_y=cosh(t)d_yt-(b_g)_y/cosh(t),
C_t+theta C=f_t/R+R d_yt+2R theta d_y,
B=f_t/rho-rho(b_g)_y.
```

All `d_yt,d_y,f_t,(b_g)_y` contributions cancel in (3). Thus

```
P[(5)]=0,        Q[(5)]=0.                                  (6)
```

Both the pure master block and the pure gauge block satisfy (6) separately.
The latter is consistent with the general identity that a pure infinitesimal
diffeomorphism has `C1=Lie_xi C0=0` when `C0=0`; the proof here is the explicit
amplitude substitution, not an assumption that every scalar perturbation is
pure gauge. No master/gauge PDE, scalar-field equation, numerical midpoint, or
error-bar cancellation is used.

## What this establishes for the candidate clock, and what remains open

If the extra marked-shell/clock data are perturbed in the same standing
quadrupole sector, the anchoring note's candidate relative directions have

```
V=-(rho/R) tau1_Y,       W=s/(pR),
b0=B-s_t/(rho p),
```

where `tau1=tau1_Y(y,t)Y` after including the effective boundary clock datum.
Although the clock transport and its boundary value determine `V`, the
same-event correction (1) vanishes for every such `V,W,b0` when the exact
registered representation (5) holds. The local cancellation therefore removes
this particular correction as a numerical unknown for that tilt class.

This does **not** construct the marked section, establish a timelike clock or
valid transport collar, choose the physical target sphere, or certify its
perturbed embedding and angular labels. It does not extend the conclusion to
arbitrary angular tilt fields, distinct background observers, general scalar
metrics, or arbitrary gauge changes. An arbitrary-tilt extension may be
investigated separately; it is unnecessary for this result and is not being
claimed here. The full covariance/gauge-transformation statement for an
augmented observer construction needs its own complete audit. The frozen
checkpoint's conservative scope remains in force.

The result concerns the given formal response and exact first-metric structure.
It supplies no new second metric, evolution theorem, detector model, physical
energy, cosmological interpretation, or actual-field curvature sign.

## Independent exact controls

`tilt-audit/check_scalar_tilt.py` verifies the angular contraction on all
18 independent symmetric-Hessian/vector basis pairs, the radial and trace
cancellations, and formal rational Laurent-polynomial substitution of (5).
Its normal and optimized reports agree.

`tilt/check_direct_tilt.py` is a separate standalone standard-library oracle.
It constructs coordinate metric two-jets at the rational sphere point
`(3/5,0,4/5,0)`, computes the inverse metric, Christoffels, Riemann, Ricci,
Schouten and Weyl with exact rational parameter dual numbers, then differentiates
the full adapted-frame projection. It imports no ADM reduction or frozen
coefficient code. Three fully populated synthetic scalar metric two-jets each
exercise radial boost, angular observer tilt, radial-screen tilt and combined
tilt, checking all nine screen components. Their `P,Q` are explicitly required
to be nonzero, so these controls reject a generic-zero shortcut. The background
Weyl vanishes in every fixture. The exact direct-curvature outputs agree with
(1) in both Python modes.

The finite fixtures independently check the signs and normalization of the
all-point analytic derivation; they are not represented as an all-point proof
by themselves. Both sets of checks use synthetic or formal data, not actual
first-field values, and constitute internal mathematical/code review rather
than external specialist review.
