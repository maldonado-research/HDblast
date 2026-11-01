# Known Gauss identity applied to the current scalar/tension shell

Self-contained application of the standard hypersurface Gauss equation.
Copied and annotated for this portable geometry component on 1 November 2026.
This is established hypersurface machinery, not a new theorem or matter model.
The model and background equations are bound by the repository-relative inputs
in `SOURCE_PINS.json`; no separate private note is a mathematical premise.

Use signature `(-++++)` and
`R^A_BCD=partial_C Gamma^A_DB-partial_D Gamma^A_CB+...`. The physical shell
is timelike, its outward unit normal n is spacelike, and
`h_AB=g_AB-n_A n_B`. The retained outward second-form convention is
`K_ab=h_a^A h_b^B nabla_A n_B`. Define the electric Weyl tensor by

```
E_ab=C(n,e_a,n,e_b)=C(e_a,n,e_b,n).
```

It is symmetric and trace-free with respect to the four-dimensional induced
Lorentzian metric h. Let `D_a phi` denote the intrinsic derivative and
`p=n(phi)`. All normal and tangent statements below use this same n.

## Direct contraction of Gauss

For a spacelike normal the Gauss equation in the retained curvature convention is

```
R4_abcd=R5_abcd+K_ac K_bd-K_ad K_bc,
Ric4_ab=Ric5_ab-R5_nanb+K K_ab-K_ac K^c_b,
R4=R5-2Ric5_nn+K²-K_cd K^cd.
```

The first R5 tensor on the right is projected into all tangent slots.
The five-dimensional Schouten tensor is
`P_AB=(Ric5_AB-R5 g_AB/8)/3`, and its decomposition gives

```
R5_nanb=E_ab+P_ab+h_ab P_nn
       =E_ab+Ric5_ab/3+h_ab Ric5_nn/3-h_ab R5/12.
```

Substitution and trace subtraction therefore yield

```
G4_ab=(2/3)Ric5_ab+[(2/3)Ric5_nn-(5/12)R5]h_ab-E_ab
       +K K_ab-K_ac K^c_b
       -(1/2)(K²-K_cd K^cd)h_ab.                       (A)
```

This also fixes the minus sign of E for its definition above; changing the
definition of an electric projection without changing (A) would reverse a term.

## Current field equations and shell laws

The retained bulk equation is
`Ric5_AB=phi_A phi_B+(2/3)V g_AB`. Set `Q=(D phi)^2=h^ab D_a phi D_b phi`.
Q is a Lorentzian contraction and is not assumed nonnegative. Then

```
Ric5_ab=D_a phi D_b phi+(2/3)V h_ab,
Ric5_nn=p²+(2/3)V,
R5=Q+p²+(10/3)V.
```

The bulk part of (A) is consequently
`(2/3)D_a phi D_b phi+[-5Q/12+p²/4-V/2]h_ab-E_ab`.
For the current pure-tension laws,

```
K_ab=(sigma/6)h_ab,   K=2sigma/3,
p=-sigma'/2,
K K_ab-K_ac K^c_b-(1/2)(K²-K_cd K^cd)h_ab
   =-(sigma²/12)h_ab.
```

The induced equation is therefore:

\[
\boxed{G^{(4)}_{ab}
 ={2\over3}D_a\phi D_b\phi
 -{5\over12}h_{ab}(D\phi)^2
 +\left({\sigma'^2\over16}-{V\over2}-{\sigma^2\over12}\right)h_{ab}
 -E_{ab}.}                                                   \tag{B}
\]

A useful trace consequence, with no electric-Weyl trace left, is

\[
R^{(4)}=(D\phi)^2-\sigma'^2/4+2V+\sigma^2/3.                 \tag{C}
\]

The corresponding known contracted Codazzi identity is
`D^a(K_ab-Kh_ab)=-D_b sigma/2=(n phi)D_b phi=Ric5_nb`, consistent
with the scalar boundary law and the signs used above.

## Background check

For the selected static-in-y background on the shell, `D_a phi0=0` and
`E0=0`. With `a=sigma/6`, `p=-sigma'/2` and the retained Hamilton identity

```
V=p²/2-6a²+6rho_b^-2,
```

the bracket in (B) becomes

```
p²/4-(p²/2-6a²+6rho_b^-2)/2-3a²=-3rho_b^-2.
```

Hence `G4_ab=-3rho_b^-2 h_ab` and `R4=12rho_b^-2`, precisely the closed
de Sitter shell `h=-d tau_b²+rho_b² cosh²(tau_b/rho_b) kappa`.
No four-dimensional matter equation of state was assumed in this check.

## Physical meaning and limits

E in (B) is an **electric projection at the shell**, with two outward normal
slots. The published H uses the **bulk geodesic observer and its angular
magnetic screen**, with a screen epsilon and one observer slot. Their slots,
locations and contractions differ. A negative selected H coefficient does
not determine E, the intrinsic Einstein tensor, a shell energy flux, or any
positive radiation density. The local shell laws alone do not close the
nonlocal electric-Weyl data; bulk evolution is still required.

Equation (B) displays induced geometry from the existing bulk scalar and
tension. Dividing it by a chosen four-dimensional gravitational constant
and calling the whole result an effective stress tensor does not identify
that tensor with actual brane matter. In particular it neither supplies a
thermal state nor adds particle degrees of freedom. A genuine matter
extension needs its specified action, state, modified junctions and matching
energy ledger; none is adopted here. Formula (C) is an available geometric
trace diagnostic for the current family, not a radiation or cosmological
origin result.

## Verification scope

The displayed calculation directly contracts Gauss, the five-dimensional
Schouten decomposition and the current scalar/tension shell laws. This is a
self-contained algebraic proof. The portable 21-identity checker additionally
checks the intrinsic four-dimensional Gaussian contraction and the downstream
trace/volume identities, but does not count this five-to-four-dimensional
application as another independently executed identity. No private checker or
private source file is needed to use the proof.
