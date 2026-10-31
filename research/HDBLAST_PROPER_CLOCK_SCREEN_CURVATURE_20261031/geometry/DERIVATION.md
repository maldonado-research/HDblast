# Proper-time Gaussian-label magnetic screen coefficient

This is a kinematic formula for the precise Gaussian construction in the
fixed reference prescription. The initial physical hypersurface, its radial/angular
labels, outward orientation, initial future normal, fixed round probability
measure, and fixed TT templates are anchor data. The observer is the unit
geodesic velocity at fixed proper time; the screen is tangent to fixed proper
time and initial radial label. Its oriented axes are the positive symmetric
inverse Gram square root of the angular connecting vectors. The result is
conditional on a common regular exact family and an invertible Gaussian map.
It supplies neither the background shooting certificate nor an evolved
numerical sign. The previous registered-clock sign is not transferred.

Use `g_e=g0+e q+e² k/2+o(e²)`, signature `(-++++)`, and
`R^a_bcd=partial_c Gamma^a_db-partial_d Gamma^a_cb+...`. Put
`R=rho cosh(t)`, `theta=tanh(t)`, `H=theta/rho`, `a=rho_y/rho`, and let
`c_chi=cosh(chi)`, `s_chi=sinh(chi)` at the background observation event.
Every derivative suffix below is an ordinary physical `t,y` derivative.
The first scalar field is the complete seven-amplitude field in the retained
screen derivation, with `Y=n0²-n1²`, `v=DY`, `h=Hess Y`.

The old time quadratic expression `F_t` is retained exactly. The new radial
expression and three first-curvature combinations are

```
eta = R(C_t-theta C)/(2rho)-R E_y/2-B/2,

F_y = ((-A+4M-L)M_y-2M L_y)/4 + C(A-L)/(2R)
      +(E/R-M_t/(2rho)) eta,

A_t = -(M_tt+theta M_t)/(2rho²)+E_t/(rho R)-(D+L)/(2R²),
A_y = -(M_yy+a M_y)/2+C_y/R-theta M_t/(2rho²)
      +theta E/(rho R)+(L-A)/(2R²),
K   = -(M_yt+theta M_y)/(2rho)+(C_t+theta C)/(2rho R)
      +E_y/(2R)-B/(2R²).
```

Let `xi^t Y,xi^y Y` be the first coordinate position variation of the
geodesic at fixed proper time and initial labels. The selected second
parameter derivative is `h_geo(Tplus-Tminus)`, where

```
h_geo = 3(c_chi U_Tx-s_chi U_Tr)/(2ell rho² cosh(t))
        +4/R {c_chi F_t+s_chi F_y
              -rho xi^t(c_chi A_t+s_chi K)
              -xi^y(s_chi A_y+c_chi K)}.                       (1)
```

Here `U_T` is the equal TT component of the exact metric second derivative,
under precisely the conditional exact-family binding retained from the exact-family source capsule.
Formula (1) uses neither a half-TT amplitude nor a first-iterate surrogate.
It requires both the time and radial TT jets. Linear scalar and vector
second sectors have zero projection by the retained SO(4) TT separation.
No first-order angular displacement coefficient occurs in (1).

## Complete baseline and radial quadratic term

At each parameter, boost the registered orthonormal time and outward spatial
radial normals through the fixed background rapidity. Keep the registered
screen, with its full positive inverse square root and physical epsilon.
The comparison contraction is linear in its normal slot, so its quadratic
term is exactly `4(c_chi F_t+s_chi F_y)/R`. The radial calculation, including
all four Codazzi, connection, epsilon, and timelike-screen terms, is in
`radial-review/RADIAL_DERIVATION.md`. It gives the full second tensor, not
just a projected term. It retains A until cancellation and does not obtain
the radial answer by simply rescaling the old time expression.

Both complete first symmetric magnetic tensors vanish for scalar first
fields: their radial Codazzi tensors contain only a trace and a Hessian,
whose symmetric curls vanish. Mixed scalar-gradient terms are isotropic
and cancel under symmetric epsilon contraction. The same is therefore true
for every radial boost. This is an identity at all points, not vanishing
only after the chosen projection. The full first Weyl tensor need not vanish.

## The first Weyl combinations

In an orthonormal frame write the Hessian coefficients of the first Weyl as

```
C_abcd = g(delta_ac h_bd+delta_bd h_ac
          -delta_ad h_bc-delta_bc h_ad) + isotropic terms,
C_a0b0 = e_t h_ab + isotropic terms,
C_ayby = e_y h_ab + isotropic terms,
C_a0by = k h_ab + isotropic terms.
```

Scalar symmetry gives `C_y0cd=0`; the isotropic parts above multiply
`Y delta` or its Kulkarni product. Define `m=M_t/(2rho)-E/R` and
`j=M_y/2-C/R`. A direct ADM/Gauss/Codazzi calculation, in the perturbed
orthonormal frame, yields the Riemann Hessian coefficients

```
g_R  = -L/(2R²)+H m-a j,
e_tR = -m_t/rho-2H m-D/(2R²)+a j,
e_yR = -j_y-2a j-A/(2R²),
k_R  = eta/R²-m_y-a m-H j.
```

For `g_R`, the intrinsic screen metric contributes `-L/(2R²)`: the
`M Hess Y` perturbation is a first angular diffeomorphism, while the L term
is a conformal perturbation. The time/radial Gauss products give `H m-a j`.
The electric terms follow by differentiating the respective second forms
and retaining lapse Hessians; the radial lapse Hessian has a minus sign.
For the mixed term use `R_a0by=D_b K_ay-D_y K_ab`; its screen derivative
and shape terms give `eta/R²-a m-H j` and its radial derivative gives
`-m_y`.

Schouten subtraction changes the angular Hessian coefficient by `-P_h`,
the time electric coefficient by `+P_h`, and the radial electric by `-P_h`.
Its mixed contribution is pure trace. Therefore

```
A_t=g+e_t=g_R+e_tR,
A_y=e_y-g=e_yR-g_R,
K=k=k_R.
```

Expanding these identities gives the three formulas above. In particular
these are invariant combinations of the linear Weyl slots, so it is safe
to use the displayed Riemann intermediate formulas only in those combinations.

## Relative orthonormal frame and all Lorentz slots

At first order the angular connecting vectors have off-screen coordinate
components `xi^t D_A Y partial_t+xi^y D_A Y partial_y`. After dividing by
R, resolve them in the boosted normal plane. The relative screen change is

```
delta e_a = p v_a U_chi+q v_a n_chi + in-plane rotation,
p=(rho c_chi xi^t-s_chi xi^y)/R,
q=(-rho s_chi xi^t+c_chi xi^y)/R.
```

First orthonormality then forces
`delta U=p v^a e_a` and `delta n=-q v^a e_a`, up to a radial boost.
The sign in `delta n` is opposite to the screen radial coefficient. The
positive Gram normalization removes the symmetric in-plane strain; its
remaining relative orthogonal rotation acts on the complete zero first H.
The first angular trajectory displacement affects that in-plane map only.

Differentiate the contraction in an oriented orthonormal frame. Its epsilon
components stay fixed under this relative Lorentz frame comparison. Every
first Weyl slot is varied:

```
delta H_ab = (1/2)epsilon_(a^{cd}[
 C(delta e_b,U,e_c,e_d)+C(e_b,delta U,e_c,e_d)
 +C(e_b,U,delta e_c,e_d)+C(e_b,U,e_c,delta e_d)].
```

The first term is zero: its U contribution has an antisymmetric repeated
normal pair and its n contribution is `C_nUcd=0`. The second term equals
`-p g Aax`, and the last two together equal `-p e_chi Aax-q k_chi Aax`,
where `Aax=epsilon_(a^{cd}v_c h_b)d`. Isotropic terms cancel under the
normalized symmetrization. Thus

```
delta H = -[p(g+e_chi)+q k_chi] Aax,
g+e_chi=c_chi² A_t+s_chi² A_y+2c_chi s_chi K,
k_chi=c_chi s_chi(A_t+A_y)+(c_chi²+s_chi²)K.
```

The retained all-point identity is `Aax=2(Tplus-Tminus)`, with fixed
probability norms `||Tplus||²=||Tminus||²=1/6` and zero cross product.
The first frame change times first Weyl contributes twice to the second
parameter derivative. Hence its scalar coefficient is
`-4[p(g+e_chi)+q k_chi]`. Substituting p,q and `c_chi²-s_chi²=1` reduces
this to the two position terms in (1). The projection normalization
`3 mean <H,Tplus-Tminus>` extracts exactly that scalar coefficient.

The epsilon/raising changes within the comparison baseline are retained
inside `F_t,F_y`; a relative Lorentz frame has constant epsilon components
only after those full metric changes have been included. Dropping the
baseline epsilon variation is a different and false operation.

## Moving events and the TT sign

The background Weyl vanishes. The first moving-event correction differentiates
the complete baseline first H, which vanishes identically, and therefore
vanishes. This justifies dropping that complete term, not isolated
`xi partial C1` terms. Second event or frame corrections multiply C0 and
vanish; radial boost and in-plane rotation corrections multiply complete H1
and vanish. The two displayed angular normal tilts are the surviving terms.

The linear TT Codazzi formula for a radial observer is
`(U_chi u) curl(T)/(2R)`. With `partial_y=-(ell rho)^-1 partial_r` and
eigenvalues `+3,-3`, the equal second TT amplitude gives exactly the first
term of (1), including its minus sign on `s_chi U_Tr`. A new bound for that
combined jet and the two first position variations is necessary before a
numerical proper-time sign can be stated.

The formula is natural when the metric, initial embedding, labels and
templates are all transported together. It still depends on those chosen
anchors, proper-time slice and screen prescription. It asserts no detector,
radiation, Big Bang, or observer-independent conclusion.
