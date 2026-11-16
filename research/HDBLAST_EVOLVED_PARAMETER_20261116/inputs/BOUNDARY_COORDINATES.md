# Coordinates, lapse and the retained boundary flow

This checkpoint quantifies the retained proper-`y` perturbation flow.
The reduced PDE remains in the conformal inward coordinate `z`. Sources
S1–S5 and their immutable public locations are listed in
[REFERENCES.md](REFERENCES.md).

## Coordinate and lapse conventions

Set

```
q_dS=-dt²+cosh(t)² kappa_S3,
g0=dy²+rho(y)² q_dS=rho(z)²(dz²+q_dS),
z=log(b/R), dy=-rho(z) dz, a=rho_y/rho,
rho_z=-a rho², shell z=0.
```

`R` is the ambient radial coordinate. The later scaled inward coordinate
`r=z/ell` is a different variable. The function `rho(z)` belongs to the fixed
background identification; `y` is not a parameter-dependent proper-distance
coordinate for the perturbed metrics. [S1 §1; S4 §§1,3]

| Quantity | Definition | Background value |
| --- | --- | --- |
| Proper-`y` radial lapse | `N_y=sqrt(g_yy)` | `1` |
| Conformal radial lapse | `Nrad=sqrt(g_zz)=rho(z) N_y` | `rho(z)` |
| ADM time lapse | `alpha` in `g=-alpha²dt²+gamma_ij(dx^i+beta^i dt)(dx^j+beta^j dt)` | `rho(z)` |
| Shell warp | `rho_b=rho(0)` | Fixed background constant |

The outward shell normal is `n_out=-Nrad^-1 partial_z` when the mixed
metric trace vanishes. We use `K=1/2 L_n h`,
`K_ab=sigma(phi)h_ab/6` and `n_out(phi)=-sigma'(phi)/2`. [S3 §4]

The normal entries of the fixed reference connection are

```
GammaBar^y_yy=0,             GammaBar^y_ab=-a rho² q_dS,ab,
GammaBar^z_zz=-a rho,        GammaBar^z_ab=+a rho q_dS,ab.
```

For comparison, `zhat=(y_b-y)/rho_b` is a scaled proper coordinate:
`GammaBar^zhat_zhat,zhat=0`. At the shell `zhat_z=1`, but
`zhat_zz=-rho_b a_b`. Its vanishing radial connection must not be assigned
to the conformal `z` coordinate.

## The literal retained perturbation flow

The twelve Neumann-block states are `q_N=(N_y,h_ab,phi)`; the four coordinate
Dirichlet states are `q_D=(g_ya)`. Subtract the known background throughout
the collar, `U=q-q0(z)`. At the shell define [S2 equations (1)–(4)]

```
f_N=(2/3)sigma(phi)N_y²-a_b N_y³ tr_h h0,b,
f_h=(1/3)sigma(phi)N_y h,
f_phi=-(1/2)sigma'(phi)N_y,
B_old(u)=-rho_b [f_b(q0,b+u)-f_b(q0,b)].
```

The generator freezes the shell background values and satisfies
`B_old(0)=0`. If `Phi_s` is its state-space flow, then

```
U_N=Phi_theta(z)(V_N), U_D=V_D
```

converts the boundary rows exactly to `partial_z V_N=0`, `V_D=0`.
It fixes the zero perturbation. Freezing this auxiliary generator leaves
the full prescribed bulk reference connection `Gamma[g0](z,t,Omega)` in
the relative wave-map gauge. All derivatives of the point transformation
remain in the transformed bulk source.

This checkpoint chooses an explicit cutoff. Let `E(s)=exp(-1/s)` for `s>0`
and `E(s)=0` otherwise, and set

```
chi(s)=E(1-s)/(E(s)+E(1-s)),
theta(z)=z chi((|z|-6ell)/(2ell)).
```

Then `theta=z` for `|z|<=6ell`, `theta=0` for `|z|>=8ell` and
`|theta|<=8ell`. The expression is smooth at zero because it equals `z`
on an open neighborhood there. This is the cutoff fixed for the present
audit; S2 prescribed only its local identity and a suitable smooth cutoff.
The claimed coordinate jets concern the identity region, including
`0<=z<=6ell`; no global cutoff-derivative budget is claimed.

## Distinction from a conformal raw-state generator

In shell-frozen conformal scaling `Nrad=rho_b nu`, `h=rho_b² M`,
`eta=phi-1`, `J=q_dS`, the physical shell rows instead give

```
A_b=-rho_b a_b, B_b=+rho_b a_b,
nu_z=-(2rho_b/3)sigma nu²+A_b nu+B_b nu³ tr(M^-1 J),
M_z=-(rho_b/3)nu sigma M,
eta_z=(rho_b/2)nu sigma'.
```

These are the same physical boundary laws in different state variables.
Their autonomous raw-state flow is a different auxiliary extension. In
particular the conformal lapse perturbation generator is
`B_new,nu=B_old,N+A_b(N_y-1)` at the shell. Its raw flow does not fix the
actual collar background without recentering. It is included here to
identify the coordinate correction; this checkpoint does not duplicate
the flow bounds for that alternative chart. Normalizing by the moving
`rho(z)` would introduce additional radial rescaling terms.

Neither auxiliary flow time `s`, cutoff value `theta(z)`, nor coordinate
time `t` is an elapsed physical clock. At background a fixed-`z` observer
has `d tau=rho(z)dt`; a nonlinear observer clock and its domain must be
constructed separately. The binding's label `Tphys` denotes coordinate
time. No proper-time cylinder or evolved finite-amplitude bound follows
from this coordinate dictionary. [S4 §§1,3]
