# A pure second harmonic gauge with smaller initial metric components

1 November 2026. This addendum changes the representative of the explicit
prepared second jet; it does not change the physical second jet. The original
five-mean map remains in `MEAN_INITIAL_GAUGE.md` and `initial_means`. The new
map is `conditioned_initial_means`. All source/parameter/frame conventions
are those of that note.

## 1. Generator and the resulting initial pair

Choose a purely second parameter generator Z, with no first generator, and
initial values

```
Z^y=-xi,  Z^t=0,  Z^i=0,
partial_t Z^y=0,  partial_t Z^t=a xi,  partial_t Z^i=0.        (C1)
```

Its second-jet pullback law is
`kbar=k+Lie_Z g0`, `wbar=w+Z(phi0)`; there are no first-generator
quadratic transformation terms. Since `Z^t|t=0` vanishes as a spatial
function, this transforms the initial section only by its radial
diffeomorphism. It gives

```
delta gamma2=-Lie_(xi partial_y)gamma0,
delta phi2=-p xi,  delta pi2=delta varpi2=0,
delta g_tt2=2a rho²xi-2rho²partial_t Z^t=0,
delta g_tI2=0.                                                (C2)
```

The vanishing momentum changes follow from natural tensor-density pullback
and `pi0=varpi0=0`. The inverse canonical map is natural under spatial
diffeomorphisms, so it commutes with this pullback exactly. Equivalently,
`delta K2=Lie_(-xi partial_y)K0=0` and `delta Pi2=0`, as the initial
section is unchanged as a set and `K0=Pi0=0`. Its future normal is geometric,
so the displayed initial time-rescaling does not alter this conclusion.

The resulting initial vector and velocity are

```
Ubar0=(f,0,-u²/3,f,-p xi),
Vbar0=ell*(
  2rho QM/(3a)+(32/3)u ut,
  rho[-f_y+a f-(13/3)u u_y-a u²/3],
  -2rho QM/(3a)+(4/3)u ut,
  -(4/3)u ut,
  0),
partial_r Ubar0=-ell rho*(f_y,0,-2u u_y/3,f_y,-p_y xi-p xi_y).
                                                               (C3)
```

Here the scalar second displacement is nonzero. Its coordinate time
velocity remains zero because `partial_t Z^y=0` and `phi0` is static.
The Hamiltonian remains satisfied: naturality gives
`DChat[Lie_Z gamma0,Lie_Z phi0,Lie_Z pi0,Lie_Z varpi0]
=Lie_Z Chat(x0)=0`. Removing only the metric Lie derivative while keeping
`phi2=0` would fail this argument and would select different physical data.

## 2. The monopole vector wave, including the shell rows

To preserve the second relative-harmonic condition throughout evolution,
solve `(Box0+Ric0)Z=0`. Derive its monopole system from the covariant warped
operator [S8, equations (2)–(3)], writing `Z=X(t,y)partial_y+b(t,y)partial_t`
and `H=tanh t`. A constant spherical harmonic has no angular-gradient vector;
there is no third angular-amplitude equation to impose. The two equations are

```
X_tt=rho²(X_yy+4aX_y)-3H X_t
     -2a rho²(b_t+3Hb)-rho²(4a_y+8a²)X,
b_tt=rho²(b_yy+6ab_y)-3H b_t+(3+3H²)b-2aX_t.                (C4)
```

The homogeneous shell-preserving coordinate rows are exactly

```
X|shell=0,  b_y|shell=0.                                    (C5)
```

The first fixes the shell as a set. For the cross metric,
`(Lie_Zg0)_yt=X_t-rho² b_y`, so (C5) preserves its Dirichlet row.
The time component b may be nonzero on the later shell; it reparametrizes
shell time and must be transported in anchored observables. The geometric
physical shell residuals have homogeneous variation under this tangential
diffeomorphism because their background values vanish. The gauge wave gives
`dC[Lie_Zg0]=0`; its boundary trace gives the normal harmonic row. Thus all
homogeneous mean boundary rows hold, with no extra scalar or vector condition.

For a direct constructive existence argument, convert to `x=t/ell,r=z/ell`.
In `Zvec=(X,b)` order,

```
Zvec_xx-Zvec_rr=Az Zvec_x+Bz Zvec_r+Dz Zvec,
Az=[[-3ell H,-2ell a rho²],[-2ell a,-3ell H]],
Bz=diag(-3ell a rho,-5ell a rho),
Dz=[[-ell²rho²(4a_y+8a²),-6ell²a rho²H],
    [0,3ell²(1+H²)]].                                      (C6)
```

Use odd reflection for X, even reflection for b, with the reflected free
wave for (C1). The retained background bounds `rho<=131`, `|a|<=14/125`,
`|a_y|<1/10000`, `|H|<=3/10^6`, and `ell=10^-6` imply, in the row-sum
max channel norm of [S9],

```
q3=3 max_i sum_j (|Az_ij|+|Bz_ij|+(3/2)|Dz_ij|)<3/250<1.
```

This fresh two-channel bound is checked by exact rational arithmetic; no
quadrupole comparator, quadrupole source, or angular-amplitude norm is used.
The same reflected Duhamel kernel therefore gives a unique C1 solution on
the retained triangle, with nonzero free initial data. All coefficients and
the compact flat initial data are smooth. The usual differentiated
characteristic/trace induction, including its finite initial-jet terms,
gives the needed smoothness; the flat corner eliminates compatibility
obstructions. This proof uses only bounded smooth coefficients and the
linear contraction, so large initial xi does not obstruct existence. It
does not provide an effective nonlinear parameter radius.

At the initial slice, (C4) gives

```
X_tt=rho²[-xi_yy-4a xi_y+(4a_y+6a²)xi],  b_tt=0.
```

Consequently

```
delta(dot k_ty)=X_tt-rho²(a xi)_y
 =rho²[-xi_yy-5a xi_y+(3a_y+6a²)xi].                         (C7)
```

Substitution of `A=f+2xi_y`, `B=f+2a xi` into the original mixed
velocity proves that (C7) is exactly the difference between the two
entries in (C3). The other metric velocities are unchanged; the time-time
change vanishes since `X_t=b_tt=0`. This checks the harmonic recompletion
against the actual homogeneous vector wave, not only the spatial pullback.

## 3. Exact-family interpretation and scope

The original prepared exact C2 family already exists under the retained
analytic premises. Let Phi_e be the flow of Z for flow parameter e²/2.
The homogeneous vector wave has a smooth solution by section 2; on the
compact common domain its flow exists for sufficiently small e. The vector
is tangent to the shell, is tangent to the initial section, and vanishes
on the initial shell collar. Hence this actual pullback preserves the shell
and its initial anchor as geometric sets. Its second derivative is

```
kbar=k+Lie_Z g0,  wbar=w+Z(phi0),
```

and its first derivative is unchanged. Naturality preserves the exact
physical Einstein/scalar/shell equations and anchored observables under
this pullback. The finite-parameter pulled-back metric is not claimed to
be in nonlinear relative-harmonic gauge. That assertion is unnecessary.
At second order, `dC[Lie_Z g0]=(Box0+Ric0)Z=0`, so the transformed second
jet solves the same forced reduced linear equations, with the unchanged
actual-first-field sources and the homogeneous boundary rows (C5).
Its full initial pair is (C3). Retained linear uniqueness identifies its
five means with the response from these conditioned initial data.

This is a statement of second-jet gauge invariance and identification.
The original family already has the required C2 parameter dependence and
lower C4 spacetime control. The second derivative at zero of the pullback
uses only `k`, `g0`, and one Lie derivative by the smooth Z, because the
first parameter generator is zero. No general lossless C2 Sobolev
composition theorem is asserted here. The retained anchored jet gauge law
applies at these finite derivative orders; an additional nonlinear harmonic
re-completion or a claim about the full conditioned finite-parameter
Sobolev solution map is not a premise of the signed coefficient.

All original first fields and quadratic bulk/shell sources remain unchanged.
The scalar component `-p xi` must now be propagated with the metric means.
The potentially large `xi_y` term has been removed from the initial metric,
but no numerical error improvement is certified by that algebra alone.
Anchored proper-time observables are unchanged through second order when
their shell clock/embedding is transported. Raw coordinate mean histories
can change. This is not a matter-production, shell-force, or sign claim.

Sources S8 and S9 are the pinned general warped-vector derivation and
mixed-kernel construction inside the retained first-gauge capsule. Their
general operators and elementary existence argument are used; their
quadrupole-specific data and numerical solution are not imported.
