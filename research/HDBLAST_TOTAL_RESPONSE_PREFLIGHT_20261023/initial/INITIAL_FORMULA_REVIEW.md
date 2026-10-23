# Independent registered initial five-mean derivation

23 October 2026. Text/formula audit only. No numerical capsule module was
imported for this review. This note does not replace the complete closure
authentication required before a numerical evaluation.

Let `ell=1/1000000`, `z=ell*r`, `ir=1/rho`, `h=U1/p`,
`E=exp(-z)`, `d=3/(rho*p)`. Use the exact registered initial first functions

```
u=C*q, v=B*(q_z+k*q), ut=E*u, vt=E*(v-d*u),
dy=-ir*dz, uy=-p*v/3-2*a*u,
py=p*(h-4*a), ay=-ir^2-p^2/3, d_y=d*(3*a-h).
```

The first identity for `uy` follows from `C=alpha/rho^2`,
`B=3*alpha/(rho^3*p)`, `k=alpha_z/alpha`, and
`rho_z=-rho^2*a`; it is valid at every point, including flat extensions.

The physical Cauchy source definitions from the retained compact initial
prescription are

```
Fmom=ir*(18*u*ut-v*vt),
Gmom=ir*(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*vy),
K=Fmom/3, Pi=(Gmom-Fmom_y)/p.
```

Here `K` is the scalar profile in `K2=(Fmom/3)*Y^2*gamma0`;
it is not the response amplitude `K0`. Substitution gives exactly the
profiles used by `N/numeric/evaluate_initial_profiles.py`:

```
K=E*ir/3*(18*u^2-v^2+d*u*v),
Pi=E*ir/p*(d*u*vy+(d*(h-ir)+10*p)*u*v
                  +(30*a-6*ir)*u^2+(2*ir-a)*v^2).
```

For the second identity, differentiate `Fmom` before substituting:

```
Gmom-Fmom_y
 =ir*(-6*u*uyt-24*uy*ut-30*a*u*ut-vt*vy+v*vyt-a*v*vt).
```

Use `E_y=E*ir`, `uyt=E*(uy+ir*u)`,
`vyt=E*(vy-d*uy-d_y*u+ir*(v-d*u))`, followed by the
displayed identities for `uy` and `d_y`. This is the all-point algebraic
proof; finite rational spot checks are supplementary only.

The evaluator's cancellation-preserving Hamiltonian source is also exactly
the original physical source:

```
A=-2*vy^2-8*p*u*vy+8*p*(h+2*a)*u*v
  +(16*p^2-192*ir^2)*u^2-(2*U2+8*p^2/3)*v^2
  -ir^2*(12*ut^2+2*vt^2),
D=ir^2*(36*u^2-2*v^2).
```

To derive it substitute
`uyy=-py*v/3-p*vy/3-2*ay*u-2*a*uy` into the original
`-24*uy^2-288*a*u*uy-24*u*uyy` expression, then use
`uy`, `py`, and `ay` above. The `a^2*u^2` terms cancel.
Set `S0=A/6+4*D/3`, `f0=rho^2*S0/6`, and
`zeta0=-3*f0_y/(2*p)`.

In mean order `(K0,L0,N0,J0,w0)`, since the normalized sphere mean of
`Y^2` is `1/6`, the registered arrays are precisely

```
U0=(0, 0, -u^2/3, f0, zeta0),

V0=ell*(
 (2*rho*K-8*u*ut)/6,
 rho*((-26*u*uy-2*a*u^2)/6-3*f0_y/2-3*a*f0),
 (-8*rho*K+44*u*ut)/6,
 (2*rho*K+4*u*ut)/6,
 (rho*Pi+2*u*vt)/6
).
```

This agrees directly with `source_adapters.initial_frame_data` followed by
`linear_operator.pack_initial`. In particular the ordinary `t` velocities
are all multiplied by `ell`, while `U0` is not. There is no additional
factor `rho` in the scaled velocity conversion.

The mixed `L0` velocity contains the `-2*a*u^2` term from the
frame-component conversion, as well as the `f0_y` and `a*f0` terms.
The nonzero lapse value `N0=-u^2/3` and scalar value/velocity are required.
Using covariant normal derivatives in place of these ordinary component
time velocities would change the problem. The mean `L0` boundary row is
Dirichlet, but its interior initial velocity need not vanish.

For a whole-cone C1 free-wave certificate, the actual radial derivative is
`d_r=-ell*rho*d_y`. Thus

```
U0_r=(0,0,-2*u*u_r/3,f0_r,zeta0_r).
```

`f0` uses initial packet derivatives through order two, `f0_y` through
three, and `zeta0_r` through four. Those derivatives are available from
the explicit packet. C1 free-wave error uses only value error in `V0`,
not a derivative of `V0`; imposing a derivative budget on `V0` would be
an unnecessary stronger requirement for the stated route.

All registered functions are flat at `r=1,2` and vanish outside `(1,2)`;
the relevant initial footprint is `[0,7/4]`. Their parity reflections are
smooth at `r=0`, and the resulting free wave vanishes when `x+r<=1`.
Any polynomial representation must preserve those properties explicitly.

The retained shell/background intervals express uncertainty about one
correlated actual background. Their widths do not prove a lower bound on
the actual error of a chosen approximant. In particular an upper error
certificate larger than `1e28` is a failed accuracy certificate for that
representation, not proof that every possible initial representation is
impossible. A claimed intrinsic obstruction needs either an information
lower bound over the retained admissible backgrounds or an exact
pointwise lower bound for the concrete representation. Spatial refinement
alone cannot shrink fixed input uncertainty.

Reviewed text sources: retained compact physical initial jet sections
1–5; initial harmonic completion equations (11)–(12); registered
`source_adapters.initial_frame_data`; numerical initial evaluator
`calculate` and initial profile proof. No full subsidiary audit or total
evolution certificate is supplied by this formula review.
