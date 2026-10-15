# Coupled physical-shell scalar activation for an explicitly registered aligned bump

15 October 2026 UTC. Independent GPT-6.1-Sol derivation. This is a conditional
analytic continuation of the frozen finite-regularity theorem, written only in
this working directory. No frozen source is changed. The new nonvanishing
argument uses smooth coefficients, a one-dimensional Duhamel formula, and an
explicit bump with both displacement and positive velocity nonzero. It does
not use a test scalar, analyticity, numerical evolution, or a claim that every
nonzero disturbance must give a scalar shell signal.

## 1. Inputs and exact coupled reduction

Fix one of the registered backgrounds and an allowed finite interval
`0 <= t <= T`, `T>0`, satisfying all hypotheses of the frozen exact nonlinear
shell evolution theorem. Its fields obey

```
g0 = dy² + rho(y)²[-dt² + cosh(t)² kappa_S3],
a = rho_y/rho,       p = phi0_y != 0,
a_y = -rho^-2-p²/3,  p_y = U'-4ap,
a_b = sigma_b/6,     p_b = -sigma'_b/2,
B_b = p_y,b/p_b + sigma''_b/2 > 0.
```

These are the original positive-detuning model and shell tension. All bounds
below concern this fixed background. In particular, no lower bound for `|p|`
uniform in the detuning is asserted. Put `dr/r=dy/rho`, `r_b=b`, and
`z=log(b/r)>=0`; then `partial_y=-rho^-1 partial_z`.

For `Y=n1²-n2²`, `Delta_S3 Y=-8Y`, write a longitudinal perturbation as

```
h_yy = -4uY,   h_ya = 0,   h_ab = 2uY (h0)_ab,
u = P/rho²,    delta phi = chi Y,
chi = -3P_y/(rho²p) = 3P_z/(rho³p).
```

The exact master expression is

```
M[P] = P_yy -(2a+2p_y/p)P_y
       +rho^-2(Box_gamma+4)P -(2/3)p²P.
```

This is a coupled Einstein–scalar reduction. The complete residual identities
are `delta E_yy=2M[P]/rho²`, `delta E_ya=0`,
`delta E_ab=-gamma_ab M[P]`, and
`delta S=-3(partial_y+2a)M[P]/(rho²p)`, with the angular factor understood.
Thus solving `M[P]=0` reconstructs both fields and all linearized bulk
equations. These identities and their regular ambient counterpart are retained
inputs; their derivation is in the 13 October regular-cone capsule.

Use identical copies and zero longitudinal first shell displacement. The full
linear Israel/scalar laws reduce to

```
u_y+2au+p chi/3=0,
chi_y+2pu+sigma'' chi/2=0.
```

The first is an identity of reconstruction. The second, together with the
bulk equation, is exactly

```
P_y = (Box_gamma+4)P/(rho_b² B_b).
```

For `l=2`, the trace-free Israel graph equation forces the even graph
displacement to zero in this longitudinal slice. All metric components,
induced scalar matching, and normal-force consistency are retained; there
is no separately chosen scalar boundary source. The complete graph proof is
in `REGULAR_SHELL_CONTACT_20261013/manuscript/SHELL_CONTACT.md`, section 5.

In inward `z`, define

```
C(z) = rho(3a+2p_y/p),       D(z) = (2/3)rho²p²,
kappa = rho_b B_b > 0.
```

The exact bulk and even shell equations are

```
P_tt-P_zz+3tanh(t)P_t-C(z)P_z+[8sech²(t)-4+D(z)]P=0,
P_tt+3tanh(t)P_t+[8sech²(t)-4]P-kappa P_z=0   at z=0.
```

Let

```
alpha(z)=exp[-(1/2)integral_0^z C(v)dv] > 0,
        =(rho(z)/rho_b)^(3/2) p(z)/p_b,
P(t,z)=alpha(z)cosh(t)^(-3/2)w(t,z).
```

Then

```
w_tt-w_zz+W(t,z)w=0,
W=(35/4)sech²(t)-25/4+C_z/2+C²/4+D,
```

and, writing `beta(t)=w(t,0)`,

```
beta''-kappa w_z(t,0)+H(t)beta=0,
H(t)=(35/4)sech²(t)-25/4+kappa C(0)/2.
```

The term `kappa C(0)/2` comes from the radial conjugation. The signs of the
normal derivatives are fixed by `z` increasing into the interior.

## 2. Explicit profile with A=B and positive velocity preserved

Choose a registered `ell>0` with `ell<T` and `2ell+c_*T<Z`, leaving
strict room for the enclosing exact-correction annulus and any fixed
transformation collar. Set `d=ell`. For example, whenever the corresponding
collar is available, `ell=T/4` is one definite choice. More generally one
may use `ell=min{T/4,(Z-c_*T)/4}` before placing the enclosing annulus.
Define the flat smooth bump

```
B(z) = exp[4-ell²/((z-d)(d+ell-z))],   d<z<d+ell,
B(z) = 0 otherwise.
```

It is positive in its support interior, has support infimum exactly `d`,
has peak one, and is increasing on `(d,d+ell/2)`. Register this bump in the *conjugated* master
variable, rather than silently treating an ambient bump as the same function.
Choose the fixed positive velocity scale `nu=1/b`, set `r=b exp(-z)`, and

```
G(r²)=r^-6 alpha(z) B(z),       K(r²)=nu G(r²)=G(r²)/b.
F|T_ambient=0 = G(|X|²)(X1²-X2²),
F_Tambient|T_ambient=0 = K(|X|²)(X1²-X2²).
```

The argument also works for any fixed `nu>0`, but the displayed registered
profile uses exactly `nu=1/b`. The symbol `T_ambient` is ambient Minkowski coordinate time; it is
different from the slab endpoint `T` and from global shell time `t`.
All powers are harmless because the support lies in a closed annulus away
from `r=0` and `r=b`. The data are smooth and zero on open neighborhoods of
both shell and vertex.

At `t=0`, `T_ambient=r sinh(t)` implies
`partial_t=r partial_Tambient` on the initial slice. Since `P=r⁴F`, the exact
mode data are

```
P(0,z)=alpha(z)B(z),       P_t(0,z)=nu r alpha(z)B(z),
w(0,z)=B(z),              w_t(0,z)=nu r B(z).
```

Both initial sectors use the same nonzero aligned matrix
`A=B_matrix=diag(1,-1,0,0)`. This construction preserves nonzero displacement,
nonzero velocity, and the positive velocity scale. It is an explicitly
specified allowed radial profile, and is not an asserted identification with
an earlier unspecified radial bump.

The regular ambient reconstruction is

```
h_AB=2F(s eta_AB-3X_A X_B),
delta phi=-3(E+4)F/(2Omega² Phi_s),   s=|X|²-T_ambient².
```

Its full residual factorization shows that the ADM data solve the linearized
constraints. They are invariant under every independent coordinate reflection:
the master polynomial and radial factors are invariant, and the reconstructed
tensor transforms equivariantly. Reflection invariance applies to metric,
scalar, metric momentum, and scalar momentum; it does not mean that every
coordinate component is an even scalar function. Consequently this seed lies
in the invariant constraint kernel used by the frozen exact localized-data
theorem. The weighted inverse gives exact canonical data tangent to it, with
the same finite Sobolev regularity and fixed annular support as in that theorem.
No scalar initial-data correction is needed by that construction. Choose the
correction annulus to enclose `[d,d+ell]` with positive gaps on both sides.
Increase the auxiliary inner collar if needed to retain its strict causal
margin. This stays inside the already stated conditional hypotheses.

## 3. An explicit causal transfer, including shell feedback

Put `j=w_t+w_z`, `k=w_t-w_z`. Then

```
j_t-j_z=-Ww,             k_t+k_z=-Ww,
j0=B'+nu r B,           k0=-B'+nu r B,
j(t,z)=j0(z+t)-integral_0^t W(s,z+t-s)w(s,z+t-s)ds.
```

For `z>=t`,

```
k(t,z)=k0(z-t)-integral_0^t W(s,z-t+s)w(s,z-t+s)ds.
```

For `0<=z<t`, replace the first term by
`2beta'(t-z)-j(t-z,0)` and the lower limit by `t-z`. Complete the system with

```
w(t,z)=B(z)+(1/2)integral_0^t[j(s,z)+k(s,z)]ds,
beta''+kappa beta'+H beta=kappa j(t,0),
beta(0)=beta'(0)=0.
```

These are causal Volterra equations, not a freely prescribed incident
waveform. They also give the required smooth master solution directly on the
positive-cone collar: the integrated boundary ODE and the characteristic
equations are a contraction on `(w,beta,beta')` on a short time step, with
Lipschitz constant at most `constant*step_length`. The constant depends on
bounded coefficient derivatives, the fixed harmonic, and `kappa`, rather
than solution amplitude. Repeated steps cover the finite interval. Repeated
differentiation preserves smoothness because all initial corner jets vanish.
The characteristic constraints `w_z=(j-k)/2` and `w(t,0)=beta(t)` propagate
from their compatible initial data. Finite propagation gives support in
`[max(0,d-t),d+ell+t]`; the strict collar margin removes any effect of an
artificial continuation deeper in the interior. This is precisely the
retained finite-angular construction specialized to this collar, and requires
no negative-cone continuation or conformal-rectangle hypothesis.

The scalar response on the physical shell is exactly

```
Phi1(t,theta)=Aphi(t)Y(theta),
Aphi(t)=3cosh(t)^(-3/2)/(rho_b³p_b)
        [j(t,0)-beta'(t)-C(0)beta(t)/2].
```

Together with `u=P/rho²`, this is a full coupled linear transfer prescription.
It includes the boundary ODE and the wave it sends back into the bulk. A
numerical waveform or closed elementary expression for this Volterra solution
is not claimed.

## 4. Smooth near-front nonvanishing with positive velocity

**Lemma.** Let `W` be bounded and smooth on a collar containing
`0<=t<=T`, `0<=z<=T`. For the above data, the coupled shell scalar
`Phi1` cannot vanish on every interval `(0,d+eta)`, for any `eta>0`
with `d+eta<=T`. Thus its activation infimum is exactly `d` in the
global clock `t`; no nonzero value or jet at `t=d` is asserted.

**Proof.** Finite propagation and the compatible zero initial shell collar
give `Phi1=0` for `0<=t<d`. Suppose it were zero on `[0,t1]`, with
`d<t1<=T`. Because `rho_b p_b !=0`, reconstruction gives
`P_z(t,0)=0` there. The physical even shell equation then gives

```
P_tt+3tanh(t)P_t+[8sech²(t)-4]P=0  at z=0,
P(0,0)=P_t(0,0)=0.
```

Ordinary differential uniqueness forces `P(t,0)=0`. Therefore
`w(t,0)=w_z(t,0)=0`. Extend `w` by zero to `z<0`, extending `W`
smoothly there in any bounded way. Zero Dirichlet and normal trace remove
all interface delta terms, so the extension satisfies the full-line wave
equation in distributions and hence its ordinary Cauchy Duhamel formula.
Classical regularity also follows by recursively using the PDE at the zero
traces. This zero extension is used *only under the contrary hypothesis*.

Let `M=max |W|` on `[0,T]x[0,T]` and `c0=nu b=1`, so
`0<=nu r<=c0`. Choose

```
0<epsilon<=min{ell/4, (t1-d)/2, 1,
               1/[2(1+M)T(3+c0)]},
t=d+epsilon,       Q=B(t)>0,
a0=M t epsilon/2.
```

Consider the backward cone `D_t={(s,y):0<=s<=t, |y|<=t-s}`.
Finite propagation from the support infimum `d` gives `w(s,y)=0`
for `y<d-s`. Every source slice of `D_t` can therefore have nonzero
width at most `epsilon`: its lower support endpoint is `d-s` and
its cone upper endpoint is `t-s`. Every backward initial interval of a
point in this cone meets the initial support only in `[d,t]`.
The monotone outer half of the registered bump gives `B(y)<=Q` there.

The free part of the full-line d'Alembert formula consequently has absolute
value at most `Q(1+c0 epsilon/2)` everywhere in this cone. If
`S=sup_Dt |w|`, Duhamel gives

```
S <= Q(1+c0 epsilon/2) + a0 S,
S <= Q(1+c0 epsilon/2)/(1-a0).
```

At `(t,0)` its free displacement term is `Q/2`; its free velocity
integral is nonnegative. The potential contribution has absolute value at
most `a0 S`. Hence

```
w(t,0) >= Q[1-a0(3+c0 epsilon)]/[2(1-a0)] > Q/4 > 0.
```

The registered epsilon bound gives
`M t epsilon(3+c0 epsilon)<1/2`, which is stronger than the needed
`M t epsilon(3+c0 epsilon)<2`. This contradicts the zero boundary trace.
The same choice can be made for every `t1>d`, so there is nonzero
scalar response somewhere in every interval `(d,d+eta)`. Smoothness makes
the set of nonzero times open. The proof uses only bounded potential,
finite propagation, positivity of the prescribed initial velocity, and the
outer monotonicity of the specified bump. It does not use analytic unique
continuation. The threshold is `d`, not `2d`. QED.

The positive lower bound above concerns the *hypothetical zero-extension
full-line solution*. It is not a positive pointwise lower bound for the
actual Wentzell solution or its scalar trace. Its role is to contradict
identically zero scalar response. In particular the argument does not
provide a numerical nonzero time or scalar amplitude.

## 5. Identification with the exact nonlinear family's first variation

The master solution is a physical longitudinal representative. To compare
with the exact family's fixed relative-harmonic reduction, use the retained
linear gauge construction from
`FULL_SECOND_ORDER_CONTACT_20261014/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`,
equation (5) in its gauge-transfer section:

```
(Box0+Ric0)xi=-C[h_long],
xi_y|b=0,       (partial_y-2a_b)xi_a|b=0,
xi|t=0=partial_t xi|t=0=0.
```

The tangential components here are covariant components. The exact coefficient
is `2a_b`, since `(L_xi g0)_ya=partial_y xi_a-2a_b xi_a` at the shell.
The finite-angular mixed wave construction solves this equation. The transformed
pair `h_wave=h_long+L_xi g0`,
`chi_wave=chi_long+xi^y p` has zero first harmonic vector and zero mixed
normal/tangential metric rows. The physical shell laws are preserved by
covariance. Its generator is tangent to the shell, and therefore
`chi_wave|b=chi_long|b`; the background shell scalar is constant, so any
remaining tangential reparameterization changes neither scalar trace.

Zero initial `xi,xi_t` preserve the canonical tangent and initial metric
values. The changed lapse/shift velocities are exactly those required by the
frozen initial harmonic completion, which solves the linear harmonic equations
uniquely for these five velocities. The linearization of the frozen reduced
evolution theorem has uniqueness for these completed initial data and boundary
rows. Its actual parameter derivative therefore equals the harmonized master
solution. This supplies the bridge to an exact constraint-compatible family;
scalar-master uniqueness alone would not supply it.

In an arbitrary gauge the shell trace is
`Phi1=delta phi+Z dot dphi0`; the second term is necessary. In this
longitudinal representative its even graph is zero. The gauge transformation
above is shell tangent, so its calculated nonvanishing is precisely that of
the physical pullback, rather than a fixed-coordinate scalar surrogate.

## 6. Consequence for the physical quadratic observable

For any registered C2 intrinsic-clock section family reducing to a background
`t=tstar`, where the lemma supplies `Aphi(tstar)!=0`, put

```
Q_epsilon=integral_section (Phi_epsilon-phi_b)² dmu_epsilon.
```

The background shell scalar is constant, so the section displacement, induced
measure variation, and angular relabeling contribute no first-order correction
to `Phi1`. Since `integral_S3 Y² dOmega=pi²/3`, the exact first-variation
coefficient is

```
Cstar=integral Phi1² dmu0
     =(pi²/3)rho_b³cosh(tstar)³ Aphi(tstar)²
     =3pi²cosh(tstar)³ P_z(tstar,0)²/(rho_b³p_b²)>0,
Q''_0=2Cstar.
```

The frozen C2 parameter theorem and a C2 section pullback imply
`Q_epsilon=Cstar epsilon²+o(epsilon²)`. Thus some parameter radius
`epsilon0>0` gives `Q_epsilon >= (Cstar/2)epsilon²` for
`0<|epsilon|<epsilon0`. This is a conditional rigorous nonlinear
nonvanishing statement at a time whose existence is established analytically.
The constant and radius are expressed through the actual linear trace and
continuity remainder, without an asserted numerical enclosure.

This linear activation infimum is not the first-arrival time of a nonzero
finite-amplitude family member. Its nonlinear constraint correction can fill
the surrounding allowed annulus at order epsilon², and its metric changes
causal cones. Neither effect changes the first derivative computed here.

The root continuation separately develops an intrinsic proper-time clock and
a fixed-window integrated observable. The present proof is enough to make its
linear squared-trace integral positive on every window ending after `d`.

## 7. What is reached and what is still missing

For the explicit registered bump above, both aligned matrices and the positive
velocity scale are retained. Conditional on the already stated background,
linear contact construction, exact initial-data, gauge-completion and nonlinear
finite-slab theorem inputs, this note proves genuine scalar activation on the
physical shell and a positive quadratic response for the exact nonlinear family.
It supplies an explicit coupled Volterra transfer, and a rigorous proof of
nonvanishing; it supplies no computed waveform, numerical first arrival for a
finite-amplitude member, explicit background-dependent numerical error bound,
global stability, cosmological origin, or external novelty conclusion.

Pure displacement with zero velocity is also a legitimate reflection-invariant
constraint tangent, but its matrix pair is `(A,0)`, not the retained `(A,A)`.
The older time-reversal/sideways-uniqueness proof applies to that alternative
and to pure velocity. It does not directly prove the mixed positive-velocity
case. The near-front Duhamel argument above is the additional mixed-profile
step, and is confined to the explicitly registered outer-monotone profile.

Finite symbolic checks accompanying this note verify the transformation,
boundary signs, angular normalization and bound algebra. They do not
computer-certify the PDE arguments or the frozen theorem.

`check_response_algebra.py` passed 17 exact identities in each of ordinary
and optimized Python; its three injected corruptions failed in both modes.
The positive outputs agree byte for byte. `CHECK_RECEIPT.json` records the
executed outcomes, interpreter and SymPy versions, and checker hash;
`SOURCE_PROVENANCE.json` records hashes of all inherited inputs. To rerun
the positive identities, use `python -I -B check_response_algebra.py`, and
repeat with `-O`. Its mutation arguments are `--inject-error potential`,
`--inject-error shell-sign`, and `--inject-error ambient-velocity`.
