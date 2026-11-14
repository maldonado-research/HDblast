# Conditional metric patch, sixteen-field principal energy and radial speed

6 November 2026. Conditional quantitative continuation of the retained 15 October
argument. All numerical conclusions in this note are conditional on the
physical metric belonging to the specified pointwise patch at every point
and on the prescribed boundary identities holding for every frozen iterate.
This note does not prove that the actual finite-amplitude nonlinear family
stays in that patch.

## 1. Actual coordinates, frame and state variables

The retained background is

```
g0 = dy² + ρ(y)²[-dt² + cosh²(t) κ_S3],
z = log(b/r),       dy = −ρ dz,       ell = 10^−6.
```

Here `t` is unscaled de Sitter coordinate time and `z` increases inward.
At the shell `ρ=ρ_b`; throughout the retained wider background collar,
`129≤ρ(z)≤131`, hence also `129≤ρ_b≤131`. The candidate slab is
`|t|≤13ell/8`, `0≤z≤6ell`. These are coordinate quantities. No proper-time
clock or observer trajectory is constructed here.

Let `(θ¹,θ²,θ³)` be a local κ-orthonormal coframe, and use the base coframe

```
ω = (dz, dt, cosh(t) θ¹, cosh(t) θ², cosh(t) θ³),
J = diag(1,−1,1,1,1),
g = ρ(z)² Γ_IJ ω^I ω^J.
```

The Euclidean matrix operator norm below is the one in this displayed
coframe. Angular frame changes are orthogonal, so the definition is global.
The scale `ρ(z)` is the fixed background warp, not an unknown evolved
field. Thus the exact background has `Γ_bg=J` throughout the collar. The
box below is a pointwise physical-metric condition; background membership
does not prove membership of the actual nonlinear family.

The direct radial component chart is

```
Γ = [[n², Dᵀ], [D, M]],
n = sqrt(g_zz)/ρ(z) = sqrt(Γ_zz) > 0.
```

`M` is the normalized tangential metric and `D` the four normalized radial
cross components, including the time cross component. This is a smooth
pointwise coordinate chart on the metric variables. It is not the temporal
ADM lapse. Away from `D=0`, it is also not the Schur-complement radial
unit-normal lapse: in general `g^{zz}≠1/g_zz`. The exact shell boundary
laws are evaluated at `D=0`, where that distinction disappears. Scalar
`φ` (or `η=φ−1`) supplies the sixteenth field.

The exact source's proper-y coordinate root is `N_y=sqrt(g_yy)`.
Changing the fixed base coordinate from `y` to `z` gives
`sqrt(g_zz)=ρ(z)N_y`, so `n=N_y` throughout this fixed coordinate chart.
The primary construction uses the literal retained proper-y perturbation
flow with its boundary coefficients frozen and with the background
subtracted throughout the collar. It does not replace that flow by a raw
z-coordinate radial-lapse flow. Reference-harmonic radial Christoffel
terms would also have to be transformed in such a replacement.

## 2. A complete physical-metric patch

Assume the real symmetric matrix satisfies

```
||Γ−J||op ≤ ε = 1/40,           129 ≤ ρ(z) ≤ 131.        (P)
```

This is a full five-by-five matrix condition, including `D`. It does not
set cross components to zero in the bulk. Since `J²=I` and `||J||op=1`,
the Neumann inverse identity gives

```
||Γ⁻¹||op ≤ 1/(1−ε) = 40/39,
||Γ⁻¹−J||op ≤ ε/(1−ε) = 1/39.                        (1)
```

For example, write `Γ=J(I+J(Γ−J))` and expand the inverse, or apply the
resolvent identity `Γ⁻¹−J=−Γ⁻¹(Γ−J)J`. Along the straight segment
`J+s(Γ−J)`, no eigenvalue crosses zero; the signature is `(−++++)`.
The same argument gives the tangential block `M` Lorentzian and invertible,
with `||M⁻¹||op≤40/39`.

Let inverse components be expressed in the dual frame
`(∂z,∂t,sech(t)E₁,sech(t)E₂,sech(t)E₃)`, where `E_A` is κ-orthonormal.
Set `a=−g^{tt}` and let `c` be the full spatial inverse block in the
`z,1,2,3` slots. Applying (1) to vectors supported in each of these blocks
gives, independently,

```
a ≥ λ* := 38/(39·131²) > 1/18000 =: λ,
c ≥ λ* I > λ I,
g^{zz} ≥ λ*,
||g⁻¹||op ≤ μ* := 40/(39·129²) < 1/16000 =: μ.        (2)
```

The spatial inverse-block lower bound is a separate conclusion; it is not
deduced just from `g^{tt}<0`. Also `|a|, ||b||₂, ||c||op≤μ*`, where
`b^i=g^{ti}`. Thus the retained coefficient-size hypothesis is satisfied
with the same safe `λ=1/18000`: its coefficient upper ceiling is
`λ^(-1)=18000`. Even the sum of all 25 inverse-component magnitudes is
less than `25μ*<1/640`. This is not a claim that that sum is bounded by λ.

One convenient sufficient subpatch is

```
|n−1|≤3/1000,   ||M−diag(−1,1,1,1)||op≤3/1000,
||D||₂≤3/1000.                                      (3)
```

Indeed the diagonal perturbation has norm at most
`max(6009/10^6,3/1000)=6009/10^6`, the off-diagonal perturbation has norm
`||D||₂`, and their sum is `9009/10^6<1/40`. Thus any independent flow
certificate whose entire physical-state image lies in (3) supplies (P).
This implication includes nonzero bulk `D`; a Neumann-only patch would be
insufficient. The flow image bound is supplied by [the boundary-flow proof](BOUNDARY_FLOW_PROOF.md).

For the primary literal flow state, distinguish the fixed background
identification from this metric subpatch. Write

```
u_n = N_y−1,       u_H = (h−h_bg(z))/ρ_b²,
D_y = g_ya/ρ_b,    u_φ = φ−φ_bg(z).
```

Tangential tensor and cross components in this formula are measured in
the same `(dt,cosh(t)θ)` coframe used to define `Γ`.

If an independent state-flow certificate gives
`|u_n|,||u_H||op,||D_y||₂≤b=101/100000` throughout the collar, the exact
physical metric reconstruction is

```
n=1+u_n,
Γ_zz=n²,
Γ_ab=J_ab+(ρ_b/ρ(z))² u_H,ab,
Γ_za=−(ρ_b/ρ(z)) D_y,a.
```

The sign uses `dy=−ρ(z)dz`. Since `ρ_b/ρ(z)≤131/129`,

```
||Γ−J||op
 ≤ max(2b+b², (131/129)²b) + (131/129)b < 1/40.       (3a)
```

The replay checks this rational inequality exactly. Formula (3a) is the
needed link from the entire literal flow state, including all four cross
components, to the physical metric patch. It remains conditional on the
state-flow image bound established in the companion proof; it does not infer it
from a finite-dimensional derivative-at-zero Cauchy estimate.

The larger componentwise flow-domain box of radius `1/100` does not by
itself imply the `1/40` metric box. To see this, let `b_out=1/100`,
`r=131/129`, `u_n=b_out`, `u_H=b_out P_A` for one spatial tangential
projector, and `D_y=b_out e_A`. The radial/A two-by-two perturbation has
diagonal entries `2b_out+b_out²`, `r²b_out` and cross entry `−rb_out`.
Its Rayleigh quotient on `(1,−1)/sqrt(2)` is
`[2b_out+b_out²+r²b_out+2rb_out]/2>1/40`. This refutes an unqualified
outer-box-to-metric-box inference; it is not a signature or PDE failure.
The inner image budget `101/100000` in (3a) avoids this defect.

## 3. Positive direct-coordinate chart and finite derivatives

Patch (P) gives `Γ_zz∈[39/40,41/40]` and, conservatively,
`n∈[39/40,41/40]`. Thus the chart and inverse are uniformly smooth.
For `s=Γ_zz=n²`, the only nonlinear part of this coordinate change is

```
|d^q sqrt(s)/ds^q|
 ≤ (2q−3)!! / 2^q · (40/39)^(2q−1),     1≤q≤16,     (4)
```

with `(-1)!!=1`. This follows from the exact derivative formula and the
displayed positive lower bound for `n`; it is not a sampled estimate.
In particular the first three bounds are each below one. The inverse
`s=n²` has first derivative at most `41/20`, second derivative `2`, and
all higher derivatives zero. Derivatives of the physical-coordinate map
`g_zz↦n` gain the factor `ρ(z)^(−2q)`; the forward physical map gains
`ρ(z)²` at each fixed base point. The remaining metric components are linear rescalings and the
scalar chart is a translation. These finite chart estimates alone do not
bound derivatives of the nonlinear boundary flow.

## 4. Explicit auxiliary energy on all sixteen fields

Let `E=E_D⊕E_N` have ranks `4,12`, a fixed positive fiber metric with the
two subbundles orthogonal, and a fixed product connection preserving them.
The sixteen reduced fields share the scalar metric wave principal symbol.
An invertible state chart and the exact boundary flow preserve this fact:
the chain rule multiplies every second derivative by the state Jacobian;
after applying its inverse, the new second-order principal coefficient is
still scalar `g^{μν}I_16`. Second derivatives of the state map multiply
products of first derivatives and enter the lower-order expression.

Use the fixed reference spatial volume

```
dμ_ref = dz dμ_κ,          ∫_S3 dμ_κ = 2π².
```

It contains no physical metric determinant. In a κ-orthonormal frame put
`W_z=D_zV`, `W_A=sech(t)D_AV`, and `Q_t=|V_t|²+|W_z|²+Σ_A|W_A|²`.
The retained principal energy density and integral are exactly

```
e[V] = 1/2 [a |V_t|² + Σ_ij c^{ij}<W_i,W_j> + |V|²],
E[V](t) = ∫ e[V] dμ_ref.
```

The moving frame only adds bounded first-order terms when rewriting the
equations with the fixed product connection; it does not change this
principal quadratic form. The index sum includes all off-diagonal spatial
terms, and the positive fiber inner product includes every one of the
sixteen fields. From (2), pointwise,

```
1/2 [λ* Q_t + |V|²] ≤ e[V] ≤ 1/2 [μ* Q_t + |V|²],
(Q_t+|V|²)/36000 ≤ e[V] ≤ (Q_t+|V|²)/2.              (5)
```

For a norm with fixed product spatial metric `dz²+κ`, use
`|t|≤T=13ell/8<1/2`. The elementary factorial bound gives
`cosh(t)≤Σ(t²/2)^k≤8/7`, hence `cosh²(t)≤64/49<2`.
Consequently `sech²(t)>1/2`, and

```
Q_fixed = |V_t|²+|D_zV|²+|D_ΩV|²+|V|²,
Q_fixed/72000 ≤ e[V] ≤ Q_fixed/2.                    (6)
```

Equations (5) and (6) integrate with the same constants. The coefficient
spatial positivity constant in fixed product coordinates can safely be
taken as `1/36000`; the extra factor `1/2` in (6) is the energy's own half.

This is an auxiliary positive analytic energy. It is not the physical
gravitational energy, a physical positivity theorem, or a stability claim.
The displayed constants are coercivity/equivalence constants, not a
Gronwall exponent, an evolved norm, a forcing/lifting bound or an error
budget. Coefficient derivatives, state-flow jets and lower-order matrices
remain necessary for the latter estimates.

## 5. Exact zero shell flux for the frozen sixteen-field problem

At `z=0`, require every allowed frozen iterate to satisfy all four exact
coordinate Dirichlet trace identities `g_za=0`, `a=t,Ω`. Since its metric
is invertible, the matrix is block diagonal there, and exactly

```
g^{zt}=0,       g^{zA}=0.                            (7)
```

The state flow leaves the four coordinate Dirichlet variables unchanged.
Thus the required structural implication is `V_D=0 ⇒ g^{za}(V)=0` at the
boundary, for every boundary state admitted by the frozen iteration. It
must not be assumed only for a putative final solution.

After the exact twelve-row state flow, the homogeneous boundary data are
`V_D=0`, `D_zV_N=0`. For the sign-normalized frozen operator

```
P V = a V_tt − 2 b^i D_iV_t − c^{ij} D_iD_jV,
```

the outward-normal shell flux is, up to its irrelevant orientation sign,

```
g^{zt}|V_t|² + Σ_j g^{zj}<D_jV,V_t>
 = g^{zz}<D_zV,V_t>
 = g^{zz}[<D_zV_D,V_D,t> + <D_zV_N,V_N,t>] = 0.      (8)
```

The first equality uses (7); the second uses the fixed positive orthogonal
fiber splitting; the last uses `V_D,t=0` and `D_zV_N=0`. No sign or matrix
symmetry of the original coupled Robin matrix is needed. Tangential/time
derivatives formed with the split-preserving product structure satisfy the
same homogeneous D/N traces. The scalar principal symbol then gives the
same zero flux for these commuted fields. A principal matrix mixing D and
N would invalidate this calculation. An artificial outer wall with all
components Dirichlet has `V_t=0` and hence zero flux as well.

No tangential harmonic boundary conditions are added. The normal harmonic
completion is already among the twelve original rows; the tangential
subsidiary Robin identities follow separately from the physical equations
and Codazzi, as in the retained gauge companion.

## 6. Radial causal speed: the conservative choice c*=2

For any nonzero causal vector, `dt` is nonzero: a vector with `dt=0` is
strictly spacelike by (P). Rescale to `dt=1` and set

```
q = (dz/dt)² + cosh²(t) |dΩ/dt|_κ².
```

Because `g(v,v)≤0`, the same covariant matrix patch gives

```
0 ≥ g(v,v)/ρ(z)² ≥ q−1−ε(q+1)
                  = (39q−41)/40.
```

Thus `q≤41/39`, including for null characteristic rays. In particular

```
|dz/dt| ≤ sqrt(41/39) < 11/10 < 2 = c*.              (9)
```

This is a coordinate causal-speed conclusion on the region where (P)
holds. It supplies the conservative radial propagation constant used in
the retained support arithmetic. With `L=5ell/2`, `T=13ell/8`, it gives
`L+2T=23ell/4<6ell`, while `T+6ell=61ell/8<8ell`.
These strict arithmetic margins do not construct a common proper-time
observer cylinder or establish the patch for the actual nonlinear family.

## 7. Two necessary qualifications with exact examples

Timelikeness of `dt` alone is insufficient for this energy. In the inverse
metric use a z slot `+1`, t/x block `[[-1,2],[2,−1]]`, and remaining slots
`+1,+1`. This has one negative eigenvalue (the t/x eigenvalues are `1,−3`),
`g^{tt}=−1`, and `g^{zz}=1`, but the spatial x quadratic coefficient is
`g^{xx}=−1`. The elementary time-multiplier energy is not positive.
The full patch (P) excludes this example.

Boundary orthogonality does not imply bulk orthogonality. For smooth
`s(z)` with `s(0)=0`, `0<|s(z)|≤1/40` away from zero, take
`Γ=J+s(e_z⊗e_t+e_t⊗e_z)`. It obeys (P) and the exact boundary trace,
but `Γ^{zt}=s/(1+s²)≠0` in the bulk. The energy and causal arguments above
allow these cross terms. Equation (8) uses their zero trace only.

## 8. Source and replay record

`SOURCE_PINS.json` identifies byte-for-byte retained public input snapshots
at an immutable GitHub commit; the original papers and capsules are unchanged. The replay checks source
hashes, retained formula markers and exact rational arithmetic, including
the counterexamples and chart bounds through order sixteen. The universal
matrix statements are proved by the norm/resolvent arguments above;
the arithmetic replay is not a numerical sampling proof of them.

Specific source locations in the snapshots:

- `ANGULAR_UNIFORM_LINEAR_CONTACT.md`: lines 15–19 (coordinates and scalar
  principal symbol), 27–60 (actual rows and positive fixed fiber splitting).
- `QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md`: lines 22–37 (both energy
  margins and every-iterate trace), 53–76 (full principal energy and flux).
- `INDEPENDENT_QUASILINEAR_MIXED_PROOF.md`: lines 25–49 and 90–115
  (independent statement of scalar principal symbol, margins and flux).
- `NONLINEAR_BOUNDARY_STRAIGHTENING.md`: lines 10–28 (direct radial chart),
  51–115 (exact homogeneous split and principal-symbol preservation).
- `INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md`: lines 194–222 (both physical
  shell laws, Codazzi and no extra tangential primal harmonic rows).
- `REGISTERED_EXACT_FAMILY_BINDING.md`: lines 42–62 (actual collar,
  `129≤ρ≤131`, unscaled coordinate convention).
- `ACTUAL_POTENTIAL_JET_BOUNDS.json`: authenticated background wider-collar
  box, including shell values and the actual fixed detuning.

The closed replay includes these metric checks together with the companion
flow checks. Follow [REPRODUCE.md](REPRODUCE.md) for the pinned command,
normal/optimized comparison and corruption controls. The arithmetic report
does not replace the analytic norm and boundary arguments above.
