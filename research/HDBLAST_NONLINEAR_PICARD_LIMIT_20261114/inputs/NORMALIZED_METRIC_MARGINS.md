# Fixed-product normalized principal margins

14 November 2026 UTC. Conditional on the retained all-state metric/flow
patch, these are algebraic consequences of the existing inverse-metric
certificate. They establish no evolved state membership or solution.

Use the actual same-state elapsed reflection, `U(s)=J_E V(-s)`, on
`0<=s<=T=13ell/8`, `ell=10^-6`, and the fixed spatial product metric
`dz^2+kappa`. The metric certificate writes the actual inverse as
`q_frame=rho(z)^(-2) Gamma^(-1)` in the moving dual frame
`(partial_z,partial_t,sech(t)E_A)`, with

```
||Gamma^(-1)-J||op<=1/39,       ||Gamma^(-1)||op<=40/39,
38/39<=-Gamma^(-1),tt<=40/39,
Gamma^(-1),spatial>=(38/39)I,   Gamma^(-1),zz>=38/39.
```

Here the time and full spatial-block lower bounds are independent
consequences of the `||Gamma-J||op<=1/40` metric patch. The mixed time-
spatial block of J is zero, so its inverse off-diagonal vector has norm
at most1/39. The flow-image certificate supplies the metric patch for
every allowed point state `|U|star<=10^-3`; background membership alone
is not substituted for all-state membership. The elapsed reflection
reverses the relevant cross signs and field slots, preserving these bounds.

Since `|t|<=T<1/2`, the retained factorial-series estimate gives
`cosh(t)<=8/7`, hence `sech^2(t)>=49/64`. For
`S=diag(1,sech(t),sech(t),sech(t))` on spatial derivative slots,

```
c_fixed=S c_frame S,             b_fixed=S b_frame,
a=-q^tt>0,
P_v u=u_ss-A^i(x,v)D_i u_s-B^ij(x,v)H_ij u,
A=2b_fixed/a,                   B=c_fixed/a.
```

Normalize at the same base point. The common positive `rho(z)^(-2)`
factor cancels in both ratios. Using global extrema129 and131 separately
would give weaker valid bounds; it is unnecessary here. Thus

```
||A||_2<=2*(1/39)/(38/39)=1/19<1,
(931/1280) I <= B <= (20/19) I,
19/20 <= B^zz <=20/19,          (B^zz)^(-1)<=20/19<2.
```

The lower spatial factor is `(38/40)*(49/64)=931/1280>1/2`.
The sharper radial lower bound uses the unchanged radial slot. No
high derivative of `(B^zz)^(-1)` is supplied here.

For a fixed positive bundle E and full covariant spatial derivative,

```
E_v[u]=1/2 integral (|u_s|^2+B^ij<D_i u,D_j u>+|u|^2),
Z1[u]^2=||u_s||L2^2+||Du||L2^2+||u||L2^2,
(931/2560) Z1[u]^2 <= E_v[u] <= (10/19) Z1[u]^2,
Z1[u]^2/4 <= E_v[u] <= Z1[u]^2.
```

All off-diagonal spatial terms are retained. The same scalar quadratic
inequalities hold on the one-copy and identical-two-copy Hilbert union.
This is an auxiliary mathematical energy, without a physical energy
interpretation.

At `z=0`, the four homogeneous Dirichlet predecessor fields imply the
retained exact physical identities `b^z=c^zA=0`. Division by positive a
preserves `A^z=B^zA=0`. The orthogonal preserved D/N split therefore
cancels the shell flux

```
(A dot n)|u_s|^2/2+B^ij n_i<D_j u,u_s>.
```

At the outer all-Dirichlet wall `u_s=0`. This flux assertion assumes the
actual homogeneous rows and uses the scalar principal coefficient on all
sixteen components. It does not apply to an arbitrary principal matrix
mixing the two boundary subbundles.

## Read-only source basis

At immutable base `d27cfd999aa2e57995c08b1228a73d3da7c973d1`:

- `research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md`,
  Sections2,4,5: the all-state metric reconstruction, inverse bounds,
  fixed/moving spatial conversion and exact flux.
- `research/HDBLAST_TEMPORAL_JETS_SEED_20261114/NONLINEAR_INITIAL_JETS.md`,
  Sections1,2: same-state time reflection and the normalized exact equation.
- `research/HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md`:
  the full-state normalized coefficient maps and fixed-product conventions.

The input byte hashes are recorded in `SOURCE_PINS.json`. Hash authentication does not prove inherited
analytic premises. State evolution, commutator production, normal
recovery, nonlinear iteration and observer conclusions are separate.
