# Effective results and verification ledger

16 November 2026 America/Los_Angeles; 17 November 2026 UTC.
All conclusions retain the original canonical family and its conditional
background, chart, constraint and construction hypotheses.

## Domain, metric and observer statements

Let `ell=10^-6`, `T=13ell/8`, `Q=25ell/16`, and
`q_star=6145ell/4096`. The exact margins are

```
T-Q=ell/16,
Q-q_star=255ell/4096.
```

Every parameter in the unchanged original interval
`I=[-10^(-10^500209),10^(-10^500209)]` has its actual intrinsic
normal-observer continuation on all `[0,Q] x S3`. This is the same
observer branch as the earlier short construction. Put
`r=rho_b q=-tau`, `V=x_q=-rho_b u`, with u the future physical
unit velocity. The normalized shell metric is `gbar=h/rho_b^2`.

The full tensor/ordered-derivative metric deviation satisfies

```
0<=delta<=1936*10^7 M<10^11 M<10^-20,
|DeltaA|<=6delta,          |D DeltaA|<=48delta.
```

These are actual field-size estimates on the original I and full
coordinate slab. They are not numerical first/second parameter
variation bounds. The effective continuation gives

```
|V^t+1|+|V^S|<=48delta q,
|t+q|<=24delta q^2,
|omega-n|_(R4)<=24delta q^2,
-T+ell/32<t(q,n)<=-7q/8<0 for q>0.
```

The Jacobi comparison, in actual normalized parallel screens, gives

```
S''-S=E S,       S(0)=Id,       S'(0)=0,
|E|_F<=2^15 delta,
|S-cosh(q)Id|_F<=2^16 delta q^2,
sigma_min(S)>3/4,
k_cut=rho_b^2 S^T S >=(rho_b^2/2)kappa.
```

All matrix/tensor factors use the declared full Frobenius norms;
the bootstrap for S uses its operator norm. The background identity
`S_0=cosh(q)Id` fixes the sign convention. Actual Gauss orthogonality,
the positive cuts and a degree-one angular covering give global
injectivity across all cuts. The physical inclusion of the shell is
separately known to be an embedding.

## Restarted causal bounds and actual observables

For `s0=19ell/40` the actual support theorem gives a background shell
and observer branch through the closed shadow interval `[0,s0]`.
With `d=(q-s0)_+`, the restarted estimates are

```
|V^t+1|+|V^S|<=48delta d,
|t+q|<=24delta d^2,
|omega-n|_(R4)<=24delta d^2,
(-t-s0)_+<=D_delay=d+24delta d^2.
```

The Jacobi forcing is zero through the shadow. Its integral equation
therefore starts at s0, while its background comparison remains
`cosh(q)Id`. This gives

```
|S-cosh(q)Id|_F<=2^16 delta d^2,
A_bg(q)=2pi^2 rho_b^3 cosh(q)^3,
|A_e(q)/A_bg(q)-1|<=2^18 delta d^2.
```

For the typed maps `varphi=iota^*phi`, `gamma=iota o x` and
`Phi=varphi o x=phi o gamma`, define positive cut density and actual
integrals by

```
dA_e=j_e dmu_kappa,
A_e=integral_S3 j_e dmu_kappa,
B_e=integral_S3 Phi_e j_e dmu_kappa,
Mobs_e=B_e/A_e.
```

The exact scalar reconstruction and the actual eleven-jet Taylor
identity give

```
|Phi_e(q,n)-phi_b|<=(10^7 M/11!) [d+24delta d^2]^11,
|Mobs_e(q)-phi_b|<=(10^7 M/11!) [d+24delta d^2]^11.
```

Positivity of the cut weights transfers the pointwise bound to the
mean without an additional denominator or Jacobian constant. Area
relative error and scalar-mean error are exactly zero through s0.
These upper envelopes do not prove a strict effect or arrival after
the shadow.

Actual C2 dependence holds through Q: joint C4 for phase, velocity,
screens and moved scalar; joint C3 for labels, cut metric, density and
moved curvature; and C3([0,Q]) for A, B and Mobs. The one-sided physical
frame includes the actual normal and obeys the accelerated confined
Fermi--Walker identities. The numerical curvature comparison is
intrinsic; no numerical physical ambient tidal estimate is asserted.

## Final finite verification

The selected immutable roster is 21 snapshots, 15 direct originating
member joins and six Git-pinned manifest anchors. The table records
the frozen evaluator/control rosters and their acceptance requirements.
Final executions against the actual externally pinned capsule are
reported by the Actions workflow and separate root receipts.

| Roster or acceptance requirement | Value |
| --- | --- |
| Declared closed-capsule file count, including manifest | `41` |
| Named checkpoint check roster in each mode | `280` |
| Own finite oracle check roster | `97` |
| Independent authentication-fault roster | `5` |
| Fresh full-checker mathematical-fault roster per mode | `47` |
| Normal/optimized checkpoint reports | Required byte-identical |
| Child optimization levels | Required 0 for normal and 1 for `-O` |
| Final execution evidence | [Actions workflow](https://github.com/maldonado-research/HDblast/actions/workflows/effective-observer-20261116.yml) and external root receipts |

The replay authenticates bytes and tests finite exact identities and
interfaces. It does not algorithmically prove the continuum mixed
embedding, geodesic continuation, Jacobi or global covering arguments.
Those remain the written mathematical proofs. Internal agent-assisted
review is distinct from external peer review or a formal certificate.
