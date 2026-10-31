# Exact reduced source and actual gauge binding

31 October 2026. This is a finite implementation/proof join on the registered
geometric coefficient map. The analytic Bianchi and subsidiary arguments
already existed. The code here supplies new exact source/operator checks; it
makes no new evolution, uniqueness, accuracy, nonlinear-family or cosmological
claim.

## Inputs and independent constructions

The registered map is `GEOMETRIC_COEFFICIENT_MAPS.json`, SHA256
`34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615`.
The full natural-field capsule is authenticated by its manifest SHA256
`5121e4304fad7ab5856977b38c6ebd9fdc205900b543bc1c3d97384d246f378f`
and by every listed file's byte count and SHA256 before and after use.

The checker uses seven arbitrary quadrupole metric profiles A,B,C,D,E,L,M and
one scalar profile s. Ordinary y/t derivatives commute. Every profile is free:
70 metric jets through order three and six scalar jets through order two.
Covariant derivatives apply the physical-frame connection to every covariant
slot, including derivative slots. Their values are therefore coherent; no
independent third covariant jets are introduced. The four angular ambient
slots are tangent to the unit S3, representing three physical angular
components; the six stored slots represent five-dimensional spacetime.

The first metric's angular block is `L Y kappa + M Hess(Y)`. Its relation to
trace/tracefree amplitudes is J=L−8M/3, H=M. Quadratic output blocks already use
`ab=[J:0,J:2,J:4,H:2,H:4,T:+,T:−]`; these receive no second conversion.

The ring is rational Laurent-polynomial algebra with exact sphere and de
Sitter quotients. Background relations are

```
rho_y=a rho, a_y=−rho^−2−p²/3, p_y=Vp−4ap,
V=6rho^−2+p²/2−6a², (Vp)_y=Vpp p, (Vpp)_y=Vppp p,
c_t=c theta, theta_t=1−theta², theta²=1−c^−2.
```

The algebraic Einstein constraint is included explicitly. Scalar and tensor
angular fields are fully reconstructed after projection, with exact
coefficient comparison modulo these relations.

From the free first metric the checker constructs

```
P_dab=(nabla_a q_bd+nabla_b q_ad−nabla_d q_ab)/2,
C1_d=g0^ab P_dab,
Ric1_ab=nabla^d P_dab−nabla_b P^d_ad,
E1=Ric1−dphi deta−deta dphi−(2/3)(Vp eta g0+V q),
S1=Box eta−q.Hess(phi)−C1.dphi−Vpp eta.
```

It checks the linear Bianchi identity and explicit first-gauge/scalar
amplitudes. The quadratic correction is then built independently as

```
deltaB_q(E1)_b=−q^ac nabla_a E1_cb
              +(1/2)q^ac nabla_b E1_ac−C1^d E1_db.
```

Neither E1/S1 nor the covector wave is obtained by rearranging a target
quadratic identity.

## Exact map joins

The coefficient checks cover every indicated output mode:

| Identity | Modes |
|---|---:|
| QC alias=QC_actual | 8 covector |
| QC_actual−QC_background=2q.C1 | 8 covector |
| gauge_gradient=sym nabla QC_actual | 20 tensor |
| PC1=2P.C1 | 20 tensor |
| QEred_general=QE−gauge_gradient+PC1 | 20 tensor |
| QEred=QE−gauge_gradient | 20 tensor |
| Qphi_general−Qphi=2p q_y{}^d C1_d−2C1^d deta_d | 3 scalar |

All contractions raise with the Lorentzian signature `(+,−,+,+,+,+)`; for
example PC1_ab=2 sum_d signature[d] P_dab C1_d. Both TT output coefficients
are tested in the all-20 tensor identities. An omitted nonzero TT PC1
coefficient produces an exact rejection; this is intentionally checked here,
since TT is invisible to B0.

The independent covector-wave producer proves the arbitrary eight-mode
operator identity

```
B0(sym nabla Z)=(1/2)L0 Z,
B0(T)_b=nabla^a T_ab−(1/2)nabla_b tr(T),
(sym nabla Z)_ab=(nabla_a Z_b+nabla_b Z_a)/2,
L0=Box0+Ric0.
```

Its rough wave is a direct contraction of two covariant derivatives. The
Ricci endomorphism is `(2/3)V Z_b+p² delta_b^y Z_y`, retaining the scalar-field
radial contribution. There is no extra temporal minus sign after raising the
contracted Ricci index. The independent operator checker has SHA256
`9d99cc8a705f42cb5441bcaafbd6225688e08c615dde2276f2f4f4be34e174c0`;
its canonical normal/optimized report has SHA256
`226117368a98a7a01255547147fcbf56d9019e1e2d29880136bd5592f8a3f549`.
The source checker imports only that producer's separately pinned
`covector_wave_rows.py`, SHA256
`cb025648deda6eec4cfd77b3bb96c6cf6e0713f2b55987b3fdfff22dbb066bbf`, and substitutes the registered
QC_actual with coherent ordinary derivatives. It verifies all eight modal
operator rows and independently evaluates full angular B0(sym nabla
QC_actual) against the composed wave.

## General reduced identity and conditional companion

Every one of the eight collected coefficients vanishes in

```
B0(QEred_general)+p Qphi_general e^y
 +2 deltaB_q(E1)+2 S1 d eta+(1/2)L0 QC_actual−B0(PC1)=0.
```

This is collected from the registered QEred_general map, independent
first-residual corrections, the independently derived wave, and B0 of the
registered PC1 map. The natural eight-row identity is also replayed within
the same free-jet calculation as a diagnostic. B0 already includes its half
trace gradient; none of its tensor inputs is trace-reversed again.

The final companion

```
B0(QEred)=−p Qphi e^y−(1/2)L0 QC_actual
```

is conditional. The first-natural master equation must hold as a
neighbourhood identity, supplying E1=0, nabla E1=0, and S1=0 through the
retained first-natural factorization and its prolongations. The established
plus-sign Lie pullback transfers those natural premises. This does not
supply the separate harmonic premise. C1=0 changes the source values;
C1=0 together with nabla C1=0 is sufficient to drop B0(PC1) and replace
B0(QEred_general) by B0(QEred). A harmonic-neighbourhood identity supplies
these prolongations. QC_actual is retained, so no unnecessary change of
lowering convention occurs inside L0. Exchanging it for QC_background
inside L0 additionally requires C1 and its covariant derivatives through
order two to vanish.

Two separately implemented coherent counterexamples make the pointwise
restrictions concrete. At rho=c=p=1, a=theta=Vp=Vpp=Vppp=0, take only
A=0, Ay=Ayy=1. Then the registered QC values agree but their L0 difference
is 3Y² e^y, with scalar modal coefficients `(1/4,3/2,3)`. A second witness
uses q=L_xi g0+d(y)Y in q_00, eta=L_xi phi0, xi=f(y)Y partial_y,
f=0, fy=1 and higher f jets zero, d=dy=dyy=0,dyyy=1. It has every E1
component zero at the point, while nabla_y E1_yy=Y/2 and
2 deltaB_q(E1)_y=−Y². Its coefficient vector is `(−1/12,−1/2,−1)`.
Both use coherent ordinary jets and direct covariant differentiation.

## Executable evidence and controls

`run_reduced_verification.py` pins the checker and all three executed algebra,
angular and wave helpers, authenticates the full old capsule and map, and
runs fresh isolated standard-library Python processes with `-I -S -B`, both
normally and with `-O`. Output paths must be absent beforehand. A successful
report is canonical; normal and optimized bytes must agree. Every corrupted
subprocess must return nonzero, emit an explicit nonzero coefficient witness,
and leave its successful output path absent. Inputs are checked before and
after each subprocess and again across the complete run.

The source controls cover premature background lowering; missing factor two
and wrong Lorentzian raising in both q.C1 and P.C1; unnormalized symmetrizing;
wrong gauge-gradient sign; dropped PC1 at value level; an omitted nonzero TT
PC1 coefficient; dropped B0(PC1); QE, QEred or principal metric forcing used
in place of general reduced QE; premature harmonic Qphi; either omitted
first-natural correction; wrong wave sign or missing one half; frozen rho
and a wrong potential chain. Independent wave controls separately cover the
connection/derivative-slot/physical-scaling/temporal-contraction and Ricci
errors. The two pointwise-premise counterexamples have their own fresh
normal/optimized corruption runs.

The formal checker implements coherent q3. This does not extend the old
numerical `actual_first_jets.py` adapter or replay a numerical-error
certificate. Any actual raw-field q3 adapter is a distinct implementation
and must carry its own higher-jet closure and chain-rule checks.
