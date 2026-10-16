# Actual angular source inventory and coefficient adapters

Work date 2026-10-16. This is a new audit in the numerical-response workspace. Frozen capsule files were read only. The implementation in `source_adapters.py` exports explicit coefficient maps for the complete scalar bulk source, the five quadratic relative-harmonic components, Einstein matter terms, physical shell and normal gauge sources, and initial spacetime values/velocities. It does **not** yet export the full geometric Einstein source coefficients.

## Provenance and convention

All source paths below are relative to `/workspace/hdblast-research-work/common-gauge-forcing-20261016/capsule-source/`:

- `first-fields/FIRST_HARMONIC_FIELDS_UNPUBLISHED.md`, equations (5), (12)–(16): actual first frame components, exact connection, natural/reduced source definitions.
- `shell-producer/COMMON_HARMONIC_PHYSICAL_SHELL_FORCING.md`, equations (2), (6)–(8): all ten physical metric rows and scalar row.
- `boundary-review/BOUNDARY_COMPATIBILITY_AUDIT.md`, equations (8), (10), (10a): quadratic normal gauge source, normal primal forcing and trace contraction.
- `inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md`, equations (6), (7), (9), (10), (12): initial physical profiles.
- `inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md`, equations (11), (12): all initial spacetime values and velocities.
- `evolution/FULL_SECOND_RESPONSE_CONSTRUCTION.md`, equations (1)–(3), (6)–(7): reduced bulk/source and boundary conversion.

Sources are second parameter derivatives, with expansion `g=g0+eps q+eps² k/2`, `phi=phi0+eps eta+eps² w/2`. Natural QE contains `−2 eta_A eta_B`; the reduced Hessian is `QEred=QE−nabla_(A QC_B)`, with parentheses including one half. The unit-principal scaled forcing is `(-2 ell² rho² QEred, +ell² rho² Qphi)`. Here `ell=10^-6`, `x=t/ell`, `r=z/ell`, and `dy=−rho dz`.

## Exact angular identities

Let `Y=n1²−n2²`, `S=n1²+n2²`, `Z=S−1/2`, `F=Y²−S/2+1/12`, and let `H_Z,H_F` denote the sphere Hessian with trace removed. In the fixed orientation `(1,2,3,4)`, use the Weyl tensor with sectional entries `C1212=C3434=−2/3`, `C1313=C1414=C2323=C2424=1/3`, all other independent entries zero. Set `C±=(C±*C)/2`, `T±_ij=C±_ia jb n_a n_b`, and `T=T++T−`.

The scalar coefficient order is `(1,Z,F)`:

| Angular scalar | Coefficients |
|---|---|
| `Y²` | `(1/6,1/2,1)` |
| `G=|DY|²` | `(4/3,2,−4)` |
| `K=|HessY|²` | `(8,4,24)` |

The vector coefficient order is `(DZ,DF)`:

`Y DY=(1/4)DZ+(1/2)DF`, and `HessY·DY=DZ−2DF`.

The tensor coefficient order is `(kappa,Z kappa,F kappa,H_Z,H_F,T)`:

| Angular tensor | Coefficients |
|---|---|
| `DY⊗DY` | `(4/9,2/3,−4/3,1/2,1/3,−2/3)` |
| `Y HessY` | `(−4/9,−4/3,−8/3,−1/4,1/6,2/3)` |
| `HessY∘HessY` | `(8/3,4/3,8,5/2,−1,−2)` |

Every coefficient of `T` is separately the same coefficient for both `T+` and `T−`. In particular `Y HessY` has a nonzero TT projection. Keeping only scalar Hessians loses the TT part of this physical shell source when `L0*d` is nonzero, where `L0=a h_s+sigma' v/3` is defined below. This angular identity does not prove that the actual evolved radial factor `L0*d` is nonzero.

Higher sphere jets reduce with

`D_c Hess_abY=−kappa_ca D_bY−kappa_cb D_aY−2 kappa_ab D_cY`,

`D_d D_c Hess_abY=−kappa_ca Hess_dbY−kappa_cb Hess_daY−2 kappa_ab Hess_dcY`.

Consequently every quadratic source contraction of the actual first pair reduces to the scalar, vector and tensor building blocks in these three tables (with isotropic scalar multiples of `kappa`). There is no independent coexact vector source. Both TT orientations must be retained even though their present source coefficients coincide.

`check_source_identities.py` verifies every tensor component by exact polynomial remainder modulo `sum(n_i²)−1`, using the standard-library `angular_exact.py` engine. It additionally verifies the third and fourth derivative tensor identities. These are all-point identities; no finite point sampling is used for these angular claims.

## Metric amplitude convention and complete scalar/gauge adapters

Write the first metric frame components as

`qyy=A Y`, `qy0=B Y`, `qya=C D_aY`, `q00=D Y`, `q0a=E D_aY`,

`qab=L Y kappa_ab+M Hess_abY`, and `eta=sY`.

For the actual first master and gauge amplitudes:

```
A = -4u+2 f_y
B = f_t/rho-rho b_y
C = f/(rho c)+rho c d_y
D = -2(u+a f+b_t)
E = c d_t-b/c
L = 2(u+a f+theta b)
M = 2d
s = v+p f
c = cosh(t), theta=tanh(t), a=rho_y/rho.
```

All derivatives accepted by the adapter are ordinary `y,t` derivatives of these frame amplitudes. They already include differentiating displayed frame factors when composed from `u,v,f,b,d`. In the scaled solver, `partial_t=ell^-1 partial_x` and `partial_y=−(ell rho)^-1 partial_r`; the factor `rho` must also be differentiated.

Define `j=1/(rho c)`, `h=theta/rho`, and `qP_i=q^{AB}P_{iAB}`. `connection_contractions` returns `qP_y,qP_0` in primitive scalar basis `(Y²,G,K)`, and `qP_a` in `(Y DY,HessY·DY)`. Every formula is explicit in the module. On the actual first harmonic solution,

`(QC)_i=−2 qP_i`.

No upper/lower time sign is implicit: `gauge_source` returns **lower** `QC_0`; the raised time component is its negative. For independent metric jets with `DC[q] != 0`, the omitted additional term is `−2 q_iD C1^D`; the registered first solution has exactly `C1=0`.

For direct scalar use put `s=eta/Y` and define three primitive coefficient triples:

```
q²:Hess(phi0) = (
  p_y(A²−B²)+a p(D²−B²+3L²−16LM),
  p_y C²+a p(C²−2E²),
  a p M²)
```

With `K_s=a s_y−theta s_t/rho²`,

```
q:Hess(eta) = (
  A s_yy−2B(s_yt−a s_t)/rho+D(s_tt/rho²−a s_y)
  +(3L−8M)K_s−8L s/(rho² c²),
  2C(s_y−a s)/(rho c)−2E(s_t−theta s)/(rho² c),
  M s/(rho² c²))
```

Thus the complete actual scalar Hessian source is

`Qphi=2[q²:Hess(phi0)−q:Hess(eta)+p qP_y]−V''' s² Y²`.

`scalar_bulk_source` composes this formula and projects it onto `(1,Z,F)`. This retains the scalar pullback `+p f`. The actual first condition `C1=0` removes the two general scalar trace terms `2q^C_D C1^D phi0_C−2C1^C eta_C`; it is a stated on-solution simplification.

`check_scalar_gauge_adapters.py` independently forms the background frame connection, all ordered first covariant metric jets, all five `qP_i`, and the covariant scalar Hessian. Eight exact generic amplitude jets agree componentwise with these adapters as regression checks. The stronger acceptance check uses independent Laurent symbols for every radial/time metric jet, scalar jet and background value. The indexed contractions agree as symbolic polynomial identities on angular witnesses of exact scalar rank three and vector rank two. Together with the proved angular primitive closure, this establishes universal radial coefficient formulas. It is not a numerical first-wave solve or a C1 estimate.

## Physical shell and normal gauge source

At the shell set `H=h_sY`, `eta=vY`, `m00=−2XY`, `m0a=V DY`, `mab=2(L_sY kappa+d HessY)` where

`h_s=−4u+2f_y`, `X=u+b_t`, `V=c d_t−b/c`, `L_s=u+theta b`.

Define

```
A0 = a h_s²/4−sigma' h_s v/6−sigma'' v²/6
L0 = a h_s+sigma' v/3
P0 = p h_s²/4+sigma'' h_s v/2+sigma''' v²/2
```

Then

```
QK00 = (−A0+2L0 X) Y²
QK0a = −L0 V Y D_aY
QKab = (A0−2L0 L_s)Y² kappa_ab−2L0 d Y Hess_abY
Qp = P0 Y².
```

The displayed tables give every coefficient explicitly. In particular the `H_Z,H_F,T+,T−` coefficients of `QKab` are `(+L0 d/2,−L0 d/3,−4L0 d/3,−4L0 d/3)`.

For the normal gauge row set

```
tau0 = 2X+6L_s−16d
A2 = (4X²+12L_s²−64L_s d)Y²−2V² G+4d² K
QC_y = [4a h_s²−3a h_s tau0−(4sigma'/3)h_s v+(sigma'/3)v tau0]Y²+2a A2
GN = [(4sigma''/3)v²+4sigma' h_s v−10a h_s²+8a h_s tau0]Y²−4a A2.
```

`shell_sources` returns all these projected values, including the Lorentzian negative sign `−2V²G`. The primal frame rows use metric forcing `−2QK`, scalar forcing `−Qp`, and normal forcing `GN` in outward `y`; converting to inward `r` multiplies each right side by `−ell rho`. The mixed `k_y0,k_ya` Dirichlet source is zero. In the frame, the homogeneous tangential trace-free rows have zero Robin coefficient.

## Complete initial projection

Denote the initial Hamiltonian source by `Q_H=A_H Y²+D_H G`, with its exact radial profiles from source equation (6). Set

`S_l=(A_H/6+4D_H/3, A_H/2+2D_H, A_H−4D_H)`,

`lambda_l=(0,8,24)`, `f_l=−rho² S_l/[2(lambda_l−3)]`, `zeta_l=−3 f_l,y/(2p)`.

Use different names `Fmom=(18u u_t−v v_t)/rho` and

`Gmom=(12u u_yt−48a u u_t−6u_y u_t−2v_t v_y)/rho`

for the physical momentum source profiles. Let `qmom=Fmom/3`, `cY=(1/6,1/2,1)`, `f=(f0,f2,f4)`, `fy=partial_y f`.

Initial frame values have only

`k00=−2u² cY`, `kab=[sum_l f_l W_l] kappa_ab`, `w=sum_l zeta_l W_l`.

Every other initial metric sector is zero. The coordinate-time derivatives of frame components at `t=0` are

```
k_yy,t = (2rho qmom−8u u_t)cY
k_y0,t = rho[ (−26u u_y−2a u²)cY−(3/2)fy−3a f ]
k_00,t = (−8rho qmom+44u u_t)cY
k_0a,t = (11u²/2−f2/2) D_aZ + (11u²−f4/2) D_aF
k_ab,t = (2rho qmom+4u u_t) [sum_l(cY)_l W_l] kappa_ab
w_t = [rho(Gmom−Fmom_y)/p+2u v_t] cY.
```

All `ya`, Hessian-TF and TT initial velocities vanish. Scaled `x` velocities are `ell` times these time velocities. `initial_frame_data` exports precisely this data. No initial frame-time correction from `c_t` survives because `c_t(0)=0`.

The mixed initial velocity is essential: its coefficient is `−26u u_y−2a u²`, after the nonlinear lapse/shift conversion. Using the second shift velocity alone would yield the wrong `−74u u_y−194a u²` coefficients. The initial `k00=−2u²Y²` is nonzero even though the second lapse value was set to zero.

## What remains for the complete numerical source projection

The unresolved piece in this audit is the geometric Einstein source

`Q_Ric_AB = nabla_C J^C_AB−nabla_B J^C_AC+2(P^C_CD P^D_AB−P^C_BD P^D_AC)`,

where `J^A_BC=−2q^A_D P^D_BC`. It requires the actual metric two-jet (hence gauge amplitudes through order three). An exact indexed evaluator must retain the connection on the derivative slot when differentiating `P`, then apply the polynomial identities above. `einstein_matter_source` already supplies every non-Ricci Einstein component; `gauge_source` supplies `QC` for the reduced correction `−nabla_(A QC_B)`. The scalar source is complete.

The capsule provides formulas and validated upper bounds, not arrays of actual evolved master/gauge values. Numerical source values additionally require actual master `W` and gauge `f,b,d` jets on the chosen mesh. Substituting independent interval endpoints or the free comparator does not instantiate the actual source. No C4 gauge derivative is required by the direct second-source construction.

## Portable replay and acceptance scope

The replay requires only Python's standard library and accepts a new output path outside the angular source tree:

```sh
python -I -S -B angular/source-audit/check_source_adapters.py --run-controls --output /tmp/source-normal.json
python -I -S -O -B angular/source-audit/check_source_adapters.py --run-controls --output /tmp/source-optimized.json
```

Run from the numerical-response capsule root, or use an absolute script path. The two output reports must be byte-identical. Each report pins the source module, ring, checkers, this note, and the shared angular engine by SHA-256 and length. Output reuse and output writes within the angular source tree are rejected. No assertions implement acceptance, so optimized mode preserves every check.

`check_shell_initial_adapters.py` uses independent symbolic first and second fields and backgrounds. It compares the physical shell adapter to the affine lapse/unit-normal Hessian before normal derivative elimination, joins its gauge source to the full bulk quadratic gauge contraction, joins `GN` to the traced physical residual, compares all Einstein matter components to direct scalar-gradient and potential products, and reconstructs initial frame data through the ADM lapse/shift conversions. The angular tensor witnesses have exact rank seven; scalar and vector rank checks are retained by the scalar/gauge checker.

Seven targeted semantic mutations must fail with their named marker and leave no output: delete the angular TT projection, flip the lower time gauge contraction sign, omit the cubic scalar-potential term, delete physical shell TT rows, flip the normal gauge source sign, omit the initial mixed metric conversion, and halve the Einstein scalar-stress Hessian. These cover the achieved source scope. No full geometric Einstein adapter is claimed by this source report.
