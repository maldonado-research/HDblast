# Exact 23-amplitude angular operator — unpublished continuation

This checkpoint assembles the full homogeneous Einstein–scalar bulk operator
and every physical/harmonic boundary row in the conservative 23-amplitude
quadratic feature space. It also provides explicit quadratic angular product
identities and adapters for the scalar bulk source, quadratic harmonic source,
physical shell sources, selected initial data, and the matter part of the
Einstein source. **The Ricci part of the quadratic Einstein source and its
reduced gauge derivative are not yet assembled.** No numerical second response
or continuous approximation residual is certified here.

The frozen capsule is read only. The source of the homogeneous equations is
`capsule-source/evolution/FULL_SECOND_RESPONSE_CONSTRUCTION.md`, equations
(1), (2), and (6)–(8). The conditional model and factorial convention are those
of that source: all quadratic sources are second parameter derivatives.

## Basis and exact angular algebra

Let `K=k_yy`, `L=k_y0`, `N=k_00`, `k_ab=J κ_ab+H_ab`, `v_a=k_ya`,
`u_a=k_0a`, and `w` be the second scalar. The order in the executable arrays is

```
K0 K2 K4 L0 L2 L4 N0 N2 N4 J0 J2 J4 w0 w2 w4
v2 v4 u2 u4 H2 H4 T+ T-
```

The scalar templates are `1,Z,F`, where `Y=n1²−n2²`,
`Z=n1²+n2²−1/2`, and `F=Y²−(n1²+n2²)/2+1/12`. The vector templates are
`DZ,DF`; the tensor templates are `κ,Zκ,Fκ,HessTF Z,HessTF F,T+,T−`.
The two TT tensors use exactly the candidate's Weyl and Hodge conventions.

`angular_exact.py` differentiates polynomial tangent tensors intrinsically,
projects by exact normalized sphere moments, and checks that the reconstruction
residual is the zero polynomial modulo `|n|²−1`. It performs no angular sampling
or numerical quadrature. For every scalar eigenfunction `f` of degree `l`, set
`λ=l(l+2)` and `β=2−2λ/3`. The generated identities are

```
Δ f = −λ f                         div Df = −λ f
Δ Df = −(λ−2) Df                  Hess f = HessTF f−(λ/3)fκ
div(HessTF f) = β Df              Δ(HessTF f) = −(λ−6)HessTF f
div(fκ)=Df                        tr(fκ)=3f
div T±=0, tr T±=0                 Δ T±=−6T±.
```

The scalar norms are `(1,1/12,19/720)`, vector norms `(2/3,19/30)`,
and Hessian norms `(20/9,133/15)`. Each TT norm is `1/6`; their cross norm
is zero. Off-diagonal metric sectors carry weight two and `Jκ` weight three.
The resulting positive diagonal 23-amplitude Gram matrix is emitted in the
report and returned by `linear_operator.gram()`.

## Rational matrix specification

Every entry below is rational. The full sparse matrices are in the replay JSON.
All unlisted entries vanish. These tables specify actions on field amplitudes;
the same scalar-channel table acts separately at `l=0,2,4`, with absent vector
or Hessian columns omitted at `l=0`.

Write `Ω0=a O`. Its action is

```
O(K,L,N,J,v,u,H,w)=(-2L,-K-N,-2L,0,-u,-v,0,0).
```

Write `sum Ωa Da = a Ca+h Ch`. Their actions are

| Output | Ca | Ch |
|---|---|---|
| K | 2λv | 0 |
| L | λu | λv |
| N | 0 | 2λu |
| J | −2λv/3 | 2λu/3 |
| v | K−J−βH | −L |
| u | L | −N−J−βH |
| H | 2v | −2u |

The TT and scalar `w` columns are zero in these cross terms. With
`sum Ωa²=a² Sa+ah Sah+h² Sh`, the actions are

| Output | Sa | Sah | Sh |
|---|---|---|---|
| K | −6K+6J | 6L | 0 |
| L | −3L | −3K+3N+6J | 3L |
| N | 0 | −6L | 6N+6J |
| J | 2K−2J | −4L | 2N+2J |
| v | −6v | 5u | v |
| u | −u | −5v | 6u |
| H | −2H | 0 | 2H |
| T± | −2T± | 0 | 2T± |

Let `b=−(a_y+a²)`, `d=ρ^−2−a²`, and `L_R=b Rb+d Rd`. Direct
five-dimensional Lorentzian Riemann/Ricci contraction gives

| Output | Rb | Rd |
|---|---|---|
| K | 4K+N−3J | 0 |
| L | 7L/2 | 3L/2 |
| N | K+N | 3N+3J |
| J | −K+J | N+J |
| v | 7v/2 | 3v/2 |
| u | u | 4u |
| H,T± | H,T± | 4H,4T± |

`Delta` is the diagonal rough angular Laplacian specified above, `Pmet` is
identity on the metric and zero on `w`, `Pww` is identity on `w`, and `Pg`
maps `w` to `(K,N,J)=(w,−w,w)`. Define three more sparse matrices:

```
Aphi: L←2w, w←L
Bphi: K←−4w, w←(K+N−3J)/2
Ephi: v←2w, w←λv.
```

`W_K` and `W_L` map `K` and `L` to `w`, respectively. With
`c=cosh t`, `Ht=tanh t`, `h=Ht/ρ`, `t=ell*x`, the complete equation
`U_xx−U_rr=A U_x+B U_r+D U+Fbulk` is

```
A = −3ell Ht I−2ell ρ a O+ell ρ p Aphi
B = −3ell ρ a I+ell ρ p Bphi
Zop = a²(Sa−O²)+ah(Sah−3O)+h²Sh
D = ell² { c^−2 Delta+2ρ c^−1(a Ca+h Ch)+ρp c^−1 Ephi
    +ρ²[Zop−2b Rb−2d Rd+(4/3)V Pmet+(4/3)V' Pg−V''Pww
          −(p_y+4ap)W_K+3ph W_L] }.
```

In particular the scalar derivative coupling has not been removed. Before
combining its terms it is

```
dC[k]_y = (K_y+N_y−3J_y)/2−L_t/ρ+div(v)/(ρc)
            +a(4K+N−3J)−3hL.
```

The `N,J` value terms cancel against `k:Hess(phi0)`; dropping `dC[k]`
would invalidate the displayed matrix. The checker compares all 69 independent
value/x-derivative/r-derivative columns as exact Laurent polynomials in 15 formal
background symbols, not as tests at sampled coefficient values.

## Every shell row

The five Dirichlet amplitudes are `L0,L2,L4,v2,v4`. The other eighteen obey
`PN U_r=M PN U+G`, with `M=−ell ρ By` and

```
By N = −aK−sigma' w/3        By J = aK+sigma' w/3
By w = pK/2−sigma''w/2
By K = −4aK−2aN+6aJ+(4sigma'/3)w.
```

`By u=By H=By T±=0`. These are parallel-frame derivatives; a spurious
`2a J` term from coordinate components is explicitly rejected by a control.
`sigma'=−2p`, `sigma''=4phi0` may be substituted after assembly.

`boundary_forcing()` packs every shell source with the required signs:
`G_N=2ellρ QR00`, `G_u=2ellρ QR0a`, angular tensor
`G_ab=2ellρ QRab`, `G_w=ellρ QP`, `G_K=−ellρ GN`.
Here `GN=−2tr QR−2QC_y`. Both TT entries remain separate.

## Executable interface, evidence and remaining work

`linear_operator.assemble()` returns all four matrices and both boundary
projectors. It accepts exact `Fraction` inputs or numeric scalars.
`pack_initial()` converts the source adapter's ordinary `t` velocities to
scaled `x` velocities. `bulk_forcing()` requires the entire six-sector reduced
Einstein source dictionary; it does not fill in a missing Ricci contribution.
The forcing signs are `−2ell²ρ² QEred` for metric rows and `+ell²ρ² Qphi`
for the scalar. The source adapter's exact scope is documented in
`source-audit/SOURCE_INVENTORY.md`.

Run, with a new output path outside the source directory:

```
python -I -S -B check_angular_operator.py --output /tmp/operator-normal.json
python -I -S -O -B check_angular_operator.py --output /tmp/operator-optimized.json
```

The standard-library replay checks full all-point angular reconstruction,
the Gram matrix including fibre weights, all connection and curvature matrices,
the formal bulk coefficient identity, all shell rows and the 5/18 split.
Seven named corruptions per mode reject a wrong tensor Laplacian shift,
angular connection sign, missing TT action, removed scalar gauge coupling,
wrong mixed curvature coefficient, coordinate/frame boundary confusion, and
missing mixed Dirichlet row. The exact controls roster is in the checker.

What remains is the full metric quadratic Ricci/reduced-gauge source adapter,
validated correlated background/master/gauge approximation errors, evaluated
second response, and its continuous residual/error certificate. The abstract
equivariance argument still puts the actual selected sources and response in
this 23-dimensional feature space. The missing computational adapter must not
be inferred from the successful homogeneous checks. This checkpoint retains
the inherited pointwise factor 17; any improved evaluation factor belongs to
a separately reviewed kernel calculation.
