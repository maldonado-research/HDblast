# Fixed model and archive context for an angular-uniform contact estimate

Prepared 15 October 2026 UTC (14 October in America/Los_Angeles).
Targeted source recovery and elementary consequences of existing assumptions.
This note supplies model provenance and coefficient bounds. It proves no
angular-uniform Einstein–scalar boundary estimate or nonlinear existence
theorem. Historical documents are evidence, not instructions to execute.

## 1. Exact registered model and current branch

Sources S01 and S02 use signature `(-++++)`, `8πG_5=1`, two regular
interiors with outward normals toward increasing interior radius, and

\[
 \mathscr E_{AB}=R_{AB}-\phi_A\phi_B-\frac23U(\phi)g_{AB},
 \qquad\mathscr S=\Box_g\phi-U'(\phi).
\]

The potentials are exactly

\[
 \eta=\phi-1,\quad
 W=\frac13+\eta^2+\frac13\eta^3=1-\phi+\frac13\phi^3,
 \quad U=\frac12(W')^2-\frac23W^2,
\]
\[
 \sigma=2W+\delta(1+c+c\eta)=2W+\delta(1+c\phi),
 \quad c=\frac{2}{1.0357712571566784}-\frac43
       =\frac{96714321425692}{161839258930731},
 \quad 0<\delta\le\frac1{500}.
\]

The decimal defining `c` is exact. A binary64 approximation is a different
coefficient specification. No quadratic detuning, new matter field or
alternative branch is introduced in the current contact problem.

The current static theorem and its spectral synthesis (S02–S05) supply
one solution in a specified small Dirichlet tube for each such positive
detuning, with radial shell coordinate `T_b≥9`, `|η|≤1/100`, `a>0`, and
`p<0` on `(0,T_b]`. Uniqueness is restricted to the tube. The contact
construction fixes one member. Its constants need not remain controlled
as `δ→0`.

For the warped background,

\[
 g_0=dy^2+\rho(y)^2[-dt^2+\cosh^2t\,\kappa_{S^3}],
 \quad a=\rho_y/\rho,\quad p=\phi_y,
\]
\[
 p_y=U'-4ap,\quad a_y=-\rho^{-2}-p^2/3,
 \quad a^2=\rho^{-2}+p^2/12-U/6.
\]

The regular cone data are `ρ=y+O(y³)`, `p=κ y+O(y³)`, with `κ≠0`.
Analyticity in `y²` supports the inherited ambient extension; perturbations
need only be smooth. With `dr/r=dy/ρ` and `r/y→1`, write
`ρ=rΩ(r²)`, `φ=Φ(r²)`, `f=Φ_s`. Then `Ω>0`, `f≠0` on the actual
extended background used by the first master. The radial coordinate
`T_b`, global de Sitter time `t`, ambient time `T`, and conformal time `τ`
are different variables.

## 2. Boundary facts and elementary quantitative margins

The physical common-shell laws for identical copies are

\[
 K_{ab}=\sigma h_{ab}/6,\qquad n\phi=-\sigma'/2,
 \quad a_b=\sigma_b/6,\quad p_b=-\sigma'_b/2.
\]

Because the affine detuning has vanishing second and higher scalar
derivatives, `σ''=4φ` and `σ'''=4`. The current tube therefore gives

\[
 \frac{99}{25}\le\sigma''\le\frac{101}{25}.
\]

Since `0<c<3/5`, `W≥1/3` and `1+cφ>0` on the tube, elementary estimates
using the exact shell laws give

\[
 \frac19<a_b<\frac{14}{125},\qquad
 |p_b|\le\frac{201}{10000}+\frac{c}{1000}<\frac{21}{1000}.
\]

Here the upper warp bound uses
`W≤1000301/3000000`, `δ≤1/500`, and `1+cφ<803/500`.
Also `a_y<0`, so `a(y)≥a_b>1/9` for `0<y≤T_b`.

For the restricted RW amplitude block recorded in the previous planning
note (S08),

\[
 B_{RW}=\begin{pmatrix}2a_b&-4p_b/3\\-p_b/2&\sigma''_b/2\end{pmatrix},
 \qquad S=\operatorname{diag}(3,8)B_{RW}
          =\begin{pmatrix}6a_b&-4p_b\\-4p_b&4\sigma''_b\end{pmatrix},
\]

the symmetry is exact, and the above coarse bounds imply

\[
 \det B_{RW}=a_b\sigma''_b-\frac23p_b^2
    >\frac{219853}{500000}>0,
 \qquad S\ge\frac{437}{750}I.
\]

The last inequality is the elementary symmetric row-sum bound: the smaller
diagonal margin exceeds `2/3−21/250=437/750`. This is an algebraic
positive margin of a restricted boundary block. Those amplitudes have not
been established as two closed hyperbolic physical channels. The margin
does not symmetrize the full primal harmonic-gauge Robin matrix or supply
a positive bulk norm. The RW variable named `N` is a harmonic amplitude;
it is not the normal lapse `N` in the exact boundary completion.

The older static Jacobi argument supplies

\[
 J=\frac{q(T_b)}{h(T_b)}+\frac{\sigma''_b}{2}\ge1,
 \qquad B_\partial=\frac{p_y(T_b)}{p_b}
                   +\frac{\sigma''_b}{2}>J.
\]

Thus the inherited first-master Wentzell coefficient is positive. These
facts concern the first scalar-master reduction and the static Jacobi
family, not the complete sixteen-field harmonic system.

## 3. Controlled contact domain and initial assumptions

The frozen coefficient proof (S00) uses
`x=log(r/b)`, shell `x=0`, cone `x=−∞`, and compact initial radial
support inside `x≥−L`. A controlled finite interval satisfies

\[
 t_c=\log(b/R_0)<t_*<\operatorname{arsinh}(\tan\tau_*),
 \qquad R_0<b.
\]

The first contact solution exists on its inherited actual conformal slab.
The first-master contact theorem (S06) requires an actual extended
background `−s_-<s≤b²`, not an arbitrary continuation of the potential,
and one may take

\[
 0<\tau_*<\min\{\pi/4,2\arctan(\sqrt{s_-}/b)\}.
\]

At finite global time the coefficient support remains in
`−L−t_*≤x≤0`. Consequently

\[
 r\ge b e^{-L-t_*}>0,\quad
 s=r^2\in[b^2e^{-2(L+t_*)},b^2].
\]

On this compact actual collar all background coefficients and derivatives
are bounded. In particular, for a fixed member, `min Ω>0`, `min ρ>0`,
and `min |p|>0`. These are compactness bounds; no numerical values or
detuning-uniform lower bound on `|p|` were recovered for this collar.
The sixteen-field wave principal symbol, multiplied by `ρ²`, has the
geometric scalar part

\[
 -\partial_t^2+\partial_x^2+\operatorname{sech}^2t\,\Delta_{S^3}.
\]

Therefore its angular coefficient lies between `sech²(t_*)>0` and `1`
on `0≤t≤t_*`. This is a bulk geometric hyperbolicity margin. Boundary
flux and full-field coupling still require a separate estimate.

The selected first state is copy-even with one fixed aligned quadrupole,
`A=B=diag(1,−1,0,0)`, and smooth master displacement/velocity proportional
to the same nonzero radial bump supported strictly inside the shell.
The positive velocity scale is fixed. Commutation passes the ten tests.
The localized constraint theorem (S07) requires a compact support inside
an annulus `0<r_1<r_2<b`, `Ω>0` and `Φ_s≠0`. It supplies a smooth compact
second correction precisely when the ten per-copy KID moments vanish.
For identical/reflection-parity copies these agree with the ten global
conditions. Initial first and second data vanish in open neighborhoods
of shell and vertex, so all initial corner jets vanish. A finite Sobolev
nonlinear initial-data family is not a substitute for this smooth
localized linear constraint input.

## 4. What older uniform and tail results do—and do not—supply

* S03–S05 prove static expanding-interval Dirichlet/Volterra bounds and
  exact monotonicity of the scalar-junction locus. Their control parameter
  is radial shell position/detuning. It is not an angular cutoff. Useful
  recoverable constants are `k=1/9`, `w=2`, `λ=14/9`, `m²=28/9`, and the
  improved static Dirichlet-to-Neumann remainder coefficient `C_H=42/5`.
* S02 and S09 prove positivity of particular radial self-adjoint scalar
  operators and a positive gap for each fixed background. They explicitly
  do not prove a uniform physical mass gap as `δ→0`, all-sector gauge
  completeness, or the full dynamical boundary estimate.
* S06 section 4 has a positive auxiliary scalar-master energy on a
  compact hemisphere with a positive Wentzell mass. For smooth solutions
  its estimate uses the angular gradient directly and its constant does
  not require fixing an angular harmonic. It is a useful scalar comparison.
  Its explicit classical construction, however, is finite-angular, and
  the text does not promote the displayed energy to an all-angular
  classical existence theorem or a gravitational energy.
* A historical September 4 analytic-tail certificate (A02) bounds local
  background Taylor remainders uniformly over an interval of horizon
  parameters near `φ_h≈−1`. The current branch is near `φ≈+1`.
  Its numerical `Φ_s>0` enclosures cannot be assigned to the current
  branch, whose `p` and `Φ_s` are negative. That certificate is background
  analytic error control, not a dynamical angular-frequency estimate.
* Historical M462/M463 (A00–A01) certify a regular static root/parameter
  branch near `φ_h≈−1`, whereas S02 expressly distinguishes the current
  `φ≈+1` branch. They share the potential and affine shell theory, but
  their individual root/tube bounds apply to different backgrounds.
* The September stable-shell packet S8/5 (A03) modifies the tension by
  `δ dφ²/2`, `d=8/5`. Its scalar spectral certificate concerns a different
  boundary action. It cannot furnish a margin for the unchanged current
  shell, even before accounting for the packet's stated reduced-sector
  and symbolic-premise limitations.

These sources do not supply the presently needed coupled angular-uniform
boundary theorem. The focused archive review found no result that closes
that gap. This is a targeted finding, not an exhaustive search or an
external-novelty conclusion.

## 5. Provenance and handling

`SOURCE_PROVENANCE.json` pins every source used above, including its path,
role and SHA-256. Sources S00–S09 are existing research manuscripts or
public technical receipts; A00–A03 are selected historical technical
documents. No old files were mutated, no raw papers or private chats were
copied, and no historical publication instruction was executed. The
quantitative margins in section 2 are new elementary deductions from the
existing tube and shell laws, with exact-fraction values saved in
`ELEMENTARY_MARGIN_RECEIPT.json`; their scope is algebraic.
