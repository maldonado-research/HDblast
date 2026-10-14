# Frozen scalar shell symbol: rank defect and a candidate wave-gauge completion

14 October 2026. Independent internal calculation for the selected copy-even
Einstein–scalar shell. This note proves frozen principal-symbol statements;
it does **not** prove the curved, forced second-order contact IBVP.

## Physical rows and the role of the coordinates

At a regular shell point choose an orthonormal frame with interior normal
coordinate \(x\leq0\), outward \(n=\partial_x\), and tangential metric
\(\eta_{ab}=\operatorname{diag}(-1,1,1,1)\). Write the linear second-order
unknowns as
\[
a=h_{nn},\qquad b_a=h_{na},\qquad c_{ab}=h_{ab},\qquad w=\phi_2.
\]
Known quadratic terms belong on the right and do not change the homogeneous
principal symbol. The first-order field and the complete quadratic source
must have been converted to the chosen wave gauge before using these rows.

For an Eulerian shell displacement \(Z\), the metric junction has principal
part
\[
\frac12\big(\partial_n c_{ab}-\partial_a b_b-\partial_b b_a\big)
 -\partial_a\partial_bZ.
\]
Thus a scalar Hessian of bending survives. Straightening the shell sets
\(Z=0\) by a boundary coordinate choice and transfers its effects into the
metric/scalar coefficients and known sources. In this **fixed-shell wave
chart**, the homogeneous principal physical rows are
\[
B_{ab}=\partial_n c_{ab}-\partial_a b_b-\partial_b b_a=0,
\qquad B_\phi=\partial_nw=0. \tag{1}
\]
The normal lapse enters the actual Robin rows through background curvature,
for example \(\delta K_{ab}\supset-K_{0,ab}a/2\) and
\(\delta p\supset-p_0a/2\). These are lower differential order in (1),
not reasons to discard lapse from the complete shell system.

For a generic \(l\geq2\) angular scalar harmonic, the seven metric amplitudes
are \(f_{yy},f_{yt},f_{tt},q_y,q_t,K,G\); including \(w\) gives eight wave
unknowns before constraint/gauge reduction. The scalar projection of the
metric Israel law has four rows: \(tt\), time-gradient, angular trace and
angular scalar trace-free Hessian. The scalar junction is the fifth physical
row. This count does not apply to exceptional \(l=0\) or \(l=1\) sectors.

## Appending all harmonic conditions leaves gauge modes

The principal harmonic vector is
\[
C_A=\partial^Bh_{AB}-\frac12\partial_Ah,
\qquad h=a+c,\quad c=\eta^{ab}c_{ab}. \tag{2}
\]
Use a decaying Laplace–Fourier mode
\(e^{st+ikz+\lambda x}\), with \(\Re s>0\),
\(\Re\lambda>0\), \(\lambda^2=s^2+k^2\), and
\(q_a=(s,ik,0,0)\). Then \(q^aq_a=-\lambda^2\).
Solving the physical metric rows gives
\[
c_{ab}=\frac{q_ab_b+q_bb_a}{\lambda},\qquad
c=\frac{2q^ab_a}{\lambda}.
\]
Substituting in (2) yields exactly
\[
\boxed{C_n=\frac\lambda2a,\qquad C_a=-\frac12q_aa.} \tag{3}
\]
Consequently all five harmonic rows together supply only **one** new
independent row after the ten metric Israel rows. The four arbitrary
\(b_a\) are residual tangential harmonic coordinate modes. They are realized
by \(\xi_n=0\), \(\xi_a=b_a/\lambda\),
\(h_{AB}=2\partial_{(A}\xi_{B)}\), with \(\Box\xi_A=0\); the shell
location is fixed while its tangential coordinates change.

| Frozen block | Incoming wave fields | Physical rank | Physical plus all harmonic rank | Remaining tangential gauge modes |
|---|---:|---:|---:|---:|
| Full metric and scalar | 16 | 11 | 12 | 4 |
| Generic angular scalar | 8 | 5 | 6 | 2 |

Thus “eleven physical laws plus five harmonic laws” is a count of rows, not
a full-rank boundary prescription. The scalar version of that same mistake
is “five physical laws plus three harmonic laws.” This defect diagnoses an
incomplete coordinate prescription; it is not a physical obstruction to the
shell theory.

## A full-rank frozen completion

One concrete candidate completion is
\[
b_a|_W=0,\qquad C_n|_W=0, \tag{4}
\]
in addition to (1). The four mixed Dirichlet rows select boundary orthogonal
coordinates. In the scalar block only \(b_t,b_L\) occur, so (4) supplies two
mixed Dirichlet rows and one normal harmonic row. These are boundary
conditions, **not** the collar identities \(g_{nn}=1,g_{na}=0\) of a material
Gaussian chart. The normal lapse remains an evolved wave unknown.

With (4), (1) gives \(\partial_nc_{ab}=0\) and (2) gives
\(\partial_na=0\) at principal level. Hence the full candidate is twelve
Neumann waves (eleven metric plus one scalar) and four mixed Dirichlet
waves. Its stable-subspace determinant is
\[
\det B_{\rm full}=\frac12\lambda^{12},\qquad
\det B_{\rm scalar}=\frac12\lambda^6, \tag{5}
\]
for the row/variable orders in the accompanying checker. The scalar block
has six Neumann and two Dirichlet waves. Both maps have full rank for
\(\Re s>0\).

The constraint boundary count also closes in the **frozen homogeneous**
problem. Reduced Einstein equations imply
\[
E_{na}=\frac12(\partial_nC_a+\partial_aC_n).
\]
The physical Codazzi identity gives \(E_{na}=0\). Thus (4) supplies Dirichlet
\(C_n\), while tangential \(C_a\) have homogeneous Neumann data. The
subsidiary system is a vector wave equation at principal level; its usual
positive auxiliary energy has zero boundary flux for these conditions.
Zero initial constraint fields and their time derivatives therefore remain
zero in this flat problem. The normal Gauss constraint remains a required
initial/compatibility condition; appending it as another independent metric
boundary row would change the count and must be justified separately.

For the actual variable tension, the needed geometric Codazzi cancellation
is exact, not a constant-tension approximation. With intrinsic dimension
four and \(\Lambda(\varphi)=\sigma(\varphi)/6\),
\[
D_bK^b{}_a-D_aK=-3\Lambda'D_a\varphi
=-\frac{\sigma'}2D_a\varphi=pD_a\varphi,
\qquad p=-\frac{\sigma'}2. \tag{6}
\]
At second order, (6) requires the **complete** material metric/scalar source
and the correct conversion to wave coordinates. Omitting the scalar shell
source spoils this cancellation. The curved subsidiary boundary condition
also has lower-order connection/gauge terms; the flat argument alone is not
its estimate.

## Background harmonic source and the glancing limitation

An absolute harmonic condition \(\Gamma^n=0\) cannot simply be imposed on
the background proper-normal Gaussian coordinate: there
\(\Gamma_0^n=-4\rho_y/\rho\). Use a wave-map gauge relative to the known
background, prescribe a compatible generalized harmonic source, or first
reparameterize the background normal coordinate.

For boundary orthogonal coordinates with normal lapse \(N=\sqrt{g_{nn}}\),
an absolute harmonic normal coordinate and exact umbilic law give in
intrinsic dimension \(d\)
\[
\partial_ng_{nn}=2d\Lambda g_{nn}^{3/2}. \tag{7}
\]
The five-dimensional coefficient is \(8\Lambda\), not the four-dimensional
vacuum coefficient \(6\Lambda\). In a relative wave-map gauge about the
proper-normal background, the corresponding *linear* normal harmonic row
is instead
\[
(C_n)_1=\frac12\big(\partial_na-h_0^{ab}\partial_nc_{ab}\big)
 +4\Lambda_0a.
\]
The difference consists of necessary lower-order terms; both have the
principal condition used in (4).

Equation (5) establishes injectivity in the open Laplace half-plane. It
does **not** establish a uniform inverse for the unweighted stable amplitude
boundary map up to glancing: with \(s=\epsilon+i\), \(k=1\),
\(\epsilon\downarrow0\), one has
\(\lambda=\sqrt{(\epsilon+i)^2+1}\to0\). Neumann rows lose a positive
lower bound in that norm. Homogeneous reflected-wave energy remains
conserved; a useful Kreiss estimate must specify its first-order variables,
boundary-data norms and any permitted derivative loss. A claim of strict
uniform Lopatinski invertibility based solely on (5) would be false.

For a **fixed finite angular harmonic space**, this all-angular glancing
limitation is not an obstruction to a coefficient construction. After
rescaling to \(x=\log(r/b)\), all components have the common principal part
\(\partial_t^2-\partial_x^2\); angular derivatives and connection terms
are bounded matrix lower-order terms on a controlled collar. A bounded,
possibly nonsymmetric Robin matrix can be handled without a positive
diagonal Robin symmetrizer: in inward coordinate \(z=-x\), let
\(p=U_t+U_z\), \(m=U_t-U_z\). A Neumann block
\(U_x=A(t)U+f(t)\) gives for its boundary trace
\[
q'=Aq+p|_0+f,\qquad m|_0=p|_0+2Aq+2f.
\]
Solve this finite-dimensional trace ODE first, then use the characteristic
Volterra equations with small-step contraction and continuation over a
finite bounded interval. Dirichlet rows give \(m|_0=2D'-p|_0\).
This route needs a proved invariant finite angular space for **both the
quadratic sources and the chosen second initial constraint correction**.
Reflection symmetry alone does not prove the latter. The full transformed
source, smooth coefficients and corner hierarchy still need to be supplied.

For the selected single quadrupole seed, an equivariant first gauge correction
keeps all first fields in the \(SO(4)\) representation \((1,1)\). A natural
symmetric quadratic map therefore takes values in
\[
\operatorname{Sym}^2(1,1)
=(0,0)\oplus(1,1)\oplus(2,2)\oplus(0,2)\oplus(2,0).
\]
The dimensions are \(1+9+25+5+5=45=9\cdot10/2\); in the spin-one product,
exchange symmetry selects pairs whose total spin sum is even. The diagonal
representations give scalar degrees \(0,2,4\), and the difference-two pair
gives degree-two TT tensors. No difference-one representation occurs, so
there is no coexact-vector forcing from this one-seed symmetric quadratic
map. This statement applies to the complete natural bulk/junction maps
and an equivariant gauge conversion, not to arbitrary independently chosen
first angular displacement and velocity seeds.

Finite angular initial corrections can then be obtained from the already
verified localized second constraint correction by projecting all its fields
onto this complete finite isotypic sum. The background linear constraint
operator \(L\) is equivariant, so
\(L(P h_2)=P Lh_2=P S=S\). This angular projection preserves radial support
and copy parity; including both \((0,2)\) and \((2,0)\) preserves reflection
symmetry. It removes no part of the quadratic source. Independent angular
TT homogeneous constraint directions can still be set to zero.

The remaining proof obligations are the actual variable-coefficient Robin
matrix, a suitable coupled energy/symmetrizer or equivalent estimate,
subsidiary lower-order boundary control, fully transformed quadratic source,
and all corner compatibility conditions on the stated contact slab.
Fournodavlos–Smulevici, *The initial boundary value problem in General
Relativity: the umbilic case*, arXiv:2104.08851v1, supplies a relevant
four-dimensional constant-tension **vacuum** method comparison (identifier
verified by the parent literature review). It does not already prove this
five-dimensional variable-tension Einstein–scalar contact problem.

## Reproduction

`principal/verify_principal_symbol.py` builds the full and scalar stable
boundary matrices, verifies (3), their ranks and (5), exhibits the gauge
kernels, and checks the Codazzi/subsidiary principal identity. It uses explicit
exceptions, accepts an absolute `--output` path, and deliberately corrupts
the harmonic trace coefficient with `--inject-error`. The normal and
optimized receipts and actual injected-error receipts are retained next to
the script: each mode passes **31 positive checks and four wrong-formula
controls**; deliberate harmonic trace-coefficient corruption exits with
failure in each mode. Symbolic checks certify these identities, not the remaining
curved existence theorem.
