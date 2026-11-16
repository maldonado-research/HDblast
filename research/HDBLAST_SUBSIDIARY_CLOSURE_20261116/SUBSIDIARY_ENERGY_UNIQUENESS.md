# Finite-regularity covector subsidiary uniqueness with the actual Robin row

16 November 2026 UTC. Original proof for the continuation based on
`6ee27c36dac67497a2454f7c1283fd7177a91888` (merged PR144). Older roots and
public capsules are references only. The theorem below closes the linear
zero-data uniqueness gate. Its application to the actual reduced solution
requires the explicitly listed reconstruction, initial-data and shell joins;
this document alone does not assert that those separate joins are complete.

There is no new amplitude interval, evaluated coefficient ceiling, copy
factor, observer construction, clock calibration or physical confirmation.
The energy is a positive auxiliary norm on five covector components.

## 1. Exact convention and statement

Let `Sigma=[0,L] x S3`, `L=6ell`, `ell=10^-6`, and let `s in [0,T]`,
`T=13ell/8`. Use the fixed background coordinate relation

```
dy = -rho(z) dz,       rho(z)>0,       rho(0)=rho_b.
```

The shell is `z=0` and the outer wall is `z=L`. In the physical proper-y
coordinate convention, on a shell with `g_ya=0`,

```
N=sqrt(g_yy)>0,       n_out=N^-1 partial_y,
n_out^flat=N dy,      partial_z=-rho_b partial_y at z=0.
```

The vector `n_out` and its covector `n_out^flat` are distinguished here.
Indices `a` run through the four shell tangent directions, including time.
The field in the theorem is the physical covector `C`, with coordinate
components `C_y,C_a`, satisfying

```
Box_g C_B + Ric_B^A C_A = 0,                           (1)
C_y=0,       partial_y C_a-(N sigma(phi)/3) C_a=0       (2)
```

on the shell. The Ricci action and the covector connection in (1) have the
retained source convention; no componentwise scalar replacement of (1) is
made. Put

```
M(s,Omega)=-rho_b N(s,0,Omega) sigma(phi(s,0,Omega))/3.
```

Then the second row in (2) is exactly

```
partial_z C_a = M C_a.                                (3)
```

For the actual row `M` is a scalar times the identity of the rank-four
tangential covector bundle. The proof permits a general real endomorphism
`M` of that bundle, without a symmetry or sign condition.

**Linear uniqueness theorem.** Suppose the following hold on this compact
slab.

1. The physical metric is Lorentzian and has the coefficient regularity
   specified in Section 2. In the fixed product coordinates put
   `a=-g^{ss}`, `b^i=g^{si}`, `c^{ij}=g^{ij}`, with spatial indices
   `i,j=z,Omega`. There are finite positive constants `a_*,c_*` such that
   `a>=a_*` and `c^{ij}xi_i xi_j>=c_*|xi|_ref^2`. The inverse coefficients
   and their first time/spatial derivatives are bounded. The spatial
   inverse-block condition is an assumption in this general theorem; it
   does not follow from `g^{ss}<0` alone.
2. At the shell the full four cross traces vanish:
   `g_ya=0` for every tangent `a`, including time. Consequently
   `g^{zs}=g^{zA}=0` there, exactly, rather than just approximately.
3. `C` solves (1), has the regularity in Section 2, and satisfies (2) as
   its actual traces. The boundary coefficient `M` has the finite
   regularity in that section. `C` vanishes on a fixed open outer collar
   for the entire slab. Thus no separate incoming subsidiary datum is
   prescribed at the outer wall.
4. On the initial slice, `C=0` and `nabla_u C=0`, where `u` is its future
   unit normal (before a possible past-time pullback).

Then `C=0` everywhere on the slab. The proof uses only finite bounded
coefficient norms for this fixed solution. It requires no favorable Robin
sign, small Robin norm or new smallness condition.

All assertions can be pulled back by `t=-s`. The covector's time slot is
pulled back with its tensor sign as well. The wave equation, zero initial
data and shell rows are invariant under that pullback; all metric cross
traces remain zero. In the rest of the proof, the notation `s` denotes the
chosen evolution coordinate after this tensor pullback when necessary.

## 2. Finite regularity sufficient for the theorem

The following explicit, stronger-than-minimal bounded-coefficient
hypotheses suffice:

```
g, g^-1, partial g, partial^2 g bounded;
g^-1 in W^(1,infinity) in time and space;
M, D_s M, D_A M, D_s^2 M, D_s D_A M, D_A D_B M bounded,             (4)
```

with continuous coefficients if a classical formulation is desired. A
fixed smooth positive density and finitely many angular coordinate/frame
charts are understood. The statements and estimates below are intrinsic
with respect to the fixed product structure and hence independent of that
choice of charts. Bounds on `g^-1` include avoidance of metric degeneracy.

For the unknown, it is enough that

```
C in C([0,T];H^2),
partial_s C in C([0,T];H^1),
partial_s^2 C in C([0,T];L^2).                         (5)
```

The equation is then an `L2` identity, and its Dirichlet, tangential-time
and normal traces in the proof exist. One may assume smoothness first and
pass to (5), or use the direct absolutely-continuous Hilbert-space
integration by parts in Section 6. No top-order weak boundary trace is
invoked. For (4) it also suffices to use essentially bounded coefficients
with the corresponding weak derivatives and strong traces; the stronger
continuous version is ample for the current application.

The current intended physical budget is much stronger:

```
partial_s^j(g,phi) in C H^(14-j),       0<=j<=14,
partial_s^j C in C H^(13-j),           0<=j<=13.        (6)
```

This budget is an application conclusion of the actual inverse-chart and
composition join, not an additional consequence claimed here from a
generic reduced equation. In spatial dimension four, Sobolev embedding
applied to each of the total orders at most three gives continuous and
bounded spacetime derivatives through order three of `g,phi`. In
particular all physical Christoffel, differentiated Christoffel and Ricci
terms needed below are bounded. The subsidiary itself has

```
C in C H13,       C_s in C H12,       C_ss in C H11.    (7)
```

The shell coefficient and its first two time jets have, respectively,
`H13.5`, `H12.5`, `H11.5` angular regularity, by the ordinary spatial trace
theorem and smooth finite composition of `sqrt(g_yy)` and `sigma(phi)`.
On the three-dimensional sphere these imply every bound on `M` in (4),
including its time/angular mixed derivatives. No estimate at a Sobolev
embedding endpoint is needed. In particular the energy argument uses a
small part of the available finite regularity, without introducing an
order-fifteen solution norm.

## 3. The complete coupled covector operator

It is useful to state precisely what is placed in the lower terms. Let
`E` be the physical covector bundle identified, using the fixed product
coordinates, with

```
E_D = span(dy),
E_N = span(ds) direct-sum T*S3,       E=E_D direct-sum E_N.         (8)
```

This is a rank-one/rank-four subsidiary split, not the rank-four/rank-twelve
split of the primal sixteen fields. The `ds` slot in (8) includes the
usual minus sign if the original physical `dt` slot was pulled back.
Use a fixed positive bundle metric `q`, equal to a positive scalar metric
on each trivial line and the round positive metric on the sphere, with
`E_D` orthogonal to `E_N`. Let `D` be its time-independent product bundle
connection. It preserves both subbundles, has `Dq=0`, and declares the
chosen `dy,ds` basis sections parallel in `s,z`. On angular covectors it
uses the round connection.

This choice does not identify `dy` with `dz`: since `dy=-rho(z)dz`, the
physical covector connection in the `dy` basis contains the corresponding
basis derivative. That derivative is included in the difference tensor
below. Alternatively one may use `C_z=-rho C_y` throughout; zero Dirichlet
trace and zero data are equivalent, but it is unnecessary to switch bases.

Distinguish this bundle connection `D` from the torsion-free product base
connection `hat nabla` of `ds^2+dz^2+kappa`. Define

```
H^D_munu Z = D_mu D_nu Z - hatGamma^lambda_munu D_lambda Z,
T^lambda_munu = Gamma[g]^lambda_munu - hatGamma^lambda_munu,
A_mu = nabla[g]_mu - D_mu  acting on E.                (9)
```

Here `A_mu` is a full five-by-five bundle endomorphism. It need not
preserve (8), need not be `q` symmetric, and is not discarded. For example,
in a coordinate covector basis its physical connection entries contain
`-Gamma[g]^A_muB`. In the `dy` frame the additional `rho` derivative is
also part of `A`. For an endomorphism-valued base covector `A`, let
`(hat D_mu A)_nu` use `D` on its endomorphism indices and `hat nabla` on
the base covector index. Expanding the physical covector wave gives the
exact formula

```
(Box_g + Ric) Z
 = g^munu H^D_munu Z + L^lambda D_lambda Z + V Z,

L^lambda = 2 g^(lambda mu) A_mu
           - g^munu T^lambda_munu I_5,

V = g^munu [ (hat D_mu A)_nu + A_mu A_nu
             - T^lambda_munu A_lambda ] + Ric_action.            (10)
```

The notation `g^(lambda mu)` here is the ordinary inverse component
`g^{lambda mu}`, not an additional symmetrization factor. Formula (10)
follows by applying `nabla[g]=D+A` twice and retaining the connection on
the second covariant derivative's base covector index. In particular its
principal coefficient is `g^munu I_5`, while every component connection
coupling and the full Ricci endomorphism remains in `L,V`.

The coefficient hypotheses imply `L,V` are bounded. This is a qualitative
finite bound for the fixed obtained solution; the previously evaluated
normalized reduced coefficient constant is not silently reused as a
numerical bound on (10).

## 4. Exact Robin conjugation, including every derivative

Choose a fixed smooth collar function `chi(z)`, supported near `z=0`, with

```
chi(0)=0,       chi'(0)=1.
```

Extend `M(s,Omega)` as the same endomorphism in `z` on this collar, using
the product trivialization, and put

```
H=chi(z) M(s,Omega),
Q=I on E_D,       Q=exp(H) on E_N,
C=Q W.                                               (11)
```

The product connection and the extension preserve the split. The matrix
exponential is invertible with inverse `exp(-H)`. Compactness and (4)
give finite bounds for `Q,Q^-1,DQ,D^2Q`, with no sign or smallness
restriction on `M`. At the shell,

```
Q=I,       D_z Q=P_N M P_N,
D_s Q=0,   D_A Q=0.                                  (12)
```

Because `D_z` differentiates the coordinate components used in (3), the
actual shell conditions are equivalent, with no discarded remainder, to

```
W_D=0,       D_z W_N=0.                               (13)
```

Here and below `P_D,P_N` are the fixed projections in (8).

For clarity, no commuting-matrix derivative identity is assumed. For all
base directions `mu,nu`, including time and angular directions,

```
D_mu Q = Dexp_H[D_mu H],
H^D_munu Q = Dexp_H[H^D_munu H]
            + D^2exp_H[D_mu H,D_nu H]                 (14)
```

on the tangential block. The first Frechet derivative is exactly

```
Dexp_H[K] = integral_0^1 exp((1-r)H) K exp(rH) dr.
```

The second is the sum of the two ordered insertion integrals over
`0<=r1<=r2<=1`:

```
integral exp((1-r2)H) K exp((r2-r1)H) J exp(r1 H) dr1 dr2
 + the same integral with K,J interchanged.            (15)
```

Thus derivatives of the actual time/angular-dependent coefficient `M`,
both orders of noncommuting factors, and the product connection's angular
terms are included. In the actual scalar case (14) simplifies, but that
simplification is not used to delete any derivatives.

Substitute (11) in (10), multiply on the left by `Q^-1`, and use the
product rule. The resulting equation is exactly

```
g^munu H^D_munu W + Ltilde^lambda D_lambda W + Vtilde W=0,

Ltilde^nu = Q^-1 L^nu Q + 2 g^munu Q^-1 (D_mu Q),

Vtilde = Q^-1 V Q + Q^-1 L^lambda (D_lambda Q)
          + g^munu Q^-1 (H^D_munu Q).                (16)
```

All orders of the displayed matrix products are literal. In particular
the terms with `D_s Q,D_A Q,D_sD_A Q,D_AD_B Q,D_s^2 Q` are present in
(16). No derivative of a boundary-dependent matrix is frozen. The
principal coefficient is still precisely `g^munu I_5`; conjugation of
the scalar identity does not create a principal block coupling. Every
new term in (16) is first or zeroth order and bounded under (4).

The energy metric used below is the fixed `q`, not the indefinite
physical inverse metric on covectors and not an unexamined conjugate of
it. `Q` need not be `q` orthogonal. This causes no problem because (16)
has a scalar principal symbol, and its non-orthogonal derivative effects
have all been retained as bounded lower terms.

Zero initial `C,C_s` implies zero initial `W,W_s`: use
`C=QW`, `C_s=QW_s+Q_s W`, and invertibility of `Q`. The time dependence
of `Q` therefore creates no initial subsidiary datum.

## 5. The exact boundary block and metric-cross requirements

In a product orthonormal spatial frame set

```
a=-g^ss,       b^i=g^si,       c^ij=g^ij.
```

For the sign-normalized principal operator

```
P W = a W_ss - 2 b^i D_i W_s - c^ij H^D_ij W,         (17)
```

the energy boundary integrand, with outward reference spatial conormal
`nu_i`, is

```
F_nu = nu_i [ b^i |W_s|_q^2
              + c^ij <W_s,D_j W>_q ].                (18)
```

All spatial cross coefficients are retained. At the shell `nu_z=-1`.
Block inversion of the actual metric at `g_ya=0`, followed by the fixed
`y,z` coordinate change, gives exactly

```
b^z=g^zs=0,       c^zA=g^zA=0,       c^zz=g^zz>0.     (19)
```

The absence of the time cross in (19) is essential; deleting only the
three angular crosses would leave `g^zs |W_s|_q^2` in (18).

The algebraic block condition is

```
P_N^* q P_D=0,
D_mu P_D=D_mu P_N=0,
[Q,P_D]=[Q,P_N]=0.                                   (20)
```

The first identity is positive-metric orthogonality of the subsidiary
split; the second permits differentiation of its traces; the third is
the block-preserving Robin transformation. With scalar principal
`g^munu I_5`, the first identity is equivalently the radial principal
block condition `P_N^* q (g^zz I_5) P_D=0` at a noncharacteristic shell.
It is satisfied by the explicit choice in Section 3, rather than
inferred from a physical Lorentzian fiber inner product.

Differentiate (13) in time in its trace space: `W_D,s=0`. Also (13)
gives `P_N D_zW=0`. Equations (19),(20) now reduce (18) exactly to

```
F_shell=-g^zz <P_N W_s,P_D D_zW>_q=0.                 (21)
```

This proves simultaneous five-component boundary cancellation even
though the bulk lower matrices in (16) couple all components. It does
not require five independent scalar boundary problems. If either
`P_N^*qP_D` or the metric cross traces were absent, (21) would in
general contain uncontrolled products of allowed boundary traces; the
argument would not apply without a different conormal boundary system.
Likewise, a general matrix-valued second-order principal operator would
require its own symmetrizer and boundary block check. That operator is
not silently substituted for the actual scalar-principal equation here.

On the outer collar, `C=0` implies `W=0` and all its first derivatives
vanish. The outer boundary contribution (18) is therefore zero as well.
No sign, size or self-adjointness of `M` was used in (21).

## 6. Positive energy and the uniqueness estimate

Use the fixed reference measure `dmu=dz dmu_kappa`. Let

```
E(s)=1/2 integral_Sigma [ a|W_s|_q^2
              + c^ij <D_iW,D_jW>_q + |W|_q^2 ] dmu.  (22)
```

By the assumed coefficient positivity there are finite constants
`kappa_-,kappa_+>0` with

```
kappa_- (||W_s||_2^2+||D_xW||_2^2+||W||_2^2)
 <= E(s)
 <= kappa_+ (||W_s||_2^2+||D_xW||_2^2+||W||_2^2).      (23)
```

No Lorentzian component norm is used in (22). The auxiliary `|W|^2`
term supplies the zero-order control, including any spatial constant
tangential component that a Neumann condition would otherwise permit.

By (16),(17),

```
P W = F,       F=Ltilde^lambda D_lambda W+Vtilde W.    (24)
```

Multiply by `W_s` in `q` and integrate with the product base and bundle
connections. Metric compatibility `Dq=0`, symmetry of `c^ij`, and the
time independence of `D,dmu` give the exact identity

```
E'(s) = integral_Sigma [
       <F+W,W_s>_q
       +(1/2) a_s |W_s|_q^2
       +(1/2) (c^ij)_s <D_iW,D_jW>_q
       -(hat nabla_i b^i)|W_s|_q^2
       -(hat nabla_i c^ij)<D_jW,W_s>_q ] dmu
       + integral_boundarySigma F_nu dS_ref.          (25)
```

For example the mixed-time term in (17) contributes
`-integral b^i D_i|W_s|^2`, so its boundary term is exactly the first
term of (18), with the sign displayed in (25). The covariant spatial
divergence in (25) includes the round connection contributions; no
angular integration boundary exists. Equation (25) also retains the
full spatial inverse block, including the `z,A` crosses in the bulk.

By Section 5 the entire boundary integral is zero. All remaining
coefficients in (25) are bounded. Cauchy--Schwarz, the elementary
inequality `2xy<=x^2+y^2`, and (23) give

```
E'(s) <= K(s) E(s) a.e.,       K in L1([0,T]).        (26)
```

One may take a finite constant `K` depending on this solution's
coercivity margins and the suprema of

```
|a_s|, |c_s|, |hat nabla b|, |hat nabla c|,
|Ltilde^s|, |Ltilde^i|, |Vtilde|, 1.
```

This is sufficient for zero-data uniqueness; there is no requirement
to evaluate it numerically. The non-symmetry and cross-component entries
of every lower matrix are bounded by their finite operator norms in
`q`, and are included in (26).

At the regularity (5), the terms in (25) are integrable, `E` is
absolutely continuous, and the boundary pairings are legitimate:
`W_s in H1` and `D_zW in H1` have `H1/2` traces, hence `L2` boundary
traces. The Dirichlet trace of `W_D` is time differentiable as an
`H1/2` trace, so its derivative is the trace of `W_D,s`, which is zero.
This justifies (21) directly. The spatial integration by parts pairs
the `L2` second derivatives of `W` with `W_s in H1`; the required
traces and Green formula hold at precisely these orders. The time
derivative of the gradient energy is justified by `W in C H2` and
`W_s in C H1`. These facts also prove (25) by time regularization if
the coefficients have only the weak bounded derivatives in (4).

Since `C=0` on the initial slice, its tangential derivatives there
vanish. Write its physical initial evolution vector as
`partial_t=alpha u+beta^i partial_i`, with `alpha>0`. The zero value of
`C` removes the connection term, and the two initial hypotheses give

```
partial_t C = alpha nabla_u C + beta^i nabla_i C = 0.
```

Tensor pullback to elapsed `s` preserves this zero. Section 4 therefore
gives `W(0)=W_s(0)=0`, and (22) gives `E(0)=0`. Applying Gronwall to
(26) proves `E(s)=0` on the whole slab. Coercivity and time continuity
give `W=0`, and invertibility of `Q` gives `C=0`. This proves the
theorem.

The proof avoids a boundary trace estimate for `W_s` against an
uncontrolled Robin trace: it removes that term exactly by (11)--(21).
Thus no impermissible estimate `||W_s||_boundary <= const ||W_s||_L2`
or favorable signed Robin energy is hidden in the argument.

## 7. The physical shell identity supplying the theorem's Robin row

This section records the exact conditional input, not a replacement
for applying the actual primal-row equivalence. Suppose the inverse
chart has supplied, on the constructed physical shell,

```
g_ya=0,
K_ab=(sigma(phi)/6) h_ab,
n_out(phi)=-sigma'(phi)/2,
C(n_out)=0,                                           (27)
```

with `K=(1/2)Lie_n h`, and suppose the reconstructed reduced equation is

```
E_AB=Ric_AB-phi_A phi_B-(2/3)U(phi)g_AB
    =nabla_(A C_B),
S=Box_g phi-U'(phi)=0.                                (28)
```

The shell dimension is four. Codazzi in the retained sign convention
and the scalar junction in (27) give

```
Ric(n_out,e_a)=D_b K^b_a-D_a tr_h K
             =-sigma'(phi)D_a phi/2
             =(n_out phi)D_a phi,
E(n_out,e_a)=0.                                      (29)
```

This includes a nonzero tangential scalar gradient; none is deleted.
Differentiating the actual scalar trace `C(n_out)=0` tangentially gives
`(nabla_a C)(n_out)=-K_a^b C_b`. The mixed component of (28),(29)
therefore gives

```
C(n_out)=0,       (nabla_n C)_a-K_a^b C_b=0.           (30)
```

At `g_ya=0`, tangential derivatives of that zero metric trace vanish,
and

```
Gamma^b_ya=N K_a^b,
(nabla_n C)_a=N^-1 partial_y C_a-K_a^b C_b
```

because `C_y=0`. Thus (30) is exactly
`C_y=0`, `partial_y C_a-2N K_a^b C_b=0`, which is (2).
The second copy of `K` is the physical covector connection term. The
coordinate Robin coefficient is `N sigma/3`, not `N sigma/6`.

The normal harmonic row is one primal row. Derivation (29)--(30) adds
no four independent tangential harmonic rows to the sixteen-field
primal problem. The normal force/conservation consequence is likewise
not an additional incoming boundary row.

At the finite budget (6), the first derivative traces in (27) and the
second derivative traces in (29),(30) have at least `H12.5` and
`H11.5` regularity, respectively. These orders are ample for each
tangential differentiation and multiplication just used. There is no
trace of an `L2` top temporal jet in this reasoning.

## 8. The exact bulk subsidiary identity at finite regularity

For the convention (28), contracted Bianchi and the scalar product rule
give the exact identity

```
nabla^A E_AB-(1/2)nabla_B tr_g E = -S phi_B.           (31)
```

Substitution of (28) and the declared Ricci commutator convention gives
one half of `(Box_g C_B+Ric_B^A C_A)` on the left. With `S=0` this is
(1). No reduced forcing, frozen-reference approximation, omitted scalar
source or constraint-independent source is allowed in this implication.
The full varying reference connection remains in the definition

```
C^A=g^BC(Gamma[g]^A_BC-Gamma[g0]^A_BC),       C_B=g_BA C^A.
```

At (6) all total third derivatives of `g` are continuous with high
spatial Sobolev reserve, so (31) is an ordinary finite-regularity
identity. Equivalently the two sides are well-defined continuous `H11`
fields by the Sobolev product calculus. Applying the smooth identity
by local finite-order approximation gives the same result. No infinite
smoothness, infinite corner sequence or arbitrary analytic approximant
satisfying physical constraints is needed.

## 9. Application ledger for the current PR144 member

The following exact joins must be established for the selected canonical
member. They are distinct from the theorem proved in Sections 1--6.

* **Inverse chart and exact reduction.** Reconstruct the physical
  `g,phi` by the literal point map and proper-y flow at the actual
  nonzero reduced solution, including all explicit position and
  reference-connection derivatives. Prove (28), the declared
  commutator convention and (6). Then (31) supplies the unforced
  equation (1), and the coefficient regularity of Section 2 follows.
* **Initial data.** Identify this solution's actual strong initial
  pair and its equation-determined accelerations with the selected
  half-action canonical family and all five harmonic completion
  equations. The exact physical Hamiltonian and momentum projections
  are projections of `H_AB=E_AB-(1/2)g_AB tr_g E`. Their vanishing,
  with initial `C=0`, yields `nabla_u C=0`; an independent prescribed
  second metric velocity is not introduced. The old doubled-action
  canonical inverse must not be substituted for the current inverse.
* **Actual shell rows and orientation.** Apply the literal inverse
  primal-row equivalence to get (27), with the physical
  `N=sqrt(g_yy)` and the full coordinate cross trace. Establish (29),
  (30), hence (2),(3), with the actual outward normal on each copy.
  The inward-z conversion uses `-rho_b`, exactly. This also gives
  all four metric inverse cross traces (19).
* **Hyperbolic patch.** Use the actual physical metric image of the
  existing low-state patch. Its full five-by-five metric condition
  supplies positivity of `a` and of the entire spatial inverse block
  `c`, separately. Bulk cross components remain present. The retained
  metric certificate proves these facts; small temporal inverse
  coefficient alone is not substituted for them.
* **Outer collar.** The actual PR144 support theorem gives
  `U=0` on `(80ell/19,6ell] x S3` for the full slab. The inverse
  chart fixes zero, so physical `g=g0,phi=phi0` there. The defining
  relative harmonic covector is then identically zero there, with
  all first derivatives needed for its outer flux. This is the
  exact datum's support theorem, not an approximant's cutoff.

Once these joins hold, all theorem hypotheses have been checked at the
actual finite-regularity solution. The present energy proof then gives
`C=0` throughout the slab. Only at that point does (28) identify the
reduced member with the unreduced Einstein--scalar system. Identical-copy
matching and any derived force identity require the actual physical
Israel/scalar laws and both outward orientations as their separate
inputs. The theorem changes neither the original amplitude interval nor
any one/identical-two-copy normalization.

The supported smooth approximants used to construct the finite reduced
solution need not have satisfied canonical constraints. This proof is
applied directly to the final solution through its exact canonical
initial traces. It does not infer initial subsidiary data by a nonexistent
constraint-preserving approximation theorem.

## 10. Bounded source receipt

The proof is original here; the following exact retained formulas were
read at base `6ee27c36dac67497a2454f7c1283fd7177a91888`:

| Immutable repository path | SHA256 | Read scope |
|---|---|---|
| `research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md` | `74bce09080d4015b2eccf8c183b33f47b80516216a9547837e3e0e34251d7b03` | Sections 1--4: conventions, initial implication, bulk/shell identities and conditional energy sketch. |
| `research/HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md` | `7b28775e6c44b8ef1213e4789c4efb3c6b82727226cd18c5953c0d20038f9f39` | Source normalization, literal reconstruction, complete reduced lower formula, reference derivatives and point-map chain rule. |
| `research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md` | `1db76f9b9bd187cf264fcbc119107c76ceb66ff356feb7a979b7db83b01b4e5d` | Literal proper-y state/generator, shell equivalence, full cross trace and zero preservation. |
| `research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md` | `e78d8b4ece66335aabe2704b4eb7cec599a3bdb5686e6704d97e006d29d29a4c` | Sections 1--5: actual coordinate/bundle distinctions, full inverse patch, positive spatial block and exact cross traces. |

The sealed prior-root `HANDOFF.md`, `SUBSIDIARY_SOURCE_MAP.md`,
`SUBSIDIARY_SOURCE_RECEIPT.json` and
`NEXT_UNREDUCED_AND_DATA_REGULARITY_PLAN.md` were also read for the current
regularity, source bindings and prospective gates. Their bytes were not
modified.

This receipt authenticates the narrow source scope. It does not certify
background hypotheses, prove any other application join, or turn finite
algebraic replay into a continuum proof.
