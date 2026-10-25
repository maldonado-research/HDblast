# Conditional component certificate for the exact-input two-feedback candidate

The finite candidate defined here has strictly negative conditional exact
`J0` intervals at the two retained shell samples `x=9/8` and `x=5/4`.
This statement concerns the five-mean coefficient response under the retained
branch, operator, source and existence premises. A separate geometric argument
is needed for an intrinsic observable, and these two samples do not establish
a sign over the whole history.

## Frozen inputs and new candidate

The source is exact main commit
`a97052a03a188e63d77cbe468b50f56e7a2e411d`. Every member of the following
complete source closures is authenticated before mathematical imports and
after evaluation, including byte counts, SHA256 hashes, path checks, absence
of symlinks and absence of undeclared regular files.

| Closure | Files including manifest | Manifest SHA256 |
| --- | ---: | --- |
| Current coupled-mean pilot | 21 | `b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de` |
| Input precision I | 55 | `c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61` |
| Common gauge forcing P | 108 | `8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6` |
| Numerical inputs N | 112 | `928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44` |
| Geometric sources G | 169 | `3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2` |
| Forced bounds C | 49 | `72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76` |
| Total preflight T | 49 | `cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72` |

The retained full response satisfies `U=B+T U`. Here `T` denotes the actual
variable-coefficient bulk wave operator and Robin lift, with all five mean
couplings. Let `S` be the same full operator construction with the constant
matrices stored in the frozen `PROFILES.json`. Let

\[
B_0=\operatorname{free}(u_{\mathrm{registered}},v_{\mathrm{P1}}),\qquad
V=B_0+SB_0+S^2B_0.
\]

Both feedback terms use the exact registered nominal displacement and its
exact derivative. The frozen continuous P1 velocity is retained exactly in
the definition of `B0`. The old displacement/derivative P1 profiles do not
define either feedback term. Their validated envelopes are used only to
bound the norms of the exact nominal profiles.

This is a new continuous candidate. Earlier frozen capsules and candidates
remain unchanged. Its moment evaluation and Simpson remainder are supplied
by `../balanced/b2/EXACT_B2_NORMAL.json`, independently reviewed in
`../balanced/b2/B2_INDEPENDENT_REVIEW.json`.

## Space, reflection and component comparison

The mean order is `(K0,L0,N0,J0,w0)`, Gram diagonal `(1,2,1,3,1)`, and reflection
is `diag(1,-1,1,1,1)`. Only `L0` is Dirichlet. The complete physical bulk source
is reflected; matrix/source reflection is not commuted. All five interior
rows and columns remain in every comparison product. Only shell traces are
projected, and the Robin matrices have their retained Neumann projection.

For `1<=a<=5/4`, put

\[
C_a=\{x,r\ge0:x+r\le a\},\quad c=a-1,\quad
\kappa=(a^2-1)/2.
\]

These triangles lie inside the retained physical cone. All response and input
differences used below vanish on the retained zero wedge `x+r<=1`.
Use the ten-state vector

\[
y(z)=\begin{pmatrix}\sup_{C_a}|z_i|\\
\max(\sup_{C_a}|z_{i,x}|,\sup_{C_a}|z_{i,r}|)\end{pmatrix}_{i=1}^5.
\]

For any operator `O` of this wave/lift form, let `A_O,B_O,D_O,M_O` be its
entrywise nonnegative coefficient bounds and put `C_O=A_O+B_O`. Then

\[
H_O(a)=
\begin{pmatrix}
\kappa D_O+cM_O&\kappa C_O\\
aD_O+M_O&aC_O
\end{pmatrix},\qquad
y(Oz)\le H_O(a)y(z).
\]

The bulk value constant `kappa` bounds the active reflected triangle area:
physical multiplicity is at most two and cancels the wave factor `1/2`.
The bulk derivative bound is `a`. The lift value and derivative constants
are `c` and `1`. No coercivity, symmetric cancellation or inverse Robin
matrix is assumed.

Actual coefficient boxes are reassembled from the certified shell Taylor
intervals and validated fifth-order remainder, with the retained potential
normalization. Nominal matrix values match the frozen profile matrices
exactly; the reconstructed whole-cone Frobenius bounds match frozen `REPORT`
exactly. Let `Delta=T-S`; its entrywise coefficient bounds use the actual
intervals minus the fixed nominal matrices, before taking magnitudes.
Write `H=H_T`, `D=H_Delta`, `Hstar=H_S`.

The maximum row sums of the actual ten-state comparisons are below
`3.08800655e-4` at `9/8` and `3.14241154e-4` at `5/4`. Each comparison used on
each integration cell is separately checked to be a positive contraction.
The exact inverse

\[
R(a)=(I-H(a))^{-1}\ge0
\]

is checked entrywise nonnegative and against both exact matrix inverse
identities. Thus it bounds the state of the actual operator resolvent.

## Input and source accounting

Let `beta` bound `y(B0)`. Its displacement and derivative norm caps come from
the frozen profile values plus their existing whole-cell envelopes; the P1
velocity norm is the maximum of its nodal magnitudes. No new profile or
initial-data computation is performed.

Let `delta0, deltar, deltav` be component error caps from the three-channel
accepted input precision certificate, summing channels within each original
cell before taking maxima. Let `epsilonv` be the frozen velocity P1 error.
The component source caps follow from

\[
F_{\rm cap}=17387513156573,\quad G_{\rm cap}=3851400,\quad
\tfrac32F_{\rm cap}+G_{\rm cap}=52162547172519/2.
\]

Each Hilbert source cap is divided by the square root of the corresponding
Gram weight using an outward rational bound. The input state bound is

\[
e_{0,a}=\delta_0+a(\delta_v+\epsilon_v)+\kappa F+cG,
\qquad
e_{0,d}=\delta_r+\delta_v+\epsilon_v+aF+G.
\]

It bounds `B-B0`, including the entire omitted physical forcing. This input
is propagated once through the full resolvent. No source term is discarded,
and the displacement error is kept separate from the much larger initial
scalar derivative error.

## Noncommutative error identity

Let `R_T=(I-T)^-1` on the retained space. Since

\[
T^2-S^2=(T-S)T+S(T-S)=T(T-S)+(T-S)S,
\]

the exact error identity is

\[
U-V=R_T(B-B_0)+T^3R_TB_0+
[\Delta+\Delta T+S\Delta]B_0.
\]

No commutation of `T` and `S` is used. The direct ten-state upper bound is

\[
R e_0+H^3R\beta+[D+DH+Hstar D]\beta.
\]

The alternate ordered coefficient bound is `[D+HD+D Hstar] beta`.
Taking the componentwise minimum of two valid complete coefficient bounds
preserves an upper bound. The full resolvent remains in the third-order
tail; replacing `H^3 R beta` by `H^3 beta` is invalid and is explicitly rejected
by the software controls.

## Sharper shell integrals

For a Neumann shell row, the reflected kernel is exactly the physical
triangle integral

\[
(Kf)(x,0)=\int_1^x\int_0^a f(s,a-s)\,ds\,da.
\]

Define the rectangular matrix

\[
\Lambda_O(a)=[aD_O+M_O,\;aC_O].
\]

Then `|Oz(x,0)| <= integral_1^x Lambda_O(a) y(z|C_a) da`.
For two possibly different operators define

\[
W_{O,P}(a;x)=a[D_O,C_O]H_P(a)
 +(x-a)M_O\Lambda_P(a).
\]

This bounds the shell value of `OPz` after integration. The second term is
the nested shell lift, obtained by Tonelli from
`integral_1^x M_O integral_1^a Lambda_P(v)y(z|C_v)dv da`.
Its weight is `x-a`, not `a-1`.

The third-order tail is bounded by

\[
\int_1^x W_{T,T}(a;x)H_T(a)R(a)\beta(a)\,da.
\]

The two complete coefficient routes are

\[
\int_1^x[\Lambda_\Delta+W_{\Delta,T}+W_{S,\Delta}]\beta\,da,
\qquad
\int_1^x[\Lambda_\Delta+W_{T,\Delta}+W_{\Delta,S}]\beta\,da.
\]

For the resolvent input, one route is the direct shell input plus
`integral Lambda_T R e0`. A second route uses
`R_T=I+T+T^2 R_T`: direct input plus its first feedback plus
`integral W_T,T R e0`. In the first feedback, the shell input is bounded
directly by

\[
\delta u_0(a)+\int_0^a(\delta v_0+\epsilon_v)\,ds
 +(a-1)\,\mathrm{SOURCE}_i.
\]

Its bulk term retains the full value/derivative bound `e0`. This avoids
substituting an initial derivative error for a shell displacement error.
The smaller of the two complete input bounds is used componentwise.

The frozen denominator-256 footprint supplies 29 intervals through `9/8`
and 61 through `5/4`. On each interval, actual coefficient and state bounds
at its right endpoint cover every smaller causal triangle. The comparison
weights are integrated as exact polynomials; the cubic tail kernel has
degree at most five. Accumulations round upward to a rational `2^-112` grid.
No floating-point quantity enters an acceptance decision.

## Measured conditional outcome

All displayed upper caps in this table are rounded outward.

| Charge to shell `J0` | `x=9/8` | `x=5/4` |
| --- | ---: | ---: |
| Complete resolvent input and omitted forcing | `<2.43300162e14` | `<5.34211165e14` |
| Full third-order response tail | `<1.26749137e14` | `<1.71520934e15` |
| Actual/constant coefficient defects, both feedback orders | `<3.56213507e14` | `<3.87841685e15` |
| Complete actual-response error | `<7.26262805e14` | `<6.12783735e15` |

The separately evaluated candidate intervals already contain arithmetic
and validated quadrature errors. Adding the actual-response error exactly
once gives the outward conditional exact-response intervals

| Shell sample | `J0` interval |
| --- | --- |
| `9/8` | `[-1.14562550e16,-1.00037137e16]` |
| `5/4` | `[-2.81184636e17,-2.68928911e17]` |

Both intervals are strictly negative. The software checks this using exact
rational endpoints rather than their decimal displays.

The simpler constant-cone ten-state bound also independently excludes zero:
its `J0` value-error caps are `<7.15287306e15` and `<2.57118587e16` at the two
samples. Thus the negative sample conclusion does not rely on the sharper
shell integrations. The final review retains both sets of exact joined
intervals.

## Reproducibility and independent controls

Normal and optimized Python produce identical error-ledger bytes and
identical review/join bytes. The candidate producer also matches under both
modes and has its separate independent 55-control review. This ledger adds
34 exact controls: independent five-node Boole integration reconstructs the
ordered mixed kernels and the full third-order polynomial on three intervals
using noncommuting positive fixtures with retained `L0` interior couplings.
Wrong shell weights, missing third operators, reversed order, a wrong
noncommutative identity and a missing resolvent are detected.
The final review also records independent analytic approval of the ordered
identity, both input routes and all nested shell kernels. The shell matrix
is time-independent, full interior couplings remain, cell-right envelopes
dominate each interval, and integration begins at the actual wedge boundary
`a=1`, as required by that review.

| Artifact | SHA256 |
| --- | --- |
| `component_ledger.py` | `8a9633409544cc52dc4436ebe322126003a2f205a79f48928afad5d572aaf460` |
| `second_feedback_ledger.py` | `49723306309358e49355812e9760d3d34cc974584d47984b22eabef4113e89cc` |
| `check_second_feedback.py` | `48b702800a4e673ceaf8a21013592330f28695c884408996ff6c3353279d9b3d` |
| `SECOND_FEEDBACK_NORMAL.json` and optimized twin | `23883937b7e5a55f1155e8b5c795644253ed5cf0abd51e0c40a09ece50220974` |
| `SECOND_FEEDBACK_FINAL_REVIEW_NORMAL.json` and optimized twin | `450c0f6e56690180e646f55377f14c1d1fdf78ff191f97a596fd61b1c0aba31f` |
| Balanced B2 candidate, normal and optimized | `3db587c30dca20f01e4964ac87f8f5ba457d55d15f7720d09eb840ed3e0e4935` |
| Balanced B2 independent review | `1da9aa9d9e56cfb91a83284a560e38c10ef644055e206301bef2f8848d9a94ea` |

Replay each producer with `python -I -S -B`, adding `-O` for the optimized
mode, and fresh `--output` paths. Both ledger scripts accept `--repo-root`.
The check script expects the authoritative normal and optimized report
filenames above and authenticates their exact pins before joining intervals.

These are conditional five-mean shell coefficient certificates. They do not
supply the remaining angular response, a pointwise Einstein residual, an
exact nonlinear finite-amplitude family, or physical evidence for a proposed
cosmological origin. Internal analytic and software review is distinct from
external peer review and proof-assistant formalization.
