# A finite analytic representation of the prescribed incoming state

7 October 2026. Primary mathematical derivation. This is an additive work
product outside the repository and every frozen packet.

The audited incoming-state prescription determines a unique continuum
linear-response target. A split into a flat-endpoint cap and eleven analytic
panels gives a finite representation whose **source representation error**
contributes less than `1.62479e-18` to the inherited finite-band direct pressure
and less than `4.09374e-21` to direct density, uniformly for `a<=t<=b` and
`K in {64,128,256}`. These are exact rational analytic upper bounds. They
require no source sampling, saved quantum-array decoding, or physical
trajectory calculation.

This result does **not** enclose the retained incoming arrays. The 1,243
Taylor coefficients per source in the finite representation have not been
evaluated here. Their arithmetic enclosure, oscillatory integration, and
comparison with the retained binary inputs remain separate tasks. The full
twelve-case certificate is UNRESOLVED, metric calibration FAIL, the actual
incoming-state error NOT_ENCLOSED, a higher-dimensional Big Bang cause
NOT_ESTABLISHED, and external mathematical novelty NOT_ASSESSED.

## 1. Repository identity and target

The source-level audit identifies the same producer bytes in:

* `research/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER/inputs/reference_code/forced_metric.py.txt`;
* `research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/independent/forced_metric.py`.

Their producer SHA256 is
`4fa625f2b534d741e2a48d6f0bfe84a89ee17bf306805906f962a8552e005906`.
The first file's lines 25–31 specify the original time domain, source family,
and settings; lines 253–256 identify the physical mode conversion; lines
368 and 387–395 give zero initial variations and the native per-panel
Duhamel update. This derivation reads source text only. The complete
source-provenance audit supplies the inventory and lineage checks.

The fixed background has `a0=L=-1/t`, `L'=L²`, physical mass squared two,
minimal coupling zero, and background modes

\[
v_0(t,k)=e^{-ikt}/\sqrt{2k},\qquad v'_0=-ikv_0.
\]

The producer uses `delta_v=v0*u`, `delta_v'=v0*(w-ik*u)`, where `u,w`
already contain the represented epsilon. Define `U=u/epsilon`,
`W=w/epsilon`, using that same nonzero represented value. Its declared
linear target is

\[
U'=W,\qquad W'=2ikW-g,\quad U(-6,k)=W(-6,k)=0,
\]
\[
g=4L^2h-2Lh'-h'',\qquad
h(s)=B(s+4)\quad\hbox{or}\quad(s+4)B(s+4),
\]
\[
B(z)=\begin{cases}e^{1-1/(1-z^2)},&|z|<1,\\0,&|z|\ge1.\end{cases}
\]

Thus the target at `a=-9/2` includes the prehistory `[-5,a]`. It is not a
new vacuum specified at `a`. WKB frequencies in the subtraction inventory
do not define this physical incoming state; no adiabatic initialization or
gap hypothesis is used in this theorem. Replacing the audited prescription
with a finite-time adiabatic state would change the target.

The exact normalization convention for epsilon and the measure constant
must be retained when comparing to binary data. The estimates below use
only the inherited `Pi>=3`; they do not substitute mathematical pi for a
represented constant. This is the same scaled `a0^4/epsilon` linear stress
target as the 7 October incoming-state theorem.

## 2. Exact initial-state identity and the removable infrared limit

For real `k`, define the entire propagator

\[
E_k(\tau)=e^{2ik\tau},\qquad
\Phi_k(\tau)=\frac{e^{2ik\tau}-1}{2ik}
 =\tau e^{ik\tau}\operatorname{sinc}(k\tau),\qquad \Phi_0(\tau)=\tau.
\]

Here `sinc(x)=sin(x)/x` with value one at zero. The unique target is

\[
W_*(a,k)=-\int_{-5}^{a}E_k(a-s)g(s)\,ds,\qquad
U_*(a,k)=-\int_{-5}^{a}\Phi_k(a-s)g(s)\,ds.\tag{1}
\]

The source is smooth and integrable on this compact interval. Therefore
both functions are entire in `k`; differentiation under the integral is
justified on every compact complex momentum set. At zero,

\[
W_*(a,0)=-\int g(s)\,ds,\qquad
U_*(a,0)=-\int(a-s)g(s)\,ds.\tag{2}
\]

Since `g` is real, `Re Phi_k=sin(2k tau)/(2k)`, and therefore

\[
c_*:=\operatorname{Re}U_*-
             \frac{\operatorname{Im}W_*}{2k}=0\qquad(k>0).\tag{3}
\]

This is an identity of the exact prescribed response, not a projection of
stored values. In the free basis its first-order Bogoliubov coefficients
obey `delta alpha=U_*-W_*/(2ik)` and
`delta beta=e^(-2ika)W_*/(2ik)`, so (3) is the first-order canonical
normalization condition. It is not an exact nonlinear normalization theorem
for a finite-epsilon truncated mode.

The coordinate `A_*=W_*/(2ik)` need not be bounded at zero. If
`G0=integral g` and `G1=integral(a-s)g`, then

\[
A_* = \frac{iG_0}{2k}-G_1+O(k).
\]

The pole is removable for `U,W`, but in general not for `A`. The correct
stress measure is `k² dk/(2 Pi²)`, so an `O(1/k)` bound for `A` is
integrable under the existing homogeneous-state theorem. Neither a uniform
`sigma` nor a condition `G0=0` is assumed.

## 3. A flat-endpoint cap with rational bounds

Set `x=s+5`, `delta=1/128`, `s_delta=-5+delta`,

\[
D=\delta(2-\delta)=\frac{255}{16384},\quad
\lambda=\frac1{5-\delta},\quad
p_\delta=\frac{2(1-\delta)}{D^2},\quad
H=\left(\frac38\right)^{63}.
\]

The bump increases on `0<x<=delta`, and

\[
B(-1+\delta)=e^{1-1/D}
 =e^{-16129/255}<e^{-63}<H.
\]

The last step follows from the elementary positive exponential series
`e>8/3`. On the cap, `|h|<=H`, `0<L<=lambda`, and the endpoint data satisfy
`|h(s_delta)|<=H`, `|h'(s_delta)|<=(p_delta+rho)H`, where `rho=0` for
`positive_B` and `rho=1` for `signed_uB`. The latter deliberately uses
the triangle bound `|(zB)'|<=|B|+|B'|`.

Let `E=E_k(a-s)` and `Phi=Phi_k(a-s)`. Since `E_s=-2ikE`,
`Phi_s=-E`, and `L_s=L²`, two integrations by parts give

\[
\int Eg\,ds=
 [-E(h'+2Lh+2ikh)]
 +\int E(4k^2-4ikL+6L^2)h\,ds,\tag{4}
\]
\[
\int \Phi g\,ds=
 [-\Phi(h'+2Lh)-Eh]
 +\int[6\Phi L^2-2(L+ik)E]h\,ds.\tag{5}
\]

Brackets denote upper-minus-lower endpoints. All lower boundary terms
vanish because the bump and its derivatives are flat at `-5`. Equations
(4)–(5) retain the full `-2Lh'` forcing term. Dropping it changes both
the boundary term and the coefficient `6L²`.

With no source evaluation, (4) proves

\[
|W_{cap}(a,k)|\le C_{cap}(k):=
H\{p_\delta+\rho+2\lambda+2k
 +\delta(4k^2+4k\lambda+6\lambda^2)\}.\tag{6}
\]

This is monotone for `k>=0`. For a separate complex initial-mode bound,
`|Phi|<=a-s<=1/2` in (5) yields

\[
|U_{cap}(a,k)|\le
H\{\tfrac12(p_\delta+\rho+2\lambda)+1
 +\delta(3\lambda^2+2\lambda+2k)\}.\tag{7}
\]

The omitted cap is one consistently omitted **real source**, so its exact
`c_cap` is zero. Its effect may be enclosed through
`|A_cap|<=C_cap(k)/(2k)`. This correlated conclusion would not apply to
independently rounded `U` and `W` endpoint values.

No analytic disk crosses `-5`. In particular, all Taylor derivatives of
the extended bump vanish at that endpoint although the bump is positive
immediately to the right; a Taylor model centered at the endpoint with a
positive holomorphic radius would be invalid.

## 4. Eleven analytic panels, with explicit source bounds

Cover `x in [delta,1/2]` by this exact rational construction:

1. Start with `left=delta`.
2. Set `right=min(3*left/2,1/2)`; use the panel `[left,right]`.
3. Replace `left` by `right`, until it is `1/2`.

There are eleven panels. For each, put

\[
x_c=(left+right)/2,\quad r=(right-left)/2,\quad R=left/2,
\]
\[
s_c=-5+x_c,\quad Q=1-x_c+R<1,\quad T=5-x_c-R>0,
\quad D_c=1-Q^2>0.
\]

The real halfwidth satisfies `r/R<=1/2`. On the complex disk
`|s-s_c|<=R`, one has `|s+4|<=Q<1`, `|s|>=T`, and
`|1-(s+4)^2|>=D_c`. For every complex `w` with `|w|<1`,

\[
\operatorname{Re}\frac1{1-w}-\frac12
   =\frac{1-|w|^2}{2|1-w|^2}>0.
\]

Apply this with `w=(s+4)^2` to obtain `|B(s+4)|<sqrt(e)<2`.
The last numerical inequality has an explicit rational series proof: for
`n>=1`, `n!>=2^(n-1)`, strictly for `n>=3`; the bound follows inductively
since the next factorial ratio is `n+1>=2`. Hence
`e=1+sum_(n>=1)1/n! < 1+sum_(j>=0)2^-j = 3 < 4`.
This proves a bounded holomorphic extension on each disk without evaluating
a physical source value.

Writing `z=s+4`, `d=1-z²`, the exact derivatives are

\[
B'/B=-2z/d^2,\quad
B''/B=-2/d^2-8z^2/d^3+4z^2/d^4=(6z^4-2)/d^4.
\]

Consequently the following **rational** disk bounds suffice:

\[
M_B=2\left[\frac4{T^2}+\frac{4Q}{TD_c^2}
 +\frac2{D_c^2}+\frac{8Q^2}{D_c^3}+\frac{4Q^2}{D_c^4}\right],\tag{8}
\]
\[
M_{zB}=Q M_B+2\left[\frac2T+\frac{4Q}{D_c^2}\right].\tag{9}
\]

Equation (9) follows from `g[zB]=z*g[B]+2B/s-2B'`. Every derivative
and geometry term is included. Let `P_j` be the **exact** Taylor polynomial
of `g` at `s_c` through degree `N=112`. Cauchy's coefficient bound gives

\[
|g(s)-P_j(s)|\le
 M_j\frac{(r/R)^{113}}{1-r/R}\le M_j 2^{-112}
 \quad\hbox{on the real panel}.\tag{10}
\]

Thus the complete interior source residual has the explicitly computable
rational L1 bound

\[
G_{int}=2^{-112}\sum_j(right_j-left_j)M_j.\tag{11}
\]

The piecewise polynomial is an auxiliary representation inside the
prehistory integrals. It does not replace the physical source or its action
contacts on `[a,b]`. Polynomial jumps at panel endpoints have measure zero
in these integrals; derivatives of a piecewise replacement are not taken.
There are `11*(112+1)=1243` exact real coefficients per source.

## 5. The proved state and stress representation errors

Define `U_T,W_T` by (1), omitting the cap and integrating each exact real
polynomial `P_j` on its panel. The representation is finite and unambiguous,
and is entire in `k` with the same removable `Phi` limit. Because the
source residual is real, `c_*-c_T=0` exactly. Equations (6), (10), and
`|E|=1` give, for `0<k<=K`,

\[
|W_*-W_T|\le C_K:=C_{cap}(K)+G_{int},\qquad
|A_*-A_T|\le\frac{C_K}{2k}.\tag{12}
\]

A full complex `U` bound follows from (7) plus `G_int/2`. Starting at
`a`, propagate both incoming directions with identical real forcing and
identical full action/subtraction contacts. Their difference is homogeneous;
the 7 October theorem applies without modifying pressure or using Ward
conservation to define it.

The exact phase norms can be bounded by

\[
N_R\le L+\frac{3L^2}{2k},\qquad
N_P\le\frac{2k}3+\frac{5L^2}{4k}.
\]

The second inequality follows by squaring: the excess is
`21L^4/(16k²)>=0`. Integrate (12) against `k²/(2 Pi²)`, then use
`L<=2/7`, `Pi>=3`, to obtain

\[
|R_*-R_T|\le C_K\left(\frac{K^2}{252}+\frac K{294}\right),\tag{13}
\]
\[
|P_*-P_T|\le C_K\left(\frac{K^3}{162}+\frac{5K}{1764}\right).\tag{14}
\]

These weighted norms retain the allowed `1/k` state amplitude rather than
postulate a uniform infrared bound. They hold for every `t in [a,b]`.
Phase-sensitive integration could reduce them, but is not needed to prove
these numbers and has not been substituted for an absolute error bound.

The following decimal values are rounded upwards from exact rational
bounds for the cap and exact-Taylor remainder together:

| Source | K | Density upper | Pressure upper |
|---|---:|---:|---:|
| positive_B | 64 | 2.03007593423e-22 | 1.99456699467e-20 |
| positive_B | 128 | 8.55585772906e-22 | 1.69228183923e-19 |
| positive_B | 256 | 4.09335125407e-21 | 1.62463522704e-18 |
| signed_uB | 64 | 2.03031565410e-22 | 1.99480252150e-20 |
| signed_uB | 128 | 8.55681027232e-22 | 1.69247024485e-19 |
| signed_uB | 256 | 4.09373100414e-21 | 1.62478594837e-18 |

These values are source-representation contributions, not achieved errors
for Taylor coefficient arithmetic, oscillatory moments, saved trajectories,
momentum quadrature, subtraction contacts, direct work, or serialization.
They neither relax nor replace the inherited `2e-8` full-integral gate.

## 6. A concrete enclosure interface and its limits

A validated implementation can enclose the 1243 exact coefficients per
source using the same source recurrences and interval backend, while keeping
coefficient errors distinct from (10). On any real panel, coefficient radii
`eta_jm` contribute at most

\[
\sum_{m=0}^{112}\eta_{jm}\frac{2r_j^{m+1}}{m+1}
\]

to the source L1 error. This formula uses a polynomial in `s-s_c`; other
coefficient normalizations require their corresponding powers explicitly.
An interval coefficient enclosure may preserve a real polynomial center,
so the source-representation center retains exact `c=0`; subsequent
independent arithmetic in `U,W` still needs its own joint enclosure.

Oscillatory integration must be separately validated. A stable small-k
route evaluates the entire `Phi`, or power-series moments, without division
by an unenclosed tiny `k`. The eleven geometric source panels are **not**
automatically admissible phase panels: a large panel can have `2Kr>4`.
If reusing an implementation with phase cap four, subdivide each source
panel deterministically into widths at most `1/64` for `Kmax=256`, using
the exact source polynomial with validated recentering. Then every phase
halfwidth satisfies `2*Kmax*r_sub<=4`. Coefficient bounds, phase tails,
propagation residuals, arithmetic and their resource budgets need a frozen
implementation and independently checked controls before physical use.

At a retained node, let `U_s,W_s` denote the exact rational values obtained
from a correctly decoded binary input divided by the same represented
epsilon. If interval evaluation supplies `[U_*],[W_*]`, then

\[
c_{err}=\operatorname{Re}U_s-\frac{\operatorname{Im}W_s}{2k},\quad
[A_{err}]=\frac{W_s-[W_*]}{2ik}.\tag{15}
\]

The zero on the target side of `c_err` is justified by (3). It does not
permit changing `U_s` or `W_s`. The exact imaginary constant in
`U_s-U_*-A_err` remains invisible to these particular linear stress
operators; its absence from a stress bound is not a full complex-mode
accuracy claim.

Equation (15), once evaluated after the required freeze, would enclose
the retained **discrete** state discrepancy at those nodes. The original
arrays do not by themselves supply a continuum interpolant, regularity
bound, or quadrature error. A continuum reconstruction target supplied by
(1) is distinct from those finite stored values. Nine rational probes can
validate nine state intervals; they do not enclose the retained continuum
state, all nodes, all twelve cases, or a momentum integral.

## 7. Degeneracies, falsifying examples, and verification

The case `k=0` uses (2) and the entire limit, never `A=W/(2ik)`.
The physical mode `v0` still has its usual `k^-1/2` normalization, and no
isolated zero mode is inserted into the continuum measure. Zero forcing
gives zero response. A zero cap contributes zero; no lower bound on a
nonzero state or on its error is inferred. The finite-band bounds do not
control `k>256` or supply a UV tail, nonlinear perturbative error, arbitrary
mass/curvature shift, or a new physical initial-state model.

Concrete failure controls are available without physical evaluation:

* For fabricated real residual `r=1` on an interval of length `D>0`,
  `W(0)=-D` and `A~iD/(2k)`, falsifying a universal uniform-`A` claim
  while leaving the weighted finite-band stress integrable.
* A smooth normalized Bogoliubov direction supported between finitely many
  probed momenta has `c=0`, agrees at every probe, and changes the direct
  stress. Neither normalization nor probe agreement bounds its continuum
  amplitude. The existing 5 October counterexample applies unchanged.
* An incorrect but homogeneous incoming direction satisfies the same Ward
  equation. Conservation therefore does not replace (15).
* Replacing the physical prehistory by zero incoming data at `a` changes
  (1); a finite-time vacuum or an adiabatic subtraction frequency cannot
  justify that replacement.
* Extending the interior analytic disks across the flat support boundary
  invalidates Cauchy's premise. Dropping `-2Lh'`, the signed-source product
  derivative, a boundary phase, or a factor of two changes explicit checked
  identities.

`verify_incoming_representation.py` uses exact SymPy algebra and standard
library rational arithmetic. It verifies both integration-by-parts
identities, bump derivatives, signed-source factors, the removable kernels,
canonical cancellation, unchanged direct Ward identity, rational cap
premises, every panel's disk and coverage constraints, the pressure
majorant, and exact integrated coefficients. It rejects five algebraic
mutations. Its 54 checks use explicit exceptions and no Python assertions.
`PRIMARY_EXACT_RECEIPT.json` contains all panel parameters and exact rational
bounds; the displayed decimals use directed upward rounding. The source
never imports the physical producer or opens saved arrays.

The analytic inequalities in the prose are ordinary mathematical proofs;
the script is a reproducible algebra/rational check, not a proof-assistant
formalization or external peer review. The next finite task is to implement
and independently validate the coefficient/phase enclosure for this exact
representation, freeze it, and only then compare authorized physical
incoming data. No claim about actual stored accuracy follows until that
comparison exists.
