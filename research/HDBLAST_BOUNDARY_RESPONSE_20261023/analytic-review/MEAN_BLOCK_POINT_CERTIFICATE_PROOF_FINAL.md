# Conditional positive mean at the shell-reaching point

23 October 2026. This separate candidate certificate sharpens the feedback
allowance for the **same** actual source integral at `(x,r)=(3/2,1/4)`.
It preserves the earlier full-state certificate unchanged. Conditional on
the retained smooth solution/operator/source premises and independent
acceptance of the pinned n16 source integral, exact arithmetic yields

```
mean W_yy(3/2,1/4)
    in [55843882.796267, 737593144520.766258].
```

These displayed decimal endpoints are rounded outward. The strictly positive
lower endpoint follows from the exact fractions in
`MEAN_BLOCK_CERTIFICATE_NORMAL.json`, not from rounded decimal subtraction.
The normal and optimized checker reports are byte-identical. The result is
one component mean of the auxiliary forced increment W; it is not the full
physical response U, a complete waveform, or an invariant force.

## Exact invariant mean block

Let P0 be SO(4) averaging of tensor sections, with probability-normalized
sphere measure. Its image in the retained second-response space consists
of the five invariant sections with amplitudes

```
z=(K_0,L_0,N_0,J_0,w_0),
```

where K is yy, L is y0, N is 00, and J is the coefficient of the sphere
metric. There is no invariant sphere vector or trace-free symmetric tensor.
The background coefficients and linear bulk/boundary operators are
SO(4)-equivariant. Thus P0 commutes with A, B, D0, M, both boundary
projectors, reflection, K, and L. The mean block is exactly closed, and
no complementary angular mode can linearly feed it. The retained matrices
independently verify this statement: every one of the 25 primitive matrices
has identically zero entries in both directions between these five indices
and the other eighteen indices.

This is a projection of the exact full response, not an arbitrary deletion
of amplitudes. Nonzero first-order harmonics still contribute to the actual
quadratic mean sources before this projection. Connection, curvature,
metric/scalar and Robin couplings within these five means remain present.

The positive H Gram matrix on this block is `diag(1,2,1,3,1)`. Only L_0
is Dirichlet; the other four means are Robin. The mean scalar basis has
means `(1,0,0)`. Inner products of the seven angular tensor basis elements
with the sphere metric are `(3,0,0,0,0,0,0)`. These facts identify J_0 as
one-third of the averaged sphere trace, and K_0 exactly as the yy mean.
They are checked by exact probability-sphere polynomial integration.

In this block the primitive matrices Delta, Ca, Ch and Ephi vanish.
No angular derivative or angular connection cross term remains. All other
primitive matrices and their coefficients are retained below.

## Complete entrywise operator majorants

All matrices in this section are the five-by-five restrictions of the
pinned `linear_operator.py` primitives. Write `|H|` for entrywise absolute
value, distinct from an operator norm. Let

```
ell=10^-6, R=131, alpha=14/125, P=7/50000,
Ht*=3/10^6, h*=Ht*/129,
b*=alpha^2+1/10000, d*=129^-2+alpha^2,
V*=1, V1*=1/1000, V2*=4, py*=1/500,
sigma1*=2P, sigma2*=4.
```

These are bounds from section 3 of the retained construction. In particular
`b=-(a_y+a^2)`, `d=rho^-2-a^2`, and `|p_y+4ap|<=py*+4alpha P`.
The actual first tension has `sigma'=-2p`, `sigma''=4phi0`; the accepted
background has `|phi0|<1`. Time/space-dependent coefficients therefore obey
the following entrywise upper matrices throughout the entire cone:

```
Abar = 3ell Ht* |I| + 2ell R alpha |O| + ell R P |Aphi|,
Bbar = 3ell R alpha |I| + ell R P |Bphi|,

Cbar = alpha^2 |Sa-O2| + alpha h* |Sah-3O| + (h*)^2 |Sh|
       +2b* |Rb| +2d* |Rd| +(4/3)V* |Pmet|
       +(4/3)V1* |Pg| +V2* |Pww|
       +(py*+4alpha P)|W_K| +3P h* |W_L|,
Dbar = ell^2 R^2 Cbar,

Mbar = ell R [alpha |By_a| +P |By_p|
              +2P |By_sigma1| +4 |By_sigma2|].
```

These expressions follow term by term from the full retained `assemble`
routine. The vanishing angular primitives account for every omitted angular
term. Terms with a common scalar coefficient are combined first, such as
`a^2(Sa-O2)`, so their exact internal cancellations are retained. The
curvature and potential terms use triangle inequalities. The matrix
inequalities `|A_ij|<=Abar_ij`, etc., hold at every point; sampled background
values are not used. The coordinate `2a` tangential boundary diagonal is
absent because the retained operator is in the parallel frame.

## Actual mean-source caps

The pinned `ACTUAL_SOURCES_NORMAL.json` supplies continuous whole-collar
mean component caps directly from its `absolute_coefficient_caps`. In the
order K,L,N,J,w the following integers round those exact caps upward:

```
F <= (5290948084880, 5290851762796, 5290755442826, 362, 577380).
```

The checker uses the original exact rational caps, not these rounded numbers.
These are actual-source enclosures including the inherited actual first
solution uncertainty. They are not comparator sources or source samples.
Their normalizations already include `-2ell^2rho^2` for metric forcing and
`+ell^2rho^2` for scalar forcing; no further factorial factor is applied.

The same report contains six complete continuous shell intervals covering
`[0,3]`. Take the maximum absolute endpoint for the constant scalar slot of
the one-copy outward QK00, QKab and Qp over those intervals, denoting the
three caps by q00, qab, qp. Let qcy be the global actual `QC/y/0` cap.
The complete boundary-source caps are

```
G_K=2ell R(q00+3qab+qcy),  G_L=0,
G_N=2ell R q00,           G_J=2ell R qab,
G_w=ell R qp.
```

Their upward-rounded integer values are
`(3607805,0,901952,901952,1353037)`. The exact source entries, intervals,
normalization, provenance hashes and arithmetic are in the checker report.
This uses the Lorentz trace `-QK00+3QKab`, the full varying-metric-lowered
QC_y, and the outward normal conversion `partial_y=-(ell rho)^(-1)partial_r`.
The absolute-value upper bound does not alter these signs in the evaluated
source integral. The already evaluated source integral retains negative L.

## Positive matrix comparison for all five coupled means

Use the same backward-closed physical cone D+ and exact zero wedge
`x+r<1` as in the original audit. For each scalar amplitude define

```
X_i=max(sup_D+ |z_i|/3, sup_D+ |partial_x z_i|,
        sup_D+ |partial_r z_i|).
```

These are scalar-component norms, not H norms, so they require no H-norm
conversion. Symmetric tensor entries are handled by the exact tensor
actions. In particular the Gram factor two for L is not incorporated into
the individual X_i; no Gram weight is removed from an H-norm argument. A positive
weighted maximum of the X_i is an equivalent comparison norm on this
finite-dimensional block.

The exact kernel estimates are `||K||<=a=3/2`, `||L||<=1` in these
component norms. Reflection is diagonal with signs of modulus one. Put

```
Sigma=Abar+Bbar+3Dbar,
C=a Sigma+3Mbar,
s=a F+G.
```

The actual C1 solution satisfies componentwise `X<=s+C X`. Every entry
of C is nonnegative and its exact largest row sum is less than `0.000853`.
Hence `(I-C)^(-1)=sum_n C^n` is entrywise nonnegative and converges. Iterating
the inequality and letting the remainder `C^n X` tend to zero gives

```
X <= b=(I-C)^(-1)s.
```

The checker solves this five-dimensional rational system exactly, checks
`b>=0` and `b=s+C b`, and verifies the contraction. For orientation, upward
integer bounds on b are

```
(7939746443463, 7937501520452, 7937707300938,
 350447890, 2878016).
```

No approximate inverse, floating eigenvalue, or numerical feedback waveform
is used. Equivalently, in the positive component weights b_i the induced
weighted maximum norm of C is less than `0.997425`.

## Point feedback allowance and signed join

Because the output K_0 is Robin, the folded kernel weights are one on D-
and one-half on D+ minus D-. The retained zero wedge removes one triangle
of area one-half from each unfused cone, leaving total kernel weight `5/8`.
Thus the full bulk feedback into K_0 has absolute value at most

```
E_bulk=(5/8)(Sigma b)_K.
```

The actual solution is C1 and z_i(1,0)=0. For every `1<=s<=d=5/4`,
`|z_i(s,0)|<=integral_1^s |partial_x z_i| <= (s-1)b_i`.
The entire segment lies in D+. Therefore the complete Robin feedback has
bound

```
E_Robin=(d-1)^2(Mbar b)_K/2=(Mbar b)_K/32.
```

These estimates need no source derivatives and no C1 numerical
reconstruction; they use the accepted actual C1 solution. They retain
all five internal couplings and every possible contribution from the full
state to the mean, since the complementary block has exactly zero coupling.
Their exact sum is enclosed by the following outward decimal brackets:

| Allowance | Lower decimal | Upper decimal |
|---|---:|---:|
| E_bulk | 509975007.954096 | 509975007.954097 |
| E_Robin | 21841336.658811 | 21841336.658812 |
| E_5 | 531816344.612907 | 531816344.612908 |

The n16 source report, SHA256
`6caf16c304fb6db940899d58d09b648066eb2356385189ccb79eb2b0a6879675`,
supplies the unchanged exact interval `[Slo,Shi]` for the complete mean
`S=KF+LG`. Its outward display is

```
S in [587660227.409174, 737061328176.153350].
```

The new checker takes precisely those pinned fractions and computes
`[Slo-E_5,Shi+E_5]`. The strict lower margin is between
`55843882.796267` and `55843882.796268`. This is a signed conclusion
about the actual auxiliary mean, conditional on the accepted inputs and
source-integral review. It does not assert that the mean feedback itself
has been evaluated. No extra factor three or angular evaluation factor
belongs here because the bounded output is already the scalar mean K_0.

## Reproduction and scope

Run `check_mean_block_certificate.py` with `python -I -S -B`, adding `-O`
for the optimized check. It accepts `--source-root`, `--compute-root` and
`--output`; output pins use stable relative role labels. It authenticates
the exact operator, angular basis, actual source report and retained
construction before import, and authenticates the normal/optimized n16
source reports before the signed join. It rereads mathematical input pins
after calculation. All acceptance decisions use Fraction arithmetic and
explicit checks that remain active under `-O`.

The checker does not independently rerun the expensive inherited PDE/source
proofs or the n16 integration. Those remain separately certified premises.
The original full-state allowance and its zero-containing interval remain
valid. This separate sharper certificate follows from exact rotational
block closure and anisotropic source bounds. W still need not separately
satisfy the full sourced physical constraints; total U retains its registered
initial data and existing conditional analytic physical completion. No
claim of a global waveform, invariant backreaction, nonlinear family,
finite-amplitude remainder, or external novelty is added.
