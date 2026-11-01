# Independent review of the weighted second initial monopole

1 November 2026. **ACCEPT, conditional on the retained analytic inputs.**
No material algebraic, normalization, support, or sign error was found in
the reviewed note. This is internal AI-assisted mathematical review, not
external specialist peer review.

The reviewed mathematical derivation is the accompanying
[weighted initial-data note](../WEIGHTED_MEAN_INITIAL_DATA.md).
The original independently reviewed note had SHA-256
`a2f9d76543aa7d16826b2331d780f70873f3ce18570e38ddb9ebc4c7596d93cf`;
this hash is review history, not a runtime dependency. Its six cited source
files were read directly and authenticated against SHA-256 pins and Git
objects at source commit `b51908d7e62e7ec80cef82438edd00923ff94723`.
The portable replay authenticates those repository-relative source bytes
using SHA-256 and byte counts, before importing SymPy and again after the
checks; it requires no Git metadata. See [source pins](SOURCE_PINS.json).
The portable note's SHA-256 is
`5958ad37ef2aa20a2b58e8554f85cfc946504c59d1caa27f7c58b9f814766514`.
Its diff was reviewed: only introductory packaging language and six public
source links changed. The scientific acceptance carries to that version;
see the [review receipt](REVIEW_RECEIPT.json).

## Findings

1. **Canonical normalization and the second derivative are correct.**
   S1 specifies action `R/2-|dphi|^2/2-V` and
   `pi=rho_gamma(K^sharp-gamma^-1 tr K)/2`. Its Hamiltonian is exactly half
   S3's physical constraint. The family convention is
   `X_e=X0+e X1+e^2 X2/2`. Consequently the second equation is
   `P X2=-D^2 C[X1,X1]`, with neither an omitted nor an additional factor
   of two. Direct expansion of the canonical kinetic density, independently
   of S3's second solution, gives the mean kinetic contribution
   `-s^2-j^2/6`.

2. **Both quadratic sources agree with a direct calculation.**
   The review independently expanded scalar curvature using
   `L=sqrt(1-4e uY)` and angular metric
   `rho^2(1+2e uY) kappa`, including angular lapse derivatives, and then
   expanded the scalar gradient, potential, and canonical kinetic terms.
   This reproduces the displayed `QH`. The probability sphere moments are
   `E[Y^2]=1/6` and `E[|DY|^2]=4/3` for the specified quadrupole.
   The exact momentum-density formula produces
   `QM=-4u s'-12a u s-2s u'+j v'/3` and angular source
   `(6us+jv) D(Y^2)`. The momentum first constraint gives the stated
   reduced radial expression.

3. **The radial weighted problem has the correct adjoint and pairing.**
   An independent linearization of the radial warped metric gives the
   stated `Pg(A,B)`. Formal integration against `rho^3 dy` gives `Jr` and
   `Jt`, with the essential tangential multiplicity three. Hence
   `A=WH Jr`, `B=WH Jt` and the Hamiltonian weak form are correct.
   For radial one-form `xi dy`, `Pm*` has components `(2xi',2a xi)`;
   its divergence gives the displayed momentum equation and weak form.
   The weight exponents `WH=d^(4alpha_w) exp(-2beta/d)` and
   `WM=d^(4alpha_w-4) exp(-2beta/d)` follow from the retained theorem.
   The physical zero-extension statement supplies the boundary control
   needed here. Potential boundary behavior must still be treated in the
   dotted weighted spaces, as the note explicitly requires.

4. **The canonical-to-physical conversion is correct.**
   Direct differentiation reproduces all displayed `K2`, `Pi2`, and
   coordinate metric velocity formulas. In particular
   `Pi2=-2ujY^2` although `varpi2=0`. The first lapse derivative cancels
   this contribution in `partial_t phi2`, so both `phi2` and its coordinate
   time derivative vanish. This cancellation uses the specified lapse and
   would not justify an arbitrary independent lapse completion.

5. **The areal-radius moment and its sign are correct.**
   Substitution of `N=rho` yields
   `Pg*(rho)=diag(-3/(2rho),-1/rho,-1/rho,-1/rho)`.
   Pairing with `(A,B,B,B)` and using `Pg(A,B)=-QH` gives exactly
   `integral rho^2(A+2B)dy=(2/3) integral rho^4 QH dy`.
   This is consistent with S4's earlier areal-radius identity after its
   coordinate and source conventions are converted.

6. **The compact metric-only kernel is valid.**
   Independent substitution gives
   `Pg(f gamma0)=-3(f''+3af')/2-Vf` and
   `Pg(Lie_(xi d_y)gamma0)=rho^-4(rho^4 p^2 xi)'`.
   Twice integrating the first formula gives
   `integral rho^4 C[f]dy=-(9/2) integral rho^2 fdy`.
   The stated zero moment condition therefore makes `xi_f` vanish beyond
   both ends of the compact support. Since `p` is nonzero throughout the
   correction annulus, division by `rho^4 p^2` is legitimate. The map
   `f -> h_f` is injective on compactly supported profiles: if `h_f=0`,
   its radial and angular components imply `xi'=a xi`, and compact support
   forces `xi=f=0`. This gives infinite-dimensional homogeneous metric
   data preserving the affine scalar sectors and the fixed moment. The
   alternative exact-family argument follows from the same conditional IFT
   because its new second source is `P h_f=0`.

7. **Weight dependence is stated at the appropriate strength.**
   The finite difference between any two admissible selected corrections
   lies in `ker P`. On a common differentiable inverse domain,
   `delta h=(I-R_W P)[(delta log W)h]` follows from differentiating
   `R_W=W P*(P W P*)^-1`. Constant multiplication of either block weight
   cancels. The note properly does not claim that every change of weight
   changes this particular correction. Source S1 does not instantiate the
   defining function or weight parameters, so these formulas do not supply
   a unique numerical radial profile.

## Independent sign proof and orientation audit

All primes are outward proper-radial derivatives. S1 and S3 explicitly give
`dy=-rho dz`, `ell=10^-6`, and first-packet support `ell<=z<=2ell`.
The increasing-y support interval runs from the point `z=2ell` to `z=ell`.
Thus its positive length is
`L=integral_ell^(2ell) rho(z) dz <= 131ell`; there is no negative measure
or sign reversal in the moment integral. Equivalently
`integral_y rho^4 QH dy = integral_ell^(2ell) rho^5 QH dz` on the support.

The independent calculation proves the pointwise identity

```
rho^4 QH = rho^4 [(4/rho^2+8p^2/3)u^2-s^2-j^2/6
                     -(v'+4pu)^2/6-(V''/6+4/(3rho^2))v^2]
            +[-2rho^4 u u'-8rho^4 a u^2]'.
```

The flux vanishes at both support endpoints. For `F=rho^2u`, the exact
first scalar reconstruction gives `F'=-rho^2pv/3`. Since `F` vanishes at
the lower endpoint, Cauchy-Schwarz gives

```
integral F^2 dy <= (L^2/2) integral F'^2 dy
                 <= (Lmax^2 pmax^2/18) integral rho^4 v^2 dy.
```

This estimate does not require a numerical solution or a chosen weight.
S6's certificate covers the entire packet, including `rho in [129,131]`,
`|p|<=7/50000`, and `|phi0-1|<=9/100000`. The exact potential polynomial
gives the rational lower bound

```
V'' >= 139954999270708393439/45000000000000000000 > 3.
```

The upper coefficient of the positive contribution is exactly bounded by

```
(4/129^2+8(7/50000)^2/3) (131/1000000)^2 (7/50000)^2/18
 = 1051339806152867/234014062500000000000000000000000000
 < 1/100.
```

Its decimal value, for orientation only, is `4.4926351644e-21`.
Dropping the other nonpositive terms proves

```
I <= -(49/100) integral rho^4 v^2 dy < 0,
integral rho^2(A+2B)dy <= -(49/150) integral rho^4 v^2 dy < 0.
```

Strictness is justified: the actual `u=alpha Q/rho^2` is nonzero because
the bump is nonzero and the retained branch has nonzero `p`. If `v` were
identically zero, the primitive relation and compact support would force
`F=u=0`. Therefore its weighted squared norm is positive.

## Validation and limits

[`checks/independent_review.py`](checks/independent_review.py) was written before inspecting the producer's
checker. It independently derives the geometric and canonical expansions,
uses exact SymPy arithmetic, and authenticates the cited source bytes.
All 29 checks pass under both `python` and `python -O`; rejection controls
detect the reversed total-derivative sign and a factor-two Hessian error.
Fresh portable replay results are [normal](results/NORMAL.json) and
[optimized](results/OPTIMIZED.json). The mathematical checks and exact
rational sign constants are unchanged by the portability edits.

Acceptance is for the derivation and the stated conditional inference.
This review does not re-prove the retained analytic weighted inverse,
background branch enclosure, nonlinear evolution theorem, or a singular
weighted boundary-value solver. It supplies no numerical radial profile,
postcontact shell-volume sign, energy interpretation, hot-matter model,
or new gauge-invariant physical observable. The initial moment depends on
the fixed first-order reconstruction gauge, exactly as the note states.
