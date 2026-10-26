# Exact sphere curl in the frozen tensor basis

This is a finite angular calculation for the registered-coordinate screen
magnetic diagnostic. It supplies neither the quadratic spacetime curvature
map nor a geometric observer prescription. The exact results below pass the
adjacent checker in ordinary and optimized Python, including the independent
polynomial calculation and deliberate fault controls.

## Frozen conventions and orientation

Use ambient indices `0,1,2,3`, unit vector `n`, metric `delta`, and
`P_ij=delta_ij-n_i n_j`. These are the source's sphere indices `1,2,3,4`
shifted down by one. The frozen ambient orientation is `EPS_0123=+1`.
The sphere orientation used here is its outward induced orientation:

```
epsilon_ijk = n^A EPS_Aijk.
```

Thus at `n=e0` the positively oriented tangent frame is `(e1,e2,e3)`.
This explicit choice is the convention needed when tying the calculation to
the spacetime orientation and ordered observer/radial pair. If those data
induce the opposite screen orientation, all curl and magnetic amplitudes in
this note acquire the same minus sign; norms do not change.

The frozen Weyl tensor has sectional entries

```
C_0101=C_2323=-2/3,
C_0202=C_0303=C_1212=C_1313=1/3,
```

with other independent entries zero. Set

```
(*C)_ij kl = (1/2) EPS_ij ab C_ab kl,
W_s=(C+s star(C))/2, T_s,ij=W_s,ia jb n_a n_b,       s=+1,-1.
```

This is exactly the `make_tt()` definition in the byte-pinned frozen
`inputs/angular_exact.py`.
No identification between the signs in these names and a curl eigenvalue
is assumed.

The source orientation and Weyl conventions are explicitly documented in
the prior repository's
`research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/angular/source-audit/SOURCE_INVENTORY.md`.
The frozen implementation copied here has SHA-256
`2c0bafec777581c79c3fd6bd9879d0d3cb8114bb899722bd33882086751ded60`.
The copied input is unchanged; no frozen repository file is edited.

## Curl calculation and complete selected-column map

With parentheses denoting normalized symmetrization, define

```
(curl t)_ij = epsilon_(i^{cd} D_c t_{j)d}.
```

The exact result in the frozen basis is

```
curl Tplus =  3 Tplus,
curl Tminus = -3 Tminus.
```

An analytic way to see the sign and factor starts at `n=e0`. The tangent
derivative is

```
D_c T_jd = W_jc d0 + W_j0 dc.
```

For `*W_s=s W_s`, contraction with `epsilon_i^{cd}` sends the first term to
`s T_ij` and the second term to `2s T_ij`, using the Weyl trace identity.
Both are already symmetric. The construction is covariant under oriented
ambient rotations, so the identity holds at every sphere point. The
polynomial checker below independently verifies every component everywhere.

For the full selected angular tensor order

```
kappa, Zkappa, Fkappa, HessTF Z, HessTF F, Tplus, Tminus
```

the curl matrix is exactly `diag(0,0,0,0,0,3,-3)`. The five zero columns are
identities of complete tensors, not merely zero TT inner products. Output
symmetry, tangency, trace, and divergence are also checked exactly.

## A quadratic angular identity needed for the complete map

Let `Y=n0^2-n1^2`, as in the registered first field. Then the complete tensor

```
A_ab = epsilon_(a^{cd} (D_cY) (D_b)D_dY)
```

satisfies

```
A = 2(Tplus-Tminus),       |A|^2 = 8S(1-S),
mean(|A|^2)=4/3,           S=n0^2+n1^2.
```

The equality is checked by full exact reconstruction and a separate raw
ambient gradient/Hessian calculation. In particular, its two TT
coefficients are `(2,-2)`, with no residual scalar, vector, or other TT
component. This identity alone says nothing about its spacetime coefficient.

## Norms, cancellation, and projection normalization

For the positive sphere Frobenius pairing and probability-normalized measure
`mean(f)=integral_S3 f/(2*pi^2)`, the exact identities are

```
|Tplus|^2=|Tminus|^2=1/6,
Tplus:Tminus=1/6-S(1-S),
mean(Tplus:Tminus)=0.
```

Consequently the two selected TT amplitudes of any square-integrable screen
tensor `B` are

```
bplus =6 mean(B:Tplus),
bminus=6 mean(B:Tminus).
```

With unnormalized sphere volume the multiplier is instead `3/pi^2`.
This pair projection is the full TT projection only after a separate proof
that no additional TT harmonic is present. Arbitrary tensors can contain
additional harmonics; the full real degree-two TT space has dimension ten.

For independent real amplitudes `a,b`,

```
curl(a Tplus+b Tminus)=3a Tplus-3b Tminus,
mean(|curl|^2)=(3/2)(a^2+b^2),
|curl|^2=(3/2)(a-b)^2+18ab S(1-S).
```

Equal amplitudes therefore do not cancel as a tensor:

```
curl(Tplus+Tminus)=3(Tplus-Tminus),
|curl(Tplus+Tminus)|^2=18S(1-S),
mean(|curl(Tplus+Tminus)|^2)=3.
```

The equal-combination curl has zeros at `S=0` and `S=1`, two coordinate
great circles. Its global angular norm is nonzero; it is not pointwise
nonzero everywhere. The orthonormal pair is `sqrt(6) Tplus,sqrt(6) Tminus`.
In the orthonormal sum/difference pair
`sqrt(3)(Tplus+Tminus),sqrt(3)(Tplus-Tminus)`, curl exchanges the two vectors
with factor three.

For the specified linear TT contribution to screen magnetic Weyl,

```
L_H[k]_TT = [1/(2 ell rho^2 cosh t)] curl(uplus_x Tplus+uminus_x Tminus),
```

the output amplitudes are
`(3uplus_x,-3uminus_x)/(2 ell rho^2 cosh t)`. For equal `ux`, its squared
normalized angular norm is `3ux^2/(4ell^2 rho^4 cosh(t)^2)`. The complete
second variation also includes the quadratic first-field/observer/screen
terms. The linear norm does not establish the full coefficient's sign or
nonvanishing.

## Finite all-point certificate and independent controls

`check_curl.py` uses exact standard-library `Fraction` arithmetic. Its first
route applies the frozen intrinsic derivative and reduces polynomial
coefficients modulo `n3^2=1-n0^2-n1^2-n2^2`. Every reconstruction residual
is the zero polynomial in this quotient. There is no point evaluation or
numerical quadrature.

The independent route has its own polynomial ring with no sphere reduction.
It reconstructs the Weyl seed from its diagonal bivector action, builds the
two frozen tensors, and takes the derivative by projecting all three slots
of the ambient ordinary derivative:

```
D_c t_jd=P_c^a P_j^b P_d^e partial_a t_be.
```

Each residual has a single degree parity. Its monomials are homogenized to
one degree by multiplication by powers of `sum n_i^2`. Exact cancellation
of all resulting coefficients proves the original residual zero on the
sphere, including chart boundaries. The independent raw curl degrees are
bounded by eight for these seven columns. The independent sphere-moment
formula is

```
mean(prod n_i^(2k_i)) = prod((2k_i)!/(4^k_i k_i!))/(sum k_i+1)!,
```

and odd monomials have zero mean. This separately checks the Gram matrix,
pointwise norm identities, and the integrated norm. The two angular
constructions share no derivative or polynomial reduction routine.

Replay from this directory, using fresh report destinations in an existing
external scratch directory:

```
python -I -S -B check_curl.py --output /tmp/CURL_NORMAL.json
python -I -S -B -O check_curl.py --output /tmp/CURL_OPTIMIZED.json
python -I -S -B run_controls.py --output /tmp/CURL_CONTROLS.json
```

Normal and optimized positive reports agree exactly. Seven deliberately
incorrect variants fail at their stated gates in both modes: reversed
orientation, a lost symmetrization factor, omitted minus column, false
equal-combination cancellation, incorrect projection normalization, reversed
independent epsilon, and an incorrect quadratic identity factor. The source
pin is checked before calculation, and bytecode writes are disabled before
imports. Both CLIs validate report destinations before importing project
inputs or doing mathematical work. They reject existing files, every symlink
path component, and destinations within their angular source directory or
an enclosing `HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026` package. Report
writes use exclusive creation. Nested control replays use isolated Python
with site initialization disabled (`-I -S -B`, plus `-O` when requested).
The APIs `curl_tensor`, `project_tt_pair`, and
`gradient_hessian_curl` return frozen ambient polynomial tensors or exact
amplitudes for use by the complete curvature calculation.
