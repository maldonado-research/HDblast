# A conditional harmonic-covector value bound from the C1 error

This calculation quantifies the second-parameter harmonic constraint of the
same integral-defined finite candidate `W1` used in PR97 and PR98. It does not
construct a new candidate or differentiate a source magnitude bound.

Write `e=W1-U`, with the actual formal target `U` from the retained physical
subsidiary completion. The prior certificate gives

```
X(e)<E=3.541e25,
X(e)=max(sup||e||H/3, sup||e_x||H, sup||e_r||H).
```

The inherited common zero wedge implies the stronger continuous value bound
`||e||H<3E/4` on `D+`, and `||e||H<E/8` at the shell time `x=9/8`.
These support premises are separate from uniform magnitude bounds on sources.

For one fixed actual first field, the quadratic gauge source `Q_C` is the
same function for `W1` and `U`. Conditional on the retained exact formal
completion having `mathfrak_C^(2)[U]=DC[U]+Q_C=0`,

```
mathfrak_C^(2)[W1]=DC[W1-U]=DC[e].
```

This cancellation is a premise-based error transfer. The calculation does
not newly evaluate `Q_C`, establish its sign, or prove the full physical
subsidiary hypotheses from a numerical reduced-system contraction.

## Full metric-to-covector map and units

The metric uses the registered radially parallel orthonormal frame, with
radial/time entries `K,L,N` and angular entries `J,v,u,H,T+,T-`. The scalar
field `w` and the transverse-traceless slots do not enter `DC`. The radial
coordinate has orientation

```
partial_y=-(ell*rho)^(-1) partial_r,
frame time derivative=(ell*rho)^(-1) partial_x,
ell=1e-6.
```

The complete harmonic covector is
`DC[k]_A=nabla^B k_AB-(1/2)nabla_A tr(k)`. Let

```
A=ell*rho*a, T=ell*tanh(ell*x), Z=ell*sech(ell*x),
lambda=l(l+2), beta=2-2lambda/3.
```

For degrees `l=0,2,4`, the scaled radial and time coefficients are

```
(ell*rho*DC)_R = -(K_r+N_r-3J_r)/2-L_x
                 +A(4K+N-3J)-3T L-lambda Z v,
(ell*rho*DC)_T = -L_r-(K_x+N_x+3J_x)/2
                 +5A L-3T(N+J)-lambda Z u.
```

The `v,u` terms are absent at degree zero. For degrees two and four, the
angular-gradient coefficient is

```
(ell*rho*DC)_A = -v_r-u_x+5A v-4T u
                 +Z[(-K+N-J)/2+beta H].
```

Thus the five frame-covector components are represented by eight modal
coefficients `R:0,2,4; T:0,2,4; A:2,4`. With scalar harmonic norms
`S0=1,S2=1/12,S4=19/720`, the positive output Gram is `S_l` for radial/time
coefficients and `lambda*S_l` for angular-gradient coefficients. The latter
follows from integration by parts on the unit sphere. It has no extra
factor two: mixed metric components have that multiplicity in the input
metric norm, while a covector does not.

The resulting norm is a positive angular `L2(S3)` comparison norm of frame
components, uniform in `(x,r)`. It is not a pointwise-in-angle maximum or a
Lorentz-invariant scalar. The physical frame-covector value and the scaled
value differ by the point-dependent positive factor `ell*rho`.

## Exact derivative Gram and continuous coefficient bound

Write the scaled map as `Rx e_x+Rr e_r+R0 e`. The checker forms all three
matrices in the full 23-dimensional input space. It proves the exact combined
derivative Gram identity

```
Gout^(1/2) [Rx H^-1 Rx^T+Rr H^-1 Rr^T] Gout^(1/2)
 =diag(7/4,7/4,7/4,7/4,7/4,7/4,1,1).
```

Consequently the derivative part has norm less than `sqrt(7/2)E` when the
two input derivative norms are each below E. The complete actual background
is interval-enclosed on the containing rectangle `0<=x<=3/2,0<=r<=7/4`.
The hyperbolic factors use valid Taylor inequalities on the complete time
interval. Per-degree weighted Frobenius bounds enclose `R0`; the largest
block bound applies because different degrees are orthogonal. This yields

```
||ell*rho*DC[e]||+ < [sqrt(7/2)+(3/4)||R0||op] E
                  < 6.625e25,
||DC[e]||+ < 5.100e29
```

uniformly over `D+`. The frame bound divides by `ell` times an authenticated
strictly positive lower bound for the actual `rho`, approximately 130.

## Shell normal cancellation and tangential sharpening

The same fixed actual M and G are used by W1 and U. Their difference obeys

```
PN e_r=M PN e,   PD e=PD e_x=0.
```

The normal row has no remaining `L_x,v,L` terms at the shell. The checker
uses the exact primitive Robin matrices to verify, with shared symbols,

```
Rr_R By_a=(4K+N-3J) row,
Rr_R By_p=Rr_R By_sigma1=Rr_R By_sigma2=0.
```

Since `M=-ell*rho[a By_a+p By_p+sigma1 By_sigma1+sigma2 By_sigma2]`,
the scaled normal row `Rr_R M+A(4K+N-3J)` cancels exactly. Hence the
normal harmonic-constraint error is **exactly zero at the shell**. This
uses shared-symbol algebra; independent interval entries would obscure the
identity. In particular scalar boundary terms cancel rather than being
discarded.

For the remaining tangential rows, the derivative operators `Rx PN` and
`Rr PD` have norms at most `sqrt(5/4)` and `sqrt(1/2)`. The Dirichlet normal
derivatives `L_r,v_r` remain present; zero Dirichlet trace does not remove
them. The connection factor A only multiplies terms already canceled or
zero under the shell projectors. The remaining tangential zero-order matrix
uses actual T and Z. At `(9/8,0)`,

```
||ell*rho*DC[e]||+ < 6.463e25,
||DC[e]||+ < 4.975e29,
DC[e]_radial=0.
```

These are absolute upper certificates under the stated assumptions. Their
numerical sizes do not establish relative smallness, a nonzero lower bound,
or a physical observable.

## Authentication and limits

`check_constraint_transfer.py` authenticates the existing coefficient driver,
profile norm convention, PR97 error report and primary radial gauge formula
before use. The driver authenticates its local dependencies and seven complete
scientific roots before importing project code; all are read back afterward.
Normal and optimized runs give byte-identical rational reports. No frozen
input or manuscript is changed.

The result bounds harmonic-covector **values** only. It does not bound their
derivatives, physical Hamiltonian/momentum residuals or the full natural
Einstein residual. It supplies no quantitative C2 error certificate, exact
finite-amplitude family, observational test, or higher-dimensional Big Bang
origin. The full first-field/background/source, initial-completion and
subsidiary-theorem premises remain conditional.
