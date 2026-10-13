# Independent analytic review of the scalar-branch continuation

Date: 11 October 2026. This is an internal mathematical review, not external
peer review, an originality certificate, or evidence for a cosmological origin.
The final source hashes will be recorded only after reviewing the completed
proof and numerical certificate.

## Scope

I read the accepted 8 October uniform-branch proof and the previously unreviewed
derivative-norm investigation. I independently derived the identities below,
checked the generalized norm estimates, and checked origin and endpoint issues.
Only a small weighted regular-cone neighborhood is considered. No global branch
uniqueness, time-dependent stability, or observational conclusion follows.

## General derivative norm

With derivative scale `s`, weight `s_T`, and norm
`max(sup|eta|/s_T,sup|eta'|/(s s_T))`, the proposed constants are consistent:

* `C_G=max(1/D,(1+m^2/D)/(4ks))` follows by dividing the unchanged physical
  derivative estimate by `s`.
* `C_F=s^2/4+M2/12` follows from the kinetic term and
  `|U(eta)-U(0)| <= M2 eta^2/2`.
* The nonlinear scalar product contributes `4 C_F(s/k)B^3 b` to `C_N` and
  `20 C_F(s/k)B^2 b^2` to `L_N`. The latter factor is `4+16`: the variation of
  scalar momentum and variation of the metric logarithmic derivative.
* `C_D=s C_G C_N` and
  `C_P=s C_G B[2 L1 B+20 C_F(s/k)B^2 b]`; the endpoint momentum derivatives
  must retain this factor of `s`.
* `C_a=4 C_F B^2/k` retains the metric scale `k`.

The generalized contraction gate already forces strict convexity of the bulk
potential on the realized tube. In fact, `2L1=M3`, `C_G >= 1/D`, and
`q=C_G L_N <= 1/2` imply

    M3 B b <= D/2,
    U''(eta) >= m^2-M3 B b >= m^2-D/2 > 0,

for `|eta| <= Bb`. This observation requires no separate numerical convexity
assumption and uses `0<D<m^2`.

## Regular-cone central variation and positivity

Write `c=eta(0)`, `h=partial_c eta`, `u=partial_c a`, and
`v=partial_c log(rho)`, so `v'=u`. At the regular center,

    eta=c+U'(c)y^2/10+o(y^2),
    p=U'(c)y/5+o(y),
    rho=y-U(c)y^3/36+o(y^3),
    h=1+U''(c)y^2/10+o(y^2),
    h'=U''(c)y/5+o(y),
    u=-U'(c)y/18+o(y),
    v=-U'(c)y^2/36+o(y^2).

Strict convexity and `U'(0)=0` give `sign U'(c)=sign c` for nonzero `c` in the
tube. The radial scalar divergence equation ensures `p` keeps this sign and
`eta` moves away from zero. For `epsilon=sign c` put `z=-epsilon u`.
The Jacobi equations can be written

    (rho^4 h')'=rho^4[U''h+4|p|z],
    (rho^2 z)'=rho^2[|p|h'/2+|U'|h/6].

The center expansions give positive `h`, `h'`, and `z` immediately to the right
of the center; the cooperative system on every interval away from zero
preserves positivity. Equivalently its regular integral kernels are positive.
Consequently `h>=1`, `h'>0`, and `z>0` for nonzero `c`, and `v` has the opposite
sign to `p`. At `c=0`, the background is exact and `h=F_lambda>0`, `h'>0`, while
`u=v=p=0`.

The damping form above avoids a potential proof gap: the alternative equation
for `-u` involving `rho^-2(-v)` has a singular positive coefficient at the
center and cannot be treated as an ordinary zero-data IVP there without its
specified quadratic center coefficient.

## Curvature-factor identity and origin term

Let `H^2(y)=rho(y)^-2` and `j=u+ph/3`. Linearizing

    a^2=H^2+p^2/12-U/6

and using `a'=-H^2-p^2/3` gives exactly

    j'+4aj=-2H^2v,
    j(T)=-2rho(T)^-4 integral_0^T rho(y)^2 v(y)dy.

The integrating-factor boundary term vanishes because `j=O(y)` and
`rho^4=O(y^4)`. Thus `j` has the same strict sign as `p` for nonzero `c`.
The relation is a consequence of the constraint and Jacobi equations, not an
additional field equation.

## Moving endpoint and the second Wronskian

Define `d=h p'-p h'`. Direct differentiation gives

    d'+4ad=4p[H^2h+p j].

Although `d(0)=U'(c)/5` is generally nonzero, `rho^4 d` tends to zero at the
center. The integrated right side has strict sign `p`, since `h>0` and
`p j>=0`. Hence `d` has strict sign `p` for nonzero `c`.

Since `h(T)>0`, the central and Dirichlet endpoint parameters are nondegenerate:

    c_beta=1/h(T),
    p_beta=h'(T)/h(T),
    partial_T p|beta=d(T)/h(T).

The accepted center IVP and finite-interval Dirichlet construction establish
`C^1` dependence before these identities are used. No implicit-function theorem
is applied at `T=infinity`.

On the scalar-junction locus, its previously imposed positive derivative
`g_beta=p_beta+W''+delta f''/2>0` gives

    beta_T=-(d/h)/g_beta.

Thus `beta_T` has the opposite sign to `p`; a new assumption `sigma''>=0` is
unnecessary. For the metric residual `F=a-W/3-delta f/6`,

    F_T=-H^2+(a_beta+p/3)(beta_T-p).

Here `a_beta+p/3=j/h(T)`, which has the same sign as `p`. Therefore

    F_T < -H^2  if c!=0,
    F_T = -H^2  if c=0.

This holds throughout the scalar-junction locus, not only at its metric roots.
Together with the existing opposite endpoint signs it proves a unique shell
inside the same tube and removes the old sufficient downward-crossing gate.

## Uniform linear Riccati comparison

For the regular linear logarithmic derivative `z_lambda=F_lambda'/F_lambda`,
put `e=lambda-z_lambda`, `d0=a0-k>0`, and `E=e/d0`. (This `d0` is unrelated to
the Wronskian `d` above.) Direct calculation gives

    E'=4lambda-2(w-3k)E+d0 E(E-3).

The regular center has `E~lambda y`, so `E` starts at zero. With
`C=2lambda/(w-3k)`, the assumption `w>4k` gives `0<C<2`, and at any putative
first upper contact with `E=C` one has `E'=d0 C(C-3)<0`. Consequently

    0 <= lambda-z_lambda
      <= [lambda/(k(w-3k))] H0^2,
    H0^2=a0^2-k^2.

The bound is valid for every positive endpoint, not only `T>=1/k`.

Automatic bulk convexity also gives `U(eta)-U(0)>=0`, hence the metric source
`F_metric=p^2/4+[U(eta)+6k^2]/6` is nonnegative. The existing positive Volterra
solution then has `R'=-(integral partial_y K F_metric R)<=0` and `R<=1`.
It follows that `H0^2=R^2 H^2<=H^2`; the nonlinear endpoint coefficient can
therefore be taken as

    C_H=lambda/[k(w-3k)].

If this positivity lemma is omitted, the conservative previous bound
`R<=4/3` instead gives the coefficient with a factor `16/9`.

## Remaining review actions

Before calling the continuation checked, read the completed statement and
proof, independently reproduce its exact constants and inequalities without
importing its producer, and record immutable source hashes. Original approved
manuscript and archive files are outside this review's edit scope.
