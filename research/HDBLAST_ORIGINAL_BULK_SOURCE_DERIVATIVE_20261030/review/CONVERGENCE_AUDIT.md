# Independent audit of the physical-side gauge C4 convergence argument

This is an internal analytic review of the retained linear first-gauge
problem. It supplies no numerical coefficient ledger, C4 approximation
error, differentiated bulk-source bound, nonlinear existence theorem or
cosmological conclusion by itself.

**Assessment: the proposed time bootstrap is valid with the explicit
stabilization offsets and regularity premises below.** Quantitative bounds
on every mixed third coefficient derivative are not needed.

## Regularity needed for the bootstrap

Write `U_n` for the Picard iterates of
`U_xx-U_rr=P U_x+Q U_r+R U+s`, with zero Cauchy data and the retained
time-independent D/N rows. The existing proof gives convergence of
`U_n`, `V_n=U_n,x` and `W_n=U_n,xx` in the weighted C1 norm and a physical-side
C3 limit. The new step requires `P,Q,R,s` to be C2 jointly on the physical
triangle and their pure third time derivatives to exist continuously there.
The explicit smooth static coefficients and temporal factors meet the
coefficient condition. The actual master C4 result makes `s_xxx` continuous,
because the source contains one master derivative. Magnitude boxes alone
would not establish these qualitative statements.

For a prescribed forcing `g` with these derivatives, three time derivatives
of the retarded physical D/N wave give

```
partial_x^3 K g = C(g_x(0))
                 + S(g_rr(0)+g_xx(0)) + K(g_xxx),
```

where `C` is the free displacement wave and `S` the free velocity wave.
Only `g_xr(0)` and `g_rr(0)`, already available from C2 regularity, enter the
free C1 estimate. This proves that the new time derivative of each finite
iterate is C1. The remaining radial fourth derivatives then follow from
the PDE, so full mixed C3 regularity of the reflected coefficients is not
an intermediate requirement. The identically zero initial shell collar
makes every required free D/N profile compatible and removes corner delta
or characteristic-join defects. It is stronger than a claim that a few
independent endpoint values happen to vanish.

## Initial stabilization must be retained

Put `J_n=partial_x^3 U_n` and

```
d=P0*s0+s_x0,
v=s_rr0+P0*(P0*s0+s_x0)+2P_x0*s0+Q0*s_r0+R0*s0+s_xx0.
```

Direct differentiation of the *finite-iterate* PDE gives:

| Iterate | Initial displacement of J_n | Initial velocity of J_n |
| --- | --- | --- |
| J_0 | 0 | 0 |
| J_1 | s_x0 | s_rr0+s_xx0 |
| J_2 | d | v-P0*P0*s0 |
| J_n, n>=3 | d | v |

In particular, the first two nonzero differentiated iterates do not have
the final common Cauchy data. The matrix product `P0*P0*s0` and the factor
two multiplying `P_x0*s0` are essential. The initial displacement derivative
is `d_r=P_r0*s0+P0*s_r0+s_xr0`. All of these profiles vanish in an open
initial shell collar, so their required odd/even free extensions are C1.

Let `L=K(P partial_x+Q partial_r+R)` and let its C1 norm be bounded by
`q<1`. For `n>=2`, the recurrence for `J_(n+1)` has the common free wave
with data `(d,v)` and the forcing

```
L J_n + K[3(P_x W_n,x+Q_x W_n,r+R_x W_n)
          +3(P_xx V_n,x+Q_xx V_n,r+R_xx V_n)
          +P_xxx U_n,x+Q_xxx U_n,r+R_xxx U_n+s_xxx].
```

Consequently, for the successive-difference convention
`delta J_n=N(J_(n+1)-J_n)`, the common-data difference recurrence begins
at `n>=3`:

```
delta J_n <= q delta J_(n-1)
 + X[3 alpha1 delta W_(n-1)
     +3 alpha2 delta V_(n-1)+alpha3 delta U_(n-1)].
```

The earlier three differences are finite exceptions. The already proved
lower differences have geometric tails times polynomials of degrees at most
0, 1 and 2. The displayed convolution therefore bounds the new differences
by a summable geometric tail times a polynomial of degree at most 3. Thus
`J_n` converges in C1, rather than merely being uniformly bounded. Derivative
identification from uniform C1 convergence yields `J=U_xxx`; hence
`U_xxxx` and `U_xxxr` are continuous on the physical triangle.

## Radial recovery and the shell

The physical forcing `T=P U_x+Q U_r+R U+s` is C2 because the existing limit
is C3 and the coefficients/source are C2. The identities

```
U_xxrr=U_xxxx-T_xx,
U_xrrr=U_xxxr-T_xr,
U_rrrr=U_xxxx-T_xx-T_rr
```

hold in distributions and have continuous right sides. Integration along
physical coordinate lines identifies the remaining classical mixed fourth
derivatives, with their one-sided traces at the shell/domain edges.
The limiting J retains its D/N rows: `J_D=0`, `J_N,r=0`.

An odd output forcing may still have a jump on an artificial negative-depth
extension. No C4 extension across the doubled shell is proved or needed.
The argument concerns the physical half-line and retains the reflection
trace terms through the retarded operator. It does not assume that positive
time Dirichlet forcing or Neumann forcing gradients vanish at the shell.

Finally the actual limit satisfies

```
N(U_xxx) <= [Nfree+X(3 alpha1 N(U_xx)
                     +3 alpha2 N(U_x)+alpha3 N(U)+S3)]/(1-q),
Nfree <= max(2||d||/X+2||v||, ||d_r||+||v||).
```

This is a valid actual-magnitude majorant once each input is authenticated
and bounded. It is not obtained by differentiating an earlier error bound
or independent interval selections.
