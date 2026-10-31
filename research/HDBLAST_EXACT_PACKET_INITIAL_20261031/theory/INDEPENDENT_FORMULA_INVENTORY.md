# Independent all23 reconstruction and regularity inventory

The independent evaluator implements the following completed table directly;
it imports neither `Context._registered` nor a new producer. Write
`c=(1/6,1/2,1)`, `f_j=fs[j]`, `g_j=dy(f_j)`, and use the registered first
profiles `u,v,ut,vt`, with `dy=-rho^-1 d_z`. Define

```
A = 2*rho*K - 8*u*ut
B_j = rho*((-26*u*dy(u)-2*a*u^2)*c_j - 3*g_j/2 - 3*a*f_j)
C = -8*rho*K + 44*u*ut
D = 2*rho*K + 4*u*ut
S = rho*Pi + 2*u*vt.
```

Here `K` and `Pi` are the registered quadratic profiles, not angular vector
labels. `m` is the highest raw packet derivative entering a reconstructed
component value. A dash denotes an identically zero component.

| Label | Completed displacement | m | Registered velocity before old P1/smoothing | m |
| --- | --- | ---: | --- | ---: |
| K:0 | 0 | — | ell*A/6 | 1 |
| K:2 | 0 | — | ell*A/2 | 1 |
| K:4 | 0 | — | ell*A | 1 |
| L:0 | 0 | — | ell*B_0 | 3 |
| L:2 | 0 | — | ell*B_1 | 3 |
| L:4 | 0 | — | ell*B_2 | 3 |
| N:0 | -u^2/3 | 0 | ell*C/6 | 1 |
| N:2 | -u^2 | 0 | ell*C/2 | 1 |
| N:4 | -2*u^2 | 0 | ell*C | 1 |
| J:0 | f_0 | 2 | ell*D/6 | 1 |
| J:2 | f_1 | 2 | ell*D/2 | 1 |
| J:4 | f_2 | 2 | ell*D | 1 |
| w:0 | -3*g_0/(2*p) | 3 | ell*S/6 | 2 |
| w:2 | -3*g_1/(2*p) | 3 | ell*S/2 | 2 |
| w:4 | -3*g_2/(2*p) | 3 | ell*S | 2 |
| v:2 | 0 | — | 0 | — |
| v:4 | 0 | — | 0 | — |
| u:2 | 0 | — | ell*(11*u^2/2-f_1/2) | 2 |
| u:4 | 0 | — | ell*(11*u^2-f_2/2) | 2 |
| H:2 | 0 | — | 0 | — |
| H:4 | 0 | — | 0 | — |
| T:+ | 0 | — | 0 | — |
| T:- | 0 | — | 0 | — |

All five mean overrides are included. Their exact formulas agree with the
listed entries because `c_0=1/6`; deleting only a redundant algebraic
assignment does not give an effective nonzero corruption. Zeroing a mean
or changing its factor does.

For every row, the retained C2 transfer theorem requires displacement
through its radial second derivative and velocity through its first radial
derivative. Thus raw packet derivatives through order five suffice for
displacement; reconstructing analytic velocity through r would need order
four. This candidate instead retains the old complete smoothed P1 velocity
as its defining function.

For each completed displacement component `d_i`, the exact interior traces
at R=7/4 are `(d_i,d_i,r,d_i,rr)(R-)`. The componentwise quadratic exterior
has exactly those same three traces; all three differences are zero for
every row above. Its third trace generally jumps, so no global C3 premise is
used. The stored P1 velocity, parity reflection, final-secant affine
exterior, and compact C2 convolution kernel give at least C4 regularity;
its value, first and second traces are continuous for every row. The exact
kernel endpoint jets and displacement trace differences are in the
independent regularity report; the old velocity node/slope values and
trace data are separately authenticated.

At r=1 all required completed jets vanish by exact packet flatness and
the independently separated nominal rho and p denominators. Both sides
of 65/64 and 66/64 use the same exact analytic q/E formula. Tail and interior
branches are enclosures of that one function. The independent tail proof
derives the derivative polynomials through five and proves their domination
by s^(2n) on s>=4096/63 using exact nonnegative shifted coefficients.
