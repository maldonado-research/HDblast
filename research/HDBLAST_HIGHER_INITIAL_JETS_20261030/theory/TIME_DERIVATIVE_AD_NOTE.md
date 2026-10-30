# Independent derivative of the original full23 assembly

This audit differentiates the authenticated original
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py`
assembly with exact rational first-order dual numbers. It is independent of
the new actual interval coefficient-derivative producer. Its four rational
manufactured backgrounds test implementation; the general identity follows
from the displayed analytic product rule, not from finite sampling.

Hold r fixed. Put `tau=tanh(ell*x)`, `s=sech(ell*x)` and `h=tau/rho`.
The actual background parameters `rho,a,p,b,d,V,V1,V2,py,sigma1,sigma2`
are r-only data, so all have zero x derivative. The only scalar x rules are

```
tau_x=ell*s^2,       s_x=-ell*s*tau,       h_x=ell*s^2/rho.
```

The original matrix coefficient is

```
D=ell^2 [s^2 Delta+2rho*a*s Ca+2rho*s*h Ch
          +rho^2 C+rho*p*s Ephi],
C=C_radial+a*h(Sah-3O)+h^2 Sh+3p*h W_L.
```

Here `C_radial` includes every time-independent curvature, potential,
normal-scalar and fibre-connection term from the original assembly. No
selected field row is omitted. Differentiating the five terms of D gives

| Original coefficient | x derivative |
| --- | --- |
| `ell^2 s^2 Delta` | `-2ell^3 s^2 tau Delta` |
| `2ell^2 rho*a*s Ca` | `-2ell^3 rho*a*s*tau Ca` |
| `2ell^2 rho*s*h Ch` | `2ell^3 s(s^2-tau^2) Ch` |
| `ell^2 rho^2 C` | `ell^3 [rho*a*s^2(Sah-3O)+2tau*s^2 Sh+3rho*p*s^2 W_L]` |
| `ell^2 rho*p*s Ephi` | `-ell^3 rho*p*s*tau Ephi` |

Thus the complete derivative is

```
D_x=ell^3 [-2s^2 tau Delta-2rho*a*s*tau Ca
           +2s(s^2-tau^2)Ch+rho*a*s^2(Sah-3O)
           +2tau*s^2 Sh+3rho*p*s^2 W_L-rho*p*s*tau Ephi].
A_x=-3ell^2 s^2 I,       B_x=0,       M_x=0.
```

The constant-shell M has no x dependence. Dividing `D_x` by ell gives the
ordinary physical-t derivative; confusing the two loses a factor ell. These
identities differentiate analytic functions before bounding them. They do
not differentiate the endpoints of a coefficient interval enclosure.

## Exact independent implementation checks

The checker first pins the original operator bytes and then imports that
standard-library module with bytecode writes disabled. It reuses the exact
primitive matrices but calls the original, unmodified assembly function.
The scalar inputs `c=1/s`, `Ht=tau`, and `h=tau/rho` are rational dual jets;
the remaining background variables have zero dual derivatives. All entries
of A, B, D and M are automatically differentiated through the original sum,
product and quotient operations and compared with the expanded derivative
above.

Four rational circle points have `tau^2+s^2=1` and `s>0`; their hyperbolic
derivative jets therefore obey the actual scalar differential identities.
The backgrounds are manufactured algebraic fixtures, not solutions of the
Einstein–scalar background equations. Seventeen check groups cover 8,464
matrix-entry equalities. Eight wrong-formula controls reject a physical-t
derivative substituted for x, omission of h_x, the wrong sech derivative
sign, and omission of five nonzero product-rule matrix terms. These checks
remain active under optimized Python, and the normal/optimized reports are
byte identical. The original operator bytes are read back unchanged.

No bound on the actual complete bulk-source derivative F_x, candidate second
jets, C2 error, natural Einstein residual or physical Big Bang mechanism
follows from this audit. The actual continuous coefficient derivative norms
are a separate interval computation using the same analytic formula.
