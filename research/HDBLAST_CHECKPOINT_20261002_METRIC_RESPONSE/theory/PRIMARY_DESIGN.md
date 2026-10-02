# Primary metric-response route

This is prospective staged work against inherited science commit
`b807a4d549a40bc78b66122cb5719987091d2ce7`. No new physical source values,
integrals or modes were evaluated while preparing it. Only symbolic algebra,
AST checks, dependency inspection and blocked-call guards were executed.

The fixed physical scalar mass and reference are `x=r=2`, with `H=1`, `xi=0`
and the fixed incoming BD state. The homogeneous variation is
`a_epsilon=a0(1+epsilon h)`, `a0=-1/eta`, `L=-1/eta`. The inherited exact
bump jet provides `h0..h5`. The canonical mode forcing is

```
g = 4 L^2 h - 2 L h' - h''.
g_n = 4 sum_j binom(n,j) (n-j+1)! L^(n-j+2) h_j
      - 2 sum_j binom(n,j) (n-j)! L^(n-j+1) h_(j+1)
      - h_(n+2), n=0..3.
```

The scalar-memory calculation applies the inherited logarithmic continuum
and finite-band history integrals to this `g`. It obtains `qmass`, its first
two derivatives, and two separately reduced direct scalar-mass stresses.
The full metric response requires additional terms:

```
q = qmass - 2 h L^2 Q0,K + Cq,K
R = Rmass - 4 h L^4 rho0,K + (h''+4 L h') L^2 Q0,K/2 + Cr,K
P = Pmass - 4 h L^4 p0,K - h'' L^2 Q0,K/2 + Cp,K.
```

`Q0,K`, `rho0,K`, `p0,K` are three separately integrated renormalized baseline
observables. They are physical values with no epsilon division. Each complete
baseline primitive is checked by differentiation against its bare integrand
minus complete subtraction, including the lower endpoint. The baseline
logarithms cancel exactly. These are evaluated with stable rational functions
of `v=K/sqrt(K^2+2L^2)` rather than subtracting large integrated powers of K.

The independent primary derivation differentiates the complete full-metric
W0/W2/W4 subtraction at fixed r. It uses `delta w=2L^2h/w`, `delta L=h'`,
`delta(-a''/a)=-h''-2Lh'`, their exact derivatives, and the full stress
operator and metric measure. Subtracting the scalar-g reduction leaves the
convergent rational inventories stored in `METRIC_CONTACT_ALGEBRA.json`.
They agree exactly with the independent theory inventory coefficient by
coefficient. Every contact moment is integrated by

```
I_n = integral_0^K k^2 (k^2+2L^2)^(-n/2) dk
    = (2L^2)^((3-n)/2) sum_(j=0)^((n-5)/2)
      (-1)^j binom((n-5)/2,j) v^(2j+3)/(2j+3).
```

An overall `1/(2 pi^2)` multiplies these moments. The continuum contact
inventory follows by setting v=1, and is independently checked. The contact
derivatives use the exact fixed-comoving-cutoff differential operator
`D=L^2 partial_L - L v(1-v^2) partial_v + sum h_(j+1) partial_h_j`.
The anomaly inventory is complementary exact theory only. It is excluded
from runtime coefficients and this numerical experiment's output quantities.

Both stresses are directly defined; conservation is a separate check.
The exact full-metric contact bridge obeys
`R' - L R + 3 L P + 3 h' L^4(rho0,K+p0,K)=0`.
The pure coordinate translation `h=L T`, whose forcing vanishes, cancels all
continuum metric responses exactly. Neither check is a definition of pressure.

The output `values` contain q, q_prime, q_second, rho, p, Q0, rho0, p0,
current. Here q and current are `a0^2 deltaQ/epsilon` and
`a0^2 delta_j/epsilon`, respectively, and are equal because phi is fixed and
`b=1,r=2`. Density and pressure are `a0^4 delta_rho,p/epsilon`. Derivatives
refer to the scaled q; `scaled_density_prime_over_epsilon` is R'.

The producer verifies every registered file, the entire experiment and gates,
and the independent manifest before any physical call. The public 40-hex
freeze commit is an explicit external provenance pin; root separately verifies
the public remote tree. Imports are guarded and perform no physical evaluation.
Fresh outputs are required outside the frozen checkpoint, with no overwrite.
Runtime is capped at 900 seconds in Python 3.12, NumPy 2.2.6, SciPy 1.15.3,
SymPy 1.14.0 and mpmath 1.3.0. Integration warnings are fatal.

Finite memory integrals report numerical quadrature estimates, and exact
contacts have floating-point arithmetic error. Only a separately supplied
directed omitted-band bound may be described as an interval enclosure. No
finite-integral certification, coupled evolution, shifted-root response,
state variation, stability or heating claim follows from this bounded task.

Development artifacts `algebra-v1/v2/v3` preserve successive additions of
baseline/Ward/coordinate and complementary anomaly checks. They all passed;
they contain no physical samples. `algebra-final-normal` and
`algebra-final-optimized` are the final identical exact inventories. No
physical failure or tolerance tuning occurred during primary development.
