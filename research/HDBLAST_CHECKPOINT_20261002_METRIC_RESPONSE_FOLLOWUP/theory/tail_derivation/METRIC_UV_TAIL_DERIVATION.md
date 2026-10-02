# Constructive UV envelopes for the prescribed metric calibration

2 October 2026. This is an analytic, prospective contribution. No physical
pulse values, source evaluations, response quadratures, modes, or roots were
computed. The executable `metric_tail_bounds` function remains deferred until
the enclosing experiment and all its source files are publicly frozen.

The result encloses only the omitted momentum band `k>K`. Its global
derivative envelopes are deliberately conservative. They do not certify
finite-k integration, mode errors, subtraction cancellation, time quadrature,
response arithmetic, a precise continuum stress, or a continuum sign.
Independently calibrated finite-K results should remain the principal
numerical claim unless additional, separately registered sharper tail
analysis supports more.

## Source and exact derivative certificates

Use the fixed H=1, r=2, minimally coupled calibration and

```
u=eta+4, d=1-u²,
B(u)=exp(1-1/d) for |u|<1, zero otherwise,
h=B or uB,
L=-1/eta,
g=4L²h-2Lh'-h''.
```

These h and g quantities have their infinitesimal amplitude epsilon removed.
Their support is `[-5,-3]`; all endpoint jets vanish. In particular `L<=1/3`
on the entire source support. Write

```
h^(n)=P_n(u)B(u)/d^(2n),
P_(n+1)=d²P_n'+[4nu d-2u]P_n,
P_0=1 for B; P_0=u for uB.
```

`certify_metric_derivative_bounds.py` extends this integer recurrence exactly
through order eight and checks it using an independent rational derivative
form. It supplies certified supremum and total-variation bounds for h
through derivative order seven. Derivative extrema are zeros of `P_(n+1)`.
Parity reduces nonzero extrema to roots in z=u². Exact polynomial root
isolation supplies disjoint rational intervals; exhaustive Sturm counts and
square-free checks exclude omitted or repeated roots. All 54 positive
critical-root intervals are preserved in the certificate. Parity supplies
the negative roots, with zero handled separately when it is critical.

Each root interval is refined to width at most `2^-60`; rational square-root
bisection and interval Horner evaluation enclose the derivative numerator
and denominator. With t=z/(1-z), the bump is `exp(-t)`. Rational positive
exponential Taylor bounds, using order 128, give

```
S_N(t)<=exp(t)<=S_N(t)+[t^(N+1)/(N+1)!]/[1-t/(N+2)]
```

for the verified t<N+2 at every interval. Reciprocals of these inequalities
bound the bump. The flat endpoint limits are exactly zero; thus maxima occur
among the exhaustive interior critical list. Summing absolute differences
between consecutive enclosed extrema and the two zero endpoints bounds
the full total variation. No floating exponential or source sampling is
used in this certificate.

The resulting integer supremum caps are:

| Derivative n | B upper | uB upper |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 3 | 2 |
| 2 | 22 | 18 |
| 3 | 507 | 448 |
| 4 | 22605 | 20593 |
| 5 | 1621069 | 1566176 |
| 6 | 221471343 | 213361523 |
| 7 | 39546260574 | 38135214987 |

Total variations are separately recorded. For m>=1,
`integral |h^m|=TV(h^(m-1))`; for m=0, the support-length bound is
`integral |h|<=2 sup|h|`.

## Global canonical-forcing norms

The exact Leibniz formula is

```
g^(n)=4 sum_(j=0)^n binom(n,j)(n-j+1)! L^(n-j+2) h^j
      -2 sum_(j=0)^n binom(n,j)(n-j)! L^(n-j+1) h^(j+1)
      -h^(n+2).
```

An independent direct differentiation check verifies this through n=5.
Supremum caps follow by absolute values and L<=1/3. L1 caps use the exact
total-variation envelopes of the corresponding h derivatives. Therefore

```
N_j(g)=|g^(j+2)(eta)|+integral_-6^eta |g^(j+3)(eta')| deta'
      <=sup|g^(j+2)|+L1|g^(j+3)|, j=0,1,2.
```

The global rational caps are:

| Profile | N0 | N1 | N2 |
|---|---:|---:|---:|
| B | 15139255/81 | 3936783602/243 | 517298697427/243 |
| uB | 14332867/81 | 3775755359/243 | 498562987423/243 |

N_j is exactly zero before the source. After the source, local g and h jets
are exactly zero while these history norms remain as conservative bounds.
The large N2 is approximately 2.13 billion for B. Its contribution alone to
the q'' tail at K=256 is of order 50 in the declared normalized units. This
shows why these rigorous global bounds cannot establish a precise continuum
pressure or its sign at the selected finite cutoffs.

## Canonical memory and scalar-proxy stress tails

Put `M²=2L²`, `v=K/sqrt(K²+M²)`, and

```
1-v=M²/[sqrt(K²+M²)(sqrt(K²+M²)+K)].
```

The stable right-hand form avoids subtracting nearly equal numbers.
For each canonical derivative g^j define

```
T_j=[3|g^j|M²/2+N_j(g)/4]/[16pi²K²].
```

The existing mass-channel integration-by-parts bound applies to the new
canonical source g. With C=1/(8pi²), its first two derivative contacts give

```
E_q=T_0,
E_q'=T_1+C L|g|(1-v³),
E_q''=T_2+C[2L|g'|(1-v³)
           +L²|g|(1-v)(3v⁴+3v³+v²+v+1)].
```

Here q=a²deltaQ_proxy/epsilon. The source-mass proxy is only an intermediate
calculation; explicit metric observable and subtraction contacts must be
added to obtain the actual metric response.

The physical baseline variance tail has the stable exact expression

```
R_Q=(1-v)²(v+2)(3v³+3v²+10v+4)/[96pi²(v+1)].
```

The normalized scalar-proxy density and pressure tails obey

```
E_rho,proxy <=(3L²E_q+L E_q')/2+L²|g|R_Q/2+E_Crho,
E_p,proxy   <=(E_q''+3L E_q'+3L²E_q)/6+L²|g|R_Q/6+E_Cp,

E_Crho=L(3L|g|+|g'|)(1-v³)/(96pi²).
```

For pressure, use the direct minimal pressure subtraction reduction:

```
E_Cp=[70L²|g|J9tail+(30L²|g|+10L|g'|)J7tail
           +(|g''|+L|g'|+9L²|g|)J5tail]/(48pi²),

J5tail=(1-v³)/6,
J7tail=(1-v)²(3v³+6v²+4v+2)/60,
J9tail=(1-v)³(15v⁴+45v³+48v²+24v+8)/840.
```

These moments retain the physical reference r=2 factors. Reusing normalized
moments without their factors would give incorrect denominators. The three
stable expressions are checked against their exact primitive differences.

## Full metric contact and baseline tails

The primary route independently supplies rational contact inventories
`sum c_n(L,h)/w^n`, `w=sqrt(k²+2L²)`, with measure `k²dk/(2pi²)`.
The inventory is faithfully copied and SHA-pinned under `reference_inputs`.
Every contact has odd n>=5; therefore

```
integral_K^infinity k²/w^n dk<=K^(3-n)/(n-3).
```

Absolute coefficient bounds multiply these positive moments. The q-contact
derivative inventories are generated and independently verified using

```
D=L²partial_L+sum h^(j+1)partial_h^j,
partial_eta[c_n/w^n]=(D c_n)/w^n-2nL³c_n/w^(n+2).
```

The `-2nL³` contribution is essential: holding the physical reference r fixed
does not keep w fixed as geometry changes.

Exact physical baseline tails `Q0(infinity)-Q0(K)`,
`rho0(infinity)-rho0(K)`, `p0(infinity)-p0(K)` are derived from the primary's
independent finite-band primitives. All factor at least `(1-v)²` and have
positive denominators proportional to `(1+v)` or `(1+v)²`. Absolute polynomial
coefficient envelopes bound their numerators. The continuum values are

```
Q0=1/(12pi²), rho0=11/(960pi²), p0=-rho0.
```

For the actual metric response, the baseline pieces are

```
q_local=-2hL²Q0+Cq,
rho_local=-4hL⁴rho0+(h''+4Lh')L²Q0/2+Crho,
p_local=-4hL⁴p0-h''L²Q0/2+Cp.
```

The corresponding tail bounds retain these factors and all reference
differences. For q' and q'', the baseline q_local tail is differentiated
with `v'=-Lv(1-v²)` before taking absolute polynomial bounds. Serialized
derivative/contact inventories and stable factored baseline tails are
checked against the original symbolic expressions. Actual metric upper
bounds add these local tail envelopes to the scalar-proxy ones. The current
tail equals q for the frozen b=1,r=2,delta_phi=0 convention.

## Deferred runtime API and checks

The self-contained runtime dependency set is:

```
metric_tail_bounds.py
DERIVATIVE_CERTIFICATE.json
METRIC_TAIL_INPUTS.json
```

All three must occupy the same directory and enter the prospective manifest.
After the public freeze, call

```
metric_tail_bounds(source, eta, cutoff)
```

for the frozen six observation times. The return contains q, q_prime,
q_second, rho, p, current, Q0, rho0 and p0 upper bounds, together with separate
canonical proxy, full metric contact and baseline metadata. A fresh
60-decimal-digit directed interval context handles geometry, momentum,
pi, all rational inputs and all bound propagation. Each nonzero binary64
upper endpoint advances toward positive infinity.

The quantities q, q derivatives and current use `a²deltaQ/epsilon` and
`a²deltaj/epsilon`; density and pressure use `a⁴delta(rho,p)/epsilon`.
Baseline bounds are physical H=1 differences, not divided by epsilon.
The only accepted pulse IDs are positive_B and signed_uB. Before the source,
response tails are exactly zero, while finite-band baseline differences
remain separately bounded.

Normal and optimized certificate runs pass 32 exact derivative/parity
identities and all 54 exhaustive rational root intervals. Normal and
optimized `verify_metric_tail_algebra.py` pass 39 serialization, derivative,
primitive and Leibniz identities and detect six synthetic wrong-formula
mutations. No deferred runtime call on physical observations has occurred.
These checks validate analytic preparation; they do not validate a later
physical implementation or certify total numerical error.
