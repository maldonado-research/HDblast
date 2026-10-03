# Adversarial audit of the primary moment candidate

Audited file:
`/workspace/hdblast-research-work/verified-integration-20261003/primary-design/verified_moments.py`

Audited current SHA256:
`de995d2edeca3944400fd90b892521f05c941b6d8a1f75c2a820740b82d995c0`

Verdict: the formulas and enclosure factors reviewed below pass. One silent
series-truncation defect was found during review and repaired by the primary
author. The repaired source callback explicitly sets/restores `ctx.cap` and
rejects an insufficient result precision. No physical source callback,
checkpoint-array decoder, trajectory, or scientific run was invoked.

## Resolved finding: global series cap

An explicit `arb_series(..., prec=27)` does not override python-flint's global
`ctx.cap`, which defaults to 10. A manufactured source h(eta)=2+3eta+5eta^2 at
c=2, outside the physical source interval, reproduced the failure: after two
derivatives, the assembled g series had `prec=8`; retrieving coefficient 24
returned exact zero even though the corresponding exact rational coefficient
was nonzero. Inflating with a degree-24 Cauchy tail would not repair that false
coefficient.

The current callback sets `ctx.cap=SOURCE_DEGREE+3` in a try/finally scope and
requires `g.prec >= SOURCE_DEGREE+1`. Inspection confirms that this resolves the
identified mechanism. The manufactured source's 25 normalized g coefficients
all contain independently computed exact rational coefficients at cap 27.
The physical callback itself was not executed.

The added pure `forcing_series_coefficients` helper was executed with that
manufactured source. All 25 g coefficients and all 25 Lg coefficients contain
their independently computed exact rational values. The helper rejects a
manufactured insufficient `ctx.cap=10` before returning any coefficients.

## Adversarial challenge to the source majorant

The M=64 proof does not replace a complex real-part estimate by an absolute
value estimate. Let w=z^2, s=|w| and x=Re(w). For s<1,

\[
 \mathop{\rm Re}\frac1{1-w}=\frac{1-x}{1-2x+s^2}.
\]

Subtracting 1/(1+s), after multiplying positive denominators, produces
(1-s)(s+x), which is nonnegative since x >= -s. Therefore

\[
 \mathop{\rm Re}\frac1{1-z^2}\ge\frac1{1+|z|^2}
 \ge\frac{64}{89}\quad\text{for }|z|\le5/8.
\]

It follows that |B| <= exp(25/89) <= 89/64. The second inequality follows from
the positive exponential series and the geometric-series bound, with no
numerical evaluation of B. The derivative and forcing bounds in
`CAUCHY_REMAINDER_PROOF.md` then give both source forcings strictly below 64.
The complex boundary and both signs of the real source coordinate are included.

## Source orientation, coefficients, and scaling

The code builds eta=c+H*x and z=eta+4. It uses the registered positive source B
and signed source z*B. Its derivatives divide by H once and twice, respectively,
so they are derivatives with respect to eta, not the normalized x variable.
The expression 4*L^2*h-2*L*h_prime-h_second matches the formula metadata and
the inherited primary's pure forcing expression.

Order `SOURCE_DEGREE+3` gives source coefficients through degree N+2. Two
derivatives leave g at precision N+1, enough to retrieve degrees 0 through N.
Arb series arithmetic and rational-to-ball construction enclose the actual
Taylor coefficients once the cap is correctly enforced; coefficient intervals
need not be independent for the resulting interval arithmetic to be safe.
The hardcoded source tail 1/(15*2^90) is exactly the M=64, R=1/8, H=1/128,
N=24 Cauchy remainder. A future change of degree or geometry must recompute
this bound rather than retain the literal unchanged.

The added work series is exactly L*g on the same normalized variable. On every
proved complex disk, |L*g| <= (8/27)*64 = 512/27 < 32. The chosen work
majorant 32 and work tail 1/(30*2^90), half the forcing tail, are consequently
valid. Multiplying a g series of precision N+1 by the L series preserves the
N+1 coefficients needed for this work series; no extra derivative is taken.

## Exact exponential majorant for the E/Q tails

A short purely rational proof supplies the required exp(8)<4096 bound:

\[
 e=\frac52+\sum_{n=3}^{\infty}\frac1{n!}
 \le\frac52+\frac{1/6}{1-1/4}=\frac{49}{18}<\frac{11}{4}.
\]

Every ratio after the 1/3! term is at most 1/4, proving the geometric tail
inequality. Hence

\[
 \exp(8)=e^8<(11/4)^8=\frac{214358881}{65536}<4096,
\]

where the final exact integer comparison is 214358881<268435456.

For |z|<=4 and x in [-1,1], u=z(1-x) has |u|<=8. After retaining m=0,...,T,
the exponential remainder is bounded by exp(8)*8^(T+1)/(T+1)!. Integration
adds a factor at most 2. This is exactly the implemented

\[
 E\text{-tail}=2\cdot4096\cdot8^{T+1}/(T+1)!.
\]

The Q series retains z^m*(1-x)^(m+1)/(m+1)! for m=0,...,T. Its remaining
pointwise terms have one extra factor (1-x)<=2 and denominator (T+2)!.
Integration adds another factor at most 2, giving exactly

\[
 Q\text{-tail}=4\cdot4096\cdot8^{T+1}/(T+2)!.
\]

Neither series tail divides by z. The code converts each proven complex disk
tail to a square by inflating both real and imaginary parts; that is a safe
outer enclosure. The exact target z is pure imaginary and satisfies the cap.
The rational-to-ball representation of z encloses that target, so evaluation
of the finite polynomial on that ball plus the exact-target tail remains safe.

## Panel factors and the additive global formula

With lambda=2ik, z=lambda*H, and r=c+H, the code's local moments are

\[
 M_0=\int_l^r g(s)\,ds,
 \quad M_{\exp}=\int_l^r g(s)e^{\lambda(r-s)}\,ds,
 \quad M_{\rm drift}=\int_l^r g(s)\phi(r-s)\,ds,
\]

where phi(d)=(exp(lambda*d)-1)/lambda is defined continuously at k=0.
Their coefficient factors H, H, and H^2 are correct. The mixed monomial
formula follows directly from y=(1-x)/2 and is exact rational arithmetic.

For a uniform source error <= epsilon_g, the first two moment errors are
at most 2H*epsilon_g because real-axis phases have modulus one. For the drift,
|phi(d)|<=d at real k, so its local error is at most
epsilon_g*integral_0^(2H) d dd = 2H^2*epsilon_g. The three implemented source
inflation factors therefore pass.

The identity

\[
 \phi(d+s)=\phi(d)+e^{\lambda d}\phi(s)
\]

proves the global weighted drift sum. The resulting response is exactly the
Duhamel solution of u'=w and w'=2ik*w-epsilon*g:

\[
 w(b)=e^{\lambda(b-a)}w(a)-\epsilon\int_a^b e^{\lambda(b-s)}g(s)\,ds,
\]
\[
 u(b)=u(a)+\phi(b-a)w(a)-\epsilon\int_a^b\phi(b-s)g(s)\,ds.
\]

The phase and forcing signs agree with the inherited primary's `phase_step`.
The reported global source-only bounds are the corresponding triangle bounds;
on a unit interval with a common epsilon_g they simplify to
|epsilon|*epsilon_g for w and |epsilon|*epsilon_g/2 for u. The output balls
also propagate coefficient, moment, phase, and incoming-state arithmetic.

## Small momentum and exported endpoints

For |2k*d|<=1, `stable_drift` uses d*sum_{m=0}^96 z^m/(m+1)!. Its remainder
is at most d*3/98!, since exp(1)<3. This includes k=0 and d=0 without a
division. In the other branch, |2k*d|>1 implies k>0 and d>0; dividing by the
exact represented nonzero frequency is safe, and Arb propagates cancellation
and division uncertainty. It is not a recurrence in inverse momentum.

`rational_endpoints` rejects nonfinite balls and exports the exact rational
values of outward lower and upper bounds. No formatted midpoint decimal is
used as a bound. Fabricated exported endpoints contained the independent
higher-precision reference enclosures in every tested component.

## Fabricated execution evidence

Command:

```
PYTHONDONTWRITEBYTECODE=1 /workspace/hdblast-cloud-setup/venv-arb/bin/python \
  /workspace/hdblast-research-work/verified-integration-20261003/primary-design/cauchy-proof/audit_fabricated.py
```

The audit uses h(eta)=2+3eta+5eta^2 at c=2 solely to challenge normalized source
jet scaling and series truncation. Its response checks use the manufactured
forcing g(eta)=2+3eta+5eta^2 on [0,1], with manufactured complex incoming data
and epsilon=7/11. They compare the candidate against independent closed
polynomial integrals evaluated at 1024 bits.

All ten response cases pass: five exact rational fabricated momenta
0, 10^-30, 10^-4, 17/3, and 256, both with zero source remainder and with a
manufactured constant source perturbation of 10^-4 covered by that remainder.
Coefficient balls also carry manufactured uncertainty. The checks cover
homogeneous evolution, source sign, phase direction, all normalized panel
factors, global additive composition, zero and very small momentum, source-only
bound factors, coefficient-ball coverage, and exact endpoint export.

This is formula and fabricated-enclosure evidence, not a scientific result or
timing benchmark. The complete physical driver, registration enforcement,
represented-input authentication, reductions, and resource limits remain the
primary author's separate work.
