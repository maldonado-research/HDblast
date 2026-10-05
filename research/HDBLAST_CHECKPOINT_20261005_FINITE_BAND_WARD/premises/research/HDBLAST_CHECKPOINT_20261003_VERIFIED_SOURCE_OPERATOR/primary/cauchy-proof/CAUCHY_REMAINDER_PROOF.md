# Analytic Cauchy remainder for the registered interior source

Status: analytic design and exact-arithmetic bound checks only. No physical
source values, trajectories, capsule arrays, or measured outputs were evaluated.
These notes are staged outside the repository and do not authorize execution.

## Formula provenance and domain

The `source_contract` in
`/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER/EXPERIMENT.json`
specifies

\[
 z=\eta+4,\quad B(z)=\exp\left(-\frac{z^2}{1-z^2}\right)
 =\exp\left(1-\frac1{1-z^2}\right),\quad L(\eta)=-\frac1\eta,
\]
\[
 h_+(\eta)=B(\eta+4),\qquad h_s(\eta)=(\eta+4)B(\eta+4),
 \qquad g=4L^2h-2Lh'-h''.
\]

The `source_jet` and `forcing_jet` definitions in the repaired primary source
file were read to confirm these formulas. They were not imported or executed.
The bump is zero outside its real support, but its interior expression is
holomorphic on any disk that excludes the singularities at eta=-5,-3. The
interval [-9/2,-7/2] lies strictly inside the support; the support cutoffs do not
enter this local analytic argument.

## One uniform complex-disk bound

Let c be any real point in [-9/2,-7/2] and let |eta-c| <= R=1/8. On this disk,

\[
 |z|\le \rho=\frac58,\quad |1-z^2|\ge d=1-\rho^2=\frac{39}{64},
 \quad |\eta|\ge\frac{27}{8},\quad |L|\le\ell=\frac8{27}.
\]

All inequalities include the closed disk. The singularities of B and L are
strictly outside it, so Cauchy's theorem applies with this radius.

For any complex w with |w| <= r < 1,

\[
 \mathop{\rm Re}\frac1{1-w}\ge \frac1{1+r}.
\]

One direct proof sets x=Re(w), s=|w|. First,

\[
 \frac{1-x}{1-2x+s^2}\ge\frac1{1+s},
\]

because after multiplication by positive denominators the difference in the
numerators is (1-s)(s+x) >= 0. Then 1/(1+s) >= 1/(1+r).
Taking w=z^2 and r=rho^2=25/64 gives

\[
 |B|\le\exp\left(\frac{25}{89}\right)\le\frac{89}{64}=b.
\]

The last bound is rational and requires no transcendental numerical estimate:
for 0 <= u < 1, the exponential power series gives exp(u) <= sum u^n =
1/(1-u), because n! >= 1.

Write D=1-z^2 and A=-2z/D^2. Differentiation gives

\[
 B'=AB,\quad A'=-2D^{-2}-8z^2D^{-3},\quad B''=(A^2+A')B.
\]

Consequently the following exact rational bounds hold:

\[
 |A|\le a=2\rho/d^2=\frac{5120}{1521},
\]
\[
 |A^2+A'|\le v=a^2+2/d^2+8\rho^2/d^3
 =\frac{70623232}{2313441}.
\]

For h_+=B, the forcing obeys

\[
 |g_+|\le b(4\ell^2+2\ell a+v)
 =\frac{951819044}{20820969}<64.
\]

For h_s=zB, product differentiation gives h_s'=B+zB' and
h_s''=2B'+zB''. Thus

\[
 |g_s|\le b\{4\ell^2\rho+2\ell(1+\rho a)+2a+\rho v\}
 =\frac{3227905133}{83283876}<64.
\]

Use the common exact bound M=64 for both source choices. The intermediate
rational bounds are approximately 45.715 and 38.758, respectively; these
decimal descriptions are not used by the proof or algorithm.

## Taylor and integrated remainders

Let a_n=g^(n)(c)/n!, P_N(t)=sum_{n=0}^N a_n t^n. Cauchy's estimate gives
|a_n| <= M/R^n. If |t| <= H < R and q=H/R, summing the coefficient tail gives

\[
 |g(c+t)-P_N(t)|\le\varepsilon_N
 =\frac{M q^{N+1}}{1-q}.
\]

Choose 64 exactly equal panels on [-9/2,-7/2]. Their centers and half width are

\[
 c_j=-\frac92+\frac{2j+1}{128},\quad j=0,\ldots,63,
 \qquad H=\frac1{128},\quad q=\frac1{16}.
\]

All these quantities are exact rationals. With M=64,

\[
 \varepsilon_N=\frac{2^{6-4N}}{15}.
\]

| Degree N | Uniform remainder upper bound |
| --- | --- |
| 20 | 2^-74 / 15, approximately 3.530e-24 |
| 24 | 2^-90 / 15, approximately 5.386e-29 |
| 28 | 2^-106 / 15, approximately 8.218e-34 |
| 32 | 2^-122 / 15, approximately 1.254e-38 |

The domain has length one. For every real phase frequency omega, multiplication
by exp(i omega eta) does not increase the modulus on the real interval.
Therefore the complete-domain phase integral error is at most epsilon_N.
A sharper bound, useful but optional, integrates each power in the same tail:

\[
 \sum_j\int_{-H}^H |g(c_j+t)-P_{N,j}(t)|\,dt
 \le\frac{\varepsilon_N}{N+2}.
\]

To see this, each term integrates to 2M H q^n/(n+1), and n+1 >= N+2 for all
n >= N+1. Summing and using 64*(2H)=1 proves the stated bound. It is below
2.072e-30 for N=24 and 3.688e-40 for N=32. Any real-axis kernel with a proved
absolute bound W multiplies this bound by W; no cancellation assumption is
needed. In particular a sine Green kernel divided by its frequency has the
continuous zero-frequency bound |sin(omega s)/omega| <= |s|.

These bounds account for analytic truncation only. Validated coefficient,
phase, kernel-moment, and accumulator arithmetic must supply separate
enclosures, or be included directly in a real/complex-ball result. A Taylor
remainder is not a bound on an unvalidated coefficient computation.

## Exact-rational coefficient construction

The proposed target is the stated analytic formula with exact rational source
constants and panel geometry. No rounded native source evaluator is involved.
Every Taylor coefficient is a rational multiple of one real exponential per
panel. This also avoids repeated transcendental source evaluations.

For z0=c+4, let Q(t)=1/[1-(z0+t)^2]=sum q_n t^n, d0=1-z0^2. Then

\[
 q_{-1}=0,\quad q_0=1/d_0,\quad
 q_n=(2z_0q_{n-1}+q_{n-2})/d_0 \quad(n\ge1).
\]

Let f0=1-q0 and f_n=-q_n for n>=1. Write
B(z0+t)=exp(f0)*sum p_n t^n. The exact rational recurrence is

\[
 p_0=1,\qquad
 p_n=\frac1n\sum_{j=1}^n j f_j p_{n-j}.
\]

Let b_n=exp(f0)*p_n and b_-1=0. Set h_n=b_n for the positive source and
h_n=z0*b_n+b_{n-1} for the signed source. The coefficients of L and L^2 are

\[
 l_n=\frac{(-1)^{n+1}}{c^{n+1}},\qquad
 m_n=\frac{(n+1)(-1)^n}{c^{n+2}}.
\]

Finally,

\[
 a_n=4\sum_{j=0}^n m_jh_{n-j}
 -2\sum_{j=0}^n l_j(n-j+1)h_{n-j+1}
 -(n+2)(n+1)h_{n+2}.
\]

For degree N, compute B through degree N+2. Compute the rational factors
exactly and enclose exp(f0) using a validated real-ball exponential. Both source
choices reuse the same bump factors. A pure ball recurrence is also possible,
but exact rational factors make dependency and coefficient provenance simpler.

## Entire phase moments and cost

With x=t/H and z=-i omega H, one panel integral is

\[
 H\exp(i\omega c)\sum_{j=0}^N a_jH^jJ_j(z),\qquad
 J_j(z)=\int_{-1}^1 x^j\exp(-zx)\,dx.
\]

All panels have the same H. The moments J_j(z) can therefore be reused for
every panel at a given exact represented frequency. They have the entire,
zero-frequency-safe series

\[
 J_j(z)=\sum_{m=0}^{\infty}\frac{(-z)^m}{m!}
 \frac{1+(-1)^{j+m}}{j+m+1}.
\]

If the registered frequencies imply |z| <= 4, truncation after m=T has the
explicit rational uniform bound

\[
 |J_j(z)-J_{j,T}(z)|
 \le \frac{2\,4^{T+1}}{(T+1)!}\frac1{1-4/(T+2)}
 \quad(T+2>4).
\]

This follows by bounding the exponential tail; successive absolute terms have
ratio at most 4/(T+2). A sharper denominator j+T+2 is available if needed.
No division by omega is used. The bound at T=48 is below 1.133e-33, at T=56
below 1.101e-42, and at T=64 below 3.514e-52. Actual frequency membership and
the |z| <= 4 claim must be checked from the registered target metadata, rather
than inferred from these analytic notes.

For orientation, if each J_j has absolute error <= delta and all other
operations are exact, |a_j H^j| <= 64*(1/16)^j implies a complete-domain moment
error <= (512/15)*delta. Thus T=48 is adequate relative to the degree-24
analytic tail, while T=64 makes this contribution negligible relative to the
degree-32 analytic tail. Phase-ball radii and all further arithmetic still need
to be propagated. A precision such as 192 bits is a reasonable starting
configuration for degree 32; achieved output radii, not the nominal precision,
determine readiness.

Degree 32 requires 35 bump coefficients at each of 64 panels, hence 64
validated real exponentials plus O(64*N^2) small exact-rational preprocessing.
Per frequency, moments require O((N+1)*(T+1)) elementary operations and all
panel reductions require O(64*(N+1)) terms per source. Center phases can be
advanced by a common enclosed phase increment with its error propagated, or
computed independently. The source preprocessing is reused across frequencies
and both source choices. This is a credible finite-cutoff design; no timing
claim is made because no physical source computation was run.

## Local checks

`check_bound_constants.py` verifies the rational arithmetic, strict M=64
inequalities, exact panel cover, remainder expressions, and phase-tail bound
constants. It contains no source evaluator and does not import project code.
The calculus and complex-analysis proof above remains necessary; these
arithmetic checks alone are not an integration certificate.
