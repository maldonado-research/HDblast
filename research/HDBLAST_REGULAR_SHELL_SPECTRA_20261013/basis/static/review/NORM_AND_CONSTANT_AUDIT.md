# Independent norm propagation and exact constants, 11 October 2026

This local continuation preserves the accepted theorem and APS package. The
calculations are for the declared polynomial boundary-value problem. Arithmetic
alone does not prove the analytical estimates, certify a floating-point shooting
trajectory, establish stability, or establish external peer review.

## Scope and principal result

For the exact-decimal model,

\[
k=1/9,\quad W(x)=1/3+x^2+x^3/3,\quad
c=2/(1.0357712571566784)-4/3,
\quad f(x)=1+c+cx,
\]

the denominator in the decimal expression is interpreted as an exact rational.
Use

\[
r=1/100,\quad \bar\delta=1/500,\quad s=14/9,
\quad B=2,\quad L=|f'(0)|/w=c/2,
\quad b=L\bar\delta.
\]

The exact producer verifies the eleven surviving construction conditions in the
companion Jacobi theorem. It also records two explicit polynomial convexity
checks and a redundant quantitative shell-slope check. Those extra checks are
sufficient for these target models; the theorem derives the required strict
bulk convexity from scalar contraction and does not require shell convexity
separately from its positive scalar-junction derivative.

With the improved linear endpoint bound, the uniform constants on this entire
interval satisfy the exact rational inequalities

\[
K_\eta<3196/1000=3.196,\qquad
K_H<254/1000=0.254,
\]

\[
q-K_H\bar\delta>423/10^6=0.000423,
\qquad q=\frac{k f'(0)^2}{24(w-2k)}.
\]

Approximate displays are
\(K_\eta=3.1951835147879217\),
\(K_H^{\rm Taylor}=0.27516506559461146\),
\(K_H^{\rm polynomial}=0.253358997308784\), and
\(q-K_H\bar\delta=0.00042328124072393033\).
All decisions in the JSON certificate use fractions, rather than these displays.

Thus, under the separately reviewed analytic estimates, every
\(0<\delta\leq0.002\) has the uniform error bounds

\[
|\beta-\alpha\delta|<3.196\delta^2,\qquad
|H^2-a_1\delta-(f_0^2/36-q)\delta^2|<0.254\delta^3,
\]

where \(\alpha=-f_1/[4(w-2k)]\) and \(a_1=kf_0/3\).
The resulting strict curvature comparison is with the separate constant-tension
model \(f(x)=f_0\). Its scalar-free solution generally fails the scalar junction
of the model with nonzero \(f_1\); it is a precisely specified comparator.

| Detuning | Readable uniform scalar error bound | Readable uniform curvature error bound |
|---|---:|---:|
| 0.0005 | \(<7.99\times10^{-7}\) | \(<3.175\times10^{-11}\) |
| 0.001 | \(<3.196\times10^{-6}\) | \(<2.54\times10^{-10}\) |
| 0.002 | \(<1.2784\times10^{-5}\) | \(<2.032\times10^{-9}\) |

These are error bounds for mathematical branch expansions. They are not error
bounds for the archived ODE output.

## Norm constants derived in physical units

Write \(\omega_T=F_\nu/F_\nu(T)\), with \(\nu=\lambda/2\),
\(\lambda=w-4k\), \(m^2=w\lambda\),
\(D=m^2-\nu(\nu+4k)>0\), and

\[
\|v\|_{T,s}=\max\left(\sup|v|/\omega_T,
                         \sup|v'|/(s\omega_T)\right).
\]

The maximum principle gives \(|G_TN|\leq\|N/\omega_T\|\omega_T/D\).
The divergence equation gives the physical derivative bound

\[
|(G_TN)'|\leq\frac{1+m^2/D}{4k}
              \|N/\omega_T\|\omega_T.
\]

Dividing only the latter by \(s\) proves
\(C_G=\max(D^{-1},(1+m^2/D)/(4ks))\).
The homogeneous solution has norm at most
\(M=\max(1,\lambda/s)\), so \(B=2M\).

For \(\|\eta\|_{T,s}\leq B|\beta|\), its metric source is
\(\mathcal F=p^2/4+[U(\eta)-U(0)]/6\), and

\[
|\mathcal F|\leq C_FB^2\beta^2\omega_T^2,
\qquad C_F=s^2/4+M_2/12.
\]

The metric kernel and its radial derivative depend on the geometric scale
\(k\), so the metric estimates retain their denominators \(k\). In particular,

\[
|a-a_0|\leq C_FB^2\beta^2\omega_T^2/k,
\]

\[
|a-\widetilde a|\leq4C_FB|\beta|d_\eta\omega_T^2/k.
\]

For clarity, put \(A_F=C_FB^2\beta^2\) and
\(h_F=2C_FB|\beta|d_\eta\). The Volterra inverse has norm at most \(4/3\),
\(2/3\leq R\leq4/3\), and
\(\|\Delta R\|\leq(16/9)h_FJ/(2k)\).
The two quotient-derivative contributions give

\[
|\Delta(a-a_0)|\leq
\left[1+4\varepsilon_V\right](h_F/k)\omega_T^2
\leq2(h_F/k)\omega_T^2,
\]

where \(\varepsilon_V=A_FJ/(2k)\leq1/4\). This independently verifies the
factor 4 in the metric Lipschitz bound.

Multiplying these metric bounds by the physical scalar bound
\(|p|\leq sB|\beta|\omega_T\) proves the essential factors \(s/k\):

\[
C_N=L_1B^2+4C_F(s/k)B^3b,
\qquad L_N=2L_1Bb+20C_F(s/k)B^2b^2.
\]

Differentiation in the endpoint parameter yields
\(\|\eta_\beta\|_{T,s}\leq B\), and the endpoint momentum bounds carry
physical factors \(s\):

\[
C_P=sC_GB[2L_1B+20C_F(s/k)B^2b],\quad
C_D=sC_GC_N,\quad C_a=4C_FB^2/k.
\]

The archived prospective formulas agree with these independent derivations.
Its executable uses an archive-relative accepted-evaluator import that is not
available at that location, so this continuation uses a self-contained producer.
That is a packaging issue rather than a counterexample to its norm estimates.

## Improvements beyond norm rescaling

The accepted scalar monotonicity estimate already gives
\(|\beta|\leq|f_1|\delta/w\). Accordingly \(L=|f_1|/w\) is valid; the
previous \((|f_1|+1)/w\) was a convenient larger choice. If \(f_1=0\),
the scalar root is exactly \(\beta=0\), and a separately positive tube radius
is chosen rather than trying to set \(b=L\bar\delta=0\).

The companion theorem proves exact monotonicity of the metric residual from
regular-cone Jacobi identities, with strict bulk convexity already implied by
\(C_GL_N\leq1/2\):

\[
M_3Bb\leq D/2,\qquad U''\geq m^2-D/2>0
\quad\text{on the realized field tube}.
\]

This removes the old absolute \(\delta^2\) slope gate. The producer retains
that gate as a failed historical comparison, not as a current requirement.
The independently derived quantitative alternative is also retained:

\[
C_v=8C_FB^2J/(3k),\quad
a_{\min}=k-C_FB^2b^2/k,\quad C_Q=C_v/a_{\min},
\]

\[
\theta=C_Q(C_{T\beta}+sB)L^2\bar\delta^2<1.
\]

The selected model has \(\theta\approx0.06821527\); this is redundant under
the exact structural monotonicity theorem. With the former loose root bound,
the refined \(a_{\min}\) estimate gives \(\theta\approx0.51809\), whereas
the coarser \(a_{\min}=k/2\) choice fails at 0.002. This demonstrates why
individual sufficient estimates should be checked without treating failures as
nonexistence of a branch.

The improved Riccati estimate follows from
\(E=(\lambda-F_\lambda'/F_\lambda)/(a_0-k)\):

\[
E'=4\lambda-2(w-3k)E+(a_0-k)E(E-3),\quad E(0+)=0.
\]

The barrier \(2\lambda/(w-3k)<2\) gives
\(\lambda-F_\lambda'/F_\lambda\leq
\lambda H_0^2/[k(w-3k)]\).
Strict bulk convexity makes the metric source nonnegative along the fixed-point
solution, hence \(R'\leq0\) and \(R\leq1\). Thus \(H_0^2=R^2H^2\leq H^2\),
and the chosen constant is \(C_H=\lambda/[k(w-3k)]=42/5\).
Candidate functions in the contraction proof still use the original absolute
metric bounds; the sign improvement is only used after solving the equations.

## Arithmetic evidence and model separation

The producer evaluates polynomial derivative bounds by
\(\sup_{|x|\leq r}|\sum a_ix^i|\leq\sum|a_i|r^i\).
The rational outer bounds for \(J\), \(\coth(1)\), and \(H_0(1/k)^2\) are
the same conservative bounds as in the accepted calculator. The sharp linear
endpoint constant needs no transcendental computation.

The polynomial curvature refinement uses
\(E=\min(L,|\alpha|+K_\eta\bar\delta)\), the exact identity
\(H^2=\delta A(\beta)+\delta^2B_0(\beta)\), and every polynomial coefficient.
There is no numerical fit or unbounded series truncation.

At \(s=k\) and the former root bound and endpoint constant, the new formulas
exactly reproduce every stored constant and polynomial in both accepted model
fixtures, including the simple polynomial model. This regression is useful
evidence of consistent propagation; it is not the independent review.

Two further certificates treat historical binary64 inputs as distinct exact
rational models. One fixes only the binary64 coupling while keeping
\(k=1/9\) and mathematical \(f_0=1+c\); the other uses the stored binary64
values of \(k,c,f_0\) as rational coefficients. Both satisfy the same range
and strict-sign requirements. Neither models the rounding of the historical
ODE solver during its trajectory, nor bounds shooting error.

## Public capsule files and replay

The portable producer is `bounds/evaluate_rescaled_bounds.py`; its frozen result
is `bounds/RESCALED_BOUND_CERTIFICATE.json`. The public regression fixture is
`bounds/public_baseline/EXACT_UNIFORM_CONSTANTS.json`. No archive checkout is a
runtime dependency.

`review/PORTABLE_EQUIVALENCE.json` records identical mathematical function ASTs
and mathematical output between the original working producer and this portable
version. Their changed provenance and review-status metadata are listed separately.
The original working producer is not needed to replay the public capsule.

Run `python -B reproduce.py --out NEW_EXTERNAL_DIRECTORY` from the capsule root.
The runner verifies all frozen hashes and runs the producer, independent checker,
and 36 symbolic checks both normally and under Python optimization. It tests wrong
proof pins and all six corrupt-certificate controls, and leaves the source capsule
unchanged. `--figures` additionally regenerates the optional visualization when
Matplotlib is installed; plotting does not affect mathematical acceptance.

Direct producer use is supported with `python -B bounds/evaluate_rescaled_bounds.py
--proof theory/JACOBI_ENDPOINT_MONOTONICITY.md --proof-sha256 THE_PROOF_SHA256
--out RESULT.json`. The proof pin checks document bytes. The separate analytical
and independent arithmetic review receipts provide the review evidence.
