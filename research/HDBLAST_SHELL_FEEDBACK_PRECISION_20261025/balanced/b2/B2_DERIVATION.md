# Exact-displacement two-feedback shell evaluation

This bounded experiment evaluates the constant-operator candidate

\[
 B_0=\operatorname{free}(u_0,v_{\rm P1}),\qquad
 Z=S B_0,\qquad V_2=B_0+S B_0+S^2B_0,
\]

where \(u_0\) is the exact registered rational displacement and \(v_{\rm P1}\)
is the frozen continuous P1 velocity. The operator \(S\) uses the frozen
constant matrices \(A_*,B_*,D_*,M_*\), the complete reflected bulk kernel,
and the Robin lift. Exact \(u_0\) and its derivative are used in both
feedback terms. The frozen P1 displacement and derivative do not occur in
this candidate.

The computation evaluates all five shell means at \(T=9/8\) and \(5/4\).
It supplies candidate intervals and quadrature errors. A separate comparison
certificate is required to relate them to the actual variable-coefficient
response, including its source, initial error, coefficient defects and full
third-iterate tail.

## Reduction to one-dimensional moments

Write \(E=\operatorname{diag}(1,-1,1,1,1)\) and \(P=(I+E)/2\). For an
initial profile, let \(U_k\) and \(V_k\) be its normalized k-fold primitives,
with \(U_0=u_0\). The velocity primitive \(V_k\) refers throughout to the
frozen P1 velocity. The reflected source of \(Z\) is
\(F(r+x)+G(r-x)\), where, for \(k\ge1\),

\[
 F_k=\tfrac12\{(A_*+B_*)(U_{k-1}+V_k)+D_*(U_k+V_{k+1})\},
\]
\[
 G_k=\tfrac12\{(B_*-A_*)(U_{k-1}-V_k)+D_*(U_k-V_{k+1})\}.
\]

The flat initial collar implies that integrating \(u_0'\) introduces no
constants. Negative arguments obey
\(U_k(-T)=(-1)^k E U_k(T)\), and likewise for \(V_k\).
The coefficient matrices act after this reflection; \(F_k,G_k\) themselves
need not have the parity of the initial profiles.

Set \(g_k=M_*P(U_k+V_{k+1})\). The first feedback at the shell is

\[
 Z(T,0)=P\{TF_1(T)-F_2(T)+\tfrac12[G_2(T)+G_2(-T)]\}-g_1(T).
\]

Define the shell, characteristic and triangle histories

\[
 \mathcal B(T)=\int_0^T Z(s,0)\,ds,\quad
 \mathcal C(T)=\int_0^T Z(s,T-s)\,ds,\quad
 \mathcal V(T)=\int_0^T\mathcal C(a)\,da.
\]

With \(H=TF_2-2F_3+[G_3(T)-G_3(-T)]/2\), direct integration of the
two reflection branches gives

\[
 \mathcal B=PH-g_2,
\]
\[
 \mathcal C=\tfrac14\{T^2F_1-TF_2+TG_2-G_3(T)-G_3(-T)+EH\}
             -\tfrac12g_2,
\]
\[
\begin{split}
 \mathcal V=\tfrac14\{&T^2F_2-3TF_3+3F_4+TG_3-2G_4(T)+G_4(-T)\\
 &+E[TF_3-3F_4+\tfrac12(G_4(T)+G_4(-T))]\}-\tfrac12g_3.
\end{split}
\]

The lift contributions follow from \(Lg(x,r)=-g_1(x-r)\) for \(x>r\)
and zero otherwise. On the causal triangle, \(Z(0,r)=0\) implies

\[
 \iint Z_x=\mathcal C,\qquad
 \iint Z_r=\mathcal C-\mathcal B,\qquad
 \iint Z=\mathcal V.
\]

Therefore the complete second feedback at the shell is

\[
 (S^2B_0)(T,0)
 =P\{(A_*+B_*)\mathcal C-B_*\mathcal B+D_*\mathcal V\}
    -M_*P\mathcal B.
\]

Projection remains outside the bulk matrix actions. The identity applies to
the frozen constant operator, and does not replace the actual coefficient
defect ledger. Only \(U_0,\ldots,U_4\) and \(V_1,\ldots,V_5\) are needed;
there is no two-dimensional numerical mesh.

## Moment evaluation and fixed budget

Only three initial displacement components are nonzero: \(N_0,J_0,w_0\).
Let \(f=u_{0,J}\) and
\(C_0=3/(2\ell\rho_{\rm nom}p_{\rm nom})\). The registered identity
\(u_{0,w}=C_0 f_r\) gives

\[
 U_{1,w}=C_0(T)f(T)-\int_0^T C_0'f,
\]
\[
 U_{k,w}=\int_0^T[p_{k-1}C_0-p_kC_0']f,\qquad
 p_k=(T-r)^{k-1}/(k-1)!,\quad k\ge2.
\]

These integrations by parts avoid requiring a seventh packet derivative.
The N and J primitives are the usual weighted integrals with weight \(p_k\).
All derivatives are derivatives of the registered nominal expressions.

The declared mesh has 256 subintervals: 16 on \([65/64,66/64]\), 112 on
\([66/64,9/8]\), and 128 on \([9/8,5/4]\). Each execution uses 257 distinct
point nodes and 128 whole-panel evaluations. All twelve displacement moment
channels for a target are batched from each profile evaluation. The second
target shares these evaluations. There is no adaptive refinement.

The unchanged nominal formulas for \(f\) and \(N_0=-u^2/3\) are extracted
by AST from the authenticated producer. Point jets have order two, with
order three at the two existing target nodes to evaluate \(f_r\) for the
scalar endpoint. Whole-panel jets have order six. For a moment integrand
\(h\), \(24\ell^4 c_4\) bounds its fourth ordinary r derivative. Composite
Simpson uses the per-panel remainder

\[
 \frac{(b-a)^5}{2880}\sup_{[a,b]}|h^{(4)}|.
\]

The cutoff is a segment endpoint; panels use one analytic branch each.
All arithmetic decisions use exact fractions, and interval operations round
outward on the \(2^{-112}\) grid. Frozen P1 velocity moments through order
five are integrated exactly, directly over their linear segments.

The report records every panel's Simpson interval, derivative enclosure and
remainder for every moment. It also records the assembled moment intervals,
all intermediate shell/history quantities, and the propagated component
quadrature bounds. The intervals in `evaluation.three_term` already include
quadrature and arithmetic uncertainty; a response ledger is added once to
those intervals.

## Candidate results

The complete J candidate intervals, rounded outward to integers, are

| T | J component of B0 + S B0 + S² B0 | Evaluation radius upper |
| --- | --- | --- |
| 9/8 | [−10729992186717241, −10729976580716336] | 7803000453 |
| 5/4 | [−275056797994296674, −275056748525388345] | 24734454164 |

The second-feedback J increments are approximately
\(+3.759247347436725\times10^{16}\) and
\(-7.770230981736475\times10^{16}\), respectively. All five means and
their complete rational intervals appear in the JSON report. The Dirichlet
L shell component is exactly zero.

These are signs of a newly evaluated candidate. This report by itself
establishes no exact-response or observable sign, nonlinear family, or
physical-origin conclusion.

## Reproduction and provenance

From this directory, pass the authenticated frozen repository explicitly:

```sh
python exact_b2.py --repo-root /path/to/frozen/repo --output NEW_NORMAL.json
python -O exact_b2.py --repo-root /path/to/frozen/repo --output NEW_OPTIMIZED.json
```

Outputs must be new files. Both recorded executions authenticate all seven
complete frozen closures before project imports and after computation. They
produce byte-identical reports:

- Producer SHA-256: `e4ea015d0526fc81a2c862cb048ecc3e5f4473bf2e0d71b55eec719e48ee2065`
- Normal and optimized report SHA-256: `3db587c30dca20f01e4964ac87f8f5ba457d55d15f7720d09eb840ed3e0e4935`

The report carries explicit hashes for the nominal producer, interval
implementation, frozen velocity/profile table, registered coefficient input,
and all capsule manifests. The reusable authenticator is the pinned
`../exact_boundary.py`. No frozen producer main routine or old numerical
campaign is executed.

Independent polynomial controls integrate the authenticated reflected kernel
exactly on both reflection branches and check the triangle reduction with
arbitrary coupled matrices. All 55 controls pass: 24 reflected-kernel
primitive cases, 18 triangle cases, nine direct producer-algebra cases and
four synthetic P1 velocity cases. Negative controls detect a missing lift
and an incorrect negative primitive argument.

The independent report checker uses a separate unrounded rational interval
implementation to reassemble every panel sum and remainder, displacement
moment (including the scalar endpoint), P1 moment through order five,
candidate/history interval, and propagated quadrature cap. It reports no
source defect. Receipts:

- `B2_EXECUTION_RECEIPT.json`: normal/optimized byte parity.
- `B2_POLYNOMIAL_CONTROLS_FINAL.json`: all 55 exact controls pass.
- `B2_INDEPENDENT_REVIEW.json`: complete independent interval reassembly
  passes; SHA-256 `1da9aa9d9e56cfb91a83284a560e38c10ef644055e206301bef2f8848d9a94ea`.

The earlier `B2_POLYNOMIAL_CONTROLS.json` records the initial 42 controls;
the final controls receipt adds the direct producer and P1 checks. Both are
retained as experiment history.
