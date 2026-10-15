# Primary-source comparison for the uniform waveform checkpoint

15 October 2026 UTC. This synthesis uses the root agent's saved results from **two paper-index searches and three targeted primary passage reads**. It adds no external search. Three sources were read in selected passages; one 2026 candidate is included only at index-metadata/abstract level. This bounded review is neither an exhaustive survey nor a novelty assessment.

## Boundary inertia and a pointwise shell trace

[Parada and Vanspranghe, *Optimal stabilization rate for the wave equation with hyperbolic boundary condition*, arXiv:2512.19167](https://arxiv.org/abs/2512.19167), study a wave equation on a smooth bounded domain of dimension at least two. Equation (1.1) couples it to a damped hyperbolic boundary equation with a boundary Laplace–Beltrami term; an optional disjoint component carries static Dirichlet data. Their reported \(1/t\) stabilization result uses a geometric control condition and high-frequency resolvent/quasimode analysis. The retrieved discussion explicitly identifies normal and tangential trace estimates as a separate part of that analysis.

This is a useful structural comparison: a boundary carrying momentum must be included in the dynamics, and derivative traces need justified estimates. Its long-time damped stabilization theorem does **not** provide the present short-time, one-dimensional pulse enclosure. We import no decay rate, trace inequality or numerical constant from it.

[*One dimensional wave equation with in-domain localized damping and Wentzell boundary conditions*, arXiv:2503.22439](https://arxiv.org/abs/2503.22439), is closer in dimension. The retrieved equation (1.1) and assumptions (A1)–(A3) specify a positive Lipschitz spatial coefficient, nonnegative damping that is positive on a nonempty interior subinterval, and endpoint oscillator variables with positive constant coefficients. The abstract reports stabilization and input-to-state methods in \(L^2\) and \(L^p\) frameworks. The returned §5 estimates, equations (5.28)–(5.30), control integrated bulk velocity/spatial-derivative quantities and endpoint state variables. The specific input-to-state theorem statement was not returned, so its exact constants and hypotheses are not reconstructed here.

The HDB model has its own time-dependent bulk potential and boundary stiffness:

\[
w_{tt}-w_{zz}+W(t,z)w=0,\qquad
-w_z(t,0)+\mu w_{tt}(t,0)+\beta(t)w(t,0)=0,
\quad\mu>0.
\]

The derivative output is \(D=w_z(t,0)-\gamma w(t,0)\). An integrated energy estimate alone is not a uniform pointwise error bound for this output. The current HDB derivations therefore use the incoming characteristic identity and the actual boundary oscillator equation. Positive \(\mu\) gives the decaying boundary kernel; integration by parts is applied to the **known bump derivative input**, retaining the boundary forcing and all other terms. Those identities and the resulting finite-time majorants are HDB proof obligations, not applications of either stabilization paper.

## Rigorous integration and the flat bump endpoints

[Johansson, *Numerical integration in arbitrary-precision ball arithmetic*, arXiv:1802.07942](https://arxiv.org/abs/1802.07942), gives a useful distinction in §2. Its algorithm can enclose an integral directly by the interval-box bound \((b-a)f([a,b])\). It uses Gaussian quadrature with a complex-ellipse remainder only when analyticity and an appropriate complex magnitude bound are verified. Section 3.2 separately discusses endpoint/improper integrals and the bounds required when truncation is used.

Our comparison profile is

\[
Q(x)=\begin{cases}
\exp\!\left(4-\frac1{(x-1)(2-x)}\right),&1<x<2,\\
0,&\text{otherwise}.
\end{cases}
\]

It is smooth and flat at its support endpoints, but it is not analytic across those endpoints. Its real analyticity inside \((1,2)\) does not justify using an analytic Gaussian remainder on a cell crossing \(1\) or \(2\).

The current HDB implementation instead uses outward integer-dyadic arithmetic, an exponential Taylor enclosure with an explicit remainder, and **real interval ranges covering every point of each cell**. Endpoint cells are split at the support boundaries. With \(v=1/[(x-1)(2-x)]\ge4\), the interior formulas are

\[
Q=e^{4-v},\qquad Q'=(3-2x)v^2e^{4-v}.
\]

The positive envelope \(v^2e^{4-v}\) decreases for \(v\ge4\) and tends to zero as \(v\to\infty\); the flat endpoints can therefore be included by real bounds without asserting a complex analytic extension. Continuous-cell integration and range bounds must retain their full outward widths. Sample values or a visually converged plot do not replace these bounds.

This describes the chosen method and its relation to the literature. It does not identify the HDB implementation with Arb or certify that code merely by citing Johansson. The implementation, elementary remainder bounds and propagation of the waveform error require their own review and replay.

## A recent comparison with limited evidence

The second search returned [*Certified and accurate computation of function space norms of deep neural networks*, arXiv:2603.06431](https://arxiv.org/abs/2603.06431). Its abstract describes local interval boxes and adaptive aggregation, contrasting them with finite point probes alone. This is a methodological comparison only: its body and authorship were not inspected, the HDB comparison profile is explicit, and this calculation uses no neural network or theorem from that work.

## Scope and provenance

The target is a conditional uniform enclosure of a specified **first-variation** waveform and its leading response coefficient. These sources do not establish a full finite-amplitude nonlinear waveform, an effective nonlinear amplitude radius, a physical time calibration, or a higher-dimensional origin of the Big Bang.

[WAVEFORM_SOURCE_LEDGER.json](WAVEFORM_SOURCE_LEDGER.json) records source-specific evidence and transfer limits. [BOUNDED_SEARCH_RECEIPT.json](BOUNDED_SEARCH_RECEIPT.json) records all five research operations. The exact invocation arguments come from a direct root-agent handoff; the saved result-body bytes have a separate SHA-256 pin. Raw paper passages remain in an internal retrieval file and are excluded from this publication-safe packet. No downloaded code was executed.
