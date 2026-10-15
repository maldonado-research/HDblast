# Exact shell-response algebra audit

The conjugation, dynamic boundary identity, characteristic system and observable normalization in the frozen source agree exactly. The physical normal trace after conjugation is

\[
R(t)=w_z(t,0)-\frac{C_b}{2}q(t)
     =\mu\bigl(q''(t)+k(t)q(t)\bigr),
\qquad k(t)=\frac{35}{4}\operatorname{sech}^2t-\frac{25}{4}.
\]

The response coefficient is therefore

\[
\boxed{\mathcal C_{T_1}=\frac{3\pi^2}{\rho_b^2p_b^2}
\int_0^{T_1}R(t)^2\,dt.}
\]

There is no remaining hyperbolic weight. These are exact algebraic conclusions, conditional on the source equations and smoothness assumptions. They supply neither a numerical background nor a numerical response or analytic proof certification.

The mathematical input is exclusively `/workspace/hdblast-research-work/shell-response-20261015/source/HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015/PHYSICAL_SHELL_RESPONSE.md`, especially equations (7), (9)–(13), and (19)–(20). The runnable checker is `check_shell_response_algebra.py`; its exact results are in `symbolic_check_results.json`.

## Conjugation and boundary normalization

Put \(h=3\tanh t\), \(D=8\operatorname{sech}^2t-4\), \(f=\cosh^{-3/2}t\), and \(V=\alpha fw\). Since \(dy/dz=-\rho\),

\[
\frac{\alpha_z}{\alpha}
=\frac32\frac{\rho_z}{\rho}+\frac{p_z}{p}
=-\rho\left(\frac32a+\frac{p_y}{p}\right)=-\frac C2.
\]

The first derivative coefficients in the conjugated equation are \(2f'/f+h=0\) and \(-2\alpha_z/\alpha-C=0\). Its zeroth order coefficient is

\[
\mathcal W=k+\frac12C_z+\frac14C^2+\frac23\rho^2p^2.
\]

At \(z=0\), \(\alpha_b=1\) and \(V_z=f(w_z-C_bq/2)=fR\). Equation (11) consequently becomes

\[
q''+kq=\mu^{-1}R,
\qquad -w_z+\mu q''+\beta q=0,
\qquad \beta=\mu k+C_b/2.
\]

The scalar trace is \(S=3fR/(\rho_b^3p_b)\). Substitution in equation (20) cancels \(\cosh^3t\,f^2=1\), giving the boxed formula. The checker also verifies \(\int_{S^3}(n_1^2-n_2^2)^2d\Omega=\pi^2/3\) using the exact unit-sphere moments, and checks that the quadratic observable coefficient depends only on the first scalar variation and background measure.

## Residual system with coefficient discrepancies

For any smooth candidate \(\widehat w\), let

\[
\widehat w_{tt}-\widehat w_{zz}+\widehat{\mathcal W}\widehat w=\widehat f,
\qquad
-\widehat w_z+\widehat\mu\,\widehat q''+\widehat\beta\,\widehat q=\widehat b.
\]

Here \(\widehat q=\widehat w(t,0)\). Define the residuals in the true equations by

\[
F=\widehat f+(\mathcal W-\widehat{\mathcal W})\widehat w,
\qquad
G=\widehat b+(\mu-\widehat\mu)\widehat q''
 +(\beta-\widehat\beta)\widehat q.
\]

For the error orientation \(e=w-\widehat w\), put \(q_e=e(t,0)\), \(v_e=q_e'\), and \(j_\pm^e=e_t\pm e_z\). Then

\[
\begin{aligned}
(\partial_t-\partial_z)j_+^e&=-\mathcal W e-F,\\
(\partial_t+\partial_z)j_-^e&=-\mathcal W e-F,\\
\mu v_e'+v_e+\beta q_e&=j_+^e(t,0)-G,\\
q_e'&=v_e,\\
j_-^e(t,0)&=2v_e-j_+^e(t,0).
\end{aligned}
\]

Initial errors are the differences of the initial data. Reconstruction uses \(e_t=(j_+^e+j_-^e)/2\), with the compatible spatial identity \(e_z=(j_+^e-j_-^e)/2\). The observable error has two equivalent forms:

\[
\boxed{R_e=j_+^e(t,0)-v_e-\frac{C_b}{2}q_e
=\mu(q_e''+kq_e)+G.}
\]

The first form permits direct control from characteristic and boundary state errors. The sign of \(G\) in the second form follows from \(-e_z+\mu q_e''+\beta q_e=-G\).

If \(\widehat R_*=\widehat w_z-C_b\widehat q/2\), and \(\widehat{\mathcal C}_*=K\int\widehat R_*^2dt\), then with \(K=3\pi^2/(\rho_b^2p_b^2)\),

\[
\mathcal C-\widehat{\mathcal C}_*
=K\int(2\widehat R_* R_e+R_e^2)dt,
\quad
|\mathcal C-\widehat{\mathcal C}_*|
\le K\bigl(2\|\widehat R_*\|_{L^2}\|R_e\|_{L^2}
+\|R_e\|_{L^2}^2\bigr).
\]

This last inequality is the usual Cauchy–Schwarz consequence of the exact difference identity; it is not itself a computed enclosure.

For a coherent normalization, write \(\gamma=C_b/2\), \(\beta=\mu k+\gamma\), and define

\[
\widehat R=\widehat w_z-\widehat\gamma\widehat q,
\qquad
d=\widehat q''+k\widehat q-\widehat R/\widehat\mu,
\quad
\delta\mu=\mu-\widehat\mu,
\quad \delta\gamma=\gamma-\widehat\gamma.
\]

Then the raw true-boundary residual and signed error forcing are

\[
G=\mu d+\frac{\delta\mu}{\widehat\mu}\widehat R
+\delta\gamma\widehat q,
\qquad
g=-G=-\mu d-\frac{\delta\mu}{\widehat\mu}\widehat R
-\delta\gamma\widehat q.
\]

In particular, the error ODE takes the convention \(\mu v_e'+v_e+\beta q_e=j_+^e+g\). The trace entering the preceding observable comparison uses exact \(\gamma\), so \(\widehat R_*=\widehat R-\delta\gamma\widehat q\). The normalization here uses exact \(k\); any independently approximated \(k\) must supply its further coefficient discrepancy as well. The checker verifies this normalized identity, coefficient conversion and sign directly.

## Auxiliary high-frequency control

The following exact family is an artificial scalar-wave control with \(\mathcal W=\beta=\gamma=k=0\), not a substitute physical background for the frozen model:

\[
u_\omega(t,z)=\left(\frac{\cos(\omega z)}{\mu\omega^2}
-\frac{\sin(\omega z)}{\omega}\right)\cos(\omega t).
\]

It satisfies \(u_{tt}-u_{zz}=0\) and \(-u_z(t,0)+\mu u_{tt}(t,0)=0\). Nevertheless,

\[
\|u_\omega\|_\infty\le \frac1{\mu\omega^2}+\frac1\omega\longrightarrow0,
\qquad R_\omega(t)=u_z(t,0)=-\cos(\omega t).
\]

For positive integer frequencies \(n\), \(\int_0^{2\pi}R_n^2dt=\pi\). Thus no fixed constant can bound the derivative trace solely by the field supremum in this auxiliary class. This tests the need to control characteristic/derivative errors. It makes no assertion that the physical model has these zero coefficients. The checker verifies the two residuals, trace, amplitude-bound limit, and exact integer-frequency integral.

## Wrong controls rejected by the checker

The checker passes 35 exact identities and rejects 9 deliberately wrong formulas with exact, nonzero witnesses:

| Wrong formula | Exact generic defect or distinguishing feature |
|---|---|
| Conjugate with \(\cosh^{+3/2}t\) | Leaves \(6\tanh t\,w_t\) |
| Use the positive \(3/2\) power in \(S\) | Leaves a \(\cosh^6t\) factor in the observable |
| Use \(q=w(t,0)\) for the normal trace | \(q-R\) is not identically zero |
| Use bare \(w_z\) for \(R\) | Misses \(C_bq/2\) |
| Omit boundary coefficient discrepancies | Misses \((\mu-\widehat\mu)\widehat q''+(\beta-\widehat\beta)\widehat q\) |
| Include the \(\beta\) discrepancy but omit the \(\mu\) discrepancy | Misses \((\mu-\widehat\mu)\widehat q''\) |
| Replace \(+\mu q''\) by \(-\mu q''\) in the trace identity | Defect \(-2\mu q''\) |
| Reverse the acceleration sign in the boundary row | Defect \(-2\mu q''\) |
| Use \(-G\) in the error trace identity | Defect \(2G\) |

Witnesses only refute universal algebraic formulas; they do not assert a physical background realization. Run with `python check_shell_response_algebra.py`. SymPy 1.14.0 was used, with rational coefficients and no floating-point arithmetic.
