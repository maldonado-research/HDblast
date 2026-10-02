# Post hoc appendix: finite-cutoff trace and UV tails

These identities were derived after the registered calculation and independently checked by two internal reviewers. They do not change its acceptance criteria. No external mathematical novelty is claimed.

For a minimally coupled real scalar on Minkowski spacetime,
\[
\ddot f_k+(k^2+x)f_k=0,\quad f_k\dot f_k^*-\dot f_kf_k^*=i.
\]
Write \(x=M^2\), \(r>0\), \(h=x-r\), \(\Omega_r=\sqrt{k^2+r}\), and \(d\mu_k=k^2dk/(2\pi^2)\).
For \(q_k=|f_k|^2\), the mode energy and pressure obey
\[
e_k=\tfrac14\ddot q_k+(k^2+x)q_k,\quad
p_k=\tfrac14\ddot q_k+\tfrac13 k^2q_k,\quad
-e_k+3p_k=\tfrac12\ddot q_k-xq_k.
\]

The direct prescription's fixed-reference integrands are
\[
q_{R,k}=q_k-\frac1{2\Omega_r}+\frac{h}{4\Omega_r^3},
\]
\[
e_{R,k}=e_k-\frac{\Omega_r}{2}-\frac{h}{4\Omega_r}+\frac{h^2}{16\Omega_r^3},
\]
\[
p_{R,k}=p_k-\frac{k^2}{6\Omega_r}+\frac{hk^2}{12\Omega_r^3}
-\frac{h^2k^2}{16\Omega_r^5}
+\frac{\ddot x(2k^2/3+r)}{16\Omega_r^5}.
\]
Their integrals from 0 to K define \(Q_K,\rho_K,p_K\). The mass-squared current is \(S_K=Q_K/2\).
The code uses rationalized versions of these expressions to reduce cancellation.

Substitution leaves the local remainder
\[
-e_{R,k}+3p_{R,k}-\tfrac12\ddot q_{R,k}+xq_{R,k}
=\frac{r[3h^2+\ddot x]}{16(k^2+r)^{5/2}}.
\]
Because
\[
\int_0^K\frac{k^2dk}{(k^2+r)^{5/2}}
=\frac1{3r}\left(\frac K{\sqrt{K^2+r}}\right)^3,
\]
the exact finite-K identity is
\[
-\rho_K+3p_K=\tfrac12\ddot Q_K-xQ_K+
v_r^3\left[\frac{h^2}{32\pi^2}+\frac{\ddot x}{96\pi^2}\right],
\quad v_r=\frac K{\sqrt{K^2+r}}.
\]
Using the infinite-cutoff offset without \(v_r^3\) leaves a known cutoff residual. Numerical differentiation adds a separate error.

For the UV expansion, the exact mode equation implies
\[
q_k'''+4(k^2+x)q_k'+2\dot xq_k=0.
\]
The declared incoming Jost state's vacuum UV expansion is
\[
q_k=\frac1{2k}-\frac{x}{4k^3}+\frac{3x^2+\ddot x}{16k^5}
-\frac{10x^3+10x\ddot x+5\dot x^2+x^{(4)}}{64k^7}+\cdots.
\]
After matching, the leading integrands are
\[
e_{R,k}=\frac{\dot x^2+2h^3}{64k^5}+O(k^{-7}),\qquad
q_{R,k}=\frac{\ddot x+3h^2}{16k^5}+O(k^{-7}),
\]
\[
p_{R,k}=\frac{2x^{(4)}+8h\ddot x+13\dot x^2-10h^3}{192k^5}+O(k^{-7}).
\]
Integrating the omitted band gives
\[
\rho_\infty-\rho_K=\frac{\dot x^2+2h^3}{256\pi^2K^2}+O(K^{-4}),
\]
\[
S_\infty-S_K=\frac{\ddot x+3h^2}{128\pi^2K^2}+O(K^{-4}),
\]
\[
p_\infty-p_K=\frac{2x^{(4)}+8h\ddot x+13\dot x^2-10h^3}{768\pi^2K^2}+O(K^{-4}).
\]
For the registered pulse at the crossing, \(x=\dot x=0,\ddot x=8,x^{(4)}=-64\). Leading tails are
\[
\rho_{\rm tail}=-\frac1{2\pi^2K^2},\quad
S_{\rm tail}=\frac7{16\pi^2K^2},\quad p_{\rm tail}=\frac1{3\pi^2K^2}.
\]

The expansion assumes a sufficiently smooth background and the stated vacuum UV behavior, including its time derivatives. This applies to the smooth pulse's specified incoming Jost state; it is not automatically valid for arbitrary finite-time instantaneous vacua. The next asymptotic order is not a certified numerical remainder bound.

## Executed diagnostics

A nine-point finite difference on the 2001 saved primary rows gives maximum absolute trace residual 1.67289355e-9 with the exact finite-K factor (1.04555847e-10 after division by m_infinity^4). Root and the adversarial reviewer independently executed the same check.

| K | Leading-tail-corrected rho at crossing | Corrected p | Corrected S |
|---|---:|---:|---:|
| 64 | -0.009677707870 | +0.006447822038 | 0.013662263071 |
| 128 | -0.009677697087 | +0.006447795345 | 0.013662256355 |
| 192 | -0.009677696510 | +0.006447793917 | 0.013662255995 |

These are post hoc asymptotic estimates from the saved first-run values, not rigorous continuum enclosures. See outputs/posthoc_tail_checks.json. The finite-Lambda PV sources retain their separate regulator dependence.

