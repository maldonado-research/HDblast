# Scalar boundary conditions and de Sitter curvature

**Research draft, 8 October 2026.** Mathematical results in a specified five-dimensional Einstein–scalar boundary problem, with internal independent reviews. External originality, APS significance, dynamical stability and a higher-dimensional origin of the Big Bang remain unestablished. This checkpoint is later GitHub work, outside the unchanged Zenodo record [23244754](https://doi.org/10.5281/zenodo.23244754).

The useful advance is an actual existence proof with a uniform error bound. The earlier scalar-curvature note assumed the required branch. Here a weighted inverse and a contraction argument construct it for sufficiently small positive tension detuning, despite the interior interval growing without bound. Uniqueness is proved within a specified weighted neighborhood, not among all possible solutions.

With two identical regular interiors, U = W′²/2 − 2W²/3, tension σ = 2W + δf, W(0) = 3k, W′(0) = 0, W″(0) = w > 4k > 0 and f(0) = f₀ > 0, the constructed branch obeys

\[
H^2=\frac{kf_0}{3}\delta+
\left[\frac{f_0^2}{36}-\frac{k f_1^2}{24(w-2k)}\right]\delta^2+R,
\qquad |R|\le K_H\delta^3.
\]

The [complete proof](branch-theorem/UNIFORM_BRANCH_THEOREM.md) supplies sufficient inequalities for the allowed interval. Its original pending-review wording is preserved; the later [22-step review](independent-math/UNIFORM_PROOF_REVIEW.md), [second adversarial review](finite-detuning-review/ADVERSARIAL_REVIEW.md), and [independent exact arithmetic checker](independent-math/INDEPENDENT_EXACT_UNIFORM_CONSTANTS.json) document completion of internal review. These are AI-assisted internal reviews, not journal peer review.

| Specified polynomial model | Guaranteed interval | Conservative curvature remainder |
|---|---|---|
| Registered exact-decimal model: k=1/9, W=1/3+η²+η³/3, f=1+c+cη, c=2/(1.0357712571566784)−4/3 interpreted exactly | 0 < δ ≤ 161839258930731/1482649651195872051200, approximately 1.0915543 × 10⁻⁷ | |R| < 335 δ³ |
| Separate simple model: k=1, W=3+5η²/2, f=1+η | 0 < δ ≤ 10⁻³ | |R| < 0.02360 δ³ using the polynomial refinement |

Exact rational values, all inequalities and parameter conventions are in [the constant certificate](numerical-audit/uniform-constants/EXACT_UNIFORM_CONSTANTS.json). Its historical pending-review field is resolved by the subsequent independently written checker cited above. Both intervals also give **strictly lower curvature than the separate constant-detuning model f=f₀**. The constant-scalar solution is an admissible comparison in that separate model; it generally fails the scalar junction of the original f₁≠0 problem. These are static comparisons, not evidence of dynamical relaxation.

The registered interval ends about 4,581 times below the smallest uploaded numerical detuning, 0.0005. Its earlier 48 integrations are not certified by this theorem. The decimal model and its binary64 numerical implementation also have distinct parameter representations. Neither distinction is hidden by the close numerical agreement.

A complementary theorem fixes any positive δ and introduces an independent weak scalar coupling ε. It proves a locally unique analytic branch and a strictly negative quadratic curvature susceptibility at every such fixed detuning. An exact w=5k example gives F(y)=cosh(ky). This weak-coupling theorem does not supply a common coupling radius down to δ=0 or an evaluated remainder at arbitrary finite coupling; see [the derivation](independent-math/INDEPENDENT_MATHEMATICAL_REVIEW.md).

## Read and reproduce

- [Focused manuscript source](manuscript/scalar-boundary-curvature.tex) and [compiled research draft](manuscript/scalar-boundary-curvature.pdf).
- [Seven-paper primary-source comparison](literature/PRIMARY_SOURCE_COMPARISON.md), including exact translations to established results.
- [Original numerical replay audit](numerical-audit/NUMERICAL_AUDIT.md): nine exact algebra checks; 24 parameter configurations at two settings, totaling 48 final integrations, all reproduced. Those are floating-point diagnostics.
- [Independent finite-detuning calculations](finite-detuning-review/FINITE_DETUNING_DIAGNOSTICS.json): 36 weak-coupling parameter configurations at two settings, totaling 72 final integrations. Stable boundary-identity estimates converge; direct geometric subtraction reaches cancellation limits in two small-response cases, retained explicitly.
- [Reproduction runner](reproduce.py), with pinned dependencies in `requirements.txt`. Run `python -I -B reproduce.py --output /tmp/hdblast-scalar-replay` for exact algebra, proof-constant checks, and arithmetic audits. Add `--numerical` to rerun the 48 original integrations and the 72 independent integrations. The output directory must be new. No network access, token or repository mutation is used by the runner.

`MANIFEST.json` binds every delivered file except itself. The `reference/` directory contains only the earlier scientific scalar note, its code and results selected from the user-supplied archive; it does not include private editorial correspondence or third-party source PDFs. Input-integrity receipts describe the original attachment inventory separately from this delivery inventory. The reproducible dependencies are Python 3.12.14, NumPy 2.2.6, SciPy 1.15.3, SymPy 1.14.0 and mpmath 1.3.0. The original numerical packet used different NumPy/SciPy versions; its inputs are preserved.

## Contribution and remaining work

Scalar junctions, curved scalar branes, regular interior caps, and effective interface actions have established precedents. The shell-curvature identity is an algebraic consequence of those established equations. The bounded comparison did not identify this exact quantified doubled-interior theorem in the selected passages; that finding does not establish priority or significance for a particular journal.

The next useful calculation is a sharper validated continuation/inverse estimate reaching the registered numerical detunings, followed by a stability analysis if making dynamical claims. An expert assessment of the theorem's originality and importance is also needed. No blast, lasting energy transfer, thermalization, observational discrimination, journal submission or acceptance is demonstrated by this checkpoint.
