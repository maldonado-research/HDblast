# Independent adversarial review of the uniform branch proof

8 October 2026. Internal AI-assisted mathematical review. This is neither external peer review nor a claim of priority or APS suitability. Reviewed proof identity is recorded in `UNIFORM_PROOF_REVIEW.json`; that JSON pins the exact file to which this review applies.

## Verdict

The argument constructs the regular branch within its specified weighted neighborhood and gives a uniform small-detuning scalar-shift and curvature remainder. The proof is stronger than the attachment's assumed-branch algebra. I found no incorrect constant or failed implication after checking each major estimate independently. This acceptance concerns the written mathematical argument. It does not turn the uploaded floating-point sweeps into certified finite-detuning numerical solutions and does not establish that their detunings satisfy the theorem's conservative sufficient conditions.

## Independently checked chain

1. The regular positive-mass mode is positive and increasing by its divergence-form equation. Its logarithmic derivative lies between zero and its asymptotic exponent by a first-crossing argument.
2. The lower-mass normalized weight dominates the normalized target-mass mode. Both the Riccati comparison and the maximum-principle alternative give the same direction of inequality.
3. The lower bound on the weight's logarithmic derivative after `2/k` follows by restricting its integral formula to the immediately preceding interval of length `1/k`. The constants `C0`, `zeta` and `J` are valid. The integral of the square of the normalized weight is uniformly bounded as the endpoint tends to infinity.
4. The Dirichlet Green inverse maximum-principle estimate uses `L_m s_T = -D s_T` with `D>0`. The derivative estimate correctly uses the fourth power of the AdS radial ratio and the factor `1/(4k)`; the derivative normalization in the Banach norm contributes the additional factor `1/k` in `C_G`.
5. The metric Volterra kernel is exactly the ratio obtained from the variation-of-parameters Green function for `rho''-k^2 rho=-F rho`. Its supremum is `1/(2k)` and its endpoint derivative is the square of the radial ratio.
6. The bound on `F` and its difference follows from `U'(0)=0`, the bounded second derivative of `U`, and the scalar derivative component of the norm. The proposed `C_F` has the correct factors `1/4` and `1/12`.
7. Under `epsilon<=1/4`, the Neumann inverse gives `|R-1|<=1/3`, hence `2/3<=R<=4/3`. No sign of `F` is assumed.
8. For the metric Lipschitz estimate, writing `L_F=2 C_F B |beta|`, I independently obtain `||Delta R|| <= (16/9)(J/(2k)) L_F d`. The `Delta R'/R` term contributes at most `(4/3)L_F d/k` and the `R' Delta R/(R R_tilde)` term at most `(2/3)L_F d/k`, after using `epsilon<=1/4`. Their sum is exactly the stated sufficient bound `2 L_F d/k` times the weight squared.
9. The scalar nonlinear remainder uses `|U'(eta)-m^2 eta|<=M3 eta^2/2`. The coefficients `4` and `20` in `C_N` and `L_N` follow respectively from one cubic metric term and its two product differences. The map preserves the weighted ball and contracts under the listed inequalities.
10. The common fixed ball of radius `Bb` gives parameter differentiability at `beta=0`; this avoids differentiating a parameter-dependent ball. Differentiating the fixed point yields the stated `C_P` and `C_a` estimates.
11. The Hamiltonian constraint derivative cancels exactly under the second-order equations. Its regular-cone value is zero, so it holds on the constructed solution.
12. The scalar logarithmic-derivative error satisfies `e'=-(w+z)e+4(a0-k)z`. The proposed integration bound `C_L exp(-2kT)` is correct. The geometric conversion to `H^2` correctly retains `R^2<=16/9`; omitting this factor without a sign theorem for `F` would be unjustified.
13. The scalar junction has derivative at least `w/2` under the stated smallness conditions. The root interval `|beta|<=|f1|delta/w` is a direct mean-value/IVT consequence, including the case `f1=0`.
14. The exact reduced curvature `K=delta A+delta^2 B_fun` is within `C_K delta^2` of `a1 delta`. The endpoint metric inequalities are in the correct direction: `9 H0^2/16 <= H^2 <= 9 H0^2/4`.
15. Metric positivity and positive tension make the signs of the metric residual and `H^2-K` coincide. Their endpoint signs prove existence without applying an implicit-function theorem at the infinite-radius limit.
16. The moving-endpoint identities at fixed boundary scalar follow by differentiating the unique regular cone IVP and eliminating the center-parameter derivative. They are exactly `p_T=U'-4ap-p p_beta` and `a_T=-H^2-p^2/3-p a_beta`.
17. The cone-to-boundary map's center-parameter derivative cannot vanish on the constructed branch: a zero derivative would make the entire regular linearized IVP variation vanish, contrary to `eta_beta(T)=1`. This supplies the required nondegeneracy.
18. The stated bounds `C_Tp`, `C_Tbeta` and `C_slope` control the endpoint derivatives. At every possible root the metric residual derivative is strictly negative. Two distinct roots would force an intermediate non-downward crossing, so uniqueness follows within the claimed tube.
19. The scalar remainder `K_eta` includes all contributions: finite-radius linear mismatch, nonlinear momentum, nonlinear `W'`, and detuned variation of `f'`.
20. Taylor expansion of the exact shell identity gives `K_H=|A'(0)|K_eta+(sup|A''|/2)L^2+sup|B_fun'|L`; the second-order coefficient reduces to the supplied susceptibility. This remainder argument uses bounded derivatives in field space and does not assume a third derivative of the branch with respect to delta.
21. The regularity hypotheses `W in C^4`, `f in C^3` supply the needed bounded `U'''`, `A''` and `B_fun'` and local parameter differentiability.
22. The sufficient inequalities are jointly feasible by choosing `b` first, including `(C_P+W3)b<w/2`, and then choosing positive `delta_bar` sufficiently small. Every problematic left-hand side tends to zero while denominators remain strictly positive for fixed `k>0`, `w>4k`, `f0>0`.

## Scope that must remain explicit

- Uniqueness is in the defined weighted small-field tube and fixed cone/proper-distance gauge. Large-field branches and different topology or normal conventions are not excluded.
- Uniformity is in the endpoint radius and sufficiently small positive detuning for fixed model constants. Uniformity as `k->0`, `w->4k`, unbounded potential derivatives or a changing function class is not established.
- The proof gives formulas for sufficient constants. Concrete numeric coverage requires certified upper/lower bounds for all model constants and inequalities; approximate decimal evaluation is not a certified threshold.
- Neither the theorem nor its review establishes dynamical stability, a physical energy scale, thermalization, observable discrimination, a Big Bang mechanism, originality, or editorial acceptance.
- A separate finite-positive-detuning weak-coupling theorem is available, but its `c->0` expansion is not interchangeable with this theorem's fixed-coupling `delta->0` statement.

## Review refinements

The initial argument's constants passed review. The author added explicit descriptions of the regular `C^1` Banach space, a common contraction ball for parameter differentiation, normalized-mode comparison, and regular-cone IVP uniqueness. A further origin-contraction clarification uses scaled variables `(rho/y-1)/y^2` and `p/y`; these avoid claiming a small Lipschitz constant in raw variables with an unscaled cross dependence. These refinements do not alter the weighted global estimates or their constants.
