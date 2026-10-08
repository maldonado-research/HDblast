# Independent adversarial mathematical review

8 October 2026. This review was performed independently from the supplied numerical producer. It is an internal AI-assisted review, not external peer review or a guarantee of novelty or APS acceptance.

## Fixed positive detuning: no substantive defect found

The fixed-positive-detuning theorem in `independent-math/INDEPENDENT_MATHEMATICAL_REVIEW.md` has a sound local proof after tightening its regular-center contraction argument. Its scope must remain a locally unique branch in the symmetric doubled-interior ansatz, with the radial gauge and cone normalization fixed. It does not assert global uniqueness, a quantified coupling neighborhood, or uniformity down to zero detuning.

Independent derivation: write `eta_b=b c+O(c^2)`. The linearized junction gives

    b = -delta g1/[2(D+w)],  y_b'(0)=0.

Expanding the exact Hamiltonian boundary identity yields

    C2 = b^2[(m^2-D^2)/12-A(2D+w)/3]
       = -delta^2 g1^2 [4A(D+w)-D']/[48(D+w)^2].

Here `D'=m^2-D^2-4AD`. This derivation agrees with the proposed coefficient. Direct differentiation also verifies that the radial constraint residual has identically zero derivative on the two evolution equations.

For the sign, `F>0,F'>0` follow from the positive-mass divergence-form equation. For `E=4A(D+w)-D'`, a first zero would have

    E' = 4[2AD^2+2ADw+6A^2D+7A^2w+k^2(2D+w)] > 0.

Since `E>0` near the cone, this excludes any zero. Thus `C2<0` for `g1!=0`. This sign is a theorem about the specified static model, not evidence about the origin of the Big Bang.

### Precise regular-center contraction

Using only `r=rho/y` and `B=phi'/y` without a weighted norm makes a blanket `O(y_*^2)` Lipschitz claim imprecise: the scalar map's dependence on `r` contains an unsuppressed factor `U'(phi_h)`. A direct correction is

    R=(rho/y-1)/y^2,  B=phi'/y,  r=1+y^2 R,
    phi(y)=phi_h+y^2 integral_0^1 u B(yu) du.

Then the fixed-point equations are

    R(y)=-integral_0^1 (1-u)u r(yu)
         [y^2 u^2 B(yu)^2/4+U(phi(yu))/6]du,
    B(y)=integral_0^1 u^4 [r(yu)/r(y)]^4 U'(phi(yu))du.

Every functional variation now has a factor `O(y_*^2)` on bounded balls with `r` separated from zero. The center values are `R(0)=-U(phi_h)/36` and `B(0)=U'(phi_h)/5`. Contraction, analytic parameter dependence and ordinary continuation on the finite positive interval give the regular solution family needed for the finite-dimensional implicit-function theorem. This correction has been communicated to, and adopted by, the proof authors.

### Explicit special case and independent representation

Let `alpha=w/k`. The regular linear solution is

    F(y)=2F1(2-alpha/2,alpha/2;5/2;-sinh(ky)^2).

For `w=5k`, this reduces exactly to `F=cosh(ky)`. Therefore

    D=k^2/A,  D'=k^2-k^4/A^2,
    C2=-delta^2 g1^2 [3k^2+20Ak+k^4/A^2]
        /[48(k^2/A+5k)^2].

This provides a closed elementary example at any fixed positive detuning. It follows from the standard hypergeometric solution and is not claimed as a new special-function identity.

## Uniform small-detuning branch: second independent review

Reviewed source: `branch-theorem/UNIFORM_BRANCH_THEOREM.md`, SHA-256

    544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51

The reviewed version includes the corrected scaled regular-center argument above. I checked the proof's substantive steps independently, including the following constants and implications:

- The positive linear comparison mode, uniform integral bound `J`, Green inverse bound `C_G`, and comparison `phi_T<=s_T`.
- The metric Volterra kernel and derivative, `2/3<=R<=4/3`, the metric logarithmic-derivative bound and its difference estimate. In particular, the latter's two contributions are bounded by `4 D_F/(3k)` and `2 D_F/(3k)`, with `D_F=2 C_F B |beta| d`, giving the stated `4 C_F B |beta| d/k`.
- The scalar contraction constants `C_N,L_N`, fixed common ball for differentiability, and the endpoint derivatives `C_P,C_a`.
- The exponential linear logarithmic-derivative estimate and conversion to the nonlinear endpoint curvature, including the retained factor `16/9` from the metric ratio.
- The monotone scalar-junction root, positive metric-junction factor, opposite signs at the finite starting radius and infinity, and existence of a shell.
- The moving-endpoint identities. Nondegeneracy follows from differentiating the endpoint identity: `eta_h(T) h_beta=1`, where `h` is the scalar cone parameter. The regular-cone IVP and Dirichlet-family differentiability make this argument legitimate.
- The estimate forcing every possible metric-junction root to be a strict downward crossing, which supplies uniqueness in the stated tube without an implicit-function theorem at an infinite endpoint.
- The scalar and curvature remainder constants and their stated powers of detuning.

**Review conclusion: no substantive mathematical defect found in the pinned proof under its stated hypotheses and sufficient inequalities.** This is a local existence and uniqueness result inside an explicitly restricted weighted neighborhood. The inequalities must be evaluated before assigning a numerical detuning range. Their conservatism is not evidence that a physical branch fails outside that sufficient range. This review does not establish external novelty, dynamical stability, reheating, or observational support.

## Independent finite-detuning numerical diagnostics

`verify_finite_detuning.py` implements the equations independently; it imports no attached or companion producer. `FINITE_DETUNING_DIAGNOSTICS.json` records all inputs, runtime versions and outputs. Six polynomial models vary `k`, `w/k`, `f0`, `delta/k`, the cubic superpotential derivative and the quadratic scalar coupling. Each has positive and negative couplings `0.04,0.02,0.01` at two solver settings: **72 final integrations**. These are distinct from the attachment's 48 configurations.

The logarithmic derivative from 80-digit hypergeometric evaluation and an independent Riccati integration differs by at most `6.62e-16`. All six predicted quadratic susceptibilities are negative. The 72 final integrations have maximum junction residual `6.22e-15`, maximum sampled relative constraint residual `1.45e-15`, maximum absolute shell-identity discrepancy `3.91e-14`, and maximum geometric curvature difference between settings `3.38e-14`.

Centered extraction cancels odd coupling terms. Direct subtraction of geometric curvatures reaches a floating-point cancellation floor in two small-response cases: those two coefficient errors do **not** decrease monotonically. This outcome is retained in the JSON. Evaluating the same exact constraint identity as a factored difference from the zero-coupling background avoids that subtraction. In that stable evaluation all six centered coefficient errors decrease; consecutive ratios are approximately `3.96` to `4.03`, consistent with a quadratic error after centering. At `|c|=0.01`, the largest relative coefficient discrepancy is `0.001291` (about `0.1291%`) in the largest-detuning model.

These are floating-point consistency checks. Neither the small residuals nor the coupling convergence certify an inverse bound, prove the branch's numerical location, evaluate its analytic coupling radius, or establish an error enclosure between samples. The proof is logically separate from these diagnostics.
