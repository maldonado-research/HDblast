# An analytic candidate using C⁴ profiles with the inherited C¹ target

This is a distinct candidate from the explicit smoothed-P1 displacement
construction. It is defined by **exact integrals** of the frozen nominal
rational displacement. Companion quadrature checks enclose one nominal
initial component and one free-wave component at selected points; they do not
evaluate the complete candidate. Its conditional C¹ transfer bound remains below the earlier
`10^29` target. The construction does not yield a quantitative C² error or a
pointwise Einstein residual, and its smoothness does not prove the existence
of an exact nonlinear family.

## Inputs and extension

Let `u0` be the expression-defined nominal initial displacement of the frozen
full23 response, and let `vP1` be its continuous P1 initial velocity. The
frozen initial report states that u0 is C¹ on `0<=r<=7/4`, with piecewise
second derivatives and matching first derivatives at all branch joins. Its
`nominal_second_derivative_H_caps[0]` is the exact H-norm upper bound

    M2 = 170416576831971060175463023195052821556359203195813205501551168471419
         / 1298074214633706907132624082305024 < 1.312842e35.

The report forms this number from the maximum componentwise absolute `U_rr`
bound over every cell of the registered footprint, then the positive
Gram-weighted H norm. At joins the second derivative may have one-sided jumps;
the global C¹ statement still makes `u0'` H-Lipschitz with constant M2 by
integrating over the finitely many pieces. This fact relies on the inherited
report and its assumptions; the checker authenticates its bytes, not a new
proof of the upstream enclosure.

Extend u0 to r>7/4 by its affine tangent at 7/4, and to r<0 by the inherited
component parity. It is zero on a neighborhood of r=0, so reflection is smooth
there. The affine continuation has zero second derivative and matches u0 and
u0' at 7/4. Consequently this full-line extension has the same H-Lipschitz
constant M2 for its first derivative. Extend vP1 by inherited parity at 0 and
affine continuation of its final P1 cell past 7/4. Both are zero through
r=65/64 on the original positive footprint.

Let the even normalized kernel be

    kδ(s) = 35/(32δ) [1-(s/δ)^2]^3  for |s|<=δ, and 0 otherwise.

Use two widths: `δu=2^-30` for the exact u0 and `δv=1/2048` for vP1. Define

    u_s = kδu*u0,       v_s = kδv*vP1,       du_s = (u_s)'.

The convolution is an exact mathematical integral over the compact kernel
support. Kernel C² regularity and u0 C¹ with bounded piecewise second
derivative imply u_s is C⁴; the corresponding P1 hinge convolution is C⁴.
Parity is preserved because k is even. Both profiles and their derivatives
still vanish for |r|<=1. Build `B0_s` by the reflected free-wave formula using
`u_s,v_s`, and build the frozen constant-matrix feedback `B1_s=K(A*B0_s,x
+B*B0_s,r+D*B0_s)+L(M P_N B0_s|shell)` with the same 5 Dirichlet/18 Robin
reflections. Here u_s and its own derivative are used consistently. The
zero initial shell collar and the inherited half-line construction are
required at corners. For odd rows the reflected coupled forcing need not be
globally smooth through r=0; the shell derivatives are understood one-sided.
Near x=r, branch-specific terms vanish with their normalized primitives because
the smoothed profiles vanish for |r-x|<1. No perturbed family or new shell law
is asserted.
The resulting nominal candidate is `V_s=B0_s+B1_s`. This is a new candidate,
not the P1-smoothed displacement used in the separate, looser note.

## Exact smoothing inequalities

Direct rational integration gives

    integral kδ = 1,     integral s*kδ(s) ds = 0,
    integral |s| kδ(s) ds = 35δ/128,
    integral s² kδ(s) ds = δ²/9.

Taylor's integral remainder and the Lipschitz bound on u0' give, globally,

    ||u_s-u0||_H <= M2 δu²/18 =: a_u,
    ||u_s'-u0'||_H <= 35 M2 δu/128 =: b_u.

This does not require u0'' to exist at the cutoff joins. For the P1 velocity,
its actual node grid has minimum spacing 1/256. Distinct δv-bands around
its interior knots do not overlap. At a single hinge, positivity and symmetry
of the kernel give `sup |kδv*z_+ - z_+|=35δv/256`. Hence

    ||v_s-vP1||_H <= (35δv/256) max_knots ||Δ slope(vP1)||_H =: a_v.

The checker computes the maximum using all 188 interior hinges, exact rational
node values, positive frozen H weights, and an outward integer square root.
It confirms the zero collar, disjoint bands, all kernel moments and provenance.

## Transfer of the existing full23 C¹ certificate

The frozen candidate `V_old=B0_old+B1_old` uses exact u0 and vP1 in B0, and
independent P1 tables p_u,p_du,vP1 in B1. Its inherited error against the
exact total formal response U is

    X(U-V_old) <= E_old ≈ 1.1787044291e28.

Let epsU and epsUr be the prior certified H sup interpolation differences
of p_u and p_du from u0 and u0', respectively, and let nA,nB,nD,nM be the
frozen H-operator upper bounds of the constant feedback matrices. Then

    ||u_s-p_u|| <= epsU+a_u,
    ||du_s-p_du|| <= epsUr+b_u,
    ||v_s-vP1|| <= a_v.

The reflected free-wave formula gives

    free = max((a_u+(3/2)a_v)/3, b_u+a_v).

For B1 use the Duhamel bound `X(K f)<=(3/2)sup||f||`, the boundary lift
bound `X(L g)<=sup||g||`, and the H-isometric parity reflection. The shell
history satisfies x<=5/4 and the profile differences vanish for x<1; hence
its velocity integral contributes at most a_v/4. Therefore

    feedback = (3/2)[(nA+nB)(epsUr+b_u+a_v)
                    +nD(epsU+a_u+(3/2)a_v)]
               +nM(epsU+a_u+a_v/4),
    X(U-V_s) <= E_old + free + feedback < 2.184e28 < 1e29.       (1)

The exact checker `checks/check_nominal_transfer.py --repo-root /path/to/HDblast`
authenticates the frozen profile, initial and response reports and reproduces
(1) in exact rational arithmetic; Python normal and optimized output agree.
The inherited first-field, background, modal operator, shell, and C¹ error
premises remain in force. This generous upper bound does not imply that the
actual approximation error is as large as the bound, and it is much too broad
to infer a nonzero physical signal.

This is a proof of an **analytic candidate and a conditional C¹ norm bound**.
Separate checks enclose one nominal initial component and one free-wave
component at selected points; the complete initial jets and feedback have not
been evaluated.
No claim is made that all higher-derivative sources are bounded or that a C²
residual is certified. A complete verified quadrature or exact integration
implementation remains separate work.
