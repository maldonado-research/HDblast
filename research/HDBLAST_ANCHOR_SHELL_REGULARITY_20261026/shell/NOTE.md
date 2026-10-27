# A C² compatibility condition for the reflected shell lift

This note concerns only the affine lift used for the formal retained
23-amplitude second-order problem. It does not establish a pointwise Einstein
residual, a full physical second-order solution, or a finite-amplitude family.
The source and solution must still satisfy the frozen bulk, initial, boundary,
and subsidiary assumptions.

## Conventions and statement

Use the frozen variables x=t/ell, r=z/ell on the half-line r≥0, with
ell=10^-6, and the source convention of second parameter derivatives (the
second metric coefficient is multiplied by epsilon²/2). Let

    (Lg)(x,r) = - integral_0^((x-r)_+) g(s) ds,    (q)_+ = max(q,0).

This is the sign and normalization in the frozen
`FULL23_SHELL_AND_SUBSIDIARY.md`, section “C1 response needs source values,
not shell-source derivatives”. Everything below applies componentwise, or in
any fixed finite-dimensional norm with the usual triangle inequalities.
Assume g∈C¹([0,T]) and g(0)=0. On x>r,

    (Lg)_x  = -g(x-r),       (Lg)_r  = +g(x-r),
    (Lg)_xx = -g'(x-r),     (Lg)_xr = +g'(x-r),
    (Lg)_rr = -g'(x-r),

and all of these derivatives vanish on x<r. The physical-side one-sided
values are intended at the domain boundary. Continuity of the first
partial derivatives across x=r follows from g(0)=0. The lift is C² across
any interior point of x=r **if and only if g'(0)=0**: the physical-side
limits of the three second derivatives are respectively -g'(0), +g'(0),
-g'(0), while the zero side has limit zero. With the two zero conditions,
the displayed derivatives extend continuously, proving the sufficiency as
well as necessity. Here C² is meant in the ordinary local interior sense;
at r=0 and x=0 use the natural one-sided derivatives of the half-domain.

For example, g(s)=s has g(0)=0 and produces the C¹ function
`-(x-r)_+²/2`, but its second derivatives jump across x=r. Thus continuity
of g, even with g(0)=0, is not enough for a C² lift. The frozen initial
packet and second corrections are supported in 1≤r≤2 with flat zero
extension, and the exact solution/source has a zero shell collar. That is a
qualitative compatibility input; an approximate source must retain the
collar or separately certify its corner defect.

## C² norm and total boundary source

For two sources g, g̃ both satisfying the above corner conditions, write
δg=g-g̃. On 0≤x≤T, r≥0,

    sup |L δg| ≤ T sup |δg|,
    max(sup |(L δg)_x|, sup |(L δg)_r|) ≤ sup |δg|,
    max(sup |(L δg)_xx|, sup |(L δg)_xr|,
        sup |(L δg)_rr|) ≤ sup |δg'|.

This is a coordinate-component bound for the *lift*, not a complete
23-field error bound; other Duhamel terms and homogeneous Robin feedback
must also be controlled. The inherited boundary equation is

    P_N U_r = M P_N U + G,       M = -ell rho B_y,

with fixed background shell coefficients. Consequently the source for the
reflected lift is the total g=M P_N U|shell+G (with the inherited packing and
sign conventions), not just its affine part G. Where the trace is
regular and M is constant,

    g_x = M P_N (U_x)|shell + G_x.

For two solutions, the corresponding derivative error includes the
trace of δU_x and δG_x, plus any separately allowed matrix approximation
error. In particular, a C⁰ bound for δG alone supplies no numerical bound
for δG_x. A genuinely C² calculation must control the differentiated bulk
source, feedback, and corner terms as well as this lift term; none is
computed numerically here.

The source-only second-jet estimate is scaled-coordinate based. For the
frozen relations partial_t=ell^-1 partial_x and
partial_y=-(ell rho)^-1 partial_r, where rho=rho(y),

    F_tt = ell^-2 F_xx,
    F_ty = -(ell² rho)^-1 F_xr,
    F_yy = (ell² rho²)^-1 F_rr + rho_y/(ell rho²) F_r.

The last term is required because rho varies in y. Evaluate the coefficients
at the relevant point; do not identify these coordinate derivatives with
orthonormal-frame components or physical energy. These are chain-rule
identities for scalar modal coefficient F, not a tensor-covariant derivative
formula.

## Why a C¹ solution bound is insufficient

An explicit source family preserves an open zero collar:

    w(q)=q²(1-q)² for 0≤q≤1, and w(q)=0 otherwise;
    g_N(s)=N^-1 w(N(s-1)),                 N=5,6,... .

Then g_N(s)=0 for s≤1, g_N∈C¹, and g_N'(0)=0. Exact rational
calculation gives

    sup |g_N| = 1/(16N),
    sup |L g_N| = 1/(30N²),
    sup |(Lg_N)_x| = sup |(Lg_N)_r| = 1/(16N),
    (Lg_N)_xx at (x,r)=(1+1/(4N)+1/16,1/16) = -3/16.

Indeed, q(1-q)≤1/4 on [0,1], with equality at q=1/2;
∫₀¹w(q)dq=1/30; and w'(1/4)=3/16.

These bounds hold on the retained domain D+={0≤x≤3/2,
max(0,x-5/4)≤r≤7/4-x}; for N≥5 each witness lies in D+, and
the entire source bump lies in the retained shell-history interval
1≤x≤5/4. The suprema refer to this domain, with physical one-sided traces.

Thus these C¹ lift norms tend to zero while a second derivative stays
nonzero. The family is an illustrative source perturbation, not a claim
that it satisfies all HDBLAST field equations. To see that qualitative smoothness and a flat collar alone do not fix the
gap, take b∈C∞ supported in (1,5/4), with b(9/8)=1, and set
g_epsilon(s)=epsilon b(s) sin((s-9/8)/epsilon). Every member is smooth
and vanishes on a fixed initial collar. As epsilon→0 the lift converges
to zero in C¹, but (Lg_epsilon)_xx=-1 at the fixed interior point
(x,r)=(5/4,1/8) of D+. This is a qualitative counterexample, not an
evaluation of any actual-field error.

An alternative without a zero initial collar is
`g_epsilon(s)=epsilon*(1-cos(s/epsilon))`: g(0)=g'(0)=0,
`sup|g|≤2 epsilon`, while `sup|g'|=1` on any sufficiently long interval.
Both examples describe why new derivative input estimates are necessary;
they do not indicate the actual error magnitude.

## Reproduction and provenance

Run `python checks.py` and `python -O checks.py`. The standard-library
program uses exact fractions and bivariate polynomial differentiation to
check the lift identities, all three second-derivative signs, continuity
conditions, the incompatible g=s example, the rational zero-collar family
and the physical-coordinate chain rule. It intentionally checks an
incorrect mixed-derivative sign and an omitted rho_y term as negative
controls. These are finite algebra checks supporting the elementary proof,
not a numerical evaluation of the retained actual first fields.

The registered equations and zero-collar context are in the public
[full 23 shell and subsidiary audit](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/shell/FULL23_SHELL_AND_SUBSIDIARY.md).
The source, smoothness, complete row and subsidiary premises stated there
remain in force. The model-dependent differentiation, quantitative high-jet
certification, and residual enclosure remain open.
