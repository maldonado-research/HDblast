# PROPOSED_UNVERIFIED_FOLLOWON: one boundary-inclusive forced-response calculation

23 October 2026. This note specifies a bounded next calculation; it contains no
new evaluated response, sign certificate, waveform, or physical-constraint
solution. It is separate from the frozen exterior result. All references below
are read-only retained files. No producer was rerun. The small cone identity and
tail arithmetic below are elementary implementation checks, without a novelty claim.

Reference roots: `NEW=/workspace/hdblast-research-work/forced-response-20261023/package-work/attempt-2/source`;
`OLD=/workspace/hdblast-research-work/geometric-bulk-source-20261016/capsule-source`.
Write `C=NEW/inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md`,
`B=NEW/inputs/FORCED_INCREMENT_BOUND.md`,
`O=OLD/projection/inputs/linear_operator.py`.
Line citations refer to these exact retained bytes.

## 1. Target and causal input region

Proposed first target: the probability-sphere mean of the covariant frame
normal-normal component, `lambda(W)=mean(W_yy)`, at

    (a,b)=(3/2,1/4),   d=a-b=5/4.

This is the retained reflected-packet location (`OLD/numeric/actual_first_jets.py:259–260`),
and its boundary history reaches beyond the flat first-contact time x=1.
It does not presume a nonzero shell signal immediately after that time.
Compute only on the backward-closed physical cone

    D+ = {0<=s<=a, u>=0, |u-b|<=a-s},
    I_boundary = [0,d],
    D- = {0<=s<=d, 0<=u<=d-s} subset D+.

The whole region satisfies s+u<=a+b=7/4<8, with s<=3/2<3.
Every bulk dependency of a point in D+ remains in D+ by the triangle inequality;
every boundary lift uses time <=s-u<=d. Thus recursive applications of the
full operator need no artificial outer boundary or uncontrolled extension.
The actual source oracle is valid here (`OLD/numeric/actual_first_jets.py:194–238`).
The existing exterior-error improvement is allowed only on entire source boxes
with u_lo>=s_hi (line 214); shell-reaching boxes require the retained global
actual-master and actual-gauge C3 errors. Do not replace them with the comparator.
The exact zero wedge s+u<1 may be imposed using C:219–223, not inferred from a
small interval output. Boundary data vanish for s<1 by that same causal premise.

The output should initially be one certified interval and its separate bulk,
reflected, material-source, and feedback accounting, with no sign requirement.
A subsequent continuous rectangle would require its own enlarged cone and
Lipschitz allowance; the exterior rectangle proof cannot simply be reused.

## 2. Exact K and L, including a useful folding identity

Let P_D select the mixed normal-tangential components k_y0,k_ya, and
P_N=I-P_D. In the selected 23-amplitude order these are five Dirichlet and
eighteen Robin amplitudes; they are not the full H dimensions 69 and 218
(C:40–46; O:10–11,109–112). Put J=P_N-P_D. For any physical source f,

    Kf(a,b) = (1/2) integral_D+ f(s,u) du ds
              + (1/2) integral_D- J f(s,u) du ds
            = (1/2) integral_(D+\D-) f + integral_D- P_N f.       (A)
    Lg(a,b) = -integral_0^d g(s) ds.                            (B)

Proof of (A): split C:211–212 at v=0 and substitute u=-v in its negative
part. The latter exists exactly for s<=a-b and 0<=u<=a-b-s. Since D- subset
D+, combine the overlap using (I+J)/2=P_N. In particular, at b=0, D+=D-
and Kf(a,0)=integral_D+ P_N f. This supplies a direct factor-of-two and
Dirichlet-cancellation check, not a modified boundary condition.

For 0<=b<=a,

    area(D-)=(a-b)^2/2,
    area(D+)=a^2-(a-b)^2/2,
    area(D+\D-)=2ab-b^2.

These follow by integrating the linear interval lengths. At the target the
areas are 25/32,47/32,11/16 respectively. The sum of the two unfused cone
areas is a^2=9/4. For a constant Robin source, (A) has weight a^2/2;
for a constant Dirichlet source it has weight ab-b^2/2. These exact constant
source tests detect omitted reflections and doubling errors.

Because yy is Robin, lambda J=lambda and lambda P_N=lambda. Its bulk value
therefore has weights 1/2 on D+\D- and 1 on D-. Other components retain their
own projection; one must not extend all 23 amplitudes evenly.

Endpoint derivatives needed in the X error ledger use no derivatives of f:

    (Kf)_x = (1/2) integral_0^x [f_ext(s,r+x-s)+f_ext(s,r-x+s)] ds,
    (Kf)_r = (1/2) integral_0^x [f_ext(s,r+x-s)-f_ext(s,r-x+s)] ds.

For x>r, (Lg)_x=-g(x-r), (Lg)_r=+g(x-r). Hence the physical-side radial
boundary derivative is +g, fixing the sign in (B). Flat corner data remove the
otherwise possible transition kink (C:208–223,243–250).

## 3. Material sources, normal projection, and frame conventions

Use the one-copy OUTWARD unit-normal affine Hessians QK (called QR in C), Qp
(QP in C), and the full varying-metric-lowered gauge Hessian QC. At r=0,

    GN = -2 tr_eta QK - 2 QC_y,
    G = -ell rho (-2 QK, -Qp, GN)   in order (m,w,K),
    G_m=2 ell rho QK,   G_w=ell rho Qp,
    G_K=2 ell rho (tr_eta QK+QC_y).                            (C)

Here ell=10^-6, K=k_yy, m=k_tangent, and outward partial_y
=-(ell rho)^(-1)partial_r. The normal source in (C) is required even when
only the yy mean is requested. It is the normal harmonic row, not an
additional physical constraint (C:161–174,190–197).

The actual shell adapter returns four distinct conventions:
`one_copy_outward_hessian`, its doubled two-copy outward value, its negated
doubled two-copy inward value, and `coordinate_physical_forcing`
(`OLD/numeric/shell_source_values.py:38–46`). Pass the FIRST into (C).
Passing the doubled source into the same formulas would add a spurious factor
two. Passing coordinate_physical_forcing and again multiplying its metric
entries by -2 would repeat the physical-row conversion.

In the registered angular basis QK_ab has seven slots: three scalar multiples
of kappa, two trace-free Hessians, and two TT tensors. Therefore for l=0,2,4,

    (tr_eta QK)_l = -QK00_l + 3 QKab_l,
    GN_l = 2 QK00_l - 6 QKab_l - 2 QCy_l,
    lambda(LG) = -2 ell rho integral_0^d
                    [-QK00_0+3 QKab_0+QCy_0] ds.              (D)

The minus in front of QK00 is Lorentzian; the factor three is the sphere
dimension, not the full shell dimension four. The last four angular tensor
slots have zero trace (`OLD/projection/angular_generic.py:71–84,121–126`).
The scalar basis has probability means (1,0,0), so no angular evaluation
factor 17 enters this mean (`NEW/compute/certify_exterior_response.py:64,94–96`).
QC_y is the `QC['y']` projected source from `source_coefficients.evaluate_sources`,
not `QC_background` (`OLD/GAUGE_COVECTOR_BOOKKEEPING.md:24–39,53–70`).

Material Robin feedback is equally explicit. With sigma'=-2p, sigma''=4phi0,

    B_y(m,w,K) = ((aK+sigma'w/3) eta_tangent,
                  pK/2-sigma''w/2,
                  -4aK+2a tr_eta(m)+(4sigma'/3)w),
    M=-ell rho B_y,
    TV = K(A V_x+B V_r+D0 V) + L(M P_N V|r=0).                 (E)

These are parallel-frame rows. In particular tangential trace-free components
have zero homogeneous Robin coefficient; a coordinate 2a diagonal must not
be inserted (C:176–188). Writing N=m_00 and J_sphere for the kappa coefficient,

    (M V)_K = ell rho [4aK+2aN-6aJ_sphere-(4sigma'/3)w],

which agrees with O:105–108,127. `J_sphere` here is unrelated to the reflection
operator J. O:177–184 implements exactly (C). The bulk forcing remains
F=(-2ell^2 rho^2 QEred,+ell^2 rho^2 Qphi), without another factorial two
(C:15–20,112–113). The retarded half in (A) is a separate wave-kernel factor.

## 4. Geometric cells and the smallest implementation task

Use rational polygon clipping, with separate lists for D+\D- and D-, plus
one-dimensional boundary intervals. In characteristic coordinates
xi=s+u, eta=s-u, the domain is

    s>=0, u>=0, xi<=a+b, eta<=a-b;
    D- adds xi<=a-b.

Split at u=0, xi=a-b=5/4, xi=1 (flat support onset), and eta=0 (exterior versus
shell-reaching error regime), then intersect a modest rational mesh with these
lines. The physical-area Jacobian is ds du=(1/2)d xi d eta; either clip in
(s,u) directly or include this factor exactly. This Jacobian is distinct from
the wave-kernel factor 1/2. The old special 20-cell rule is invalid on this cone.

For tighter input ranges, optionally split the inherited bump branches
xi and -eta at +/-1, +/-2 and their retained 1/64 tail transitions
(`OLD/numeric/actual_first_jets.py:96–155`); these are accuracy improvements,
not necessary validity assumptions. All enclosing (s,u) source rectangles
must stay in the oracle domain. This target has ample margin to x+r=8.
Boundary partition includes 0,1,5/4 and evaluates actual shell sources at r=0;
never normally differentiate shell-eliminated physical identities
(`OLD/numeric/shell_source_values.py:4–7`). Evaluate the QC source on that same
actual boundary jet. Preserve source correlations where algebraically available.

Exact area times whole-cell source range encloses a value integral without
a separate floating quadrature error, just as in the frozen exterior evaluator.
For X control, the characteristic endpoint integrals above also require whole
segment enclosures. The three subsequent operator applications need values,
both first derivatives, and the boundary trace of the preceding 23-amplitude
field. A table of point intervals or cellwise constant centers supplies no
certified C1 reconstruction. Use a compatible C1 reconstruction with rigorous
whole-cell and segment defects, including derivative and boundary traces.

## 5. Error ledger for four terms, or one full residual

Let S=KF+LG and W4=S+TS+T^2S+T^3S. The retained global constants are

    q=2109/1000000,
    S_cap=3*17387513156573+3851400=52162543321119.

They remain valid on this backward-closed domain. Exact rational arithmetic gives

    ||W-W4||_X <= q^4 S_cap/(1-q)
      =1031965259713929363466110159/997891000000000000000000
      =1034.146274206230... <1035.                              (F)

This tail includes all bulk couplings and Robin feedback; it excludes input,
quadrature, arithmetic and reconstruction errors (B:11–24,35–47,74–89).
The value-only mean error is THREE times an X error, so (F) alone gives <3105
for mean W_yy, not <1035. First x/r mean derivative errors have factor one.

One usable non-double-counting ledger: construct compatible C1 fields V0,...,V3
with certified defects

    ||V0-S||_X <= e0,
    ||Vj-T V(j-1)||_X <= ej, j=1,2,3.

Then V=V0+V1+V2+V3 satisfies

    ||W-V||_X <= tail4
       +(1+q+q^2+q^3)e0+(1+q+q^2)e1+(1+q)e2+e3.               (G)

Proof: compare Vj with T^jS and use d_j<=e_j+q d_(j-1), then sum the four
inequalities. The constants in (G) are exact bounds, not measured errors.

For example allocate e0<=3 deltaF+deltaG+Q0+R0. Here deltaF/deltaG are
continuous H-norm source discrepancies, Q0 is the certified X error of
integration (values AND first derivatives), and R0 is the C1 reconstruction
defect. For j>0 allocate ej<=3 deltaHj+deltaGj+Qj+Rj, where

    Hj=A(Vj-1)_x+B(Vj-1)_r+D0 Vj-1,
    Gj=M P_N Vj-1|boundary.

Coefficient approximation errors are included through, for example,
deltaHj<=deltaA sup|V_x|+deltaB sup|V_r|+deltaD0 sup|V|, and
deltaGj<=deltaM sup|P_N V|boundary, with arithmetic and evaluation defects
also enclosed. All norm conversions use the probability-sphere Gram matrices
and symmetric off-diagonal multiplicity two (O:130–142;
`NEW/compute/certify_exterior_response.py:42–52`). Inherited actual-first-field
uncertainty belongs in deltaF/deltaG or directly in source enclosures, once.
If interval area integration already includes within-cell variation and input
uncertainty, do not add that same width again as quadrature error.

An alternative is a single compatible reconstructed field V and rigorous
full residual R=V-S-TV. Then ||W-V||_X<=||R||_X/(1-q) (C:262–265).
This already accounts for all omitted iterations: do NOT add (F) again.
Neither route can certify a feature unless the resulting mean interval excludes
zero. No such interval has been computed in this follow-on.

## 6. Genuine remaining work and scope boundary

The retained material-source algebra, normal projection, actual jet oracle,
and full operator recipe exist. What is missing from the frozen exterior
evaluator is a boundary-inclusive clipped integration implementation, joined
actual QK/Qp/QC boundary enclosures for this history, and certified C1 fields
or residuals for the feedback terms. The evaluator currently requests only
bulk metric/scalar maps and sets shell terms to zero by exterior causality
(`NEW/compute/certify_exterior_response.py:81–106,119–123`); changing its apex
alone cannot perform this calculation. No claim is made that inherited global
post-contact uncertainty is small enough to resolve a feature. Finer spatial
cells do not automatically shrink the fixed inherited solution-error caps.

Stop the next implementation at one honest interval plus complete error
accounting; if it contains zero, report that without enlarging into a global
campaign. The retained construction already gives a conditional analytic
series for total U and a conditional physical completion in C:310–418. What
remains unevaluated is a certified full numerical waveform; this note does not
claim that its original boundary formula or analytic completion is missing.
A numerical physical total-response computation must retain the registered
second initial prescription and subsidiary constraints. W alone,
although driven by the actual sources and full reduced operator, is an
auxiliary summand and is not a standalone sourced physical-constraint solution
(B:62–67; `NEW/compute/EXTERIOR_RESPONSE_PROOF.md:9–14,138–147`).

## Read-only input pins

SHA256 of the principal retained inputs:

* C: `5de45b21236cfe8c6d99fb403bacc4624069c8d96d6cf0ad2ec89b8e2d3617e1`
* B: `317262b4581037f00e6f65cba6771ed49a0bc8c89bf5d6b5bf7759b758eb6d7f`
* NEW/compute/EXTERIOR_RESPONSE_PROOF.md: `a50bd56d2322bc0dbb7e8161e39bd8980b941575066d3d98e4cde8e76f737111`
* NEW/compute/certify_exterior_response.py: `58303378240715b415d0ad602b6953f666b8bf066a423276687edaf6a9c3f50a`
* OLD/numeric/shell_source_values.py: `ed9c6e084eae776509481fee4f490cea349eb5c1665139e286bd0a58be483d7f`
* OLD/numeric/actual_first_jets.py: `a68fa8025d2170efca6c3e85d443910a796d0c30a9b7dd17ed988cd9e00bbb00`
* O: `a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189`
* OLD/projection/angular_generic.py: `fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9`
* OLD/GAUGE_COVECTOR_BOOKKEEPING.md: `bc563ac6fbae1f14b133b9d8ccf5c61d57de2128090dfa8af0e7e19465d8a526`
