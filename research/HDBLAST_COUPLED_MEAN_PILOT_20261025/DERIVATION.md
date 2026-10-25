# Exact reflected ridge integration and its bounded response use

## Retained space

Use mean order (K0,L0,N0,J0,w0), Gram diagonal (1,2,1,3,1), reflection
J=diag(1,-1,1,1,1), PN=(I+J)/2 and PD=(I-J)/2.
Only L0 is Dirichlet. The norm is

X(z)=max(sup |z|H/3, sup |z_x|H, sup |z_r|H).

The physical cone is
D+={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x}.
The exact zero wedge x+r<=1 is essential to the improved contraction.
It gives |z|H<=(3/4)X(z), and |z(x,0)|H<=(1/4)X(z).

These are retained premises, not new PDE-existence conclusions.
See the pinned parent
[LOCAL_TOTAL_REVIEW.md](../HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/LOCAL_TOTAL_REVIEW.md)
and
[FULL_SECOND_RESPONSE_CONSTRUCTION.md](../HDBLAST_FORCED_RESPONSE_BOUNDS_20261023/inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md).

## Half-line kernel

Let continuous F,G have finite piecewise-polynomial descriptions and set

f(s,r)=F(r+s)+G(r-s), r>=0;
f_ext(s,-r)=J f(s,r).

The matrices used to form f need not commute with J. Reflecting the complete
source is necessary.

Let PF(y)=integral_0^y F, QF(y)=integral_0^y PF, and similarly PG,QG.
All four primitives vanish at zero. For a=r+x and q=r-x, q>=0,

Kf=(1/2){x PF(a)-(QF(a)-QF(q))/2
              +(QG(a)-QG(q))/2-x PG(q)}.

For q<0, put tau=x-r. Then

Kf=(1/2){
 x PF(a)-(QF(a)+QF(tau))/2
 +(QG(a)+QG(-tau))/2-r PG(-tau)
 +J[tau PF(tau)-QF(tau)+(QG(tau)+QG(-tau))/2]
}.

Both formulas equal the retained integral

Kf=(1/2) integral_0^x integral_(r-x+s)^(r+x-s) f_ext(s,v) dv ds.

Proof: for x<=r the diamond stays on the physical half-line, so integrate the
two ridge terms directly. For x>r, split the inner integral at v=0. The negative
part occurs exactly for 0<=s<=x-r. Reflect it with v -> -v and multiply its
complete vector by J. Integrate each primitive with the substitutions
v+s and v-s. Normalization at zero removes all constants. No quadrature
approximation is used.

For p=(Kf)_x+(Kf)_r, in both branches,

p=x F(a)+(PG(a)-PG(q))/2.

For q>=0,

m=(Kf)_x-(Kf)_r=(PF(a)-PF(q))/2+x G(q).

For q<0,

m=(PF(a)-PF(tau))/2+r G(-tau)
    +J[tau F(tau)+(PG(tau)-PG(-tau))/2].

Thus (Kf)_x=(p+m)/2 and (Kf)_r=(p-m)/2. These agree at x=r, where
m=PF(2x)/2+x G(0). Continuity of F,G suffices for Kf to be C1, including
polynomial joins. No derivatives of F or G are needed. A possible odd source
jump at the reflected shell does not invalidate this kernel calculation.

At r=0, m=Jp, while the value formula lies in PN. Consequently
PD Kf=PD(Kf)_x=0 and PN(Kf)_r=0. Initial values and first derivatives are zero.
For finite polynomial pieces of degree<=d the integrated pieces have
degree<=d+2 and exact rational evaluations at rational arguments.

## Shell lift

For continuous g=PN g with g(0)=0,

Lg(x,r)=-integral_0^((x-r)+) g(s) ds.

For x>r, (Lg)_x=-g(x-r) and (Lg)_r=+g(x-r); both vanish for x<r.
The corner condition gives C1 across x=r. Its shell radial derivative is g.
The sign is fixed by the retained boundary law PN U_r=M PN U+G.

Sufficient zero-wedge conditions are F=0 on[0,1], G=0 on[-1,1], and g=0
on[0,1]. They imply Kf+Lg=0 when x+r<=1.

The shell matrix is singular: M(0,0,3,1,0)^T=0. The construction uses no
inverse, coercivity or symmetric-energy cancellation.

## Full coupled finite forcing

Keep the exact registered rational displacement u0 in B0. Replace only its
registered velocity by a single continuous rational P1 profile v, with
certified error epsilon_v. Its primitive V is integrated exactly. B0 is the
reflected free wave of u0,v.

For feedback only, let u,d be independent continuous P1 approximations to
u0 and its radial derivative. They need not satisfy d=u'. They are not
replacement C1 displacement data.

Extend u and v by J, d by -J. Then V(-y)=-J V(y).
With full constant five-by-five matrices A*,B*,D*,M*, define

F=1/2[(A*+B*)d+(A*+B*)v+D*u+D*V],
G=1/2[(B*-A*)d+(A*-B*)v+D*u-D*V],
g=M* PN[u+V].

This is exactly the approximate physical forcing formed from the free-wave
value and derivative profiles. No internal scalar/metric coupling is dropped.
F,G are piecewise quadratic; Kf is piecewise quartic; Lg is piecewise cubic.
Their primitives reduce to the first three exact primitives of P1 profiles.

The stored profile nodes are shared, include the cutoff joins65/64,66/64,
and are exactly zero through65/64. The initial footprint is[0,7/4].
Each P1 error is bounded by endpoint quantization plus h^2/8 times a whole-cell
second-derivative cap. For d this requires the third derivative of u0;
the native producer therefore uses order-six packet jets on each smooth
cutoff segment, with the correct factorial and ell factors.

The actual coefficient enclosure uses the retained shell Taylor coefficients
and the validated fifth-order remainder, not nominal candidate coefficients.
The shifted scalar is eta=phi0-1. In the retained operator,
V,V',V'' are UP[0],UP[1],UP[2] without an extra half,
sigma'=-2p, sigma''=4(1+eta).
The initial-profile quantity UP1/p is distinct from the time connection
h=tanh(ell*x)/rho.

## Complete conditional error identity

Let S=KF_exact+LG_exact denote the retained inhomogeneous forcing contribution.
Its norm is bounded by S_cap=52162547172519/2. It is charged once, even when
not numerically integrated in this pilot.

Let q=3409/8000000, U_cap=1.03e33, e_free the pinned accepted initial error,
and

e0=e_free+epsilon_v+S_cap;
X0=(1+q)U_cap+e0 >= X(B0).

The last inequality follows from B0-U=(B0-U_free)-TU-S.

For Gram-induced operator error bounds epsilon_A,...,epsilon_M, and measured
P1 errors epsilon_u,epsilon_d,

delta_f <= (epsilon_A+epsilon_B+(3/4)epsilon_D)X0
              +(norm(A*)+norm(B*))epsilon_d+norm(D*)epsilon_u;
delta_g <= (1/4)epsilon_M X0+norm(M*)epsilon_u.

Use a Gram-weighted Frobenius upper bound to enclose each matrix operator norm.
No symmetry or cancellation is assumed. Then

e1=(3/2)delta_f+delta_g >= X(TB0-[Kf+Lg]).

For U=B+TU and Uhat=B0+B1, B1=Kf+Lg, the exact identity is

U-Uhat=T^2 U+(I+T)(B-B0)+(TB0-B1).

Therefore

X(U-Uhat) <= q^2 U_cap+(1+q)e0+e1.

The native pilot measures this complete expression. It uses the fixed mesh
sequence64,128,256 and stops at its finite limit. A proposed suballocation is
not a measured bound. An insufficient upper bound proves no lower bound on
actual error.

A finite Uhat generally has a nonzero remaining Robin defect; the error
identity covers it. Do not replace the full boundary law by homogeneous
Neumann or report the finite candidate as an exact solution.

## Verification and limits

The independent JavaScript checks use exact rational polygon clipping and
Green monomial integration of the diamond, separate from the primitive
formula, plus exact endpoint-integral derivatives. They check both reflection
signs, polynomial and continuous hinge profiles, shell/initial/characteristic
joins, and wrong-sign/wrong-half controls. Python is compared against
independently generated exact fixtures and direct piecewise integration.

The numerical source closures are authenticated before imports and after
evaluation. Normal and optimized runs must agree. Source pins and numerical
outcomes are separate from software test success.

A passing result is a conditional C1 five-mean approximation. It is not a
full angular second-order spacetime, a pointwise Einstein residual certificate,
an invariant shell-volume value or sign, a nonlinear finite-amplitude family,
or physical confirmation. This standard kernel identity establishes no
originality or priority claim, and internal AI-assisted review is not external
specialist peer review.

## Measured bounded outcome

The first actual runner execution used the declared meshes64,128,256.
The total upper caps were approximately3.470126e30,1.188711e29,
and1.609032296423508e28 respectively. Only the last meets1e29.
All displayed scientific upper bounds use outward rounding; the exact fraction
and all per-cell quantities are in results/REPORT.json.

At the accepted level the main ledger is:

| Contribution | Safe outward upper cap |
| --- | ---: |
| Retained initial error plus new velocity reconstruction and omitted forcing, e0 | 8.518098e27 |
| Complete first-feedback approximation error, e1 | 7.381566e27 |
| Omitted two-term tail q^2 U_cap | 1.870300e26 |
| Total, including (1+q)e0 | 1.6090323e28 |

These are X-norm bounds in the retained scaled coordinates, not physical
measurements. The source contributions and initial uncertainty are included
once. Sample intervals enclose arithmetic evaluation of the specified
candidate; the separate global X bound encloses its difference from the
conditional exact response.

The kernel, full profile coefficients and exact matrix coefficients define
the approximation continuously over D+. Six stored evaluations are examples,
not the definition or the extent of the certificate. The full angular
second-order spacetime and intrinsic-observable sign remain unevaluated.
