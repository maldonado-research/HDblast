# Full23 shell source, nullspace and physical completion audit

26 October 2026. New exact algebra; frozen inputs remain read only. This note
distinguishes the inherited formal second coefficient from a newly evaluated
numerical reconstruction. It asserts neither an exact nonlinear family nor a
physical origin interpretation.

## Authenticated scope

`check_full23_shell.py` authenticates complete manifests of the seven closures
named in the latest `balanced/exact_boundary.py`, and the latest shell precision
closure itself, before importing any scientific code. It verifies the helper's
SHA256 and compares its independent complete-closure authentication with the
new checker. It imports only the already authenticated G linear operator and
source adapter and the N exact polynomial ring. It executes no old producer,
mesh campaign or source certificate. Every closure is reauthenticated at exit.

The read-only checkout HEAD `4d532de8e698494f8fc7e2e31d7c7e69ede3e7a4` is the
verification-only descendant of the registered main
`9a45ee0deaf3496fee93491c2803cb168bb36c92`. Scientific manifest pins, rather
than a working-tree branch name, select this audit's inputs.

## Complete primal boundary rows

Use the parallel background frame, `ell=10^-6`, `x=t/ell`, inward `r=z/ell`
and `partial_y=-(ell*rho)^-1 partial_r`. Source conventions are second parameter
derivatives, with the second coefficient multiplied by `epsilon^2/2`.

The 23 ordered amplitudes comprise `K,L,N,J,w` at scalar degrees 0,2,4;
`v,u` at degrees 2,4; `H` at degrees 2,4; and both `T+`, `T-`.
There are five homogeneous Dirichlet rows `L0,L2,L4,v2,v4` and eighteen
Robin rows

`PN U_r = M PN U + G`, `M=-ell*rho*By`.

For each scalar degree the nonzero homogeneous rows are

```
By K = -4*a*K-2*a*N+6*a*J+(4*sigma1/3)*w
By N = -a*K-(sigma1/3)*w
By J = +a*K+(sigma1/3)*w
By w = (p/2)*K-(sigma2/2)*w.
```

`By u=By H=By T+=By T-=0`. These are parallel-frame derivatives.
Adding a coordinate-component `2*a*J` term changes the problem.

The physical one-copy affine Hessians are `QR`, `QP`, and the normal gauge
source is `GN=-2*tr(QR)-2*QC_y`. Their solver signs are

```
G_K=-ell*rho*GN
G_N=2*ell*rho*QR00
G_u=2*ell*rho*QR0a
G_ab=2*ell*rho*QRab
G_w=ell*rho*QP.
```

The scalar and metric rows share one common tension, with
`sigma1=-2*p`, `sigma2=4*phi_b`, `sigma3=4`. Four mixed metric Dirichlet
rows and twelve Robin rows refer to physical fibre components; their modal
counts are five and eighteen. Tangential subsidiary constraints are not extra
independently prescribed primal rows.

## New exact homogeneous boundary simplification

At each scalar degree put `s=(-N+3*J)/2`, `q=N+J`. Then

```
N=-s/2+3*q/4, J=s/2+q/4,

By [s,K,w]^T = B3 [s,K,w]^T,
B3 = [[0,2*a,2*sigma1/3],
      [4*a,-4*a,4*sigma1/3],
      [0,p/2,-sigma2/2]],
By q=0.
```

Thus the eighteen-dimensional Robin matrix is similar to three identical
`-ell*rho*B3` blocks and nine exact zeros. Its nine universal right-null
vectors are `K=w=0,N=3*J` at each scalar degree and the six unit slots
`u2,u4,H2,H4,T+,T-`. Its corresponding left-null modes are `N+J` and
the same six unit slots. The exact checker proves both null identities
symbolically and exhibits a rational background parameter witness of rank
nine. Exceptional parameters may reduce rank further.

The transformation is orthogonal in the positive comparison Hilbert norm:
`N^2+3*J^2=s^2+3*q^2/4`. Within every active block `s,K,w` have the
same scalar harmonic weight. Consequently `||M||_H` equals the ordinary
Euclidean operator norm of `-ell*rho*B3`, with the useful bound

```
||M||_H^2 <= (ell*rho)^2 *
 (36*a^2+20*sigma1^2/9+p^2/4+sigma2^2/4)
 = (ell*rho)^2*(36*a^2+329*p^2/36+sigma2^2/4)
```

in the registered `sigma1=-2*p` convention. This is an upper norm bound,
not coercivity. No inverse of M is used. Null rows can still have nonzero
inhomogeneous forcing and therefore are retained.

## New eight-profile source factorization

At the shell write the actual first fields as
`H=h*Y`, `eta=v*Y`, `m00=-2*X*Y`, `m0a=V*D_aY`,
`mab=2*(L*Y*kappa_ab+d*Hess_ab(Y))`. Define

```
alpha=a*h^2/4-sigma1*h*v/6-sigma2*v^2/6
beta=a*h+sigma1*v/3
pi=p*h^2/4+sigma2*h*v/2+sigma3*v^2/2
tau=2*X+6*L-16*d.
```

The primitive `GN` coefficients in the angular basis
`(Y^2,|DY|^2,|Hess Y|^2)` are

```
g0=4*sigma2*v^2/3+4*sigma1*h*v-10*a*h^2+8*a*h*tau
   -4*a*(4*X^2+12*L^2-64*L*d)
g1=8*a*V^2
g2=-16*a*d^2.
```

The complete source is exactly `G=B*c` for the fixed rational 23-by-8 matrix
stored in the report and profiles

```
c = [-ell*rho*g0,-ell*rho*g1,-ell*rho*g2,
     2*ell*rho*(-alpha+2*beta*X),ell*rho*pi,
     -2*ell*rho*beta*V,
     2*ell*rho*(alpha-2*beta*L),-4*ell*rho*beta*d].
```

Its eight columns are respectively the K-slot projections of
`Y^2,|DY|^2,|Hess Y|^2`, the N-slot `Y^2`, w-slot `Y^2`, u-slot
`Y*DY`, angular `Y^2*kappa`, and angular `Y*Hess Y`. The source's universal
angular rank is exactly eight, established by the rational coefficient rank
over independent symbolic background/field amplitudes. The report supplies
`B^T H B`, giving the exact error-norm identity
`||B*delta c||_H^2=delta c^T*(B^T H B)*delta c`.

In particular

```
G_T+=G_T-=-(8/3)*ell*rho*beta*d
G_u4=2*G_u2
G_H4=-(2/3)*G_H2.
```

The TT source may be nonzero despite homogeneous Neumann TT rows and zero
initial TT data. This rank-eight affine source factorization is not an
eight-dimensional invariant response subspace: bulk couplings and the Robin
operator act on the full retained section multiplicities. Orientation parity
can identify the two TT responses only after its separate all-operator and
data check.

## C1 response needs source values, not shell-source derivatives

For any continuous Robin source `g` with `g(0)=0`, use the zero extension
for negative time and

`(Lg)(x,r)=-integral_0^max(x-r,0) g(s) ds`.

The physical-side derivatives are exactly
`(Lg)_x=-g(x-r)`, `(Lg)_r=g(x-r)` on `x>=r`, and zero on `x<r`.
They agree continuously on `x=r`. Hence on `0<=x<=3`, with
`||U||_X=max(sup|U|/3,sup|U_x|,sup|U_r|)`,
`||Lg||_X<=sup|g|`. The same result holds for an error `delta G`.
No derivative of a nonmean source, and no radial derivative of a shell-only
source, enters this bound. The reflected bulk Duhamel map likewise has
`||K f||_X<=3*sup|f|` for continuous f, including the physical one-sided
Dirichlet reflection. These estimates apply to all 23 retained amplitudes in
the correct Hilbert norm.

The Robin feedback is `L(M PN U|b)`. It uses only the boundary trace of U,
and `||U|b||<=3||U||_X`. A fixed-point or residual certificate can therefore
use C0 bulk and shell approximation errors for its C1 response target.
Differentiating the boundary row itself to request a numerical G_x estimate
is unnecessary for this target. A C2 response target would be different.

## Initial/corner compatibility

The inherited physical second initial prescription and harmonic initial
one-jet are not replaced by zero scalar data. They include
`g2_tt=-2*rho^2*(u*Y)^2`, the mixed metric velocity conversion, nonzero
scalar correction and scalar normal-velocity correction. The ordinary
t-velocity is converted to the solver's x-velocity by multiplication by ell.
The first gauge wave has `xi=xi_t=0` initially and its gauge acceleration
matches the inherited first spacetime one-jet used in this completion.

The initial packet and every second correction have support in
`1<=r<=2` with flat endpoint zero extension. The shell has an open zero
initial collar. Therefore U0, V0 and the actual quadratic shell source have
zero corner jets; reflected free data, the affine boundary lift and the
Robin feedback respect the fixed zero wedge `x+r<1`. The differentiated
corner recursion is compatible at every qualitative smoothness order.
Flatness supplies these exact identities, not new quantitative high-jet caps.
A concrete numerical reconstruction still must preserve or bound its actual
initial, boundary and corner defects.

## Inherited formal subsidiary theorem and numerical scope

The inherited C `inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md`, sections 5
and 7, already defines the actual smooth full finite-angular second
coefficient by a convergent Volterra/Picard series. It supplies the C1
contraction and then a qualitative regularity induction. The construction is
conditional on the registered smooth background and actual first solution.
It is not a new numerical full23 waveform or an exact nonlinear existence
result. The same construction is retained in P's evolution source.

Selected-space closure alone does not prove any physical constraints. The
formal physical completion additionally uses all of the following premises:

1. The background solves its Einstein/scalar equations. The same actual
   first metric/scalar representative solves the first natural equations,
   its physical shell rows and `DC[q]=0` on a neighborhood. Boundary
   orthogonality is `q_y_mu=0`, including its tangential derivatives.
2. The complete natural quadratic Einstein Hessian QE and scalar Hessian
   Qphi obey differentiated Bianchi. The reduced source is precisely
   `Qred=QE-nabla_(A QC_B)`, with the registered one-half symmetrization
   and `QC=D^2C[q,q]`. An omitted Ricci, gauge-gradient or first scalar
   term would break this premise. The G exact geometric projector supplies
   the previously missing full Ricci/reduced-gauge coefficient adapter;
   a numerical integration must connect that adapter to the actual first
   field enclosures, rather than treating approximate first fields as exact
   harmonic solutions.
3. The same physical affine shell Hessians obey the natural Codazzi join
   `QE_ya=D^b QR_ba-D_a tr(QR)`. The scalar junction is essential to the
   cancellation, even though the background tangential scalar gradient is
   zero and no independent QP term remains in that displayed identity.
4. The normal primal source is exactly
   `GN=-2*tr(QR)-2*QC_y`, verified here symbolically at all three scalar
   degrees. It gives the normal component of `C2=DC[k]+QC` zero at the
   shell. This is component C2_y=0, not the equation partial_y C2=0.
5. The completed initial harmonic one-jet gives C2=0 initially. The inherited
   physical Hamiltonian and four momentum constraints, together with the
   reduced equations' accelerations, give partial_t C2=0 initially.
6. Qualitative smoothness justifies classical subsidiary traces. It follows
   from the retained smooth/flat master and the characteristic trace-ODE
   induction for the first gauge and second coefficient, without extracting
   a new numerical fourth gauge derivative from the C3 certificate.

With these premises, Bianchi gives `(Box0+Ric0) C2=0`. The normal
component has Dirichlet data. Codazzi and the reduced equations give
`partial_y(C2)_mu-2*a*(C2)_mu=0` for tangential coordinate covector
components, or `partial_y(C2)_mu-a*(C2)_mu=0` in the parallel frame.
The finite angular subsidiary system is normally hyperbolic. Its bounded
coefficient characteristic trace ODE gives local uniqueness; continuation
gives C2=0 throughout the causal domain. This step requires no favorable
Robin sign, energy coercivity or inversion of a singular boundary matrix.
Identical copies provide matching; the inherited force/conservation identities
supply the remaining shell identities.

A new continuous numerical full23 certificate can approximate that exact
formal coefficient and inherit its physical interpretation if it certifies
actual initial, full bulk, shell and operator approximation errors with the
same conventions and the complete inverse bound. The reconstructed candidate
need not satisfy the physical constraints exactly at each finite iterate.
Its numerical proximity to the exact physically completed coefficient must
be established. A low-dimensional selected-space declaration or a mean5
sample certificate does not provide that proximity by itself.

## Primary frozen sources

- G `projection/inputs/linear_operator.py`: all 23 homogeneous bulk/boundary
  amplitudes, source packing and initial velocity scaling.
- G `projection/inputs/source_adapters.py`: actual first-gauge shell Hessians.
- N `angular/EXACT_ANGULAR_OPERATOR.md`, every-shell-row section, and
  `angular/SELECTED_SPACE_CLOSURE.md`: angular conventions and closure scope.
- P `boundary-review/BOUNDARY_COMPATIBILITY_AUDIT.md`, sections 2-6:
  normalized physical/gauge Hessians, Codazzi, initial/corner compatibility.
- C `inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md`, sections 5-7:
  the C0-to-C1 integral inverse, qualitative smoothness and subsidiary uniqueness.
- Latest shell precision `balanced/exact_boundary.py`: seven complete upstream
  closure pins and their authentication implementation.
