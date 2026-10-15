# Independent angular-uniform Robin estimate for the full wave bundle

15 October 2026. Independent derivation. This is a linear estimate and a
smooth linear existence construction on the registered positive-cone collar.
It does not construct an exact nonlinear initial-data family or assert a
nonlinear contact theorem. No external publication has been made.

The conclusion is positive: the nonsymmetric full Robin matrix does not
prevent an angular-uniform linear estimate. A local bundle conjugation
converts its homogeneous boundary conditions into ordinary mixed
Dirichlet/Neumann conditions. The positive norm is a control norm on all
fifteen metric components and the scalar, rather than gravitational energy.
Bare L2 boundary-source control is not part of this conclusion.

## 1. The complete background principal operator

Write x=log(r/b), z=-x, and c(t)=sech(t). Then

    g0 = rho(-z)^2 [dz^2 - dt^2 + cosh(t)^2 kappa_S3].

For scalars, direct calculation gives

    rho^2 Box_0 = -partial_tt + partial_zz + c^2 Delta_kappa
                  - 3 tanh(t) partial_t - 3 a rho partial_z.

The reduced metric equations have principal part -Box_0 k/2 and the scalar
equation has principal part Box_0 w. Multiplication of the metric equations
by -2 rho^2 and the scalar equation by -rho^2 therefore makes the common
principal operator

    P0 = partial_tt - partial_zz - c(t)^2 Delta_D.

Here D is a fixed angular connection on the full bundle

    E = Sym^2(R dy + R dt + T*S3) + R,

using the product trivialization in the radial and time variables. The
difference between the actual background tensor connection and this fixed
angular connection contributes only first-order and zeroth-order bundle
operators. Thus, on the controlled collar, the complete sixteen-component
linear system is

    P0 U = A^t U_t + A^z U_z + A^i D_i U + A^0 U + F.       (L1)

All coefficients and their required derivatives are bounded on a finite
time slab and the controlled radial support. Constants may depend on those
bounds and the slab, but will not depend on an angular representation label.
No scalar-master decomposition or division by an angular eigenvalue is used.

Choose a positive product fiber metric H using dt^2, dy^2, and kappa, with
the induced positive metric on symmetric tensors and a positive scalar
weight. This does not use the indefinite Lorentzian contraction of k.
The fixed angular connection preserves H. Changing the smooth product
trivialization or the positive fiber metric only changes the constants.

## 2. The actual twelve by twelve matrix

Let s=sigma'_0, b2=sigma''_0, a=a_b and H0=h0|boundary. Order the Neumann
bundle as U_N=(r,K,w), where r=k_yy, K=k_ab is a symmetric four by four
tangential covariant tensor, and w is the scalar. Hold H0 fixed in the
normal extension of this boundary trivialization. The exact proper-y
homogeneous linear boundary law is

    partial_y U_N + B_y U_N = 0,

    B_y(r,K,w) = ( 4a r - 2a tr(H0^{-1}K) - (4s/3)w,
                   -a H0 r - 2a K - (s/3)H0 w,
                   (s/4)r + (b2/2)w ).                     (L2)

There are ten independent K components. An off-diagonal K_ij contributes
2(H0^{-1})^{ij}K_ij to the trace. The four remaining metric components,+k_ya, have Dirichlet boundary data. Couplings to specified Dirichlet traces
can be included in the known boundary source.

Since partial_y=-rho_b^{-1}partial_z, the inward form is

    U_D=d,                 -U_N,z + M(t) U_N = g,
    M(t)=rho_b B_y(t).                                      (L3)

Equation (L2) follows independently by differentiating the exact laws:

    F_yy = (4 sigma/3) N^3 - 2a N^4 tr_h(H0),
    F_h  = (sigma/3) N h,
    F_phi= -(sigma'/2) N,
    N=sqrt(g_yy).

Their derivative at the background is -B_y. The accompanying checker also
checks every factorial order-two source in equation (5) of the boundary
completion note, including k_yy=2n2+2n1^2 and the 4n1 n1,y term.

The nine-dimensional trace-free subspace has B_y=-2a I. On (r,tau,w),
where tau=tr(H0^{-1}K) with H0 held fixed in the normal extension, the matrix
is

              [ 4a   -2a   -4s/3 ]
    B_trace = [ -4a  -2a   -4s/3 ].                         (L4)
              [ s/4   0     b2/2 ]

For example, a=b2=0, s=1 gives eigenvalues 0 and +/- i/sqrt(3).
This example is an algebraic check of the generic matrix structure, not an
assertion that these are registered physical background parameters. It
shows why a positive symmetrizer of B_y must not be assumed. The derivative
of tau would instead contain the derivative of H0 if a y-dependent frame
were used; (L4) explicitly uses the fixed normal extension.

H0(t)=-rho_b^2 dt^2+rho_b^2 cosh(t)^2 kappa is angular-parallel. Therefore
B_y, M, the Dirichlet/Neumann projections, and all transformations below are
angular-parallel and SO(4)-equivariant on the complete bundle.

## 3. An exact collar conjugation

Take a smooth real function chi(z) supported in a small collar, with
chi(0)=0 and chi'(0)=1. Set

    T_D=I,             T_N(t,z)=exp[chi(z) M(t)],
    U=T V.                                                  (L5)

Then T is invertible everywhere, is the identity at the boundary and away
from the collar, preserves the two subbundles, and satisfies T_N,z=M T_N
at z=0. Consequently (L3) becomes exactly

    V_D=d,                        -V_N,z=g.                  (L6)

This calculation does not require M to be symmetric, normal, diagonalizable,
or time-independent. The inequalities

    ||T||, ||T^{-1}|| <= exp(|chi| ||M||)

and smoothness of matrix exponentiation give bounded coefficient derivatives
on the finite slab. When M and M_t do not commute, the correct derivative is

    T_t = integral_0^1 exp[(1-r)chi M] (chi M_t)
                            exp[r chi M] dr.                (L7)

It is generally not chi M_t T. The checker tests a noncommuting example.

Because the principal part is a scalar differential operator times the
identity on the bundle, T^{-1}P0(TV) has precisely the same principal part.
Angular-parallelness even gives D_i T=0 and [Delta_D,T]=0. Thus the conjugated
system has the form (L1), with different bounded first-order coefficients
and source T^{-1}F, and ordinary mixed boundary conditions (L6).

## 4. Positive energy, independent of all angular cutoffs

First take d=g=0. On z>=0 define

    E1[V](t) = (1/2) integral ( |V_t|_H^2 + |V_z|_H^2
                      + c(t)^2 |D_Omega V|_H^2 + |V|_H^2 ) dz dvol_kappa.

Since c(t) stays strictly positive on a fixed finite time slab, this is
uniformly equivalent to the H1 x L2 phase norm there. Integration by parts
in z produces the flux -<V_t,V_z>_H at z=0. Its Dirichlet contribution
vanishes because V_D,t=0; its Neumann contribution vanishes because
V_N,z=0. Angular integration by parts has no boundary. The derivative of
c^2, the term <V,V_t>, and each bounded first-order bulk term are bounded
by C E1. In particular the angular lower-order coupling is controlled by
||D_Omega V|| ||V_t||, not by a representation-dependent coefficient norm.
Hence

    dE1/dt <= C E1 + sqrt(2 E1) ||T^{-1}F||_L2,

and, after the usual regularization at E1=0,

    sqrt(E1(t)) <= exp(Ct) [ sqrt(E1(0))
                            + C integral_0^t ||F(s)||_L2 ds ].           (L8)

The constants depend on collar coefficients, T, the positive norm and the
time slab. They are independent of any SO(4) angular cutoff. Since T and
its first derivatives are bounded, (L8) is also an estimate for the original
full field U in a uniformly equivalent positive control norm.

This proves that the lack of a symmetric Robin-flux cancellation in the
original variables is not a failure of the linear estimate. A nonsymmetric
Robin boundary matrix may lead to finite low-frequency growth; (L8) allows
that. It rules out unbounded growth rates as angular frequency increases.

## 5. Smooth existence and higher derivatives

Let Pi_N be the compact-group isotypic projection onto finitely many SO(4)
representations, applied to the full bundle, rather than a componentwise
scalar-Laplacian cutoff. All background operators, boundary projections,
T and the chosen lift below commute with Pi_N. Each projected problem is
a finite collection of radial vector wave equations, with smooth bounded
coefficients for that fixed N. The characteristic/trace-ODE construction
in the 14 October finite-angular note gives its smooth compatible solution.
Applying (L8) to the difference of two such solutions bounds the difference
by the tails of the initial data and source. Therefore these solutions
converge in the phase norm to a unique full-angular solution. This is an
existence proof independent of an unproved general Robin bundle theorem.

For smooth regularity, let J_alpha be the infinitesimal SO(4) actions on
complete tensor sections. They are first-order angular operators, generate
angular Sobolev norms on S3, and commute with the full background system
and its boundary operators. Time derivatives also preserve the homogeneous
boundary conditions after (L5). Define

    Em[V] = sum_{j+|alpha|<=m} E1[partial_t^j J^alpha V].      (L9)

Time commutation with c(t)^2 Delta_D creates terms such as
(partial_t c^2) Delta_D partial_t^{j-1}J^alpha V. They are controlled by
Em: their derivative order is j+|alpha|+1, with one fewer time derivative
and one more available angular derivative. Other commutators are first
order with bounded coefficient derivatives. Thus, for compatible smooth
data,

    sqrt(Em(t)) <= C_m exp(C_m t) [ sqrt(Em(0))
          + integral_0^t sum_{j+|alpha|<=m}
                     ||partial_t^j J^alpha F(s)||_L2 ds ].   (L10)

The constants remain independent of N. The equation solves for V_zz.
Repeated differentiation then recovers higher normal derivatives by
induction on their number: two normal derivatives are exchanged for two
tangential/time derivatives and lower-order terms. Source derivatives up
to total order m-1 enter this recovery. This yields the usual full smooth
regularity and spatial H^{m+1} x H^m control, with the explicitly stated
time-jet/source norms in (L9)-(L10). This note does not claim an optimal
anisotropic boundary trace norm.

Initial time jets in Em(0) are determined recursively by the wave equation.
They must satisfy V_j,D|boundary=0 and partial_z V_j,N|boundary=0 for every
required order. In the 14 October construction, every perturbative datum
and source vanishes on an open initial shell neighborhood, so these corner
conditions follow immediately from the background solution and recursion.
For general data they are hypotheses, not automatic consequences of the
lowest-order boundary conditions.

Radial finite propagation follows either from the finite-angular
characteristic construction or from the energy argument outside moving
surfaces z=L+t. Angular propagation cannot increase the maximum radial
speed above one. Thus an arbitrary smooth extension of coefficients beyond
the controlled support cannot affect the constructed solution there.

## 6. Inhomogeneous boundary sources: a safe, explicit norm

Choose compactly supported radial functions eta0 and eta1 such that
eta0(0)=1, eta1(0)=0, eta1'(0)=1. In the conjugated variables take

    L_D(t,z,Omega)=eta0(z)d(t,Omega),
    L_N(t,z,Omega)=-eta1(z)g(t,Omega).

Then W=V-L has homogeneous mixed boundary conditions. Its source is

    F_W=T^{-1}F-P_T L,

where P_T is the full conjugated differential operator, including its
bounded first-order terms. The exact source estimate is (L10) with F_W,
plus sqrt(Em[L](t)) on the left-hand reconstruction and the lifted initial
jets W_j(0) on the initial side. For example at energy order one,

    ||P_T L(t)||_L2 <= C [ ||d_tt||_L2 + ||g_tt||_L2
                     + ||d||_H2(S3) + ||g||_H2(S3)
                     + ||d_t||_L2 + ||g_t||_L2 ].            (L11)

The analogous differentiated expression controls (L10). Pointwise phase
norms of L use only d,g and their first time/angular derivatives. This is
a deliberately sufficient, nonoptimal source norm and is independent of
the cutoff. All actual smooth Taylor boundary sources satisfy it. A claim
of a sharper source norm needs an additional trace/lifting theorem.

## 7. Correct Neumann normalization at glancing

Freeze the flat scalar principal operator at the boundary. At Laplace
frequency s_L=gamma+i tau, gamma>0, and angular covector eta, let

    kappa=sqrt(s_L^2+c^2|eta|^2), Re(kappa)>0,
    Lambda=sqrt(1+|s_L|^2+c^2|eta|^2).

The decaying mode is exp(-kappa z)a. The raw Neumann boundary symbol is
kappa, whereas the boundary symbol acting on unit phase amplitude is
kappa/Lambda. This normalized symbol tends to zero in the glancing regime.
That is not a proof of homogeneous ill-posedness: ordinary homogeneous
Neumann waves already have this symbol and satisfy the exact zero-flux
energy estimate.

At tau=c|eta| and fixed gamma, |kappa| and Re(kappa) are of order
sqrt(gamma Lambda). For unit prescribed Neumann data a=1/kappa, the
integrated squared H1 phase norm has the scale

    Lambda^2 / ( |kappa|^2 Re(kappa) )
                  approximately Lambda^(1/2) gamma^(-3/2).   (L12)

Thus a cutoff-independent estimate with only bare L2 boundary data is
false even for the scalar flat Neumann problem. The explicit smooth lift
in section 6 avoids that false claim. A fixed bounded M changes the raw
Robin symbol to kappa I+M. At high-frequency glancing and fixed gamma,
the ratio ||M||/|kappa| tends to zero. In any case the exact conjugation
and (L8) give the finite-time uniform estimate without an inference from
an unweighted determinant. No physical instability is deduced here.

## 8. Subsidiary system and scope

The background subsidiary equation (Box_0+Ric_0)C=0 and the first gauge-vector
equation have the same scalar principal structure. Their boundary has one
normal Dirichlet row and four tangential Robin rows. The same conjugation,
positive norm and cutoff-independent uniqueness apply. Consequently zero
subsidiary initial value and velocity imply C=0 in the full angular bundle,
provided the exact boundary/gauge completion and physical initial constraints
have supplied those data as in the earlier checkpoint.

This establishes the immediate angular-uniform linear target for smooth
compatible data on the controlled background collar, including the full
coupled metric/scalar and subsidiary fields. It does not prove an exact
nonlinear family. Nonlinear boundary straightening and the additional
quasilinear obligations are recorded separately.

The companion executable is verify_uniform_robin.py. Its saved results are
verification_results.json. Those checks establish algebra and symbol
identities, including deliberate corruption controls; they do not replace
the energy and limiting arguments above.
