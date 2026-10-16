# Full second coefficient on the retained causal collar

Conditional mathematical construction, 16 October 2026. This note defines an
actual second coefficient by a convergent integral series, conditional on the
retained smooth regular background and actual first master. It uses the actual
first relative-harmonic metric/scalar pair and the previously constructed
physical/harmonic second initial data. It is a formal second parameter
coefficient, not an exact nonlinear family or a finite-amplitude remainder.
The numerical joins use authenticated accepted source certificates. Independent
reviews here are internal AI-assisted reviews, not external specialist review.

## 1. Objects, angular space and norm

Set ell=10^-6, x=t/ell, r=z/ell, dy=-rho dz and
D={0<=x<=3, r>=0, x+r<=8}. The unknown is U=(k,w), with factorial
convention g=g0+epsilon q+epsilon²k/2 and phi=phi0+epsilon eta+epsilon²w/2.
Here (q,eta) is the *actual transformed first pair*, including all gauge
metric/scalar terms. Let Qred and Qphi denote its direct affine Hessians of
Ered=Ric-dphi dphi-(2/3)Vg-nabla_(A C_B) and S=Box phi-V'. The parentheses
contain one half. No covariance expression requiring a separately bounded
fourth gauge derivative is used.

Use the radially parallel orthonormal background frame
e_y=partial_y, e_0=rho^-1 partial_t, e_a=(rho cosh t)^-1 E_a.
The positive fibre norm is |k|_F²=sum_(A,B)|k_AB|², counting symmetric
off-diagonal components twice, plus |w|². Raising Lorentzian indices changes
signs and is an isometry for this positive norm. Integrate against probability
measure dvol_kappa/vol(S3); the resulting angular Hilbert norm is denoted |.|_H.
This is a fixed positive comparison norm, not a gravitational energy.

The first state is the SO(4) representation (1,1). The elementary identity
Sym²(V tensor W)=Sym²V tensor Sym²W plus wedge²V tensor wedge²W, with
Sym²(spin1)=spin0 plus spin2 and wedge²(spin1)=spin1, gives

  R={(0,0),(1,1),(2,2),(2,0),(0,2)}.

Its representation dimensions add to45. This number is NOT the dimension of
the full tensor-section space. The full metric/scalar space H consists of:

* five scalar copies: k_yy,k_y0,k_00,angular metric trace,w, each with l=0,2,4;
* two angular one-form copies k_ya,k_0a, with scalar gradients l=2,4;
* one angular trace-free metric copy, with scalar Hessian-TF l=2,4;
* both transverse-traceless l=2 tensor orientations.

Thus dim H=5(1+9+25)+2(9+25)+(9+25)+10=287. The Dirichlet subspace k_y mu
has dimension35+34=69, and its orthogonal Robin complement dimension218.
All these are sections of the tensor bundles, not scalar coordinate-component
cutoffs. The absence of coexact-vector representations follows from R.

For completeness the homogeneous-bundle rough Laplacian has eigenvalue
2[jL(jL+1)+jR(jR+1)]-s(s+1) on isotropy spin s. For the scalar, one-form,
and trace-free tensor bundles s=0,1,2 respectively. The relevant eigenvalues
are l(l+2), l(l+2)-2, l(l+2)-6, and6 for the two TT l=2 sectors.
Consequently ||Delta_kappa||<=24 and ||D_kappa||<=sqrt24<5 on H, with
the derivative norm summed over its angular slot. These facts also follow
directly from scalar/gradient/Hessian commutation on the unit three-sphere.
The scalar Hessian-TF l=0 sector is absent, so no exceptional inversion occurs.

Background coefficients and every linear operator below are SO(4)-equivariant.
The first solution and first gauge are equivariant linear images of the same
fixed quadrupole. Their symmetric quadratic sources, including all boundary
sources, lie in R. The inherited explicit second initial profiles have only
l=0,2,4 scalars and their gradients. Hence the entire problem preserves H.
No TT sector is discarded merely because its initial data vanish.

## 2. Exact reduced linear operator

Define

  L_R(k)_AB = (Ric_A^C k_CB + Ric_B^C k_CA)/2 - R_ACBD k^CD,
  dC[k]_A = nabla^B k_AB - (1/2)nabla_A tr(k).

With the curvature convention of the retained Ricci formula, direct variation
of Ricci and subtraction of nabla_(A dC_B) give

  Box k = 2 L_R(k) - 2(dphi0 tensor dw + dw tensor dphi0)
          -(4/3)V k -(4/3)V' w g0 + 2Qred,
  Box w = V'' w + k^{AB}(nabla²phi0)_AB + dC[k]^A phi0_A - Qphi.       (1)

These equations define every coupling; they do not assume a scalar-only metric
ansatz. In particular the scalar equation retains the derivative coupling dC[k].

Here is an explicit frame recipe for the unit-principal matrices, avoiding any
choice of poorly conditioned harmonic basis. Write c=cosh t, Ht=tanh t,
h=Ht/rho. Let Omega_0 and Omega_a be the connection actions on covariant
two-tensors induced by

  nabla_0 e_y=a e_0, nabla_0 e_0=a e_y,
  nabla_a e_y=a e_a, nabla_a e_0=h e_a,
  nabla_a e_b=h delta_ab e_0-a delta_ab e_y

after the intrinsic sphere connection has been assigned to D_kappa. They act
as zero on w. Their positive operator norms obey
||Omega_0||<=2|a| and ||Omega_a||<=2(|a|+|h|). There is no Omega_y.
The frame wave operator, multiplied by rho², is exactly

  rho² Box = -partial_t²+partial_z²
    -(3Ht I+2rho Omega_0)partial_t -3a rho partial_z
    +c^-2 Delta_kappa +(2rho/c)sum_a Omega_a D_a +rho² Z,
  Z=-Omega_0²-3h Omega_0+sum_a Omega_a².                              (2)

The term partial_t Omega_0 vanishes because a depends only on y. The traced
frame connection supplies +4a partial_y-3h nabla_0 before changing y to z;
partial_y²+4a partial_y=rho^-2 partial_z²-3a/rho partial_z.
In the angular square the structural Omega_a is sphere-parallel, with its
angular label included, so no missing derivative-of-Omega term occurs.

Insert (1) in (2), use partial_t=ell^-1 partial_x and
partial_y=-(ell rho)^-1 partial_r, and collect the derivatives and values.
This defines the exact, finite-dimensional operators A,B,D0 and source F:

  U_xx-U_rr=A U_x+B U_r+D0 U+F,
  F=(-2ell²rho² Qred, ell²rho² Qphi).                                (3)

Equations (1),(2), including the definition of dC and the explicit frame
connection, are a complete algebraic recipe for A,B,D0. A separate287-square
table is unnecessary: all actions are fixed covariant tensor contractions and
bounded operators on the stated H.

## 3. Actual coefficient bounds

The authenticated whole-collar first-exit report gives rho in[129,131],
|a|<=14/125, |a_y|<1/10000, |p|<=7/50000,
-9/100000<=phi0-1<=-8/100000 and
|Ht|<=3/1000000. These bounds cover the entire8ell collar, rather than the
older shorter prose interval. The registered potential polynomial is

  V(1+v)=-2/27+(14/9)v²+(50/27)v³-v⁴/6-(4/9)v⁵-(2/27)v⁶.

Direct rational polynomial majorants give |V|<1, |V'|<1, |V''|<4;
indeed |V'|<1/1000. The scalar equation gives |p_y|<1/500, and thus
||nabla²phi0||_F<=|p_y|+2|ap|<1. Warped-product curvature has sectional
coefficients -(a_y+a²) and rho^-2-a², so each frame curvature component
is bounded by R*=a²+|a_y|+129^-2<13/1000. Its Ricci operator is bounded
by4R*, and its contraction on symmetric tensors by25R*. Hence ||L_R||<1.

The following fully explicit positive majorants suffice. Let R=131,
alpha=14/125, P=7/50000, h*=3/(1000000*129), L=ell, and set

  A*=L(3*3/1000000+4R alpha)+7LRP,
  B*=3LR alpha+7LRP,
  Z*=(2alpha)²+3h*(2alpha)+3[2(alpha+h*)]²,
  D*=L²[24+40R(alpha+h*)+R² Z*
          +R²(2+4/3+4+4+1+P(8alpha+6h*))+30RP].                  (4)

The metric derivative coupling contributes at most4LRP per x or r derivative;
the scalar dC coupling at most3LRP. The angular part of these couplings is
bounded by30RP L². The dC zeroth-order connection contribution is at most
P(8alpha+6h*) before rho²L² scaling. The other terms in D* are, in order,
the angular Laplacian, angular connection cross term, Z, twice L_R, potential
metric term, potential metric/scalar cross term (||g0||_F=sqrt5<3), scalar
potential, scalar Hessian contraction, and the dC connection contribution.
Using sqrt3<2 and sqrt24<5 gives the factor40 in the angular cross term.
The use of sums of block norms is conservative and applies to the full Hilbert
direct sum. Exact rational evaluation verifies

  ||A||<=A*<10^-4, ||B||<=B*<10^-4, ||D0||<=D*<10^-6.              (5)

## 4. All primal rows and their frame conversion

Let gamma be the background induced shell metric, K=k_yy and
T=tr_gamma(k_tangent). The independent boundary audit establishes the
copy-even physical Hessians QR,QP and normal gauge source GN. In coordinate
components the rows are

  k_y mu=0,
  partial_y k_mu nu-2a k_mu nu=(aK+sigma'w/3)gamma_mu nu-2QR_mu nu,
  w_y=pK/2-sigma''w/2-QP,
  K_y=-4aK+2aT+(4sigma'/3)w+GN.                                  (6)

Their sources are evaluated from the actual first pair; GN is exactly
-2tr_gamma QR-2(Q_C)_y. Thus the last row is C2_y=0, with
C2=dC[k]+Q_C, not an additional physical condition. The normal is outward
partial_y, sigma'= -2p, sigma''=4phi0.

In the parallel frame m=k_tangent the first derivative in the second row
of(6) becomes precisely m_y: differentiating the frame factors subtracts
2a k_tangent. Therefore the homogeneous frame boundary operator is

  m_y=(aK+sigma'w/3)eta_tangent,
  w_y=pK/2-sigma''w/2,
  K_y=-4aK+2a tr_eta(m)+(4sigma'/3)w.                              (7)

In particular the tangential trace-free block has coefficient zero. Writing
s=tr_eta(m)/2, its scalar trace block is
[[0,2a,2sigma'/3],[4a,-4a,4sigma'/3],[0,p/2,-sigma''/2]]
on(s,K,w). The matrix with2a in the first diagonal entry belongs to unconverted
coordinate components and cannot be inserted into this frame evolution.

Let P_D,P_N be the orthogonal projections of H defined above. On r=0 the
complete problem is

  P_D U=0,  P_N U_r=M P_N U+G,
  M=-ell rho B_y, G=-ell rho(-2QR,-QP,GN)                          (8)

in the order(m,w,K) on the Robin complement. Every coupling is a value
coupling; there is no tangential or time derivative of a second unknown.
The maps eta_tangent and tr_eta have norms2. Thus

  ||M||<=ell R(10alpha+2|sigma'|+P/2+|sigma''|/2)
         <=ell R(10alpha+4P+P/2+2)<1/2000.                       (9)

The four independent Dirichlet fibre components and twelve Robin fibre
components are all retained. There are no extra tangent harmonic primal rows.

## 5. Closed convergent series on the physical causal triangle

Define J=P_N-P_D. For any angular-H-valued field f on r>=0 extend it by
f_ext(x,-r)=J f(x,r). For a source F define the reflected Duhamel operator

  (K F)(x,r)=1/2 integral_0^x integral_(r-x+s)^(r+x-s)
                         F_ext(s,v) dv ds.                       (10)

For Robin data g with values in P_N H, define

  (L g)(x,r)= - integral_0^max(x-r,0) g(s) ds.                      (11)

Then P_D K F|r=0=0, P_N partial_r K F|r=0=0, while partial_r Lg|r=0=g.
The zero extension in(11) has continuous first derivatives because every
boundary source here is zero near the corner. Arbitrary iterates are taken
in the closed C1 subspace which vanishes on the fixed region x+r<1; the maps preserve
it by causality. Formula(10) remains valid for a reflected Dirichlet source
whose values jump at r=0: its integral has continuous first derivatives.

Let U_free be the ordinary d'Alembert solution from the inherited initial
values U0 and scaled velocity V0, reflected with J. Their flat zero corner
makes this free solution smooth. Explicitly

  U_free=[U0_ext(r+x)+U0_ext(r-x)]/2
           +(1/2)integral_(r-x)^(r+x) V0_ext(v)dv.

The exact coefficient is the fixed point

  U=U_free+K(A U_x+B U_r+D0 U+F)+L(M U|r=0+G).                    (12)

No arbitrary selections from interval coefficient/source boxes occur: all
functions in(12) are the actual fields defined in the earlier constructions.
The boxes only bound the norm of this fixed operator. Every characteristic
point used for a point(x,r) in D satisfies s+|v|<=x+r<=8. The boundary lift
uses s<=x-r. Consequently no outer-boundary extension beyond the controlled
collar is needed or affects the solution.

Give C1(D;H) the norm

  ||U||_X=max(sup_D |U|_H/3,sup_D |U_x|_H,sup_D |U_r|_H).

For continuous source f, differentiating the integration endpoints yields
||Kf||_X<=3 sup|f|_H. The value estimate is x² sup|f|/2 and each first
derivative estimate x sup|f|; J is an isometry, so no component-count factor
is hidden. For continuous corner-compatible g, ||Lg||_X<=sup|g|_H.
The boundary trace of U has norm at most3||U||_X. Thus the linear part T of
the right side of(12) satisfies

  ||T||<=3(10^-4+10^-4+3*10^-6+1/2000)
         =2109/1000000=:q<1.                                    (13)

Put b=U_free+KF+LG. The explicit complete second coefficient is

  U=sum_(j=0)^infinity T^j b,
  ||U-sum_(j=0)^(N-1) T^j b||_X<=q^N ||b||_X/(1-q).             (14)

This is a constructive integral definition of the actual response, with a
continuous error bound for every truncation. A residual for any compatible
candidate V satisfies ||U-V||_X<=||V-b-TV||_X/(1-q). The complete exact
series has zero integral residual; no floating-point U2 waveform is asserted.
Uniqueness follows in the same norm. The formula uses physical half-line
reflections and the explicit boundary lift, not an even extension across all
boundary rows or a presumed symmetric Robin energy cancellation.

## 6. Initial/source norm joins and a coarse actual response cap

The quantitative initial checkpoint gives adapted-frame component caps
|k0|<1.06e24, |nabla k0|<6.4e29, |w0|<5.13e33, |dw0|<3.06e39.
The first spacetime initial one-jet used there agrees with the actual first
gauge wave: xi=xi_t=0 and its radial initial acceleration equals the completed
shift velocity. In particular the inherited mixed second metric velocity and
g2_tt=-2rho²(uY)² are retained. The second scalar data are not set to zero.

Positive Frobenius norm gives |k|_H<=5 times an individual component cap.
Since the frame is radially parallel,
U0_r=-ell rho nabla_y U0 and
(k0)_x=ell rho(nabla_0 k0-Omega_0 k0). The scalar has no Omega_0.
These exact conversions imply

  sup|U0|_H<6e33, sup|U0_r|_H<5e35, sup|V0|_H<5e35,
  ||U_free||_X<10^36.                                           (15)

The authenticated accepted bulk and shell reports imply the deliberately rounded envelopes
|Qred_AB|<10^21, |Qphi|<4e15, |QR_mu nu|<1.5e9, |QP|<4.4e9 and
|(Q_C)_y|<6000. The trace join gives |GN|<1.3e10. Hence

  sup|F|_H<=ell²R²(10*10^21+4e15)<2e14,
  sup|G|_H<=ell R(8*1.5e9+4.4e9+1.3e10)<4e6.                    (16)

Combining(13)-(16) gives ||U||_X<1.02e36. The bound is intentionally coarse
and dominated by the previously certified initial scalar one-jet. It is an
upper bound, not the evaluated amplitude or a finite-amplitude validity range.

For an angular pointwise cap, the reproducing kernel of H has constant trace
dim H=287 under the transitive SO(4) action and probability normalization.
Its evaluation norm is at most sqrt287<17. Therefore throughout D times S3,

  every adapted-frame component |U|<6e37,
  every scaled frame-component derivative |U_x|,|U_r|<2e37.        (17)

These are value and first-derivative estimates only. Formula(14) gains the
same factor17 (and3 for values) when used as a pointwise tail estimate.
It does not bound second time/normal derivatives or a high-Sobolev norm.

## 7. Qualitative smoothness, subsidiary equations and physical completion

C0 forcing alone proves the C1 fixed point in(12), not a classical C2 field.
The retained master has a smooth/flat initial packet and smooth background
coefficients. Its boundary is dynamic, not an ordinary Robin row. To justify
its qualitative smoothness, use the exact retained scaled relations
H_master=W_r-ell gamma y, y=W(x,0),
y''=(ell/mu)H_master-A_master(x)y, with mu>0 and
A_master=ell²[5/2-(35/4)tanh²(ell x)]. Set z=y' and let
p_master=W_x+W_r be the outgoing characteristic at the shell. Then the
boundary is the finite smooth trace ODE

  y'=z,
  z'=(ell/mu)(p_master-z-ell gamma y)-A_master y,
  m_master=2z-p_master.                                         (18a)

Together with the scalar bulk characteristic equations, this is a local
Volterra system: outgoing p is its initial value plus a time integral;
the solved two-state trace ODE supplies the incoming m. On a short interval
the source-integral Lipschitz constant is O(delta), and the trace-state
difference gains another time integral. Thus it is a contraction for small
delta, with no need to assume a favorable sign for gamma or A_master. Smooth
coefficient/data differentiation and fixed flat corner give the same transport/
trace-ODE regularity induction detailed below. This proves smoothness of the
existing actual master by uniqueness and continuation; its numerical C3 box
is not the reason for the qualitative conclusion. The exact dynamic equations
are in the retained interior-master `jets/HIGHER_JETS.md`.
The already constructed first gauge solves
a finite system of the form(3) with homogeneous mixed boundary rows and smooth
forcing from that master. The following elementary regularity argument applies
to it and then to(12), giving qualitative smoothness without asserting a new
numerical bound for higher gauge/source derivatives.

Set p=U_x+U_r, m=U_x-U_r. The equations become
p_x-p_r=f, m_x+m_r=f and U_x=(p+m)/2, where
f=(A+B)p/2+(A-B)m/2+D0 U+F. The outgoing p trace is its initial value plus
the source integral along the characteristic ending at r=0. For Robin trace
v=U_N(x,0), v'=p_N(x,0)-Mv-G and incoming m_N=p_N-2Mv-2G.
The Dirichlet incoming trace is m_D=-p_D. These identities are equivalent to
(10)-(12). Their characteristic integral system has a local contraction on
continuous(p,m,U,v), obtained by first solving the finite trace ODE; bounded
coefficients suffice after choosing a short time interval.

Differentiate this first-order characteristic system in x. The highest
derivatives satisfy the same transport/trace system, and differentiated
coefficients and data multiply only lower derivatives. Their initial values
are determined by the differential equations from smooth initial radial jets.
The boundary corner identities are exactly the same differentiated equations;
they hold to every order because all data and forcing vanish on an open initial
shell collar. Difference quotients in these integral equations converge to the
unique differentiated solution: subtracting the differentiated candidate gives
the same contraction bound, with an inhomogeneity tending uniformly to zero.
This justifies, rather than presupposes, each time differentiation. Repeating
gives all time derivatives continuous. The transport equations recover all
radial derivatives, successively. Repeating on finitely many bounded time
subintervals covers D. Angular smoothness follows from finite dimensionality.
Thus the actual first gauge and actual second coefficient are smooth. Constants
for this qualitative bootstrap may depend on higher derivatives which the C3
numeric certificate has not bounded. No such constants enter(13)-(17).

One can verify this induction directly on the Picard series, including the
finite initial-jet offsets. For physical smooth f,
partial_x Kf=S(f(0,.))+K(partial_x f), where S is the mixed free wave with
zero displacement and indicated velocity. This follows by integration by parts
in source time, not by differentiating a reflected coefficient in r. Since g
vanishes near x=0, partial_x Lg=L(partial_x g). At differentiation order m,
the highest time derivative V_n=partial_x^m U_n therefore satisfies
V_(n+1)=T V_n+R_n, with the SAME operator T. R_n consists of K and L applied
to coefficient derivatives times lower-level fields partial_x^j U_n, j<m,
and their first x or r derivatives, plus smooth free waves. For example the
term m A_x partial_x^m U_n is controlled as the x derivative of the lower-level
field partial_x^(m-1)U_n; it is not claimed to have strictly smaller total
time order. Those fields converge in C1 by induction, so the K inputs converge
in C0, as do the L trace inputs, and hence R_n converges in C1. The free-wave initial
jets are eventually constant in n: the iteration has the fixed U0,V0 initially,
and each new time jet is obtained from U_xx=U_rr+A U_x+B U_r+D0 U+F using
already fixed lower time jets. One must discard this finite offset, rather than
claim every early iterate has every final initial acceleration. Thus R_n converges
in C1. The elementary recurrence with ||T||<1 then makes V_n converge in C1
(subtract its limiting fixed point, and split the geometric convolution of
R_n-R into an initial finite part and a uniformly small tail). Uniform convergence
of functions and first derivatives identifies the limits with actual derivatives.
Recovering repeated r derivatives from the PDE gives all mixed orders. Every
iterate is smooth on the physical side: the only reflection transition originates
at the corner, where the fixed zero wedge eliminates all jump jets. This also
explains explicitly why a mere C0 source bound would not prove this bootstrap.

Now all classical subsidiary traces are justified. The exact trace-reversed
Bianchi identity, differentiated twice using E0=E1=S0=S1=C0=C1=0, yields

  (Box0+Ric0) C2=0, C2=dC[k]+Q_C.                                (18)

The normal row in(6) gives C2_y=0 at the shell. The physical metric/scalar
rows and Codazzi give E2_na=0; using the reduced equation and C2_y=0 yields
partial_y(C2)_mu-2a(C2)_mu=0 in covector coordinate components. In the
radially parallel frame this is partial_y(C2)_mu-a(C2)_mu=0. The inherited
initial gauge one-jet gives C2=0 initially; the physical Hamiltonian and four
momentum constraints, using the reduced accelerations from(3), give
partial_t C2=0 initially. These accelerations are not extra Cauchy data.

The finite angular covector subsidiary system is another smooth normally
hyperbolic system with one normal Dirichlet block and four Robin fibre rows.
The same characteristic/trace-ODE contraction on short intervals proves
uniqueness for its zero data; it does not require a favorable Robin sign.
Therefore C2=0. Equations(1) then solve all natural second Einstein/scalar
rows with their actual Hessian sources. Equations(6) supply every physical
copy-even shell row. Taking identical copies supplies all induced matches;
the inherited exact force/Codazzi identities give the shell force and
conservation rows. No independent incoming constraint row is appended.

Finite propagation is explicit in(10)-(12): sources and first/second initial
data are zero beyond r=2+x, so U is zero there. The retained8ell causal
triangle covers the entire construction through x=3 and through first contact.

## 8. Scope and verification

The new mathematical content here is the full finite-angular integral response
and its explicit inverse bound for the actual common-gauge data, not a numerical
time stepping approximation. The accompanying rational checker verifies the
angular multiplicity/eigenvalue counts, actual inherited coefficient joins,
operator majorants, initial/frame conversions and source-to-response arithmetic.
It does not replace the displayed all-point operator/Volterra/regularity proof,
nor independently replay every prior master/gauge/initial-data certificate.

Primary inputs: the first-gauge COMPLETED_HANDOFF and FIRST_HARMONIC_FIELDS
note; its GLOBAL_GAUGE_COEFFICIENTS_NORMAL report; the explicit physical
second-initial-jet and initial-harmonic-completion notes; the second-initial
norms handoff; and the current bulk, shell and boundary audit reports. All are
read only. No old frozen capsule or manuscript is changed by this construction.

The coefficient result is conditional on the retained background/first-solution
model hypotheses. It makes no exact nonlinear existence, Taylor-remainder,
all-angular-frequency estimate, physical-origin, external-priority or journal
readiness claim.
