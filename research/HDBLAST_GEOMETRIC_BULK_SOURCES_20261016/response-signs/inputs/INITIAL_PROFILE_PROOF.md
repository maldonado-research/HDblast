# Actual initial profiles, complete packet cells and a shared background approximation

The initial fields are those already selected by the accepted physical
constraint and local relative-harmonic completion. We evaluate them on
t=0, ell<=z<=2ell, ell=10^-6, times the full sphere. All coefficients are
second parameter derivatives (epsilon^2/2 convention). Earlier frozen files
are unchanged.

## Correlated actual functions and numerical arithmetic

At every r-cell, X(ell r) belongs to the authenticated shell interval plus
ell*r times the authenticated actual ODE image, X=(rho,eta,p,a). The original
whole-collar first-exit proof supplies this enclosure. We check rho>0,p<0.
For k=rho(5a/2-U1/p), its retained interval is strictly negative; alpha(0)=1
and alpha=exp(integral_0^z k) is enclosed by exp(z*kbox). All exponentials
use alternating rational Taylor bounds of degree28 at argument<=1 after
halving, then outward interval squaring. This does not introduce binary
floating values into mathematical decisions. The dyadic arithmetic grid is
2^112 and every primitive rounds outward.

The actual ODE is differentiated by normalized Taylor recurrence through
order four at each possible cell point. C=alpha/rho^2, B=3alpha/(rho^3 p),
h=U1/p and k are formed from the same state expressions. Thus these are
bounds on one actual smooth background, never choices of mutually unrelated
coefficient tables. The conditional shell uncertainty remains present.

For the numerical initial packet use s=1/(1/4-(r-3/2)^2), v=3-2r.
The exact derivative polynomials Q^(n)=Pn(s) exp(4-s) times v for odd n
are evaluated through n=4. Centering each polynomial at its cell midpoint
reduces cancellation while preserving an exact polynomial identity. Division
by ell^n*n! converts these to physical normalized q jets. At both endpoint
tails of width1/64, s>=4096/63 and
|Q^(n)|<=s^(2n)exp(4-s), n<=4. This bound decreases because s>=2n.
The compact profiles have their exact smooth zero extension outside (1,2).
The 496 interior cells plus two tails form a checked gap-free cover of
[1,2]. These are whole-cell enclosures, not a sampled maximum.

## Physical and harmonic initial formulas retained

The program evaluates the accepted cancellation-preserving formulas
u=Cq, v=B(q_z+kq), ut=exp(-z)u,
vt=exp(-z)[v-3u/(rho p)]. Let dy=-rho^-1 dz, vy=dy(v).
The scalar Hamiltonian source A, D, f_l, zeta_l, pure trace K profile and
Pi profile are the exact reduced formulas in the retained physical initial
producer. In particular zeta_l=-3(f_l)_y/(2p), and its spatial derivative
uses the packet's fourth derivative. No scalar initial value is set to zero.
The common quadratic momentum cancellation and the -2vy^2 term are retained.

The full metric values include k00=-2u^2Y^2 and kab=f delta_ab.
The five normal derivative classes and scalar normal derivative are precisely
those in the retained initial spacetime frame proof. In particular

  (nabla0 k)0y=-26 U Uy-(3/2)f_y-3af,
  partial0(k0y)=(nabla0 k)0y-2a U^2,
  nabla0(phi2)=Pi2+2U Pi1.

This restores the mixed shift/lapse velocity terms when converting to actual
scaled component partial_x velocities. Other component classes have no
additional connection term because k0y=kyy=0 initially. The radial frame is
parallel, so partial_r=-ell*rho partial_y. The reported radial partial bounds
may use the larger full spatial covariant bound and remain valid.

Angular amplitude evaluation uses |W0|<=1, |W2|<=1/2, |W4|<=7/12.
Their first sphere derivative factors are 0,1,5/2. For the last factor,
W4=Y^2-W0/6-W2/2, |DY^2|<=2 and |Ds|<=1 for
s=n1^2+n2^2, giving |DW4|<=5/2. The radial-angular projector connection
contributes |a f| to the spatial metric bound. These accepted pointwise
factors apply over the whole sphere, not just at chosen angles.

The output records component upper bounds and continuous scalar mean,
l=2/l=4 amplitude cells. The sphere probability mean is zeta0 because W0=1
and the other harmonics have mean zero. The cell records can be plotted as
bands; their midpoints are display values and are not a solved background.
The program certifies no interpolant approximation error for initial data.

## A separate concrete background approximator

evaluate_background_approximation.py defines one shared affine rational
state Xhat(z)=mid(shell X)+z F(mid(shell X)), and
alphahat=1+z k(mid(shell X)). All derived coefficients C,B,h,k and the
static part of the master potential are rational functions of this same
state. This is a concrete correlated approximation, not independent
selections inside old interval boxes. Every actual-minus-hat derivative
through scaled radial order four is enclosed over a complete 64-cell collar
cover [0,8ell]. Actual and hat are differentiated independently before
interval subtraction. All temporal potential terms cancel exactly because
the same -35 tanh^2(t)/4 is retained in both.

Its error is dominated by inherited shell uncertainty and the deliberately
low affine order. It has not been used to tighten the master/gauge errors
or to evaluate the quadratic Einstein-scalar forcing. It establishes an
explicit input representation with honest absolute derivative errors.
A validated regular-branch shell solve and higher-order correlated Taylor
approximation would materially improve its accuracy.
