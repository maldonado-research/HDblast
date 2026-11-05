# The registered first master through past shell contact

4 November 2026, US Pacific (5 November UTC). Conditional mathematical
continuation of the *linear first master*, under the retained regular background
and reconstruction identities. It constructs that selected linear packet
through past shell contact, gives an estimating-energy comparison and a uniform
upper bound for its first physical scalar coefficient, and proves qualitative
activation. It does not complete second constraints, a nonlinear past family,
or a physical incoming-state law.

## 1. Actual dynamic boundary and verified signs

For the registered copy-even quadrupole, the retained reference fixes

```
w_tt-w_zz+mathcalW(t,z)w=0,
-w_z+mu_b*w_tt+beta(t)w=0 at z=0,
mu_b=1/(rho_b B_b)>0,
B_b=p_y,b/p_b+sigma''_b/2,
beta(t)=C_b/2+mu_b*(35 sech(t)^2/4-25/4).
```

Here C is the longitudinal master coefficient, not the small matrix coefficient
in the different five-mean problem. In the registered model, sigma''=4phi and
the static ODE gives p_y/p=h-4a, h=V_bulk'/p. Therefore

```
B_b=h_b-4a_b+2phi_b,
C_b/2=rho_b*(h_b-5a_b/2).
```

The companion arithmetic forms these quantities directly from the raw
analytic shell intervals for rho_b,p_b,V_bulk',a_b,eta_b=phi_b-1. It checks
the identity sigma''=4(1+eta) against the stored second-tension interval.
This reconstructs mu_b from its definition rather than assigning it by name.
The raw enclosure gives approximately .00205977<mu_b<.00227167; safe exact
bounds are

```
1/1000<mu_b<1/400,
189<C_b/2<300-1/160,
180<beta(t)<300 for every real t.
```

The last bound uses 0<=sech(t)^2<=1, so its temporal summand lies between
-(25/4)mu_b and (5/2)mu_b. Positivity of the spring coefficient is decisive.
None of these is a sign assumption about an arbitrary Wentzell boundary.

The fresh precontact calculation already supplies the actual static collar
0<=z<=8ell, ell=10^-6, with rho in[129,131] and

```
mathcalW(t,z)=rho(z)^2 Qbg(z)-(35/4)tanh(t)^2,
1/3<Qbg<6.
```

Using |tanh(t)|<=1 gives the stronger all-real-time bounds on that same collar:

```
mathcalW>129^2/3-35/4=22153/4>5500,
mathcalW<131^2*6=102966<103000.
```

The static branch and validity of its analytic shell box remain retained
premises. This is explicit time algebra, not a transfer of evolved fields or
the future five-mean operator caps to negative time.

## 2. Correct scaling and the positive mass/spring energy

Set x=t/ell,r=z/ell,U(x,r)=w(ell*x,ell*r), q(x)=U(x,0). Then

```
U_xx-U_rr+V(x,r)U=0,
-U_r+Lambda*q_xx+b(x)q=0,
Lambda=mu_b/ell, b(x)=ell*beta(ell*x),
1000<Lambda<2500, 180ell<b<300ell,
5500ell^2<V<103000ell^2.
```

The derivative U_r points into the collar. For a regular finite-energy solution
on the half-line define

```
E(x)=1/2 integral_0^infty[U_x^2+U_r^2+V U^2]dr
       +Lambda*q_x^2/2+b*q^2/2.                           (1)
```

All three bulk and both boundary contributions are nonnegative. Integration by
parts in the bulk produces the lower-endpoint flux -q_x U_r(0). The dynamic row
U_r(0)=Lambda*q_xx+bq cancels it exactly against the derivative of the two
boundary summands. Consequently

```
E'(x)=1/2 integral V_x U^2 dr+b_x q^2/2.                  (2)
```

Neither boundary energy may be discarded. Changing the sign of the normal row
would change this cancellation and must not be hidden in a boundary convention.

## 3. Separable time dependence and the existing endpoint-mass realization

The positive endpoint-mass Hilbert space and closed-form construction are
established mathematics, already used in the project's
`HDBLAST_REGULAR_SHELL_SPECTRA_20261013` domain note, sections 2--3 and 6.
After substituting `p_y/p=h-4a`, that note has Liouville bulk potential
`V_L=rho^2 Qbg-5/2`,
endpoint coefficient `s(0)=C_b/2`, and mass `d_b=mu_b`. Before scaling, the
static expression below is its expression shifted by `(5/2)I`; this also
shifts the boundary spring to `C_b/2+(5/2)mu_b`.
The snapshot `inputs/PRIOR_SPECTRAL_DOMAIN.md` records the earlier construction.
The contribution here is its application to this actual packet and short
past-contact domain, with quantitative first-scalar and energy bounds. It is
not a new selfadjointness or general energy principle. The auxiliary spatial
extension below is not identified with the earlier full regular-cone operator
or its global spectrum.

Put

```
omega(x)=-(35/4)ell^2*tanh(ell*x)^2,
V(x,r)=V0(r)+omega(x), V0=ell^2*rho^2 Qbg,
b(x)=b0+Lambda*omega(x),
b0=ell*(C_b/2+(5/2)mu_b).
```

Thus b_x=Lambda*V_x exactly; V_x is independent of r. A useful rigorous
construction is in the real Hilbert space

```
H=L2(0,infty) direct-sum R,
||(u,q)||_H^2=integral u^2 dr+Lambda*q^2.
```

The form domain is `{(u,u(0)):u in H1(0,infty)}`. Its components are
trace-linked on this domain although general Hilbert vectors have independent
components. Density follows
by approximating a desired independent boundary component with a thin smooth
boundary layer, whose bulk L2 norm tends to zero. The static positive form is

```
a0(u,u)=integral[u_r^2+V0u^2]dr+b0*u(0)^2.
```

For closedness, a sequence Cauchy in `a0+||.||_H^2` converges in bulk H1,
and its trace therefore converges to the trace of that same limit. This norm
is equivalent to the H1/trace graph norm. The positive-form representation
gives a selfadjoint A0 on
`D(A0)={(u,u(0)):u in H2(0,infty)}`. It acts

```
A0(u,q)=(-u_rr+V0u, [-u_r(0)+b0q]/Lambda), q=u(0).
```

Compact bulk tests give `-u_rr+V0u in L2`, hence u in H2, and tests with
nonzero trace identify the endpoint component. Conversely integration by
parts for H2 functions gives the representation identity. The +u_r(0)v(0)
from the bulk cancels the displayed negative boundary derivative; endpoint
products at infinity vanish by H2/H1 integrability. There is no separate
static Robin condition imposed on D(A0). In D(A0^2), the trace identity for
A0Z additionally reads
`-u_rr(0)+V0(0)q=[-u_r(0)+b0q]/Lambda`; this is the required higher
compatibility, not an extra datum at D(A0) level. The wave equation
with dynamic boundary is precisely

```
Z_xx+[A0+omega(x)I]Z=0, Z=(U,q).                         (3)
```

For s=-x, the data are Z(0)=(Q,0), Z_s(0)=(-ell exp(-ell r)Q,0), with the same
flat compact Q supported in[1,2]. They belong to every power domain of A0:
all their boundary jets vanish, and static differentiation preserves support
away from the boundary. Every initial A0 power has zero endpoint component,
so all initial differentiated compatibility conditions hold. The scalar
bounded omega commutes with A0. The wave
group formula therefore provides the backward Volterra construction. Here
omega is even, so `omega(-v)=omega(v)` after x=-s:

```
Z(s)=cos(s sqrt(A0))Z0+A0^(-1/2)sin(s sqrt(A0))Z1
 -integral_0^s A0^(-1/2)sin((s-v)sqrt(A0))*omega(v)*Z(v)dv. (4)
```

Since the sine propagator has H-operator norm at most s-v, the integral
operator on 0<=s<=3 has norm at most

```
9*sup_|x|<=3|omega(x)|/2
 <=(2835/8)ell^4<1.                                    (5)
```

The estimate uses |tanh(ell*x)|<=ell|x|. A convergent Neumann series solves
(4) uniquely. Since the sine/cosine operators and scalar omega commute with
A0, the same construction works in every power-domain graph norm with the
same bound. The equation supplies two time derivatives in each lower power
domain; differentiation and smooth omega iterate this assertion. Bounded
smooth V0 and one-dimensional elliptic/Sobolev regularity yield joint
smoothness up to r=0. This supplies the classical traces and dynamic
boundary, rather than inferring smoothness from energy alone. It is an
application of standard linear form/wave machinery, not a numerical PDE
evolution or a nonlinear theorem.

For this argument extend only the *auxiliary static coefficient* V0 beyond the
verified collar by a smooth convex interpolation to 6000ell^2: keep it exact
for r<=6, and constant for r>=7. This preserves
5547ell^2<V0<102966ell^2 and all time identities. No physical background is
asserted outside the actual collar. The following finite-propagation argument
makes that extension irrelevant on this slab.

## 4. Support, shell compatibility and the postcontact domain

For backward s, V_s<=0 and b_s<=0. The identities can first be integrated
with compact radial cutoffs and then passed to infinity by H2/H1
integrability from the graph construction. In the outer region r>=2+s, the ordinary
moving-cut energy identity has boundary contribution

```
-U_s U_r-e=-[U_s+U_r]^2/2-VU^2/2<=0.
```

Its initial energy is zero, so the constructed solution vanishes there.
Before s=1, the inner-region energy on 0<=r<=1-s includes the same mass and
spring at r=0. Its physical boundary flux cancels, while its moving outer
contribution is

```
U_s U_r-e=-[U_s-U_r]^2/2-VU^2/2<=0.
```

It too starts at zero and hence remains zero. Consequently

```
support U(-s,.) subset[1-s,2+s] for0<=s<1,
support U(-s,.) subset[0,2+s] for1<=s<=3.                 (6)
```

For the entire closed slab the outer support is at most r=5, strictly inside
the actual depth r=8. The boundary law is now retained through contact; no
assertion of an open shell gap is made for s>1. All differentiated initial
corner conditions hold because the data are flat in a shell neighborhood.
The retained full first reconstruction and dynamic-row equivalence apply
where the solution meets the shell. A common-harmonic identification still
retains its full initial one-jet and gauge/uniqueness premises.

This is a first-linear support theorem. Larger nonlinear constraint-correction
annuli, all second data, and their earlier contacts remain outside the claim.

## 5. Energy is monotone toward the present and almost constant

The exact temporal derivative is

```
d(x)=V_x=-(35/2)ell^3*sech(ell*x)^2*tanh(ell*x),
b_x=Lambda*d(x).
```

Thus d>=0 for x<=0, and (2) strengthens to

```
E'(x)=d(x)*||Z(x)||_H^2/2>=0,
E(x)<=E(0), -3<=x<=0.                                  (7)
```

No upper growth factor is needed for the backward energy bound. Positivity
also gives a comparison from below. Indeed

```
E(x)>=(5500ell^2/2)||Z(x)||_H^2,
```

because `b/Lambda>180ell/2500>5500ell^2`. Combining this with
`E(v)<=E(0)` for v in[x,0]
and integrating d exactly proves

```
0<=E(0)-E(x)
 <=[-omega(x)/(5500ell^2)]E(0)
 <=[7 tanh(ell*x)^2/4400]E(0),
(1-63ell^2/4400)E(0)<=E(x)<=E(0), -3<=x<=0.              (8)
```

The fractional energy change is at most 63/(4400*10^12), about 1.44e-14.
This is a mathematical estimating energy with time-dependent weights, not
a physical conserved Noether or injected energy.

At x=0 the boundary contribution is zero. If D=||Q_r||_2, endpoint Cauchy gives

```
D^2/2<=E(0)
 <=[1+(103001ell^2)/4]D^2/2.                             (9)
```

The supplied first velocity is retained. Write `q_p(s)=U(-s,0)` for the
past-time boundary trace, with s=-x. The boundary mass controls its velocity,
and its zero history for 0<=s<=1 gives

```
|q_p'(s)|^2<=E(0)/500,
|q_p(s)|^2<=[(s-1)_+]^2 E(0)/500<=E(0)/125.             (10)
```

This controls the master trace. Energy (1) alone does not give a pointwise
boundary acceleration or physical scalar-response amplitude. The next
graph estimate supplies the missing control for this smooth registered packet.

## 6. A quantitative upper bound for the actual first scalar trace

The unit packet has explicit first two derivative bounds. On 1<r<2 put
`y=1/[(r-1)(2-r)]>=4` and `a=3-2r`, so |a|<=1 and

```
Q=exp(4-y),
Q_r=a*y^2*Q,
Q_rr=[a^2*y^4-2a^2*y^3-2y^2]*Q.
```

For n=2,3,4, the function `y^n exp(4-y)` is nonincreasing on y>=4,
because its derivative has sign n-y. Endpoint flatness extends the bounds
over the whole line. The unit support length then gives

```
||Q||_2<=1, ||Q_r||_2<=16,
||Q_rr||_2<=256+128+32=416.                              (11)
```

Let `Y=A0^(1/2)Z`. It need not have compact spatial support. However the
power-domain regularity above makes it a smooth abstract solution of
`Y_xx+(A0+omega I)Y=0`; scalar omega commutes with A0. Its complete
mass-space/form energy obeys the same identity as (7), by abstract
selfadjointness, so `E_Y(x)<=E_Y(0)` backwards. This does not apply a
finite-propagation assertion to the nonlocal square root.

Set `kappa0=5500ell^2`, `epsilon=63ell^2/4400`. Static coercivity gives
`A0>=kappa0 I`, and on -3<=x<=0, `|omega|/kappa0<=epsilon`. Thus

```
E_Y(x)=1/2||Y_x||_H^2+1/2||A0 Z||_H^2
                         +omega(x)||A0^(1/2)Z||_H^2/2
       >=(1-epsilon)||A0 Z||_H^2/2.                    (12)
```

At x=0, omega=0. Write `f=-ell exp(-ell r)Q`, the backward initial
velocity. Because Z0 and f are flat at the shell,

```
A0Z0=(-Q_rr+V0Q,0),
2E_Y(0)=||A0 Z0||_H^2+a0(f,f)
 <=[416+103000ell^2]^2
      +ell^2[(16+ell)^2+103000ell^2].                    (13)
```

The form identity `||A0^(1/2)f||_H^2=a0(f,f)` retains every endpoint
contribution of the nonlocal vector. No boundary component of Y is set to
zero separately. The exact rational comparison of (12)--(13) proves

```
||A0 Z(x)||_H<417, -3<=x<=0.                            (14)
```

Also (11) in the basic initial energy gives `E(0)<129`; (10) then
gives `|q(x)|<2` on the same interval.

The physical scalar coefficient is not q. In the retained copy-even
longitudinal representative it is

```
S_phys(t)=3*V_z(t,0)/(rho_b^3 p_b),
V=alpha(z)cosh(t)^(-3/2)w,
alpha_b=1, alpha_z,b=-C_b/2.
```

The normal embedding term is included by the retained straight-shell/gauge
identification. Using the actual dynamic row and the endpoint component
of A0 gives the exact cancellation, with t=ell*x,

```
S_phys(ell*x)=3mu_b cosh(ell*x)^(-3/2)/(rho_b^3 p_b)
               *[-(A0Z)_b/ell^2+(5/2)q].               (15)
```

For example the endpoint row is `U_r=b0q-Lambda(A0Z)_b`, so
`U_r/ell-(C_b/2)q=mu_b[-(A0Z)_b/ell^2+(5/2)q]`.
The time-dependent omega cancels exactly. From the Hilbert mass norm,

```
mu_b |(A0Z)_b|<=sqrt(mu_b ell)||A0Z||_H,
sqrt(mu_b ell)<1/20000.
```

The retained shell box also gives rho_b>=129 and |p_b|>=3/25000.
Since cosh(t)^(-3/2)<=1 for real t, (14)--(15) imply the actual bound

```
|S_phys(t)|
 < 3/[129^3*(3/25000)]
       *[417/(20000ell^2)+(5/2)*(1/400)*2]
 < 250000000, -3ell<=t<=0.                              (16)
```

The physical first field is
`Phi1(t,Omega)=S_phys(t)Y_sphere(Omega)`, where
`Y_sphere=n1^2-n2^2` has absolute value at most one and angular mean zero.
Thus (16) is also a bound on this actual first field coefficient. It is
not a scalar mean, monopole volume coefficient, radiation density or
heating estimate. At the background proper shell clock `tau=rho_b*t`,
the shell scalar is constant, so the first clock pullback adds no
tangential shift term. No SI calibration is supplied. This is a bound
on the selected linear coefficient, not a finite-amplitude nonlinear
family or an effective common perturbation window. It is an upper bound
rather than a computed response profile and is distinct from the
auxiliary contradiction estimate used next.

## 7. Qualitative first physical-scalar activation

Suppose S_phys vanished on a past interval (-T1,0), ell<T1<=3ell.
The original V boundary equation, with its zero initial value/velocity,
then forces V(t,0)=0, hence w=w_z=0 there. Extend w by zero to z<0 and extend
the smooth potential boundedly as an auxiliary coefficient. There is no
distributional interface source. Cauchy uniqueness identifies this with
the whole-line solution with backward velocity -exp(-z)q_ell(z).

Use unscaled backward elapsed time S=-t and A(S,z)=w(-S,z), with initial
velocity `-exp(-z)q_ell(z)`. Choose
`0<eta<=min(ell/4,(T1-ell)/2)` and set
`S_obs=ell+eta<T1`, `Qfront=q_ell(S_obs)>0`.
On the backward diamond of (S_obs,0), initial samples lie in the increasing
edge [ell,S_obs], and the nonzero part of each source slice has width at
most eta: for source time v<=ell the positive-depth overlap is contained
in `[ell-v,S_obs-v]`, and for v>=ell it is in `[0,S_obs-v]`.
Under the zero-trace hypothesis the negative-depth part vanishes.
Thus, with M=103000 and beta=M*S_obs*eta/2,

```
Ustar<=Qfront*(1+eta/2)/(1-beta),
free observed value>=Qfront*(1-eta)/2,
observed value>=Qfront*[1-eta-3beta]/[2(1-beta)]>0.        (17)
```

The decisive scalar guard eta+3beta<1 holds for every eta<=ell/4; the worst
case is beta=5Mell^2/32. This contradicts the proposed zero trace. Therefore
S_phys is nonzero somewhere in every interval `(-T1,-ell)` with
`ell<T1<=3ell`. It is zero on `[-ell,0]` by (6) and continuity; smoothness
and the open zero interval give zero jets at t=-ell itself.

As in the retained activation proof, the positive lower value in (17)
belongs to the auxiliary whole-line solution under the contradictory
zero-trace hypothesis. It is not a numerical pointwise lower bound for
the actual dynamic-boundary response, nor a nonlinear first-arrival theorem.

## 8. Unresolved physical preparation

This gives a local first master development with its actual homogeneous
dynamic shell boundary through the first possible contact. It does not
construct second-order bulk/junction sources on that past domain, a full
nonlinear compatible family, a remotely prepared incoming state, a matter
sector or an observational signal. The principal-channel shares from the
precontact theorem must not be continued through reflection by energy
conservation alone. The earlier characteristic-precursor capsule and all
other frozen work remain unchanged.
