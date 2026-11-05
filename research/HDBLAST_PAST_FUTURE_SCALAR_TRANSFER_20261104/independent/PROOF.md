# Independent odd-part proof and retained future-certificate join

4 November 2026, America/Los_Angeles. Independent internal mathematical
review for a possible continuation; no earlier frozen file is changed.
All conclusions below retain the analytic branch, shell enclosure, full
linear reconstruction and uniqueness premises stated in the selected
contact and physical-shell references. The older actual future waveform
enclosure theorem is accepted as a prior conditional certificate, rather
than rerun or proved by hashing.

## 1. Isolate the odd initial-data part

Use the normalized coordinates x=t/ell, r=z/ell, ell=10^-6 and the common
positive endpoint-mass realization

```
H = L2(0,infinity) direct-sum R,
||(u,q)||_H^2 = integral u^2 dr + Lambda*q^2,
Lambda = mu_b/ell > 1000,
A0(u,q) = (-u_rr+V0*u, (-u_r(0)+b0*q)/Lambda),
V0 = ell^2*rho^2*Qbg,
b0 = ell*(gamma+(5/2)*mu_b), gamma=C_b/2,
Z_xx+(A0+omega(x)I)Z=0,
omega(x)=-(35/4)*ell^2*tanh(ell*x)^2.
```

The form domain links q=u(0); the operator domain has no additional static
Robin condition. Both endpoint kinetic and spring terms are included in
the form/wave construction. The smooth compatible packet data belong to
every power domain because they are flat near the shell. Scalar omega
commutes with A0. These facts supply smooth unique developments for both
signs of the initial velocity, not merely an energy-space boundary value.
The static auxiliary extension is common to both developments and agrees
with the actual coefficients on their support r<=5 for |x|<=3; no claim
about a physical background beyond the collar is required.

Let Zplus solve the forward equation with displacement Q and velocity
g=ell*exp(-ell*r)*Q; let Zminus have displacement Q and velocity -g.
The complete bulk equation and dynamical boundary law are even in x.
Consequently uniqueness gives

```
Zcurrent(+s)=Zplus(s), Zcurrent(-s)=Zminus(s), 0<=s<=3.
```

Now define O=(Zplus-Zminus)/2. This independently useful decomposition has
O(0)=0 and O_s(0)=(g,0). The desired past-minus-future difference is
Delta=Zminus-Zplus=-2O. Keeping the factor two here avoids introducing it
again in the physical observable.

For positive s, omega_s<=0. The complete energy of any such solution is

```
E[X] = (||X_s||_H^2 + ||A0^(1/2)X||_H^2
                         + omega*||X||_H^2)/2,
E'[X] = omega_s*||X||_H^2/2 <= 0.
```

The available bounds Q<=1, support length1 and exp(-ell*r)<=1 imply
E[O](s)<=E[O](0)<=ell^2/2. Finite propagation applies to O itself, so its
shell trace o(s) vanishes for 0<=s<=1. Since the endpoint kinetic energy
is at least 500*|o_s|^2, integration over 1<=s<=3 yields

```
|o(s)|^2 <= 4*E[O](0)/500 <= ell^2/250,
|Delta_q(s)|^2=4|o(s)|^2 <= 2ell^2/125 < ell^2.
```

This is a bound on the local master trace, not on its normal derivative.

## 2. Graph energy and physical derivative control

The higher vector A0^(1/2)O is generally nonlocal. Apply the abstract
selfadjoint form identity to it; do not apply a compact-support argument
to this vector. Its energy is nonincreasing because omega is scalar, and
static coercivity A0>=5500ell^2*I gives

```
2 E[A0^(1/2)O](s) >= (1-epsilon)*||A0 O(s)||_H^2,
epsilon = 63ell^2/4400 < 1.
```

Indeed |omega|<= (35/4)*ell^2*(3ell)^2. Since O(0)=0, its initial
displacement contribution vanishes exactly. The initial form is

```
2 E[A0^(1/2)O](0)=a0(g,g)
 <= ell^2*((16+ell)^2 + 103000ell^2).
```

Here g_r=ell*exp(-ell*r)*(Q_r-ell*Q), the initial endpoint trace is zero,
and ||Q_r||_2<=16. The latter follows analytically by writing
Q=exp(4-y), y=1/[(r-1)(2-r)]>=4:
|Q_r|<=y^2*exp(4-y)<=16. No bump sampling is used. The identity a0(g,g)
retains all endpoint contributions of the nonlocal square-root vector.
Therefore

```
||A0 Delta(s)||_H^2
 <= 4ell^2*((16+ell)^2+103000ell^2)/(1-epsilon)
 < (33ell)^2.
```

The embedding-included first physical scalar coefficient is a linear
functional of the fixed first master. The retained boundary cancellation
has the exact form

```
S_X(ell*s)=3mu_b*cosh(ell*s)^(-3/2)/(rho_b^3*p_b)
                 *[-(A0X)_b/ell^2+(5/2)*X_b].
```

For clarity, U_r=b0*q-Lambda*(A0Z)_b implies
U_r/ell-gamma*q=mu_b[-(A0Z)_b/ell^2+(5/2)q]. The temporal omega cancels.
The spatial normal is not reversed by x->-x, and cosh is even. Thus this
functional gives S_Delta(t)=S_current(-t)-S_current(t) with no extra sign.
The constant background shell scalar also removes first-order clock/label
pullback terms; the paired proper times are -rho_b*t and +rho_b*t at the
same anchor and quadrupole normalization.

The H norm includes Lambda times the squared endpoint component, hence

```
mu_b*|(A0Delta)_b| <= sqrt(mu_b*ell)*||A0Delta||_H,
sqrt(mu_b*ell)<1/20000.
```

Using rho_b>=129, |p_b|>=3/25000, mu_b<1/400 and |Delta_q|<ell gives

```
|S_current(-t)-S_current(t)|
 < (25000/129^3)*(33/(20000ell)+ell/160)
 = 264000000001/13738809600 < 20,
                            0<=t<=3ell.
```

This is an absolute bound on the first linear coefficient in model units.
It gives neither an individual scalar amplitude below20 nor a relative
error near zeros, physical energy law or finite nonlinear amplitude window.

## 3. The older certificate concerns exactly this observable

The selected October15 INERTIAL_SHELL_WAVEFORM primary proof specifies
ell=10^-6, the same delta and model parameter, Q, initial velocity exp(-z)Q,
quadrupole n1^2-n2^2, conjugation alpha*cosh(t)^(-3/2) and complete boundary
law beta=gamma+mu_b*(35sech(t)^2/4-25/4). Its actual report pins the same
raw shell and wide-background JSON bytes as the contact construction.
Different safe potential ceilings 102000 and103000 do not change the field.
Under the retained full reconstruction and linear uniqueness premises,
this identifies the future certificate with S_current(t), rather than
with a zero-velocity state or a free comparison waveform.

The observation is explicitly embedding-included at tau=rho_b*t:

```
H(x)=ell*(w_z(ell*x,0)-gamma*w(ell*x,0)),
ell*S(ell*x)=[3/(rho_b^3*p_b)]*cosh(ell*x)^(-3/2)*H(x).
```

The report's ell_S intervals therefore must be divided by ell before
applying the difference bound. Its H intervals include the proved
actual-minus-comparator error. Comparator H0 alone is insufficient. The
independent arithmetic additionally checks that the displayed ell_S
intervals contain the full gain/cosh/actual-H interval products; this
does not replace the older continuous enclosure theorem.

If [L,U] encloses a retained future value of S, then [L-B,U+B] encloses
the matched past value, with B=264000000001/13738809600. Exact endpoints
are in the independent arithmetic report. In particular,

```
S_current(-5ell/4) < -37376,
S_current(-7ell/4) > 37379.
```

The smooth actual scalar coefficient has at least one zero in
(-7ell/4,-5ell/4) by the intermediate value theorem. This says neither
uniqueness nor the total number of zeros. Because the angular quadrupole
changes sign, these are signs of S, not a globally signed angular field.

The four older central/late narrow-bracket witnesses and the continuous
negative future tail have margins below B. Their transported intervals
straddle zero. This method consequently does not transfer those narrow
zero brackets or tail sign. A quadratic response integral also needs a
separate norm/cross-term argument; interval endpoints cannot simply be
widened by B as if it were pointwise.

## 4. Status and scope

Time reflection, positive graph-energy estimates and interval transport
are established principles. This is a conditional application to the
registered packet and its physical scalar trace, not a new general law
or an established novelty/significance claim. Internal independent review
is not external peer review or proof-assistant certification.

The earlier future-contact project already constructs a conditional full
finite-angular second Taylor coefficient under its stated correction,
gauge, corner and support hypotheses. The separately pinned scope
addendum preserves that formal existence result. The quantitative chosen
past coefficient/observable, a nonlinear defect-to-error map, effective
parameter interval, physical incoming mechanism, radiation and
cosmological-origin evidence remain separate open tasks.

The checker is new standard-library rational arithmetic over seven
selected inputs. It imports no project program and executes no legacy
producer, numerical PDE solve or background/source/bump oracle. Analytic
arguments and retained certificate premises remain explicit.
