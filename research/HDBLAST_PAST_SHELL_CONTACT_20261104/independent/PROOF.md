# Independent past-contact energy and physical scalar bounds

4 November 2026, America/Los_Angeles. This companion derivation concerns the
registered first linear master on `-3<=x<=0`, including shell contact. It supplies
an independent coefficient route and a bound for the actual first scalar
coefficient. The smooth background, master reconstruction and gauge identification
remain retained premises. This is neither physical Noether energy nor a nonlinear
past preparation theorem.

## 1. Raw coefficients and the boundary normal

The [retained master normalization](inputs/PHYSICAL_SHELL_RESPONSE.md) gives

```
w_tt-w_zz+W(t,z)w=0,
-w_z+mu_b*w_tt+beta(t)w=0 at z=0,
mu_b=1/(rho_b B_b), B_b=h_b-4a_b+2(1+eta_b),
C_b/2=rho_b(h_b-5a_b/2), h=U_potential'(eta)/p,
beta=C_b/2+mu_b[(35/4)sech(t)^2-25/4].
```

Here z increases into the collar, `dy=-rho dz`, so the minus normal derivative
must retain its displayed sign. For `x=t/ell,r=z/ell`, `ell=10^-6`, the row is

```
-U_r+Lambda*q_xx+b(x)q=0,
Lambda=mu_b/ell, b=ell*beta(ell*x), q=U(x,0).
```

The checker evaluates the exact model polynomial derivative at the **raw eta
endpoints**, instead of using the stored derived Uprime/h/B/mu/C intervals.
On the outer eta interval, `U3>0` proves U2 increasing, and `U2>31/10` proves
U1 increasing. Endpoint division by the negative raw p interval therefore
encloses h. The identity `sigma''=4(1+eta)` is checked against the raw tension
interval. Interval multiplication then reconstructs B, C/2 and `1/(rho B)`.
The outward safe consequences are

```
1/1000<mu_b<1/400, 189<C_b/2<300-1/160,
180<beta(t)<300 for all real t,
1000<Lambda<2500, 180ell<b<300ell.
```

For the bulk, four fresh first-exit images use the raw shell box and the vector
field on the selected `8ell` outer rectangle; stored flow tubes are not premises.
The joint h/a monotonicity of
`Qbg=2h^2-10ah+35a^2/4-U2-p^2/6` gives `1<Qbg<16/3`.
Since `W=rho^2 Qbg-(35/4)tanh(t)^2`, this independently implies the shared
all-real-time safe bounds `5500<W<103000` on the actual collar. These are
explicit time identities, not transferred future evolution enclosures.

## 2. Positive form, construction and support through contact

Put

```
omega(x)=-(35/4)ell^2*tanh(ell*x)^2,
V=V0(r)+omega(x), V0=ell^2*rho^2 Qbg,
b=b0+Lambda*omega, b0=ell*(C_b/2+(5/2)mu_b).
```

Extend only the auxiliary static V0, keeping it exact for r<=6 and smoothly
interpolating to `6000ell^2` by r=7. The interpolation preserves the shared
positive bounds; no physical background outside the actual collar is asserted.
Use

```
H=L2(0,infinity) direct-sum R,
||(u,q)||_H^2=integral u^2+Lambda*q^2,
D(a0)={(u,u(0)):u in H1},
a0(u,u)=integral(u_r^2+V0*u^2)+b0*u(0)^2.
```

This is the standard positive endpoint-mass form, already used in the selected
[earlier domain construction](inputs/SCALAR_DOMAIN.md). Trace linkage holds on
the form domain, while general H components are independent. Thin boundary
layers prove density in H; the H1 trace theorem proves closedness. Its operator is

```
D(A0)={(u,u(0)):u in H2},
A0(u,q)=(-u_rr+V0*u,[-u_r(0)+b0*q]/Lambda).
```

Integration by parts fixes the minus sign in the endpoint component. There is
no additional static Robin condition on D(A0). The new local expression is the
earlier Liouville expression shifted by 5/2 before scaling; this shifts its
spring by `(5/2)mu_b` too. The auxiliary half-line extension is not the earlier
global regular-cone spectral operator.

The equation is `Z_xx+(A0+omega I)Z=0`, with `Z=(U,q)`. In backward time s=-x,
the data are `(Q,0)` and `(-ell exp(-ell*r)Q,0)`. Flat compact support in [1,2]
puts both in every A0 power domain, including all endpoint compatibility
conditions. Scalar omega commutes with A0. The wave-kernel Volterra operator
on 0<=s<=3 has norm at most `(2835/8)ell^4<1`, so its Neumann series constructs
the solution in every power-domain graph norm. The equation and one-dimensional
regularity give smooth classical bulk and boundary fields.

Moving-cut energy identities give outer flux
`-(U_s+U_r)^2/2-VU^2/2` at r=2+s. Before s=1, the inner wedge has flux
`-(U_s-U_r)^2/2-VU^2/2` at r=1-s and includes the boundary mass and spring.
Coefficient derivatives in backward time are nonpositive. Starting from zero
truncated energy, these identities prove support inside [1-s,2+s] before
contact and [0,2+s] afterward. Cutoffs and graph regularity justify passage to
infinity. Thus the full support stays within r<=5<8; the auxiliary extension
cannot affect this physical slab. No inner gap is asserted after s=1.

## 3. Direct energy comparison and the master trace

For

```
E=1/2 integral(U_x^2+U_r^2+V U^2)
    +Lambda*q_x^2/2+b*q^2/2,
```

the bulk lower-end flux is `-q_x U_r(0)`. Both boundary-energy derivatives
cancel it using `U_r=Lambda*q_xx+bq`. Hence, with d=omega_x>=0 for x<=0,

```
E_x=d||Z||_H^2/2>=0, E(x)<=E(0).
```

The bulk and boundary coefficients give
`E>=5500ell^2||Z||_H^2/2`. Integrating d exactly, then using
`|tanh(ell*x)|<=3ell`, gives

```
(1-epsilon)E(0)<=E(x)<=E(0),
epsilon=63ell^2/4400.
```

No exponential growth estimate is needed. At zero, endpoint Cauchy yields
`E(0)<=D^2[1+103001ell^2/4]/2`, where `D=||Q_r||_2`.
The boundary mass gives `|q_s|^2<=E(0)/500`; integrating from the exact zero
history s<=1 gives `|q|^2<=E(0)/125` on 0<=s<=3.

## 4. Higher-graph control

The square-root state `G=A0^(1/2)Z` is generally nonlocal. No support assertion
for G is used. Its abstract energy identity follows on the propagated power
domains from selfadjointness and the scalar commuting omega:

```
Ehalf=||A0^(1/2)Z_x||_H^2/2+||A0Z||_H^2/2
           +omega||A0^(1/2)Z||_H^2/2,
Ehalf_x=omega_x||A0^(1/2)Z||_H^2/2>=0.
```

The spectral lower bound `A0>=5500ell^2 I` implies
`Ehalf>=(1-epsilon)||A0Z||_H^2/2`. At zero, with
`f=±ell exp(-ell*r)Q`,

```
2Ehalf(0)=||-Q_rr+V0Q||_2^2+a0(f,f).
```

The boundary coordinate of A0Z0 and the trace of f vanish because the data are
flat near the shell. This does not set the boundary coordinate of G to zero:
the form identity retains all its contributions.

For `y=1/[(r-1)(2-r)]>=4`, `a=3-2r`, the exact bump derivatives are

```
Q_r=a*y^2*Q,
Q_rr=(a^2*y^4-2a^2*y^3-2y^2)Q.
```

The functions `y^n exp(4-y)` decrease for y>=4 and n<=4, so unit support gives
`||Q||_2<=1`, `||Q_r||_2<=16`, `||Q_rr||_2<=416`. Consequently

```
2Ehalf(0)<=(416+103000ell^2)^2
              +ell^2[(16+ell)^2+103000ell^2],
||A0Z||_H<417.
```

These same envelopes give `E(0)<129` and therefore `|q|<2`.

## 5. An upper bound for the actual first physical scalar

The retained embedding/gauge bridge identifies
`Phi1(t,Omega)=S(t)Y_sphere(Omega)`, where `Y_sphere=n1^2-n2^2`,
`|Y_sphere|<=1`, and

```
S=3 Vphys_z/(rho_b^3 p_b),
Vphys=alpha(z)cosh(t)^(-3/2)w,
alpha_b=1, alpha_z,b=-C_b/2.
```

The dynamic row gives `w_z=mu_b*q_xx/ell^2+beta*q`. Also
`q_xx=-(A0Z)_b-omega*q` and
`beta-C_b/2=mu_b(5/2+omega/ell^2)`. The omega terms cancel exactly:

```
S(ell*x)=3mu_b cosh(ell*x)^(-3/2)/(rho_b^3 p_b)
                  *[-(A0Z)_b/ell^2+(5/2)q].
```

The boundary coordinate has weight Lambda, so
`mu_b|(A0Z)_b|<=sqrt(mu_b ell)||A0Z||_H`, with
`sqrt(mu_b ell)<1/20000`. Together with rho_b>=129, |p_b|>=3/25000,
cosh(t)^(-3/2)<=1, the graph and master-trace bounds give

```
sup_{-3ell<=t<=0,Omega}|Phi1(t,Omega)|
 < 25000/129^3 *[20850000000+1/80]
 = 1042500000000625/4293378
 < 250000000.
```

This is a deliberately coarse upper bound for the **actual selected first
coefficient**, in fixed model units. It is not the auxiliary activation value,
an observed signal, a scalar mean, radiation density or a finite nonlinear
amplitude. The background shell scalar is constant, so the first proper-clock
pullback introduces no tangential shift term. No common nonlinear epsilon
window, SI calibration or physical Noether-energy interpretation follows.

## 6. Activation and computational scope

For the separate qualitative statement, suppose the physical scalar trace
vanishes on (-T1,0), ell<T1<=3ell. Its original boundary ODE then forces
Vphys=w=w_z=0 at the shell. Zero extension across the shell creates no interface
source. Whole-line Cauchy uniqueness applies under this contradictory hypothesis.
Use backward elapsed time `s_obs=-t=ell+eta`, with
`0<eta<=min(ell/4,(T1-ell)/2)`. The increasing front of Q and the source-strip
width eta give the positive auxiliary factor

```
(1-eta-3beta)/(2(1-beta)), beta=103000*s_obs*eta/2,
eta+3beta<=ell/4+3*(5*103000ell^2/32)<1.
```

This contradicts zero trace and proves some activation in every (-T1,-ell).
It gives no actual pointwise lower response or first-arrival amplitude. The
front itself has zero jets by smooth flatness. Neither the upper-bound nor
activation argument extends precontact principal-weight shares through reflection.

`verify.py --output ABS_FRESH_OUTPUT` recomputes the interval and exact scalar
gates from five authenticated local snapshots; use `-I -S -B`, optionally `-O`.
The normal/O reports coincide. Finite gates do not replace the continuous
closed-form, reconstruction, smoothness or uniqueness arguments above. No old
checker, background/bump/source oracle or numerical PDE evolution is executed.
No second constrained preparation, nonlinear past family, remote incoming law,
matter production or cosmological-origin result is supplied.
