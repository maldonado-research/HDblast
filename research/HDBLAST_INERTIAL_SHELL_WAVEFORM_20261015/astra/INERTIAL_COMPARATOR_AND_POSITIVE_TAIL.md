# Inertial comparison resolves a signed first-variation tail

15 October 2026 UTC. Independent analytic derivation and exact-rational
checks. This note concerns the conditional coupled **linear first variation**
at the newly specified width ell=10^-6. It does not establish the localized
nonlinear inverse, gauge bridge, support control, reflection hypotheses,
finite-regularity evolution theorem, or an effective nonlinear remainder.

The retained near-phi=+1 branch, exact parameters and analytic background
theorem are inputs. Its outer rectangle continues to Z=8ell by the strict
first-exit argument independently replayed here. No narrow-width response
coefficient is transferred. The calculation proves a useful fact beyond an
absolute tail bound: the actual normalized first-variation trace is positive
throughout 2<=t/ell<=3. Since the retained physical gain is negative, its
physical quadrupole coefficient is negative there.

## 1. State, comparator and exact boundary identity

Set ell=1/1000000, T=3ell, Z=8ell. Write q(z)=Q(z/ell) and
h(z)=exp(-z)q(z), with the retained smooth nonnegative bump Q supported on
[1,2] and bounded above by one. The actual smooth compatible solution obeys

    w_tt-w_zz+W(t,z)w=0,
    -w_z+mu w_tt+beta(t)w=0 at z=0,
    w(0,z)=q(z), w_t(0,z)=h(z).

Here gamma=C_b/2, mu=(rho_b B_b)^(-1), and

    beta0=gamma+5mu/2,
    delta_beta(t)=beta(t)-beta0=-(35mu/4)tanh(t)^2.

Let y(t)=w(t,0), j0(t)=q'(t)+h(t), and

    I(t)=integral_0^t W(s,t-s)w(s,t-s) ds.

Integrating (partial_t-partial_z)(w_t+w_z)=-Ww on its incoming
characteristic gives the exact relations

    w_t(t,0)+w_z(t,0)=j0(t)-I(t),
    mu y''+y'+beta0 y=j0-I-delta_beta*y,
    y(0)=y'(0)=0.

Define the free inertial comparator u by u_tt-u_zz=0 with the same data and
the constant dynamic condition -u_z+mu u_tt+beta0*u=0. Its boundary value
y0 solves mu y0''+y0'+beta0*y0=j0, with zero initial boundary data. If G
is the scalar impulse response,

    mu G''+G'+beta0 G=0, G(0)=0, G'(0)=1/mu,

then y0=G*j0. This definition applies whether the characteristic roots are
real or complex; no root-discriminant division is needed in this proof.
The normalized outputs are

    H(t/ell)=ell*(j0-I-y'-gamma*y),
    H0(t/ell)=ell*(j0-y0'-gamma*y0).

For f=I+delta_beta*y, the exact error identity is therefore

    H-H0=ell*(-I+(G'+gamma*G)*f).                 (1)

Here convolution always integrates from 0 to t. The entire preceding
history of I appears in (1), including the main pulse; replacing it by a
post-pulse bound inside that full convolution is not justified.

## 2. Independently checked coefficients and majorants

Throughout the retained rectangle and the wider collar the checker proves

    m=1/500 < mu < 1/400,
    200 < gamma < Gamma=250,
    0 < beta0 < B=250, |beta(t)|<B,
    |W(t,z)|<M=102000, |W_z(t,z)|<Kz=30000000.

The independent interval for W_z is approximately
[-21392210.894804,28190745.284942]; these decimals are display only.
For completeness, let U be the background scalar potential and h_bg=U_eta/p.
The retained derivative-free bulk coefficient is

    W(t,z)=V(z)-(35/4)tanh(t)^2,
    V=rho^2 F,
    F=2h_bg^2-10a h_bg+35a^2/4-U_etaeta-p^2/6.

The ODE gives

    rho_z=-rho^2 a, eta_z=-rho p,
    p_z=-rho(U_eta-4ap), a_z=1/rho+rho p^2/3,
    (h_bg)_z=rho*(h_bg^2-4a h_bg-U_etaeta),
    F_z=(4h_bg-10a)(h_bg)_z+(-10h_bg+35a/2)a_z
        -U_etaetaeta*eta_z-p*p_z/3,
    W_z=2rho*rho_z*F+rho^2*F_z.

The checker differentiates the independently expanded exact polynomial

    U=-2/27+(14/9)eta^2+(50/27)eta^3-eta^4/6
      -(4/9)eta^5-(2/27)eta^6

and evaluates those formulas by rational interval arithmetic. No finite
difference or sampled background trajectory supplies the derivative bound.

Use the retained positive matrix majorant, independently checked here:

    U0=1+T, r=M*T^2/2, s=B*T^2/(2m),
    d=(1-r)(1-s)-T*r/m>0,
    A=U0*(1-s+T/m)/d, Vb=U0/(m*d).

It supplies sup|w|<=A on the causal triangle and sup|y'|<=Vb. In particular
|y(t)|<=t*Vb. Finite propagation and the data gap imply y(t)=0 for t<=ell,
so also |y(t)|<=(t-ell)_+*Vb. The strict support and wider-background
first-exit hypotheses are the same as in the preceding linear theorem.

Since beta0>0, the kernel energy satisfies

    d/dt [mu*(G')^2+beta0*G^2]=-2*(G')^2.

Its initial value is 1/mu. Consequently |G'|<=1/mu and |G(t)|<=t/mu,
and hence |G'+gamma G|<=(1+Gamma*T)/mu for 0<=t<=T.

## 3. Uniform comparison error

The characteristic integral has the global bound |I(t)|<=M*t*A.
Also |delta_beta(t)|<=(35mu/4)t^2. Insert these and |y(t)|<=t*Vb
into (1). The exact mu factor cancels in the beta variation contribution.
For all 0<=t<=T,

    |H-H0| <= E_global,
    E_global=ell*[M*T*A+(1+Gamma*T)*M*A*T^2/(2m)
                         +(35/16)*(1+Gamma*T)*Vb*T^4]
            < 3.067e-7.                              (2)

The exact rational calculation is approximately 3.066900779405216e-7.
Before the first front, H=H0=0 for 0<=t<=ell. Equation (2) measures the
actual-minus-inertial-comparator error, rather than the much larger
actual-minus-Dirichlet error of the preceding feasibility calculation.

## 4. Odd-reflection cancellation and the sharper tail error

Let uD be the free Dirichlet solution with the same initial data. The
halfline Duhamel representation is

    w(s,z)=uD(s,z)+1_(s>=z)*y(s-z)+D(-Ww)(s,z),

where the Dirichlet wave kernel has nonnegative mass at most the elapsed
time. Thus |D(-Ww)(s,z)|<=M*A*s^2/2. Decompose uD into its q and h parts.
On the incoming ray, for t>=2ell,

    uD_q(s,t-s)=q_odd(t-2s)/2,
    integral_0^t uD_q(s,t-s) ds=0.

More precisely, pairing positive and negative arguments gives

    integral_0^t W(s,t-s)*uD_q(s,t-s) ds
      =1/4 integral_ell^(2ell) [W((t-r)/2,(t+r)/2)
                            -W((t+r)/2,(t-r)/2)]*q(r) dr.

The paired spatial coordinates differ by r, as do the paired times.
Because |d(tanh^2 s)/ds|<=2T for s in [0,T], the bracket has absolute
value at most (Kz+35T/2)*r. The elementary bound
integral r*q(r)dr <=2ell^2 therefore yields

    |weighted q contribution| <= (Kz+35T/2)*ell^2/2.       (3)

The h component is nonnegative on this ray: its defining odd-extension
integral runs from t-2s to t, with upper endpoint beyond the positive
support. Direct interchange of integration gives

    integral_0^t uD_h(s,t-s) ds
       =1/2 integral_ell^(2ell) r*exp(-r)*q(r) dr <=ell^2.

Thus its weighted contribution has absolute value at most M*ell^2.
The delayed lift and Duhamel remainder give respectively

    integral_0^t |y(2s-t)|*1_(2s>=t) ds
        <=Vb*(t-ell)_+^2/4,
    integral_0^t M*|D(-Ww)(s,t-s)| ds <=M^2*A*t^3/6.

Combining these observations for 2ell<=t<=T gives |I(t)|<=I_tail, where

    I_tail=(Kz+35T/2)*ell^2/2+M*ell^2
           +M*Vb*(T-ell)^2/4+M^2*A*T^3/6.

Use this only for the instantaneous -I term in (1), retaining the global
history estimate in the convolution. The resulting uniform tail error is

    |H-H0| <= E_tail,
    E_tail=ell*[I_tail+(1+Gamma*T)*M*A*T^2/(2m)
                       +(35/16)*(1+Gamma*T)*Vb*T^4]
          < 2.97e-10,       2<=t/ell<=3.               (4)

Its exact rational value is approximately 2.961665226886478e-10. The
cancellation is essential for this proof: substituting the global error
(2) into the tail-sign argument below fails, as the checker verifies.
No Born iteration or pointwise trace inferred from an L2 estimate is used.

## 5. A separate coarse proof of a strictly positive actual tail

Integration by parts in y0=G*(q'+h), using q(0)=0, gives

    y0(t)=integral_0^t [G'(t-s)+exp(-s)G(t-s)]q(s) ds,
    y0'(t)=q(t)/mu
            +integral_0^t [G''(t-s)+exp(-s)G'(t-s)]q(s) ds.

For t>=2ell, j0(t)=q(t)=0. Therefore

    H0(t/ell)=ell integral_ell^(2ell) L(t-s,s)q(s) ds,
    L(r,s)=-G''(r)-gamma G'(r)
                    -exp(-s)[G'(r)+gamma G(r)].

At r=0,

    L(0,s)=1/mu^2-(gamma+exp(-s))/mu
          >=1/mu^2-(Gamma+1)/mu >=59600.

The last expression decreases with mu on [1/500,1/400], since
(Gamma+1)mu<2. Its minimum is the endpoint value at mu=1/400.
This deliberately loose sign proof needs no correlation between mu and
gamma; correlated intervals can make a much tighter displayed waveform.

The kernel equation and energy estimates give

    C2=(1+B*T)/m^2 >=sup|G''|,
    C3=(C2+B/m)/m >=sup|G'''|,
    CL=C3+(Gamma+1)*C2+Gamma/m >=sup|partial_r L|.

Hence L(r,s)>=59600-T*CL>0 for the complete convolution domain. On
[4/3,5/3], Q>=exp(-1/2)>1/2, so integral_1^2 Q(x)dx>1/6. Consequently

    H0(x)>ell^2*(59600-T*CL)/6
         =9.808075427083333...e-9,    2<=x<=3.

With (4), the actual solution satisfies the strict uniform conclusion

    H(x)>9.5e-9,                    2<=x<=3.            (5)

This establishes a nonzero signed post-pulse tail for the conditional
linear problem. It does not establish a nonlinear tail at a specified
finite perturbation amplitude. Its physical sign follows from the retained
negative gain, and refers to the quadrupole coefficient; it is not a
statement that the angular quadrupole harmonic has one sign everywhere.

## 6. Reproduction and limits of the certificate

Run the independent checker against either the frozen previous capsule or
a relocation containing its two pinned input files:

```sh
python -I -B check_inertial_error.py --source-root /path/to/prior/capsule --output INERTIAL_ERROR_NORMAL.json
python -I -O -B check_inertial_error.py --source-root /path/to/prior/capsule --output INERTIAL_ERROR_OPTIMIZED.json
```

All numerical decisions use exact fractions. The checker imports no
producer, uses no Python assertions, authenticates both mathematical input
files, reconstructs the wider first-exit margins and coefficient intervals,
checks both majorant rows, and rejects five deliberate invalid cases.
The derivative formula and the analytic energy, propagation, characteristic
and cancellation arguments are written mathematical proofs, not machine
verified theorems. Independent internal derivation and finite arithmetic
are distinct from external peer review or formal proof verification.
