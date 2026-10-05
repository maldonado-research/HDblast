# Independent action, source-work and finite-band Ward derivation

Prepared 5 October 2026 UTC. This route derives the identities from the
minimal-scalar action and independently transcribes the fixed-reference
subtraction inventory. It imports no other route's numerical or algebra helper.
The derivation and first exact receipt preceded reading the new primary
derivation. No real pulse callback, saved physical mode, archive array,
momentum quadrature or remote operation was used.

The result is an exact algebraic conservation/error framework in the specified
linear model, and a continuous-momentum bound on its **analytic source
truncation component**. It is not a completed state, pressure/contact, cutoff
or cosmological certificate. The historical metric failures remain unchanged.

## 1. Direct action and the physical mass source

Use signature (-+++), minimal coupling, prescribed real mass squared x and
the matter action

\[
\Gamma=-\frac12\int\sqrt{-g}
  [g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi+x\phi^2].
\]

Let v be the canonical mode, L=a'/a, D=v'-Lv and
`kappa=k^2+a^2 x`. The equation of motion gives

\[
v'=D+Lv,\qquad D'=-LD-\kappa v.
\]

Three real bilinears suffice: U=|v|^2, A=|D|^2, C=Re(Dv*).
Their independently derived derivatives are

\[
U'=2C+2LU,\quad A'=-2LA-2\kappa C,\quad C'=A-\kappa U.
\]

Direct metric and mass variation defines

\[
\rho_k=\frac{A+\kappa U}{2a^4},\quad
p_k=\frac{A-(k^2/3+a^2x)U}{2a^4},\quad Q_k=\frac{U}{a^2}.
\]

Substitution, without defining pressure through conservation, proves

\[
\rho_k'+3L(\rho_k+p_k)=\frac{x'Q_k}{2}.             \tag{I1}
\]

Equivalently `(a^3 rho_k)'+p_k(a^3)'=a^3 x' Q_k/2`.
If an approximate mode has residual
`F=v''+[k^2+a^2 x-a''/a]v`, the additional exact defect is
`Re(D*F)/a^4`. This term vanishes on shell; a small integrated cancellation
does not independently bound its absolute magnitude.

The fixed positive-reference W0/W2/W4 inventory is also checked independently
for arbitrary symbolic FRW and mass jets. In numerator normalization
`S_R=a^4 rho_sub`, `S_P=a^4 p_sub`, it obeys

\[
S_{R,0}'-L(S_{R,0}-3S_{P,0})=0,
\]
\[
S_{R,2}'-L(S_{R,2}-3S_{P,2})=a^4x'Q_0/2,
\]
\[
S_{R,4}'-L(S_{R,4}-3S_{P,4})=a^4x'Q_2/2.
\]

The verifier treats `a^2`, `L,L',...`, and `D_mass=a^2(x-r)` and its
derivatives as arbitrary jets, with constant positive r and
`w^2=k^2+a^2r`. It checks the published density and **direct pressure**
expressions separately. Thus the finite-band subtracted identity follows
before taking any ultraviolet limit. An added common local action must shift
stress and Q together, as specified in the inherited common-action convention.
Ward consistency does not choose the finite renormalization coefficients.

## 2. The registered metric pulse has fixed physical mass

At `a0=L=-1/eta`, `L'=L^2`, `x=r=2`, write
`a=a0(1+epsilon h)` and

\[
g=4L^2h-2Lh'-h''.
\]

The canonical background mode is `v0=exp(-ik eta)/sqrt(2k)` for k>0.
Raw perturbations include epsilon:
`delta v=v0 u`, `delta v'=v0(w-ik u)`.
Put `X=Re u`, `Z=Re w`, `T=Im w`. Direct stress variation yields

\[
R_m=\frac{(2k^2+3L^2)X-kT-LZ}{2k\epsilon},\quad
P_m=\frac{(2k^2/3-L^2)X-kT-LZ}{2k\epsilon}.
\]

For arbitrary real forcing f, `u'=w`, `w'=2ik w-epsilon f` gives

\[
c_R=X-T/(2k),\quad c_R'=0,
\]
\[
R_m'=L(R_m-3P_m)+\frac{Lf}{2k}.                     \tag{I2}
\]

The bare metric contacts are obtained from the D operator, mass term and
physical-volume factor separately:

\[
C_{R,b}=Lh'/(2k)-2kh-2L^2h/k,\quad
C_{P,b}=Lh'/(2k)-2kh/3,
\]

with `R0,b=(2k^2+3L^2)/(4k)` and
`P0,b=(2k^2/3-L^2)/(4k)`. They satisfy

\[
C_{R,b}'=L(C_{R,b}-3C_{P,b})-3h'(R_{0,b}+P_{0,b})
         -\frac{Lg}{2k}.
\]

Linearizing the independently checked subtraction Ward identities gives
the same relation for the full renormalized contacts and matching baseline.
Consequently, for `R=R_m+C_R`, `P=P_m+C_P`, `B0=R0+P0`,

\[
R'=F+\frac{L(f-g)}{2k},\qquad
F=L(R-3P)-3h'B0.                                   \tag{I3}
\]

With f=g the fixed-mass physical Ward source is zero. The two apparent
`Lg/(2k)` terms in the raw-mode and contact equations cancel. The certified
primitive `integral Lg` is one component of this bookkeeping, not the full
action pressure/contact work and not a physical heating measurement.

## 3. Canonical source work, retarded state and instantaneous contacts

After the canonical field redefinition and its declared boundary convention,
`omega^2=k^2+a^2x-a''/a`. The canonical Hamiltonian obeys

\[
H_{can}=(|v'|^2+\omega^2|v|^2)/2,\quad
H_{can}'=(\omega^2)'|v|^2/2.
\]

At first order its normalized variation is

\[
\delta H_{can}/\epsilon=k c_R/\epsilon+g/(4k),
\qquad (\delta H_{can}/\epsilon)'=g'/(4k).
\]

The explicit `g/(4k)` is an instantaneous Hamiltonian source potential.
It is distinct from the physical fixed-mass stress and its pressure work.
On a complete compact pulse with g vanishing at both ends it gives zero
net first-order canonical energy change. Positive particle occupation energy
would begin at higher perturbative order; this identity cannot certify it.

The state contribution has retarded support. With zero incoming perturbation,

\[
w(t)=-\epsilon\int_a^t e^{2ik(t-s)}f(s)ds,\quad
u(t)=-\epsilon\int_a^t\Phi_k(t-s)f(s)ds,
\]

where `Phi_k(d)=(exp(2ikd)-1)/(2ik)`, `Phi_0(d)=d`.
For nonzero incoming perturbations a separate homogeneous solution is added.
Neither the retarded source integrals nor their error enclosures bound an
unverified incoming quantum state.

The metric/subtraction contacts are instantaneous functions of h, its
derivatives, geometry and the fixed reference; they have no state history.
For the raw variance `q_m=X/(k epsilon)`, direct pressure reduction proves

\[
R_m=kc_R/\epsilon+(3L^2q_m-Lq_m')/2,
\]
\[
P_m=kc_R/(3\epsilon)
 +(q_m''-3Lq_m'-3L^2q_m)/6+f/(6k).                 \tag{I4}
\]

The positive instantaneous `f/(6k)` is required when converting the raw
pressure into variance derivatives. It is not an independently selectable
energy-transfer law. At a matched continuous cutoff its integrated value
is `M_K f/3`, so the corresponding analytic pressure contact contains
`-M_K g/3`. Omitting either term changes the direct pressure target and
spoils the contact ledger.

## 4. Full approximate-state defect, including canonical drift

For any differentiable raw state define

\[
e_u=u'-w,\qquad e_w=w'-2ik w+\epsilon f.
\]

The direct mode density has the additional residual

\[
E_R=\frac{(2k^2+3L^2)\operatorname{Re}e_u
          -k\operatorname{Im}e_w-L\operatorname{Re}e_w}{2k\epsilon}.
                                                               \tag{I5}
\]

Define an arbitrary independently assessed contact/subtraction defect

\[
\chi_C=C_R'-L(C_R-3C_P)+3h'B0+Lg/(2k).
\]

Then the exact off-shell full ledger is

\[
R'-F=E_R+L(f-g)/(2k)+\chi_C.                       \tag{I6}
\]

The commonly used primitive `G=(3L^2X-LZ)/(2k epsilon)` has only
`E_G=[3L^2 Re(e_u)-L Re(e_w)]/(2k epsilon)`.
The missing term is

\[
E_R-E_G=(kc_R/\epsilon)'
       =k\operatorname{Re}e_u/\epsilon
          -\operatorname{Im}e_w/(2\epsilon).
\]

A G-only error budget can therefore miss a density drift even when its own
residual is small. A Wronskian projection cannot remove this term as an
unrecorded numerical correction.

For normalized amplitudes `U=u/epsilon`, `W=w/epsilon`, the residuals are
`r_u=e_u/epsilon`, `r_w=e_w/epsilon` and formula (I5) has no epsilon
denominator. The forcing/contact terms in (I6) retain the displayed scaling.
The exact verifier checks both the raw equation and this normalization.

## 5. Finite-band integration and mixed targets

Use the fixed measure `dmu=k^2 dk/(2 Pi^2)`. Pi denotes one declared positive
constant throughout the calculation. It may be mathematical pi or the
registered represented constant
`14488038916154245685/4611686018427387904`; substituting one for the other
mid-ledger changes the target. All rational bounds below use only Pi>=3.

\[
M_K=\int_0^K d\mu/(2k)=K^2/(8\mathrm{Pi}^2).
\]

The exact finite-band identity requires suitable local integrability and
time-differentiability, not the existence of the K-to-infinity renormalized
limit. On a compact band bounded residual envelopes make the `1/k` factors
integrable because `dmu/k=O(k) dk`. The k=0 plane-wave normalization itself
is not evaluated. More general initial states require explicit infrared
integrability assumptions. A fixed time-independent discrete measure also
obeys the same identity when all operators use that same measure.

If a discrete mode measure has `M_d=sum weights*mu/(2k)` and contacts
are instead the analytic finite band with `M_A=M_K`, the direct hybrid
integrand obeys, for f=g and exact mode flow,

\[
\int_a^b F_A
=\Delta(G_d+C_R^A)+(M_A-M_d)\int_a^b Lg.            \tag{I7}
\]

For general f the additional term is
`-M_d integral Lf+M_A integral Lg`. The pressure shift
`C_P^mix=C_P^A+(M_A-M_d)g/3` removes the displayed moment correction
by declaring a changed pressure target. It cannot silently replace the
original analytic pressure. This one moment correction does not bound the
remaining continuum momentum error.

With uniform normalized complex residual bounds
`|r_u|<=q_u`, `|r_w|<=q_w`, a simple sufficient band bound is

\[
\left|\int_0^K E_R d\mu\right|
\le\frac{K^4+3L^2K^2}{8\mathrm{Pi}^2}q_u
 +\frac{K^3/3+|L|K^2/2}{4\mathrm{Pi}^2}q_w.         \tag{I8}
\]

Using joint complex-modulus control can improve the second coefficient to
`[(K^2+L^2)^(3/2)-|L|^3]/(12 Pi^2)` by Cauchy--Schwarz.
For raw residual bounds the entire right side is divided by |epsilon|.
The forcing mismatch adds `|L|M_K |f-g|`; independently bounded contacts,
temporal quadrature, endpoint rounding and initial-state errors require
their own terms.

A time-dependent cutoff or measure adds boundary or weight-derivative terms.
Fixed-comoving-K conservation is not automatically the identity at a
fixed physical momentum cutoff.

## 6. Independent closed retarded stress kernels

Consider a real forcing change delta=f-g with **identical incoming states
and unchanged contacts**. Raw-mode differences have zero c_R difference.
Let d=t-s and define

\[
A_0(d)=\int_0^K\sin(2kd)dk=(1-\cos(2Kd))/(2d),
\]
\[
A_1(d)=\int_0^Kk\cos(2kd)dk
=K\sin(2Kd)/(2d)+(\cos(2Kd)-1)/(4d^2),
\]
\[
A_2(d)=\int_0^Kk^2\sin(2kd)dk
=-K^2\cos(2Kd)/(2d)+K\sin(2Kd)/(2d^2)
  +(\cos(2Kd)-1)/(4d^3).
\]

These are entire functions by their compact integral definitions, with
`A0(0)=0`, `A1(0)=K^2/2`, `A2(0)=0`. Their poles in the displayed quotient
formulas are removable. An actual implementation should use their entire
series near zero, not subtract large nearly equal terms.

Independent substitution into direct R_m and P_m gives

\[
\mathcal K_R=L A_1/(4\mathrm{Pi}^2)-3L^2 A_0/(8\mathrm{Pi}^2),
\]
\[
\mathcal K_P=A_2/(6\mathrm{Pi}^2)
             +L A_1/(4\mathrm{Pi}^2)+L^2 A_0/(8\mathrm{Pi}^2).
\]

The band stress change is `delta R_K=integral K_R delta` and
`delta P_K=integral K_P delta`. The verifier proves all three upper-cutoff
derivatives and zero lower limits, all zero-lag limits, `A0'=2A1`,
`A1'=-2A2`, and

\[
(\partial_d+L^2\partial_L)\mathcal K_R
 =L(\mathcal K_R-3\mathcal K_P),\qquad
\mathcal K_R(L,0)=LM_K.
\]

Thus differentiation of the retarded integral returns exactly the source
mismatch in (I6). This pressure kernel was obtained from the direct pressure
operator before applying the Ward check.

## 7. Cauchy source-tail consequence; explicit limited scope

The inherited degree-24 exact Taylor polynomials on 64 panels have
radius 1/8, halfwidth 1/128 and complex bound |g|<=64. Their analytic
pointwise tail is `1/(15*2^90)`. Integrating the coefficient-tail envelope
gives

\[
\alpha=\int_a^b |g-P_{24}|ds\le1/(390\,2^{90}).
\]

This is a theorem about **exact Taylor coefficients**. It is not an error
bound on rounded coefficients, phase computation, the numerical state or
the complete nine-probe output. It holds on the registered interior interval
only; it says nothing about the nonanalytic bump support endpoints.

Since `|A0|<=K`, `|A1|<=K^2/2`, `|A2|<=K^3/3`,

\[
|\delta R_K|\le\alpha\frac{|L|K^2+3L^2K}{8\mathrm{Pi}^2},
\]
\[
|\delta P_K|\le\alpha\left[\frac{K^3}{18\mathrm{Pi}^2}
       +\frac{|L|K^2}{8\mathrm{Pi}^2}+\frac{L^2K}{8\mathrm{Pi}^2}\right],
\]
\[
\left|M_K\int L(P_{24}-g)ds\right|
\le\alpha M_K\sup|L|.
\]

All prefixes of the unit interval have the same upper bounds by positivity.
Using `|L|<=2/7`, `Pi>=3`, `K<=256`, exact rational arithmetic proves

| Component | Exact K=256 bound with Pi>=3 | Strict decimal upper |
| --- | --- | --- |
| Raw density source truncation | 899/1663385213724160574061721681920 | 5.407e-28 |
| Raw pressure source truncation | 3219337/14970466923517445166555495137280 | 2.152e-25 |
| Integrated source mismatch | 1/1856456711745714926408171520 | 5.389e-28 |

These cover the continuous momentum band for this component; they do not
reinterpret nine computed momentum probes as a complete continuum result.
The pressure component exceeds the separate `1e-26` source-moment gate.
That gate has not been changed or assigned to a different observable.

Remaining requirements include true incoming-state bounds, validated
coefficient and propagation arithmetic, normalized ODE defects over the
continuous band, full action/contact evaluation, geometry, temporal
quadrature, endpoint serialization and cutoff tails. No full twelve-case
pressure/contact gate can be declared passed from the above component alone.

## 8. Verification and limitations

`verify_independent_ward.py` with SymPy 1.14.0 passes **42 exact identities**
and rejects **10 symbolic omission controls**, in both normal and optimized
Python. No check relies on Python assertions. Independent source-free,
active-source, matched-stress and source/operator documents were read as
model specifications; their numerical producers were not run.

The initial receipt and the pre-contact extension receipt are retained as
development history, together with their exact verifier sources. The final
normal and optimized receipts name the final verifier's SHA256. Stable
content agrees after removing only the explicit optimization-mode flag.

This is internal independent algebra review, not external peer review or
proof-assistant formalization. The scalar action and its conservation law
are established mathematics; external novelty is not assessed. Fixed mass,
linearization, exact model background and declared finite contacts remain
premises. No new evidence for thermalization, a hot Big Bang, nonlinear
Einstein evolution or a higher-dimensional cause is claimed.
