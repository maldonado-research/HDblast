# Exact theory inputs for the bounded causal-response calibration

Prepared 2 October 2026. This document and its companion verifier contain
analytic inputs, not a completed numerical experiment. The enclosing public
registration must freeze the actual code, source identifiers, mode methods,
cutoff and time quadrature ladders, tolerances, mutation rules, and stopping
gates before any source response is evaluated. The inherited normalization
and finite prescription are those of `PROSPECTIVE_CAUSAL_RESPONSE_PROTOCOL.md`.

## Geometry, state, pulse, and dimensions

Choose de Sitter units with H=1, r=2, a(eta)=-1/eta, eta_i=-6,
eta_f=-3/2. H denotes the unit inverse length; epsilon=10^-4 is an amplitude
in units of H^2, the center is -4/H, and the half-width is 1/H.
Define u=eta+4 and the normalized smooth bump

    B(u)=exp(1-1/(1-u^2)) for |u|<1, and zero otherwise.
    s_positive(eta)=epsilon B(u),
    s_signed(eta)=epsilon u B(u).

Thus support is (-5,-3). The prescribed observations are
(-11/2,-9/2,-4,-7/2,-5/2,-3/2). They include a strictly earlier observation,
three interior observations, and two later observations. B(0)=1; this is a
factor e times the unnormalized example in the inherited proposal. Every
source derivative vanishes at both support endpoints and near eta_i.

The physical mass perturbation is delta x=s/a^2=eta^2 s. In this interval
|delta x|<=25 epsilon; hence 2+delta x>0 without source sampling. The initial
state is the original incoming BD state, held fixed. eta_i is an integration
boundary in a source-free past, not a reset vacuum condition.

In physical dimensions, [eta]=-1, [k]=[M]=[H]=1, [a]=0,
[x]=[r]=[s]=[epsilon]=[Q]=2, [rho]=[p]=4. The inherited phi and b are
dimensionless, so [j]=4. Then [A_eta]=4 and the omitted-tail expression has
dimension two. All logarithms below have dimensionless arguments. Under
eta->L eta, k->k/L, H->H/L, x,r,s->(x,r,s)/L^2, a is invariant and
delta Q->delta Q/L^2. Derivatives and integration measures must be rescaled
with the same L.

## Independent normalization and matched contact

The canonical interaction is H_I=(s/2) integral v^2 and the unperturbed
mode is exp(-ik eta)/sqrt(2k). Wick's connected commutator has spatial
integral `-i/(2 pi^2) integral dk sin(2k Delta)`; multiplying by the Kubo
factor -i/2 gives the negative coefficient -1/(4 pi^2). Independently,

    delta v_k''+k^2 delta v_k=-s v_k,
    delta v_k(eta_i)=delta v_k'(eta_i)=0,
    delta |v_k|^2= -1/(2k^2) integral_{eta_i}^eta s(t)sin[2k(eta-t)]dt.

The momentum measure k^2 dk/(2 pi^2) gives the same result. Starting from
the physical chi interaction gives a(t)^4 delta x(t) times
[a(eta)a(t)]^-2, exactly a(eta)^-2 s(t).

At fixed geometry and fixed r, the order-two subtraction varies by
`delta |v|_sub^2=-s/(4 omega_r^3)`, omega_r=sqrt(k^2+a^2 r).
No derivative of the source enters omega_r. The combined renormalized
momentum response is

    delta Q = 1/(8 pi^2 a^2) integral_0^infinity dk [
        s(eta) k^2/(k^2+M^2)^(3/2)
        -2 integral_{eta_i}^eta s(t)sin(2k(eta-t))dt ],
    M=a sqrt(r).

The terms are combined before the cutoff is removed. Writing T=eta-eta_i,
the exactly equivalent finite-cutoff response is

    delta Q_K = 1/(8 pi^2 a^2) [
      s(eta)(asinh(K/M)-K/sqrt(K^2+M^2))
      -integral_0^T s(eta-tau)(1-cos(2Ktau))/tau dtau ].

The removed-cutoff finite part is

    delta Q = -1/(8 pi^2 a^2) [
       integral_0^T (s(eta-tau)-s(eta))/tau dtau
       +s(eta)(ln(MT)+gamma_E+1) ],

or, using s(eta_i)=0,

    delta Q = -1/(8 pi^2 a^2) integral_{eta_i}^eta
       s'(t)[ln(M(eta-t))+gamma_E+1]dt.

The finite `+1` and Euler constant are inherited prescription terms. Their
deletion changes the model. A useful stable implementation splits at support
boundaries and observation time. For a remaining logarithmic endpoint,
subtract its endpoint coefficient analytically or apply a frozen endpoint
transformation; the quadrature implementation must declare which it uses.

Before eta=-5 both terms vanish exactly. At post-pulse observations,

    delta Q = -1/(8 pi^2 a^2) integral_{-5}^{-3} s(t)/(eta-t)dt.

It is strictly negative for the positive pulse. For the signed pulse it is
also strictly negative: pair u and -u on (0,1), obtaining
`u B(u)[1/(d-u)-1/(d+u)]>0`, with d=eta+4>1. This signed-source sign is
an analytic consequence of this chosen odd pulse, not a universal sign for
signed sources.

## Pulse derivatives and constructive ultraviolet tail

For |u|<1 all derivatives below are with respect to u; outside, including
the endpoints, they are zero:

    B'   = -2u B/(1-u^2)^2,
    B''  = 2(3u^4-1)B/(1-u^2)^4,
    B''' = -4u(6u^6+3u^4-10u^2+3)B/(1-u^2)^6,

    (uB)'   = (u^4-4u^2+1)B/(1-u^2)^2,
    (uB)''  = 2u(u^4+4u^2-3)B/(1-u^2)^4,
    (uB)''' = -2(3u^8+24u^6-26u^4+3)B/(1-u^2)^6.

For arbitrary half-width sigma, s^(n)=epsilon f^(n)(u)/sigma^n.
Consequently the derivative norm and A scale as |epsilon|/sigma^2.
Three integrations by parts, with zero initial derivatives, give

    | integral_0^T s(eta-tau)sin(2ktau)dtau -s(eta)/(2k) |
        <= A_eta/(8k^3),
    A_eta=|s''(eta)|+integral_{eta_i}^eta |s'''(t)|dt.

Combining this with the elementary subtraction inequality gives

    |delta Q-delta Q_K|
      <= [3|s(eta)|M^2/2+A_eta/4]/[16 pi^2 a^2 K^2].

The independent analytic certificate in
`review/INDEPENDENT_ANALYTIC_AUDIT_20261002.md` proves
`integral |B'''|<104`, `sup |B''|<23`,
`integral |(uB)'''|<102`, and `sup |(uB)''|<19`.
Thus a simple safe bound is A<=127 epsilon for B and A<=121 epsilon for uB;
after the pulse it improves to 104 epsilon and 102 epsilon. These are proved
envelopes, not numerical norm estimates.

The companion `tail_bounds.py` uses a tighter analytic total variation.
For q=1/(1-u^2), define

    b(q)=e^(1-q)(4q^4-12q^3+6q^2),
    f(q)=sqrt(1-1/q)e^(1-q)(4q^4-12q^3+2q^2).

Let alpha_B,beta_B be the two q>1 roots of
`2q^3-14q^2+21q-6`, bracketed in (3/2,13/8) and (5,81/16).
Let alpha_F,beta_F be the two q>1 roots of
`4q^4-32q^3+64q^2-36q+3`, bracketed in (7/4,15/8) and (5,6).
There are exactly two roots above one for each polynomial, as proved by
Descartes' rule after q=1+z and the sign-changing brackets. Put
D_B=b(beta_B)-b(alpha_B), D_F=f(beta_F)-f(alpha_F). Exact A/epsilon is

| eta | positive_B | signed_uB |
| --- | --- | --- |
| -11/2 | 0 | 0 |
| -9/2 | 2D_B | 2D_F |
| -4 | 2D_B | 2D_F |
| -7/2 | 2D_B-4+(832/81)e^(-1/3) | 2D_F+(992/81)e^(-1/3) |
| -5/2, -3/2 | 4D_B-4 | 4D_F |

The interval implementation freezes 256 exact rational root bisections and
60 decimal digits of directed interval arithmetic. It encloses the primitive
values, performs the tail arithmetic including pi and a in intervals, and
advances the final binary64 upper endpoint toward positive infinity.
Its normalized bound directly covers y=a^2 delta Q/epsilon. Importing the
module performs no numerical source evaluation. Its source functions are
to be called only after registration. Static independent review checked the
closed formulas and enclosure construction; runtime validation belongs to
the subsequent registered run.

The analytic tail covers only k>K. Mode error, finite-k quadrature error,
time quadrature error and floating-point cancellation require separate
evidence. A valid frozen acceptance gate may compare the actual route
disagreement to this tail plus fixed declared numerical allowances. A
finite ladder must stop at its declared maximum; exceeding a budget is a
diagnostic failure, not permission to change contact terms, pulse data,
state, or tolerances after seeing outputs.

## Static derivative, current contact, state, and Ward controls

For the separate past-infinite source s=delta x/(H^2 eta^2), the matched
susceptibility is

    Q_x=-(2 gamma_E+ln(r/H^2))/(16 pi^2),
    Q_x at r=2H^2=-(2 gamma_E+ln2)/(16 pi^2).

The common-action identity Psi(2)=1-2 gamma_E gives the same value.
Letting r track x would add Q_r=1/(16 pi^2); it is a wrong-reference control.
The stationary source is not the registered compact-pulse history.

For x(phi)=r[1+b(phi-phi0)/2]^2,

    x_phi,0=br, x_phiphi,0=b^2r/2,
    delta j=(br/2)delta Q+(b^2r Q0/4)delta phi,
    Q0=H^2/(12 pi^2).

The second term is a separate mass-law contact. The source amplitude gamma
multiplies the complete translated current once. At b=0 its response is zero.

For an occupation-diagonal Gaussian state with n(k), Wick contractions give
(1+n)^2-n^2=1+2n. The extra retarded kernel is

    delta K_state=-theta(Delta)/(2 pi^2 a^2) integral dk n(k)sin(2kDelta).

For n>=0 supported below k_max, it is negative when
0<Delta<pi/(2k_max). A general squeezed Gaussian state also has anomalous
correlations; the stated diagnostic specifically uses the occupation state.
For fixed zero source and a small initial negative-frequency amplitude beta,

    delta Q_initial=1/(2 pi^2 a^2) integral dk k Re[beta(k)e^(2ik(eta-eta_i))].

The initial Wronskian changes only at second order in beta. This finite
boundary variation is absent in a source-only response and must not be
silently reset during a run. Moving eta_i within the same zero-source past
preserves the source response when the same incoming state is maintained.

The linear common-action Ward identity remains

    delta rho_dot+3H(delta rho+delta p)=Q0 delta x_dot/2.

In conformal time it is
`delta rho'+3(a'/a)(delta rho+delta p)=Q0 delta x'/2`.
The product delta Q delta x_dot first enters at second order. This calibration
has not supplied stress kernels or demonstrated that the physical stress
response obeys this identity; the symbolic check verifies its order and
dimensions only.

## Analytic spectral companion and scope boundary

For a real compact source define q_M=a^2 delta Q and choose a constant
auxiliary conformal scale mu. Then

    q_M=q_mu-ln(M/mu)s/(8 pi^2).

The canonical free oscillator modes after the pulse have leading excitation
energy per comoving volume

    E_can^(2)=1/(32 pi^2) integral_0^infinity omega |S(omega)|^2 domega,
    S(omega)=integral s(eta)e^(-i omega eta)deta.

Fourier analysis of the retarded finite-part kernel, or the finite-cutoff
identity, gives

    (1/2) integral s' q_M deta
      =E_can^(2)+1/(32 pi^2) integral (M'/M)s^2 deta.

The drift sign is positive after integrating the contact contribution by
parts; compact endpoints eliminate the boundary term. Here M'/M=-1/eta>0.
The finite-cutoff local primitive satisfies
`A_K'=-(M'/M)K^3/(K^2+M^2)^(3/2)`, confirming the same sign and limit.
Subtracting this known reference drift is necessary before equating that
matched source work to the canonical excitation energy. These are analytic
identities; no additional spectral grid or numerical test is introduced.

Canonical oscillator energy is not the complete minimally coupled physical
stress tensor. This scalar-response calibration does not establish stress
or metric kernels, quantum-corrected shell stability, a cosmological initial
state, physical radiation production, heating, thermalization, or discovery.
