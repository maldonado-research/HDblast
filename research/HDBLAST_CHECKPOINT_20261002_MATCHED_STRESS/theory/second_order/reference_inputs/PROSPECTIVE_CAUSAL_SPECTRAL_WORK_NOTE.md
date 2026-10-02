# Compact-pulse spectral energy and the matched work contact

Prepared 2 October 2026. **Analytic companion only; no numerical response,
mode evolution, or new test has been evaluated for this note.** The identities
below use established linear-response and Bogoliubov methods. They are not a
novelty or priority claim and do not establish shell dissipation, thermalization,
or quantum stability.

## 1. Setting and Fourier convention

Use exactly the fixed-geometry, fixed-state reference in
`PROSPECTIVE_CAUSAL_RESPONSE_PROTOCOL.md`: minimal coupling, incoming BD state,
`x0=r=2H²`, `a=-1/(H eta)`, and `v=a chi`. In the registered units `H=1`,
`r=2`. Let the real canonical mass source

    s(eta)=a(eta)² delta x(eta)

be smooth and compactly supported strictly before `eta=0`. Extend `s` by zero
to the real line when taking Fourier transforms. This is a mathematical
extension of the source and free canonical equation; it does not extend the
physical de Sitter chart through `eta=0`. Our convention is

    S(omega)=integral_R s(eta) exp(-i omega eta) deta,
    s(eta)=integral_R [domega/(2pi)] S(omega) exp(+i omega eta).

Write `s=epsilon f` with fixed smooth compact `f`. Every work identity below
is for the coefficient quadratic in the source. Higher-order finite-amplitude
particle production is not computed here.

## 2. Positive canonical excitation energy

The forced equation and the post-pulse coefficients are

    v_k''+(k²+s)v_k=0,
    v_k=(alpha_k exp(-ik eta)+beta_k exp(+ik eta))/sqrt(2k),
    beta_k^(1)=i S(2k)/(2k).

The canonical Hamiltonian after the pulse is the free plane-wave Hamiltonian.
Its excitation energy per comoving volume has quadratic term

    Ecan^(2)[s]
      =1/(2pi²) integral_0^infinity dk k³ |beta_k^(1)|²
      =1/(8pi²) integral_0^infinity dk k |S(2k)|²
      =1/(32pi²) integral_0^infinity domega omega |S(omega)|².

Smooth compact support makes the ultraviolet integral convergent; bounded
`S(omega)` makes its infrared endpoint convergent. For real, nonzero `s`,
`Ecan^(2)>0`: vanishing weighted spectral integral would force `S=0` away
from zero, and continuity and Fourier uniqueness would then force `s=0`.
The exact finite-amplitude excitation energy equals this quadratic term plus
`O(epsilon³)`. Only the quadratic term is claimed in the present calculation.

In particular this functional is proportional to a fractional Sobolev
seminorm. Define its convention explicitly:

    ||s||_(Hdot^(1/2))²
      =integral_R [domega/(2pi)] |omega| |S(omega)|².

Then

    Ecan^(2)=||s||_(Hdot^(1/2))²/(32pi)
      =1/(64pi²) integral_R deta integral_R du
                    [s(eta)-s(u)]²/(eta-u)².

The last factor follows from Parseval and
`integral_R [1-cos(omega z)]/z² dz=pi|omega|`. It expresses nonlocal memory
through a positive quadratic form, without assigning a thermal interpretation.

## 3. The inherited composite has a background-dependent contact

Let `q_M=a² delta Q_ren` be the inherited matched first-order response and
`M(eta)=a(eta)sqrt(r)`. The protocol gives

    q_M(eta)=-1/(8pi²) integral_-infinity^eta du s'(u)
                    {ln[M(eta)(eta-u)]+gamma_E+1}.

For any constant conformal mass scale `mu>0`, introduce the auxiliary
decomposition

    q_mu(eta)=-1/(8pi²) integral_-infinity^eta du s'(u)
                    {ln[mu(eta-u)]+gamma_E+1},
    q_M=q_mu-ln(M/mu)s/(8pi²).

This defines an analytic bookkeeping device, not a change of the inherited
finite physical action. The `mu` dependence cancels in `q_M`.

For the Fourier convention in section 1, the retarded multiplier of `q_mu` is

    chi_mu(omega)
      =[ln(|omega|/mu)-1]/(8pi²)+i sgn(omega)/(16pi).

For example, insert a convergence factor `exp(-zeta tau)`, `zeta>0`, in the
logarithmic kernel and use

    integral_0^infinity exp[-(zeta+i omega)tau]
       [ln(mu tau)+gamma_E+1] dtau
      =[1+ln(mu)-ln(zeta+i omega)]/(zeta+i omega).

Multiplication by `-i omega/(8pi²)` followed by `zeta -> 0+` gives the stated
multiplier as a distribution. With forward transform `exp(+i omega eta)`,
its imaginary part instead has a minus sign. No physical sign changes.

Parseval now gives the source-work identity

    W_mu^(2)=(1/2) integral_R s'(eta) q_mu(eta) deta
      =(1/2) integral_R [domega/(2pi)]
                     omega Im[chi_mu(omega)] |S(omega)|²
      =Ecan^(2).

Constant real contact terms make zero compact-pulse work because
`integral s's=0`. The inherited time-dependent contact does not:

    W_M^(2) := (1/2) integral_R s' q_M deta
      =Ecan^(2)+1/(32pi²) integral_R (M'/M) s² deta.

On this de Sitter geometry, `M'/M=a'/a=-1/eta>0`. Therefore the exact
canonical-energy identity in the inherited convention is

    Ecan^(2)=W_M^(2)-1/(32pi²) integral_R (a'/a) s² deta.

Both terms on the right are required. The extra term is independent of the
constant positive reference `r`, although the response itself depends on `r`.
It is explicit background dependence of a local primitive:

    F_local(eta,s)=-ln(M/mu)s²/(32pi²),
    partial F_local/partial s=(q_M-q_mu)/2.

Its endpoints vanish for a compact pulse, while its explicit time derivative
does not. Consequently its source work is
`-integral partial_eta F_local deta=integral (M'/M)s²/(32pi²) deta`.
Calling this term a harmless total derivative would miss the background drift.

## 4. Independent regulated time-domain proof

An independent proof uses precisely the protocol's combined hard cutoff. Set

    A_K(M)=asinh(K/M)-K/sqrt(K²+M²),
    J_k(eta)=integral_-infinity^eta du s(u) sin[2k(eta-u)],
    q_M,K=[s A_K(M)-2 integral_0^K J_k dk]/(8pi²).

Integrating the bare contribution to `(1/2)integral s' q_M,K` by parts gives

    -1/(8pi²) integral_0^K dk integral_R deta s' J_k
      =1/(8pi²) integral_0^K dk k |S(2k)|²
      =Ecan,K^(2).

The step uses `J_k'=2k integral_-infinity^eta s(u)cos[2k(eta-u)]du` and
symmetrizes the triangular integration domain. Compact support removes the
boundary term. Direct differentiation of the local subtraction gives

    A_K'=-[M'/M] K³/(K²+M²)^(3/2),

and hence the exact finite-cutoff identity

    W_M,K^(2)=Ecan,K^(2)
      +1/(32pi²) integral_R deta s² (M'/M)
                                  K³/(K²+M²)^(3/2).

Taking `K -> infinity` proves the same matched identity without using the
Fourier transform of a singular kernel. Smooth compact support and bounded
positive `M` on the pulse support justify the removed-cutoff limit.

## 5. Why this does not determine physical minimal-stress energy transfer

Canonical excitation energy is not the full minimally coupled stress energy.
Let `h=a'/a`, and compare the evolved state with the original BD state after
the source has switched off. For these state differences, local counterterms
cancel because the geometry and mass have returned to the same values.
With `Delta q=Delta<v²>` including its coherent part, the exact operator
relation on this reference gives

    a⁴ Delta rho
      =Delta Ecan-(h/2)(Delta q)'+(3h²/2)Delta q.

Indeed `a⁴ rho` is the expectation of
`[(v'-hv)²+(grad v)²+a²x0 v²]/2`, and `a²x0=2h²` here.
The last two terms contain polarization and coherence; they can already
occur at first order in the pulse. Thus positive canonical occupation energy
alone is not a positive physical `Delta rho`, nor a complete redshift law.

The distinction also appears directly in the induced source-work terms:

    (1/2) integral a⁴ delta x' delta Q_ren deta
      =W_M^(2)-integral h s q_M deta.

The physical Ward identity further includes expansion work,
`(a⁴rho)'=h a⁴(rho-3p)+(a⁴Q/2)x'`. Stress response and its matching contact
terms are required to apply it. The present note neither computes those
kernels nor assigns the prescribed external pulse to a dynamical shell.

The Bogoliubov route, regulated time-domain work route, and spectral route
agree analytically. A separate theory agent independently checked the energy
coefficient, Fourier-sign convention, and background-contact drift. Any
future numerical test of these identities requires its own public choices
and error controls before evaluation.
