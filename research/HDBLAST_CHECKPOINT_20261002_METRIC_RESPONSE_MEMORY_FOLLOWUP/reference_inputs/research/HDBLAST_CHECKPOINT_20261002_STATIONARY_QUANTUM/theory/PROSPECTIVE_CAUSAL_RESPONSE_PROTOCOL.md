# Prospective fixed-geometry causal-response calibration

Prepared 2 October 2026. **Analytic proposal only; no new numerical experiment
is registered or executed by this document.** The next bounded task is to
calibrate the homogeneous retarded variance response on one prescribed de
Sitter geometry, at the exact mass reference already checked independently
by modes and by the Euclidean action. This supplies one missing causal
ingredient. It does not calculate quantum-corrected stability, solve an
initial-value shell problem, or establish energy transfer or heating.

The stationary project has selected a finite common-action prescription and
found static numerical roots. Its sphere-restricted action does not determine
the time-dependent retarded Hessian. The construction below derives one
component directly in real time and fixes its local term by the inherited
positive-reference subtraction. It is suitable as the next automation
priority before any proposed full quantum stability or initial-data study.

## 1. Exact reference, state, and perturbation class

Use physical units with a dimensionless scale factor, conformal time eta of
dimension length, and comoving momentum k of dimension mass:

    ds²=a(eta)²[-deta²+dX²],
    a(eta)=-1/(H eta), eta<0, H>0,
    x0=r=2H²,       xi=0,
    chi=v/a,        v_k(eta)=exp(-ik eta)/sqrt(2k).

This is a minimally coupled massive scalar, not a change to the conformally
improved stress tensor. Its canonical frequency happens to be exactly k:

    v_k''+[k²+a²x0-a''/a]v_k=0,
    a²x0=a''/a=2/eta².

The exact plane waves define the same Euclidean/Bunch-Davies state used in the
inherited de Sitter checkpoint. The subtraction remains that of the original
minimal action. In particular the renormalized reference values are

    Q0=H²/(12pi²),    rho0=11H^4/(960pi²),    p0=-rho0.

They need not vanish just because the canonical modes are plane waves.

Keep geometry and r fixed. Prescribe an infinitesimal physical mass-squared
perturbation delta x(eta), and define

    s(eta)=a(eta)² delta x(eta).

Work first on a compact observation interval eta_i<=eta<=eta_f<0. Take s to
be smooth and identically zero in a neighborhood of eta_i; a smooth compact
pulse entirely inside this interval is the simplest test. The initial density
matrix is the original BD state and is held fixed as the source is varied.
The background exists before eta_i; eta_i is a convenient lower integration
limit, not a new finite-time vacuum prescription. Require x0+delta x>0 when
a finite source is used only to check the linear limit.

The first calibration treats s as an externally prescribed source. It does
not impose a bulk or shell equation on it. The same formulas can be expressed
in inherited length units by rescaling eta,k,H,x,r consistently with L.

## 2. Retarded commutator and its independent mode check

At fixed geometry the canonical interaction Hamiltonian is

    H_I(eta)=s(eta)/2 integral d³X v(eta,X)².

With Wightman function

    G^>(eta,X;eta',0)
      =integral d³k/(2pi)³ exp[-ik(eta-eta')+ik.X]/(2k),

Wick's theorem gives the connected commutator

    <[v²(eta,0),v²(eta',X)]>
      =2[(G^>)²-(G^<)²].

Disconnected products cancel. Writing Delta=eta-eta'>0, the spatial integral
is

    integral d³X (G^>)²
      =1/(8pi²) integral_0^infinity dk exp(-2ik Delta),

    integral d³X <[v²(eta,0),v²(eta',X)]>
      =-i/(2pi²) integral_0^infinity dk sin(2k Delta).

The Kubo sign follows from

    delta<O(eta)>=i integral <[H_I(eta'),O(eta)]> deta'
                   =-i/2 integral s(eta')<[O(eta),v²(eta')]> deta'.

Consequently

    delta Q_bare(eta)
      =-1/[4pi² a(eta)²] integral_eta_i^eta deta' s(eta')
                                      integral_0^infinity dk sin(2k Delta).

Equivalently, starting with the physical chi interaction introduces
a(eta')^4 delta x and the chi commutator contributes
[a(eta)a(eta')]^-2. Their product is exactly a(eta)^-2 s(eta'); this checks
both scale-factor powers in the displayed response.

An independent forced-mode derivation uses

    delta v_k''+k² delta v_k=-s(eta)v_k,
    delta v_k(eta_i)=delta v_k'(eta_i)=0,

    delta v_k(eta)=-integral_eta_i^eta
                    sin[k(eta-eta')]/k s(eta') v_k(eta') deta'.

It gives exactly

    delta|v_k(eta)|²=-1/(2k²) integral_eta_i^eta
                              s(eta') sin[2k(eta-eta')] deta'.

Inserting the momentum measure k²dk/(2pi²) reproduces the Kubo result.
These two routes independently fix the minus sign, Wick factor two, and
factor 2k in the memory kernel.

## 3. A regulated kernel, and why its separated-time limit is incomplete

An Abel factor exp(-epsilon k), epsilon>0, gives the ordinary retarded kernel
acting on s:

    K_epsilon(eta,eta')
      =-theta(Delta) Delta/[2pi² a(eta)²(epsilon²+4Delta²)].

For separated times its limit is

    K_nonlocal=-theta(Delta)/[8pi² a(eta)² Delta].

This is not an integrable kernel at coincidence. A bare statement of
"PV(1/Delta)" is insufficient: the response is one-sided, and there is no
symmetric principal-value cancellation across Delta=0. A retarded finite-part
extension and its matched local term must be specified. Multiplying unrelated
singular distributions by theta does not define that extension.

The regulator should multiply the bare response and its subtraction before
taking the limit. A hard comoving momentum cutoff used consistently in both
pieces gives an equally useful finite regulator below. Neither cutoff is a
new physical EFT cutoff. The combined removed-cutoff answer is the target.

## 4. Fixed-reference subtraction fixes the local contact term

The inherited order-two variance subtraction uses the positive constant r,
with

    omega_r=sqrt(k²+a²r),
    U2=[a²(x-r)-a''/a]/(2omega_r)
         -omega_r''/(4omega_r²)+3omega_r'^2/(8omega_r³),
    |v|²_sub=1/(2omega_r)-U2/(2omega_r²).

Since the metric and r are fixed, only the first term of U2 varies:

    delta |v|²_sub=-s(eta)/(4omega_r³).

No source derivative is silently inserted into omega_r. The perturbation
x-r has adiabatic order two in this prescription; derivatives of it occur
above the required subtraction order for Q. The matched linear response is
therefore the convergent combined momentum integral

    delta Q_ren(eta)=1/[8pi² a(eta)²] integral_0^infinity dk {
       s(eta) k²/[k²+M(eta)²]^(3/2)
       -2 integral_eta_i^eta s(eta')sin[2k(eta-eta')] deta' },

    M(eta)=a(eta)sqrt(r).

This equation is already a concrete, mathematically well-defined target on
the stated source class. The two terms must be combined before removing the
momentum cutoff. Their separate integrals diverge logarithmically.

The combined integrand is integrable at k=0 and is O(k^-3) at infinity when
the source and initial derivatives vanish as stipulated. This result uses
the inherited matching. Adding a finite local (delta x)² action coefficient
would add a local response and would define a different finite model unless
the other matched couplings were compensated.

## 5. Exact finite-part memory formula

Let T=eta-eta_i>0, s0=s(eta), and M=a(eta)sqrt(r). At a finite hard cutoff K,
the preceding expression is exactly

    delta Q_K=1/[8pi² a²] {
       s0[asinh(K/M)-K/sqrt(K²+M²)]
       -integral_0^T s(eta-tau)[1-cos(2K tau)]/tau dtau }.

For finite K the integrand at tau=0 has a regular removable limit. Split
s(eta-tau)=s0+[s(eta-tau)-s0]. The exact identity

    integral_0^T [1-cos(2K tau)]/tau dtau
      =gamma_E+ln(2KT)-Ci(2KT)

and the large-K asymptotic of the local subtraction give

    delta Q_ren(eta)=-1/[8pi² a(eta)²] {
       integral_0^T [s(eta-tau)-s(eta)]/tau dtau
       +s(eta)[ln(MT)+gamma_E+1] }.

The first integrand has limit -s'(eta) as tau->0. Every logarithm has a
dimensionless argument. The constants gamma_E+1 belong to this explicitly
matched finite prescription and cannot be discarded as an unspecified
contact term.

An independent Abel-regulator derivation uses

    integral_0^infinity exp(-epsilon k)
               k²/(k²+M²)^(3/2) dk
      =-ln(epsilon M/2)-gamma_E-1+o(1),

and gives exactly the same finite term. This equality is analytic; no new
numerical regulator comparison has been run.

Integrating by parts gives a particularly convenient logarithmic memory form:

    delta Q_ren(eta)=-1/[8pi² a(eta)²] integral_eta_i^eta s'(eta')
          {ln[a(eta)sqrt(r)(eta-eta')]+gamma_E+1} deta'.

The logarithmic endpoint singularity is integrable. This equivalent form
uses s(eta_i)=0. Resetting an instantaneous vacuum, inserting an abrupt
nonzero source at eta_i, or differentiating an initial density matrix adds
boundary/state terms and is outside the displayed derivation. For the
declared compact smooth pulse, moving eta_i further into its zero-source past
leaves the response invariant.

One may describe the nonlocal part as a one-sided finite-part distribution
with an arbitrary auxiliary conformal scale mu, provided its local term is
simultaneously changed by ln(M/mu). The sum above has no such auxiliary-scale
ambiguity. The physical positive reference r remains fixed.

## 6. Constructive momentum-tail bound for future numerical registration

Let

    I_k=integral_0^T s(eta-tau)sin(2k tau)dtau,
    A_eta=|s''(eta)|+integral_eta_i^eta |s'''(eta')|deta'.

Repeated integration by parts, using the vanishing initial derivatives,
gives for every k>0

    |I_k-s(eta)/(2k)| <= A_eta/(8k³).

Also

    |k²/(k²+M²)^(3/2)-1/k| <= 3M²/(2k³).

Thus the omitted tail in the COMBINED momentum response obeys

    |delta Q_ren-delta Q_K|
      <= [3|s(eta)|M²/2+A_eta/4]/[16pi² a(eta)² K²].

The bound is explicit once a pulse and derivative norm are frozen. It covers
the ultraviolet tail, not finite-k quadrature, time quadrature, derivative
evaluation, or floating-point error. The compact source and eta_f<0 keep
all required scale-factor and derivative norms bounded on the test interval.
An implementation must evaluate the actual bound rather than cite the
asymptotic order alone.

## 7. Exact stationary susceptibility check, with its state-limit caveat

A useful independent analytic check allows constant physical delta x present
from the infinite conformal past, so

    s(eta')=delta x/(H² eta'²).

This source is not the finite compact-pulse experiment. It defines a separate
analytic past-infinite limit with the unchanged incoming BD condition; s
decays in that past. Set A=-eta>0. The finite-part memory integral satisfies

    integral_0^T [s(eta-tau)-s(eta)]/tau dtau
      =s(eta)[-ln(1+T/A)-T/(A+T)].

Taking T->infinity in the full, already subtracted expression gives

    delta Q_ren/delta x
      =-[2gamma_E+ln(r/H²)]/(16pi²).

At r=2H² this is exactly

    Q_x|_(x=r=2H²)=-(2gamma_E+ln2)/(16pi²).

The inherited common-action expression independently yields this value:
Psi(2)=psi(2)+psi(1)=1-2gamma_E and
Q_x=[Psi(2)-ln(r/H²)-1]/(16pi²). Agreement tests the finite local constant,
the Kubo normalization, and fixed-r differentiation simultaneously. It does
not justify replacing a finite-time retarded history by an equilibrium
susceptibility. Different switching/state limits must be analyzed explicitly.

At the same point the inherited action has Q_r=1/(16pi²). Allowing r to
track x would instead add this term to the static derivative. This is an
exact wrong-reference control, not an allowed change to the causal model.

## 8. Source translation and the missing metric response

For the final stationary model's quadratic mass law at phi_b0,

    x_phi,0=b r,       x_phiphi,0=b²r/2,
    delta x=b r delta phi,
    delta j=(b r/2)delta Q+(b²r Q0/4)delta phi.

The second term is the instantaneous mass-law contact term. It is separate
from the renormalization contact term already present in delta Q. Both are
required. For b=0 the scalar-coupling perturbation vanishes. The loop/source
amplitude gamma multiplies the complete translated shell response once.

These formulas calibrate the response at the exact source reference
x=2H². A finite-amplitude stationary root generally has a slightly different
x/H² ratio. Its full causal response requires the propagator at that actual
root, or a separately justified expansion away from the reference; the exact
plane-wave kernel is not automatically its exact Hessian.

Even at fixed metric, delta Q alone does not give delta rho and delta p. The
linearized common-action Ward identity is

    d(delta rho)/dt+3H(delta rho+delta p)=(Q0/2)d(delta x)/dt.

It requires stress-current response kernels and explicit stress contact terms
derived from the same action. The product delta Q*d(delta x)/dt first enters
the energy-exchange equation at second order. A first-order variance memory
response cannot by itself demonstrate net particle energy, dissipation,
radiation, or heating.

For genuine shell stability one also needs metric-metric and metric-scalar
retarded kernels, their contact terms and Ward constraints, the bulk scalar
and gravitational perturbations with their boundary conditions, and a
specified initial state. Those ingredients are not supplied by the present
fixed-geometry calibration. Small gamma and a passed stationary dimensional
screen do not remove these missing equations or prove perturbative stability.

## 9. Exact controls to freeze before any numerical run

The following tests are mathematically distinct and should not be replaced
by several comparisons of the same implementation.

1. **Kubo/mode normalization.** The composite commutator and forced-mode
   derivations must give the same negative coefficient and sin(2kDelta).
   Deliberately reverse the Kubo sign, omit the Wick factor two, and replace
   2k by k; each mutation must fail the exact reference checks.

2. **Strict causality.** A source supported after the observation time gives
   identically zero response. An advanced or time-symmetric replacement must
   fail this test. Before a pulse begins, both renormalized formulas vanish.

3. **A signed memory test without local ambiguity.** For a nonnegative,
   nonzero smooth pulse entirely BEFORE observation, s(eta)=0 and

       delta Q=-1/(8pi²a²) integral_pulse s(eta')/(eta-eta') deta' < 0.

   This exact sign survives every local contact convention because the source
   vanishes at observation. A flipped sign, instantaneous-only response, or
   memory window excluding the earlier pulse must fail.

4. **Matching/contact calibration.** The analytic past-infinite constant
   delta x check must give -(2gamma_E+ln2)/(16pi²). Dropping the '+1', dropping
   gamma_E, or treating r as a moving reference changes that value. The
   finite-K and finite-part formulas must agree within their separately
   controlled errors, with the actual combined-tail bound recorded.

5. **Fixed state.** As a deliberately different Hadamard Gaussian state, use
   a smooth nonzero occupation n(k) with compact momentum support. Its
   commutator response contains the factor 1+2n(k), so the extra kernel is

       delta K_state=-theta(Delta)/(2pi²a²)
                       integral_0^infinity n(k)sin(2kDelta)dk.

   This term is ultraviolet finite and cannot be removed by state-independent
   local counterterms. For support k<=k_max and 0<Delta<pi/(2k_max), n>=0
   implies a strictly negative extra kernel. A code that silently returns the
   vacuum kernel for this changed state must fail. This mutation is a state
   diagnostic, not an alternative state in the registered primary experiment.

6. **Initial-state boundary term.** With zero source, add a small initial
   negative-frequency component beta(k), smooth with compact momentum
   support. At first order it produces

       delta Q_initial=1/(2pi²a²) integral dk k
                          Re[beta(k)exp(2ik(eta-eta_i))],

   using plane waves phased at eta_i. A source-only Kubo response is zero;
   this explicit nonzero term detects an unregistered reset or variation of
   the initial state. Wronskian normalization changes only at second order
   for this infinitesimal Bogoliubov variation.

7. **Conformal factors and dimensions.** Derive the response once from the
   physical chi interaction and once from v. Both must reproduce the source
   a(eta')² delta x and output factor a(eta)^-2. Omitting either factor must
   fail this equality and a consistent rescaling-of-coordinates check.

8. **Mass-law contact.** Translate to delta j using BOTH terms in section 8.
   Dropping x_phiphi Q0 delta phi/2 must fail the exact quadratic-law chain
   derivative. No scalar-current test replaces the absent stress Ward test.

Every automated positive and negative control must use explicit exceptions
or equivalent checks that remain active under Python optimization. Preserve
failed runs, registration hashes, source hashes, and output provenance.

## 10. Concrete future experiment and stopping boundary

The next registration should select one smooth compact pulse and one signed
pulse, for example a fixed multiple of

    B(t)=exp[-1/(1-t²)] for |t|<1, B(t)=0 otherwise,

with all amplitude, center, duration, and observation times frozen in declared
de Sitter units. Fix eta_i before both pulses and eta_f<0. Include observations
before, during, and after each pulse. Actual values, precision, cutoff ladder,
time quadrature, allowed errors, and failure/stopping gates must be public
before numerical evaluation. They are deliberately not selected from new
outputs in this analytic proposal.

Use three independently implemented routes where feasible: the forced
linearized canonical modes with the fixed initial state; the combined
momentum integral with its explicit tail bound; and the logarithmic
time-domain convolution. The latter two test quadrature and subtraction;
the first additionally tests the real-time mode response. A separately frozen
small-amplitude central difference of full prescribed-background modes may
test the linear limit, but it is optional and is not a coupled bulk solve.

The evidence package should distinguish analytic exact identities, numerical
quadrature/integration error estimates, actual tail bounds, initial-state
assumptions, and every wrong-formula control. Check the reference Q_x and
the post-pulse sign/memory identities before extending the domain or adding
variables. Failure stops at a bounded diagnostic; it must not trigger an
unregistered state change, local counterterm adjustment, or mass-law tuning.

Before implementation, the pulse derivatives and their tail-bound norms must
be derived/evaluated, a stable quadrature for the logarithmic endpoint fixed,
the finite-mode initial conditions and Wronskian checks frozen, and the
independent implementations reviewed. The fixed-geometry kernel and its
local contact coefficient above have been analytically derived and checked
by two agents; no numerical accuracy or implementation readiness claim
follows from that algebra alone.

Completing this calibration would establish one matched retarded scalar
susceptibility around the exact inherited reference. A subsequent separately
registered task could derive stress and metric response before assessing the
physically screened small-gamma stationary case. There is no present claim
of quantum stability, causal relaxation, a cosmological initial state,
particle-production yield, radiation transfer, thermalization, observation,
or scientific discovery.
