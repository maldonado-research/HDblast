# Independent analytic check of the prospective matched stress proposal

Scope: read-only analytic review of the inherited
`PROSPECTIVE_MATCHED_STRESS_RESPONSE.md`. No source samples, mode evolution,
stress integrals, response quadratures, or physical numerical tests were
evaluated. The expressions below were reduced algebraically from the
proposal's definitions. This is a targeted check, not a full derivation of
all fourth-order counterterms or a numerical implementation signoff.

Use constant H, x0=r=2H², a=−1/(H eta), L=a'/a, and a fixed comoving
cutoff K. Let d=delta x, s=a²d, and reserve z for the dimensionless cutoff
ratio so that it is not confused with a mode function:

    M²=a²r,
    z=K/sqrt(K²+M²),
    dot(z)=−H z(1−z²).

Dots in this note are cosmic-time derivatives. The state is the unchanged
incoming BD state, with a smooth source zero near the initial boundary.

## Minimal stress and explicit operator contacts

Write A_k=delta|v_k|². The forced-mode equation gives

    A_k''+4k²A_k=−s/k,
    delta|v_k'|²=−k²A_k.

These identities reproduce the proposal's direct bare density and pressure
and the explicit terms +s/(2k) and −s/(2k) in their respective brackets.
In particular the scale-factor bilinears survive although the reference
canonical frequency is k.

For the pressure reduction, delta Q_k=A_k/a² gives

    (delta Q_k''+L delta Q_k'−3L² delta Q_k)/(6a²)
      = [−(4k²/3+L²)A_k−L A_k'−s/(3k)]/(2a⁴).

The difference between this expression and the physical pressure is
−s/(12ka⁴)=−Q0,bare,k d/6. This checks the coefficient in the closed
pressure formula and distinguishes it from the explicit mass variation in
the original operator. The density reduction similarly gives
(L² delta Q_k−L delta Q_k')/(2a²)+Q0,bare,k d/2.

## Exact finite-band baseline

For w=sqrt(k²+M²), the reference second-order frequency term reduces to

    U=−L²/w−3L²M²/(4w³)+5L²M⁴/(8w⁵).

Insert this in the proposal's definition of Q0,K. The logarithms in the
first two integrated terms cancel because M²=2a²H². The remaining exact
expression is

    Q0,K=H²/(32pi²) [8z²/(1+z)−z³/3−z⁵].

It vanishes as K tends to zero and tends to H²/(12pi²) as K tends to
infinity. At finite fixed comoving K it is time dependent. Direct
differentiation gives

    dot(Q0,K)/2+H Q0,K
      =H³/(192pi²) [25z³+6z⁵−15z⁷].

The derivative of z is essential. A fixed physical cutoff would instead
move the comoving integration boundary, and requires the corresponding
boundary-flux terms; that is a different finite-band identity.

## Finite-K contact cancellation in the Ward identity

The density contact in cosmic time is

    C_K=z³(H²d−H dot(d))/(96pi²).

Specializing the proposal's finite trace remainder, using its explicit
J5,J7,J9 polynomials, gives

    delta A_K = [2z³ ddot(d)+6H z⁵ dot(d)
                 +H²d(−27z³−12z⁵+15z⁷)]/(192pi²).

The proposal's finite-band formulas are therefore

    delta rho_K=(H² delta Q_K−H dot(delta Q_K))/2+Q0,K d/2+C_K,

    delta T_K=−2H² delta Q_K−Q0,K d
              +(ddot(delta Q_K)+3H dot(delta Q_K))/2+delta A_K,

    delta p_K=(delta T_K+delta rho_K)/3.

This last equation is used only for the present algebra check, not as an
independent numerical definition of pressure. In
dot(delta rho_K)+H(4 delta rho_K+delta T_K), all derivatives of delta Q_K
cancel. The remaining contact combination obeys

    dot(C_K)+4H C_K+H delta A_K
      =−H³d(25z³+6z⁵−15z⁷)/(192pi²)
      =−[dot(Q0,K)/2+H Q0,K]d.

The complete identity is thus exactly

    dot(delta rho_K)+3H(delta rho_K+delta p_K)=Q0,K dot(d)/2.

It is not Q0 dot(d)/2 at finite K. Replacing either the baseline or the
finite contact/trace terms by their continuum limits inside a purported
exact finite-K comparison breaks the displayed cancellation. This check
supports the proposal's finite-band formulas; it supplies no quadrature or
mode-integration accuracy claim.

## Continuum and state checks

At z=1 the trace contact above becomes

    delta A_r=−H²d/(8pi²)+[ddot(d)+3H dot(d)]/(96pi²),

which agrees with the conformal-time expression in the proposal. The
direct continuum contacts satisfy

    −C_rho+3C_p
      =[d''+2Ld'−12L²d]/(96pi²a²)=delta A_r.

For the separate stationary past-infinite source, dot(d)=0 and delta Q=Q_x d.
The density formula then gives

    rho_x=H²Q_x/2+Q0/2+H²/(96pi²)
         =H²[5/3−2gamma_E−ln2]/(32pi²),
    p_x=−rho_x.

This tests the declared finite contacts and fixed-r differentiation.
It does not identify a compact pulse with an equilibrium susceptibility.

The Ward-plus-trace reconstruction still permits a homogeneous C/a⁴
solution. Unchanged initial state and a source-free initial interval fix
C=0; a reset state would introduce independent data. A passed Ward test
alone would not fix the individual stress contacts. Pressure derivatives,
stress tail bounds, and direct mode/subtraction matching must therefore
remain independent checks in the new frozen numerical experiment.
