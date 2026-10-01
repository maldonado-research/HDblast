# Action-derived stress, current and state differences

Analytical formulation, 1 October 2026. These are standard field-theory identities applied to the declared HDBLAST spectator extension. External mathematical novelty and a completed renormalized source are not claimed.

## One action gives both stress and current

For S_chi=-1/2 integral sqrt(-h)[(partial chi)²+M²(phi)chi²] on a spatially flat FRW shell, M²=m0²+g²(phi-phi_star)² and signature (-+++), normalized modes f=chi obey

fddot+3H fdot+(k²/a²+M²)f=0,
a³(f fdot*−fdot f*)=i.

With P=k/a, omega²=P²+M², the mode contributions are

e=(|fdot|²+omega²|f|²)/2,
p=(|fdot|²−(P²/3+M²)|f|²)/2,
j=(M²_phi/2)|f|².

Integrate each against k²dk/(2pi²). Then edot+3H(e+p)=phidot*j, or covariantly nabla_mu T^mu_nu=−J partial_nu phi. The scalar/geometry sector must carry the opposite work.

The covariance equations qdot=2r, rdot=s−3Hr−omega²q, sdot=−6Hs−2omega²r (q=|f|², r=Re(f*fdot), s=|fdot|²) independently give this identity.

## The numerical work ledger

A fixed comoving band and fixed quadrature weights preserve the same identity. Check it using

D(t)=[a³rho]_initial^t + integral 3Ha³p dt − integral a³phidot J dt.

This must be integrated independently alongside the modes or from adequately resolved outputs. Merely substituting the mode equation into an analytical rhodot expression is not a numerical evolution check. Normalize with absolute energy/work contributions and a declared floor.

For moving endpoints k_L,U(t), add boundary flux B=[k² kdot e]lower^upper/(2pi²). A fixed physical momentum cutoff k=a Lambda has B=H[k³e]lower^upper/(2pi²). An energy cutoff omega<Lambda instead uses K=a sqrt(Lambda²−M²), with Kdot/K=H+(Lambda Lambdadot−M²dot/2)/(Lambda²−M²). Spectral-boundary flux must not be attributed to scalar particle production.

An a⁴rho ledger has the additional trace term:
d(a⁴rho)/dt=a⁴[H(rho−3p)+phidot J].
Dropping it assumes radiation behavior before it is established.

## Clock mapping and physical energy

X=a^(3/2)f, u=af and d eta=dt/a give

Xddot+[omega²−(3/2)Hdot−(9/4)H²]X=0,
u''+[k²+a²M²−a''/a]u=0.

X=sqrt(a)u and Xdot=(H/2)X+u'/sqrt(a) map the state and derivative. Reconstruct physical energy from

e=[|Xdot−3HX/2|²+omega²|X|²]/(2a³)
 =[|u'−aHu|²+(k²+a²M²)|u|²]/(2a⁴).

Choosing independent instantaneous vacua in the two clocks changes the state. A negative canonical X-frequency squared does not alone establish a physical tachyon or a meaningful particle number.

## Occupation and coherent contributions

At one time define the physical Hamiltonian phase-space basis b=1/sqrt(2a³omega), bdot=−iomega b. For f=alpha b+beta b*, let n=|beta|² and C=Re(alpha beta*). This is an instantaneous phase-space basis, not one exact reference solution throughout time.

Then
e=omega(n+1/2)/a³,
p=P²(n+1/2)/(3a³omega)−(2P²/3+M²)C/(a³omega),
j=M²_phi(n+1/2+C)/(2a³omega).

Thus physical energy can exactly equal an occupation integral in this particular basis, while pressure and force still require the coherent terms. The exact identity ndot=(3H+omegadot/omega)C makes an occupation-only gas residual equal to integral k² omega ndot dk/(2pi²a³). This is a measurable failure of an uncompleted gas closure during production, not extra energy supplied for free.

When occupation instead diagonalizes the canonical X oscillator, physical energy also includes geometry/coherence corrections. The choice of particle basis must always be stated.

## Exact in-minus-out pair

Two exact modes A and B on the same metric and mass history have constant Bogoliubov coefficients

alpha=i a³(B* Adot−Bdot* A),
beta=−i a³(B Adot−Bdot A).

Their drift is a numerical error diagnostic. The differences Delta rho, Delta p and Delta J obey the same exchange law and comoving work ledger. At the declared final physical-Hamiltonian out data, Delta rho_f=integral k² omega_f |beta|² dk/(2pi²a_f³). This is a specified finite-time excitation comparison, not an invariant particle population independent of the chosen state and endpoint.

For two Hadamard states of the same operator, state-independent local counterterms cancel in the difference. On a sufficiently smooth background, choosing identical Hadamard modes outside a finite band and regular exact modes inside it gives a well-defined finite relative observable. Compact momentum modifications have smooth two-point differences. The finite differentiability of our numerical spline means this is a conditional continuum construction, not a demonstrated Hadamard completion of the executed surrogate.

Changing the band, initial state, initial slice or out slice changes the state pair. Those tests measure sensitivity; they are not automatically regulator removal for one fixed state. Subtracting g=0 or a different background is a comparison between different operators and lacks this cancellation argument.

Delta J is a difference of forces. It cannot replace absolute J_A in the shell junction without accounting for the reference vacuum polarization/current and common finite matching terms.

## Why conservation does not certify renormalization

The instantaneous zero-point quantities e0=omega/(2a³), p0=P²/(6a³omega), j0=M²_phi/(4a³omega) themselves obey the finite-band exchange identity algebraically. Their subtraction can therefore pass a Ward check without being one evolved reference state or a complete covariant renormalization. Higher ultraviolet divergences and finite matching terms remain.

Fourth adiabatic order is relevant to stress subtraction in four dimensions, while field-square/current has lower superficial divergence. One must derive the common prescription, ordering of the mass background and finite counterterms. Finite adiabatic order is not automatically the Hadamard condition.

## Conformal null control

At M=0, xi=1/6, u=exp(-ik eta)/sqrt(2k) has zero conformal-basis production on flat FRW. If stress is checked, include the nonminimal terms
e_xi=xi(3H²q+3H qdot),
p_xi=xi[−(2Hdot+3H²)q−qddot−2Hqdot].
Then e=k/(2a⁴), p=e/3 mode by mode. This does not establish zero absolute renormalized stress; the anomaly and local terms remain.

The primary numerical run uses xi=0. The conformal null test checks geometry and state mapping for a different coupling.
