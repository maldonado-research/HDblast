# What the spectral-work identity adds

For a real smooth compact canonical mass pulse, let
`s=a^2 delta x`, `q_M=a^2 delta Q`, and `M=a sqrt(r)`. The matched response
can be split using any constant auxiliary conformal scale mu:

    q_M=q_mu-ln(M/mu)s/(8 pi^2).

The arbitrary mu cancels from the full response. This is a way to display
the local contact already fixed by the inherited subtraction, not an added
physical parameter or a changed prescription.

The canonical plane-wave modes after the pulse have a positive quadratic
excitation-energy functional

    E_can^(2)=1/(32 pi^2) integral_0^infinity omega |S(omega)|^2 domega,
    S(omega)=integral s(eta)exp(-i omega eta)deta.

Its equality to source work requires the matched background contact:

    E_can^(2)=(1/2) integral s' q_M deta
               -1/(32 pi^2) integral (a'/a)s^2 deta.

Here a'/a=-1/eta>0. A constant contact makes no net compact-pulse work because
the integral of s's vanishes. The inherited contact depends on time through
M, so integration by parts leaves the displayed nonzero drift term. Omitting
it would misidentify matched source work as canonical excitation energy.
The subtraction is fixed analytically; it cannot be fitted after a work
comparison. It is independent of the constant r even though q_M itself
depends on r.

There are two independent analytic checks. The retarded Fourier multiplier
has imaginary part `sgn(omega)/(16 pi)` for the stated Fourier convention,
which gives the positive spectral functional. The exact finite-cutoff
subtraction `A_K=asinh(K/M)-K/sqrt(K^2+M^2)` obeys

    A_K'=-(M'/M)K^3/(K^2+M^2)^(3/2),

giving the same contact sign and coefficient before taking the cutoff limit.
These identities would support a useful separately registered numerical
work/spectral comparison. No such energy integral was evaluated in the
completed scalar-response run.

Canonical energy also differs from physical minimally coupled stress.
After the pulse, state differences obey

    a^4 Delta rho=Delta E_can-(h/2)(Delta q)' +(3h^2/2)Delta q,
    h=a'/a,  Delta q=Delta<v^2>.

Coherence and polarization appear in the last two terms and can be linear
in the pulse; canonical occupation energy starts at quadratic order. Their
presence explains why spectral positivity does not establish a positive
physical energy change, radiation, or heating. The physical Ward identity
also contains expansion work and needs the matched rho and p responses.
