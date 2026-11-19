# Can a change of particle basis supply a passive decay law?

Changing which oscillator defines an empty state also changes what an amplitude-damping bath tries to prepare. In the HDblast shell model, these two choices matter separately: which particles are counted, and which Hamiltonian measures their energy. The calculation below gives the energy that a bath or external drive must supply when those choices disagree. We work at one time with one real oscillator, use $\hbar=1$, and suppress the wavevector. The field Hamiltonian is positive with $E>0$. The same calculation applies separately to independent modes.

### Two definitions of power

Let $A$ annihilate the physical ground of $H=E(A^\dagger A+1/2)$. Suppose a zero-temperature Markovian channel damps a different normalized oscillator $B$ at rate $\Gamma\geq 0$. Write $A=\alpha B+\beta B^\dagger$, $|\alpha|^2-|\beta|^2=1$, and $N=|\beta|^2$. Let $n=\langle B^\dagger B\rangle$ and $z=\langle B^2\rangle$. The oscillator attenuation equation returns the covariance to the $B$ ground. At a fixed instant the power **removed from the physical system Hamiltonian** is

$$
P_H=-\operatorname{Tr}(H\,{\cal D}\rho)
=\Gamma E(\langle A^\dagger A\rangle-N)
=\Gamma E\big[(1+2N)n+2\operatorname{Re}(\alpha\beta^*z)\big].
$$

The expectation of the changing Hamiltonian has the separate term $\langle\partial_s H\rangle$. A state-independent subtraction in $H$ drops out of $P_H$, since the channel preserves trace. If each counted quantum is instead deposited into massless daughters with energy $E$, their assigned power is $P_B=\Gamma E n$. A model with both this deposit and the specified channel must account for a further external power $P_B-P_H$. In a time-dependent shell, attributing it to the shell requires a corresponding term in its dynamical energy equation.

Because $|z|^2\leq n(n+1)$, the sharp interval for $P_H/(\Gamma E)$ at fixed $n$ is

$$
(1+2N)n\ \pm\ 2\sqrt{N(N+1)n(n+1)} .
$$

Pure squeezed vacua with opposite phases reach the endpoints. For $n>0$ the entire interval is nonnegative precisely when $n\geq4N(N+1)$. At $n=0$ both powers vanish. At the physical $A$ vacuum, $n=N$ and $P_H=-\Gamma EN$ while $P_B=+\Gamma EN$. If $N>0$, using $P_H$ as the sole source for an initially empty daughter fluid would drive its density below zero over an infinitesimal step with nonzero rate. This is a diagnostic of the proposed accounting.

### A condition covering multiple linear channels

At fixed positive $E$, take any collection of completely positive Markovian Lindblad terms with nonnegative rates $\gamma_j$ and jump operators $L_j=u_j A+v_j A^\dagger$. For the physical ground $|0_A\rangle$,

$$
\left.\frac{d\langle H\rangle}{ds}\right|_{\rm diss}
=\sum_j\gamma_j\langle 0_A|L_j^\dagger(H-E/2)L_j|0_A\rangle
=E\sum_j\gamma_j|v_j|^2 .
$$

A Hamiltonian commutator makes zero contribution in an energy eigenstate. Thus every active jump must have $v_j=0$ for these linear dissipators to remove nonnegative physical power on **all** states. Conversely, jumps proportional to $A$ give $P_H=E\sum_j\gamma_j|u_j|^2\langle A^\dagger A\rangle\geq0$ at that instant. This all-states energy condition concerns a bath acting with fixed couplings. Equilibrium passivity under arbitrary cyclic unitaries is a distinct definition. The test leaves open microscopic interactions with memory, physical interaction energy, driving, and restrictions to a particular class of prepared states.

### Applying the condition to an expanding shell

In dimensionless shell time $s=H_0\tau$ write $h=d\ln a/ds$, $E^2=k^2/a^2+m^2$, $\Omega=aE$, and conformal time $d\eta=ds/a$. For the rescaled coordinate $u=a\chi$ choose the leading adiabatic mode data $u'/u=-i\Omega-\Omega'/(2\Omega)$, $|u|^2=1/(2\Omega)$, where prime means $d/d\eta$. This defines the proposed instantaneous Gaussian target even when an adiabatic approximation to the true evolution fails. With $\Pi=a^3 d\chi/ds$ set

$$
\ell=h+\tfrac12 d\ln\Omega/ds=h+C/E^2,\qquad
C=\tfrac14[d(m^2)/ds+2hm^2] .
$$

The target has $\operatorname{Var}\chi=1/(2a^3E)$, $\operatorname{Cov}(\chi,\Pi)=-\ell/(2E)$, and $\operatorname{Var}\Pi=a^3(E^2+\ell^2)/(2E)$. Its number in the physical basis is $N=\ell^2/(4E^2)$. The physical ground has the same number in the target basis. At finite mass and fixed smooth background $N=(ah)^2/(4k^2)+O(k^{-4})$. A positive $B$-count deposit applied to the physical ground therefore has the opposite sign to the signed physical power, mode by mode.

For the published rate $\Gamma=y^2m^2/(8\pi E)$, real-field mode measure and cutoff $K$, the hypothetical deposit per unit physical volume at such a ground-state slice is exactly

$$
Q_B=-Q_H={y^2m^2\over64\pi^3}
\int_0^{K/a}\frac{p^2}{p^2+m^2}\left(h+\frac{C}{p^2+m^2}\right)^2\,dp .
$$

The difference between the two rates is $2Q_B$. In the formal large-$K$ limit on a prescribed history the coefficient of $K$ in $Q_B$ is $y^2m^2h^2/(64\pi^3a)$. This is a test of a specific state and specifies no physical continuation of the shell theory beyond its frequency screen.

### A finite unitary realization

A vacuum auxiliary oscillator $D$ with independently specified Hamiltonian $E(D^\dagger D+1/2)$ makes the missing work measurable. Switch a beam-splitter interaction on and off so that $B_{\rm out}=\sqrt{\eta}B+\sqrt{1-\eta}D$, $D_{\rm out}=-\sqrt{1-\eta}B+\sqrt{\eta}D$. Free Hamiltonians at the endpoints are equal. Tracing over the auxiliary system gives exactly the attenuation channel. In the physical system vacuum,

$$
\Delta H_{\rm sys}=(1-\eta)EN,\qquad
\Delta H_{\rm aux}=(1-\eta)EN,\qquad
W_{\rm ext}=2(1-\eta)EN .
$$

The explicit map is one possible driven realization. Its bare energy cost does not identify a physical driver for the cosmological shell. An independent two-mode calculation exponentiates the interaction in the physical number basis and tests both equal-count input phases and the target ground. [Collision experiment](collision/COLLISION.md)
