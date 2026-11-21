# Relativistic drag for two equal counterstreams

Consider two dilute populations of particles of mass $m$ on a locally homogeneous patch of wall. Give each population the same proper number density $n$ and rest temperature $T$. In the wall frame their mean velocities are $+v$ and $-v$ along $z$. The calculation assumes isotropic Maxwell–Jüttner distributions in the two rest frames, binary elastic scattering, negligible changes of the background fields during an encounter, and an angular-even differential cross section in the center-of-momentum frame. Units are $c=\hbar=k_B=1$.

## The momentum transferred in one encounter

Write $P=p_1+p_2$, $D=p_1-p_2$, $s=P^2$, $P\cdot D=0$, $\beta=(1-4m^2/s)^{1/2}$. Each outgoing particle has mean momentum $P/2$ when the angular distribution is even. The mean $z$ momentum lost by stream 1 is $D_z/2$. For the opposite stream the sign reverses. The two proper distributions are

$$
f_\pm(p) = A \exp[-\gamma(E\mp vp_z)/T],\qquad
A=\frac{n}{4\pi m^2 T K_2(m/T)}.
$$

With $p_1$ drawn from the positive stream, the wall-frame momentum loss per unit volume and unit time is

$$
F_z = \int\frac{d^3p_1}{E_1}\frac{d^3p_2}{E_2}
 f_+(p_1)f_-(p_2)\sqrt{(p_1\cdot p_2)^2-m^4}\,
 \sigma_T(s)\frac{D_z}{2}.
$$

For the selected even elastic channel $\sigma_T=\sigma_{\rm tot}$. A spin-independent event-counting rate generally has different momentum and flux weights. Energy transfer vanishes *on average* in the wall frame by symmetry of the beams. Directed motion can still become random internal energy.

## Reduction to one invariant-mass integral

The invariant two-body measure at fixed $P$ is $\beta\,d\Omega^*/2$, and the invariant flux is $s\beta/2$. In the center-of-momentum frame orient the polar axis along the projection of the wall's $z$ axis into $P^\perp$. Its cosine $u$ is uniform and $D_z=\beta\sqrt{s+P_z^2}\,u$. The angular average of $u e^{xu}$ is $S'(x)$ with $S(x)=\sinh(x)/x$. Integration over transverse $P$ gives

$$
\int \frac{d^2P_\perp}{2P^0}e^{-\gamma P^0/T}
=\frac{\pi T}{\gamma}\exp[-\gamma\sqrt{s+P_z^2}/T].
$$

Use $P_z=\sqrt{s}\sinh y$, $E=\sqrt{s}\cosh y$, $a=\gamma/T$, $b=\gamma v\beta/T$. The force becomes

$$
F_z = \frac{A^2\pi^2T}{\gamma}
 \int_{4m^2}^{\infty} ds\;s\beta^3\sigma_T(s)\, I(a,b,s),\qquad
I=\int_0^\infty dy\, E^2 e^{-aE}S'(bE).
$$

The rapidity integral can be done in modified Bessel functions. With $x_\pm=(a\pm b)\sqrt{s}$,

$$
I=\frac{\sqrt{s}}{2b}[K_1(x_-)+K_1(x_+)]
  -\frac{K_0(x_-)-K_0(x_+)}{2b^2}.
$$

For numerical evaluation measure energies in units of $m$. Put $\Theta=T/m$, $r=\mathrm{atanh}|v|$, and $\sqrt{s}/m=2\cosh\eta$. Then $x_\pm=2\cosh(\eta\pm r)/\Theta$. The closure momentum density of stream 1 is $n m h\gamma^2 v$. This yields

$$
C(\Theta,r)=\frac{m^2\nu}{n}
=\frac{2}{\Theta\gamma^3 h K_2(1/\Theta)^2}
 \int_0^\infty d\eta\,\sinh^4\eta\,(m^2\sigma_T)\,\frac{I}{v}.
$$

Here $I$ is dimensionless and $\nu=F_z/(nmh\gamma^2 v)$. The live fluid uses the Taub–Mathews specific enthalpy $h=5\Theta/2+\sqrt{1+9\Theta^2/4}$. The collision integral itself uses the exact Maxwell–Jüttner distribution. Replacing $h$ by $K_3(1/\Theta)/K_2(1/\Theta)$ gives the rate per unit exact Maxwell–Jüttner momentum density. Using the solver’s $h$ makes the drag source equal the integrated particle force for the specified $n,T,v$. Each of the six sign-and-axis populations has proper number density $nL^3=S_0 D_c/(6 a^3\mu\gamma)$, with $\mu=mL$ and $D_c/a^3$ the total rest-mass density in the solver's tension units. Thus $\nu L=(nL^3)C/\mu^2$.

## Vector pole and controls

Use one unpolarized spin-one elastic channel with $M/m=8300$ and $\epsilon=g^2/(12\pi)=0.30$. In units $m=1$ its isolated massive running-width cross section is

$$
\widehat\sigma(\eta)
=\frac{3\pi}{\cosh^2\eta}
\frac{[\epsilon(s+2)]^2}
 {(s-M^2)^2+[\epsilon(s+2)\tanh\eta]^2},\qquad
 s=4\cosh^2\eta.
$$

This formula chooses a single elastic decay channel and the angular-even squared pole amplitude. Other particle species, crossed amplitudes and medium corrections could change the drag. The calculation also leaves the microscopic origin of the assumed local rest distributions open.

At $\Theta=0$ each stream is monoenergetic. The exact result is $C=2\sinh r\,\widehat\sigma(r)$, including the vanishing $r=0$ limit. At $r=0$, $\Theta\ll1$, $C\sim16\widehat\sigma(0)\sqrt\Theta/(3\sqrt\pi)$. For $\Theta\gg1$, $r=0$, a constant cross section gives $C\to 5\widehat\sigma/4$. To see the last coefficient, independent massless bath momenta have $\langle p\rangle=3T$, $\langle p^2\rangle=12T^2$ and angular average $\langle(1-\cos\psi)(p_{1z}-p_{2z})^2\rangle=10T^2$. A fully saturated single spinless s-wave, $\widehat\sigma=4\pi/\sinh^2\eta$, gives $C\to\pi/\Theta^2$.

## Numerical implementation

The quadrature splits in $\eta$ at the thermal peak and around the pole. The cutoff obeys $4\sinh^2[(\eta-r)/2]/\Theta=60$ on either accessible tail. This controls the extracted exponential. The reported numerical checks below measure the resulting integral. The integrator combines scaled Bessel functions with $\exp[-4\sinh^2((\eta\pm r)/2)/\Theta]$ to retain cold, rapidly moving beams. Around $v=0$, the odd-power series of $S'$ evaluates $I/v$ directly, also at zero speed.

The lookup stores $\ln C$ on 89 uniformly spaced rows in $\log_{10}\Theta\in[-5,6]$ and 401 columns in $y=\mathrm{asinh}(|\gamma v|/\sqrt{\Theta/(1+\Theta)})\in[0,20]$. Four-point polynomial interpolation evaluates four $y$ rows at the requested physical momentum before interpolating in temperature. Positive out-of-range high temperature or rapidity raises an error. Below the grid the reader uses $C\simeq2\sqrt\Theta\,\widehat\sigma(r) g(z)$, $z=|\sinh r|/\sqrt\Theta$, where

$$
g(z)=\left(z+z^{-1}-\frac{1}{4z^3}\right)\mathrm{erf}(z)
 +\left(1+\frac{1}{2z^2}\right)\frac{e^{-z^2}}{\sqrt\pi}.
$$

The nonrelativistic relative-velocity Gaussian has component variance $2\Theta$ and mean $2v$ along $z$. Its momentum-transfer moment gives this function. Using $\sinh r$ and $\widehat\sigma(r)$ extends the small-$\Theta$ expression to its exact cold-beam limit at each finite $r$. The finite-$\Theta$ extension is approximate.

The 636 off-grid, limit and resonance-focused probes in [BOOST_CHECK.json](BOOST_CHECK.json) compare the production reader against direct quadrature. Doubling the quadrature order changed the result by at most $1.4\times10^{-11}$ relative on these probes. The largest sampled lookup difference was 0.83%, the 95th percentile was 0.068%, and the largest sampled cold continuation difference was 0.030%. The independent calculation in [sampling.py](sampling.py) draws particles directly in the two rest frames and boosts both four-momenta. Its scrambled Sobol estimates and estimated sampling errors are recorded in [SAMPLING.json](SAMPLING.json). These tests concern the numerical evaluation of this specified effective collision model.

At the broadest discrepancy in the first six sampling probes ($\Theta=10,r=4$), twelve fresh scrambles of $2^{21}$ pairs brought the direct estimate within 0.0093% of the one-dimensional kernel. The between-scramble standard error was 0.0099%. [Sampling convergence](../scattering/CONVERGENCE.json)

Halving both table spacings and increasing the integration order from 16 to 20 changed the hot-interval radius factor by at most 0.00051% across four complete wall runs. The largest sampled interpolation discrepancy for the denser table was 0.054%. [Whole-history refinement](REFINEMENT.json)
