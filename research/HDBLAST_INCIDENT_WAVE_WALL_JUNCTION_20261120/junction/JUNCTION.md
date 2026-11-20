# How much radiation reaches a self-gravitating wall, and can it start expansion?

We replace the scalar-supported HDblast geometry by a negative-cosmological-constant spacetime containing a spherical pulse of null dust. We join two retained interiors to a closed wall. The incoming radiation may differ between the two sides. The wall radius follows its Einstein junction equations. The source and perfect absorption of any incident radial stream are supplied as boundary conditions. The radiation runs assume that the deposited energy is isotropic wall matter with pressure $p=\rho/3$. The scalar runs below assign the incident work to a coherent field. Both treatments determine the resulting geometry and energy budgets for their specified absorption laws. A microscopic coupling and a thermal distribution remain open.

Five-dimensional Vaidya brane cosmologies and transport across irradiated branes have a substantial prior literature. Langlois, Sorbo and Rodríguez-Martínez studied a brane *emitting* gravitons. Jennings et al. derived mass and flux equations for unequal sides. Keresztes and Képíró studied a closed brane irradiated from one side, allowing transmission. Pressure imbalance can then compete with self-gravity. Here we ask whether the radius can have a local minimum on the ordinary two-interior branch with fully absorbed incident radiation. We derive the transfer law for equal streams and then prove an acceleration bound for independent streams, including irradiation from just one side. [Emitting brane](https://arxiv.org/abs/hep-th/0206146) · [Asymmetric transport](https://arxiv.org/abs/hep-th/0412281) · [Irradiated closed branes](https://arxiv.org/abs/hep-th/0610112)

## Geometry, source and clock

Use AdS length $\ell$ as the unit of length and time, and $\ell^2$ as the unit of the mass parameter $\mu$. The following dimensionless equations set $\ell=1$. The tension is tuned to $\lambda=6/(\kappa_5^2\ell)$. Write $e=\rho/\lambda$, $H=\dot a/a$, and $v=\dot a$, with dots denoting wall proper time. The bulk cosmological constant is $-6/\ell^2$. On each retained side choose retarded time $u$, increasing toward the future on the wall. The outgoing Vaidya metric and null stress are

$$ ds^2=-f\,du^2-2\,du\,dr+r^2d\Omega_3^2,\quad
f=1+r^2-\frac{\mu(u)}{r^2},\quad
T_{uu}=-\frac{3}{2\kappa_5^2 r^3}\frac{d\mu}{du}.$$

Here the null stress expression uses $\ell=1$. Constant-$u$ null rays travel to *larger* $r$, from the retained bulk toward the wall. Positive energy requires $\mu_u\le0$. The test places the emitting world tube at $r_s=1.6$ and the wall initially at $a=2$. The initial apparent horizons lie below the tube, $f>0$ throughout the calculated source-to-wall annulus, and no packet passes out through a future black-hole horizon. Data on the tube specify a smooth finite packet with finite total energy on each three-sphere. We leave the internal structure and powering of this source open.

Set $X=\sqrt{f_b+v^2}$. Proper-time normalization and the covector normal to the wall, oriented to increasing radius in each retained copy, give $\dot u=1/(X+v)>0$, $n_A=(-v,\dot u)$, $n^u=-\dot u$ and $n^r=X$. The positive energy current into the wall on one side is $-T_{AB}u_{\rm wall}^A n^B=T_{uu}\dot u^2$. For two identical interiors and identical incident streams, the angular extrinsic curvature is $K^i{}_j=(X/a)\delta^i_j=(1+e)\delta^i_j$. These choices fix the sign before squaring the angular equation. The Friedmann and two-sided energy equations follow,

$$ H^2=2e+e^2+\frac{\mu}{a^4}-\frac{1}{a^2},\qquad
\dot e+4He=-\frac{\mu_u\dot u^2}{2a^3}.$$

A decreasing $\mu$ sends positive energy into wall matter. Choosing an increasing mass in this retarded chart would reverse the bulk energy condition. The quantity $\mu/a^4$ in the Friedmann equation is the bulk gravitational mass term and is distinct from the locally deposited radiation $e$.

## A transfer law and a no-bounce result

Let $F=a^4e$ and $B=1+e+H=(a\dot u)^{-1}$. The exact radiation and acceleration equations imply

$$ dF=-\frac{d\mu}{2B},\qquad
\frac{\ddot a}{a}=-2e-3e^2-\frac{\mu}{a^4}+
       \frac{\dot\mu}{2a^4B},\qquad
\dot B=-\frac{1}{a^2}-2(B-1)^2.$$

For example, the acceleration also follows directly from $K^\tau{}_\tau=1-3e$ and $K^\tau{}_\tau=(\ddot a+f_r/2+f_u\dot u^2/2)/X$. Differentiating $B$ and using the two equations above and Friedmann cancels the explicit null flux. These equations extend through $H=0$. As long as $B>0$, $e\ge0$, $\mu\ge0$ and $\mu_u\le0$, the radius is concave in proper time. A wall initially at rest cannot begin expanding under this irradiation. Its turning point at positive radius is a maximum. For unequal bulks the transfer equations change. We treat the acceleration on this branch below. Positive detuning, different wall stresses, reflection, transmission and the full scalar-supported HDblast junction require their own equations.

Since $B$ strictly decreases on a closed trajectory, a smooth packet delivering $\Delta\mu_{\rm lost}>0$ obeys

$$ \frac{1}{B_i}\ \leq\ \frac{2(F_f-F_i)}{\Delta\mu_{\rm lost}}\ \leq\ \frac{1}{B_f} .$$

As packet duration tends to zero at fixed bounded initial state and fixed mass loss, its duration in wall proper time tends to zero too. The equations for $a$ and $B$ have finite right-hand sides throughout the short-pulse limit on this untrapped branch, so their jumps vanish. The limiting positive radiation kick and equal decrease in expansion rate are

$$e_f-e_i=H_i-H_f=\frac{\mu_i-\mu_f}{2a_i^4 B_i}.$$

Thus a sufficiently short incoming pulse reverses a slow initial expansion if $\Delta\mu_{\rm lost}>2a_i^4 B_iH_i$ and enough positive mass is available. This statement takes a limit of regular smooth pulses. It does not assume matching rules for a distributional shell.

## Finite-width pulses

We start with no wall radiation on the two-sided symmetric geometry. Within a pulse let $s=u/\Delta u$ run from zero to one and take $\mu_i-\mu(s)=(\mu_i-\mu_f)(35s^4-84s^5+70s^6-20s^7)$. The endpoint derivatives through third order match the two constant-mass regions. Each run keeps $\mu_f\ge0$.

| Initial $H$ | Retarded width | Lost fraction of initial $\mu$ on each side | $H$ at pulse end | Wall $e$ at pulse end | First maximum radius |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.35355 | 0.20 | 0 | 0.16735 | 0 | 2.44949 |
| 0.35355 | 0.20 | 0.50 | 0.12059 | 0.04652 | 2.33635 |
| 0.35355 | 0.20 | 1.00 | 0.06825 | 0.09861 | 2.24669 |
| 0.04000 | 0.025 | 0.15 | 0.00881 | 0.01817 | after pulse |
| 0.04000 | 0.025 | 0.50 | −0.03385 | 0.06082 | during pulse |
| 0 | 0.05 | 0.50 | −0.08916 | 0.06440 | 2 at launch |

For $H_i=0.04$, the half-mass case turns around 52.17% of the way through the finite retarded-time pulse. The zero-width critical mass fraction is 0.330684. The full smooth calculation finds critical fractions 0.243931 at $\Delta u=0.02$, 0.308828 at 0.005 and 0.329807 at 0.0002. A slower pulse allows the positively curved universe to approach its maximum before the pulse finishes. The threshold quoted here is the mass fraction for $H_f=0$, and it depends on the initial conditions, width and profile. Comparing runs with different mass loss and identical retarded width compares the same source protocol. The elapsed wall times differ slightly.

## What energy was deposited?

The ratio $2\Delta F/\Delta\mu_{\rm lost}$ measures comoving radiation, not absorption probability. All null dust incident on the wall is absorbed in this model. Every ray emitted at fixed $u$ reaches the wall with frequency relative to the wall observer divided by frequency at the stationary emitting surface equal to $\sqrt{f_s}\dot u$. The proper emitter and wall energies on both sides, divided by $\lambda\,{\rm Vol}(S^3_{\rm unit})\ell^3$, satisfy

$$Q_{\rm emitter}=-\int\frac{d\mu}{2\sqrt{f_s}},\quad
Q_{\rm wall}=-\int\frac{\dot u}{2}\,d\mu,\quad
Q_{\rm wall}=\Delta(a^3e)+\int e a^2\,da.$$

The last equality is the fluid first law including the work of expansion. For the faster half-mass run of width 0.2 the proper energy emitted is 1.13285, wall heat is 0.55850, and the internal-energy increase and expansion work are 0.53410 and 0.02440. The wall registers 49.30% of the proper energy registered at the emitting tube, while $2\Delta F/\Delta\mu_{\rm lost}=0.80326$. The source-to-wall redshift varies during the pulse. Neither number identifies a measured photon efficiency or a temperature. The assumed equation of state by itself cannot distinguish photons, another relativistic species or an out-of-equilibrium isotropic distribution.

The primary solver uses $(a,B,\tau)$ in retarded time with $e=[(B-1)^2+a^{-2}-\mu a^{-4}]/(2B)$. A separate integrator uses $(u,a,v,e)$ in proper time with the unsquared temporal Israel condition and an independently accumulated positive wall flux. We also integrate the source-to-wall energy flux by scalar quadrature, check the fluid work against the internal-energy change, reproduce the closed dark-radiation solution without a pulse and solve for the first maximum on the final radiation plateau. The independent endpoint and energy checks differ by less than $1.1\times10^{-11}$ in the stated dimensionless units for the tabulated finite runs. [Calculation](junction.py) · [Numerical entries and controls](JUNCTION.json)

## Can irradiation from just one side start expansion?

Allow independent mass profiles $\mu_i(u_i)$ and retarded clocks on the two sides, $i=1,2$. Use a normal toward larger radius on each side and write $b_i=\sqrt{f_i+\dot a^2}/a>0$. The angular and temporal Israel conditions on the retained interiors are

$$\frac{b_1+b_2}{2}=1+e,\qquad
\frac12\sum_i\frac{\ddot a+f_{i,r}/2+f_{i,u_i}\dot u_i^2/2}{a b_i}=1-3e.$$

Each incident current contributes with the same sign to Codazzi. Squaring the angular condition gives a positive term when the two local bulk mass parameters differ,

$$ \dot e+4He=-\sum_i\frac{\mu_{i,u_i}\dot u_i^2}{4a^3},\qquad
H^2=e(2+e)-\frac1{a^2}+\frac{\mu_1+\mu_2}{2a^4}
 +\left[\frac{\mu_1-\mu_2}{4a^4(1+e)}\right]^2 .$$

The sign of the acceleration follows from the **unsquared temporal equation**. Define $h_m=2/(1/b_1+1/b_2)$, $w_i=(1/b_i)/(1/b_1+1/b_2)$ and

$$ q_i=\frac{\mu_i}{a^4}-\frac{\mu_{i,u_i}\dot u_i^2}{2a^3}\ge0,
\qquad \bar q=\sum_i w_iq_i,\qquad
\frac{\ddot a}{a}=(1-3e)h_m-1-\bar q.$$

The harmonic mean of two positive numbers never exceeds their arithmetic mean, so $0<h_m\le 1+e$. When $0\le e\le 1/3$, this gives $\ddot a/a\le-2e-3e^2-\bar q$. For $e>1/3$, the coefficient of $h_m$ is negative and $\ddot a/a<-1-\bar q$. Equality with zero would require $e=0$ and both masses zero. The spatially closed Friedmann equation then gives $H^2=-1/a^2$, which is impossible at finite real radius. The radius is therefore strictly concave while these assumptions hold. This excludes a local minimum and rules out starting expansion from rest, even for one-sided irradiation. It also applies to initially unequal nonnegative bulk masses. A relative change of expansion between histories and the sign of the acceleration within either history answer different questions.

We tested this result by starting with the same mass on each side and varying how the two profiles subsequently change. The examples below specify the loss on each side as a fraction of that side's initial mass. The calculation runs until the final supplied stream ends in its own retarded time. The end proper times can differ between examples.

| Initial $H$ and $e$ | Width | Loss left / right | End wall time | End $H$ | End $e$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.35355 and 0 | 0.20 | 50% / 50% | 0.53629 | 0.12059 | 0.04652 |
| 0.35355 and 0 | 0.20 | 75% / 25% | 0.54295 | 0.11941 | 0.04605 |
| 0.35355 and 0 | 0.20 | 100% / 0% | 0.54977 | 0.11884 | 0.04495 |
| 0.04000 and 0 | 0.025 | 100% / 0% | 0.05323 | −0.03093 | 0.05907 |
| 0 and 0 | 0.05 | 50% / 0% | 0.10021 | −0.05625 | 0.03151 |
| 1 and 0.4 | 0.05 | 50% / 0% | 0.22751 | 0.61969 | 0.20507 |

The first three runs lose the same total bulk mass parameter. At the common wall time $\tau=0.600$, after all three pulses have ended, the equal split gives $a=2.272278$ and $H=0.106258$. The one-sided pulse gives $a=2.273805$ and $H=0.107320$. Both are decelerating, with $\ddot a/a=-0.209053$ and $-0.214088$, respectively. The pulse-free history from the same initial state gives $a=2.310207$ and $H=0.152559$ at this time. In this example the one-sided history has a slightly greater expansion rate than the equal split at a common time, and both have lower expansion rates than the pulse-free history.

In the third run the asymmetry contributes $+0.003018$ to $H^2$ at the end of the pulse. The fourth run reverses its expansion at wall time 0.02769. The last run crosses $e=1/3$ and tests both parts of the bound. Exchanging left and right gives the same wall dynamics. A proper-time integrator evolves the temporal Israel equation through reversal and checks the angular equation, Friedmann equation, integrated comoving flux and acceleration bound. For cases expanding through the end of the pulse, an independent integrator uses the expanding Friedmann branch and Codazzi directly. The largest absolute discrepancy across these checks is below $1.8\times10^{-11}$. [Asymmetric calculation](asymmetric.py) · [Numerical entries and controls](JUNCTION.json)

An independent generalization derives the acceleration equation using spherical curvature, bulk mass and the pressure normal to the wall. It also gives an exactly solvable positive-vacuum-energy wall for which the closed geometry does bounce. [Normal pressure and a comparison with vacuum energy](NORMAL_PRESSURE.md)

## An independent geometric check on scalar absorption

The companion scalar example assigns incident energy to a coherent scalar field on the wall. During its short interaction we reconstruct the two identical bulk mass profiles $\mu(u)$ from the scalar work and the retarded clocks. We then fix that mass profile and re-integrate the motion with the **unsquared temporal Israel condition**, the canonical scalar equation, and the lapse of the outgoing metric. In this second integration the bulk stress determines the incident work independently of the scalar force. The prescribed input force and plateau potential remain common to both calculations.

On the plateau, $P=a^3\dot\psi$ is constant and the squared speed is

$$G(a)=\dot a^2=V(2+V)a^2-1+\frac{\mu}{a^2}
       +\frac{(1+V)P^2}{a^4}+\frac{P^4}{4a^{10}}.$$

For positive $V$ and nonnegative $\mu$, this function is strictly convex for $a>0$. In the example its two positive roots are 1.260483444 and 1.899038171. The smaller lies below the specified emitting surface. The contracting wall arrives on the outer allowed branch and reverses at radius 1.899038171 and proper time 0.794581514. With $V=0$ and otherwise nonnegative mass the derivative of $G$ is nonpositive, which excludes such an outer turning point in the kinetic-only case.

The stationary tube emits 1.147261 units of local energy on the two sides. The wall receives 0.940845 units as scalar work. Pressure work during absorption is −0.001215. The independently measured wall flux, fluid first law, temporal evolution and spatial junction agree to better than $5\times10^{-9}$ in these units. This replay tests the geometry of the reconstructed pulse. It leaves the conversion of incoming radiation into a coherent field as an assumed coupling. [Independent replay](scalar_validation.py) · [Scalar model and controls](SCALAR_CAPTURE.md) · [Numerical results](JUNCTION.json)

A separate forward calculation fixes the emitted mass profile in retarded time before solving the wall trajectory. On the chosen scalar potential the 94% mass-loss pulse of retarded width 0.005 reaches the plateau and produces a bounce at radius 1.899614. It subsequently reaches an expanding interval in which the radiation fluid exceeds 90% of wall matter. [Forward pulse and controls](PRESCRIBED_SOURCE.md)

The force-driven scalar model also admits a finite downhill exit. Its nominal run expands for 0.319917 e-folds above the 90% wall-material threshold before turning around again. At the first crossing, the linear radiation term supplies 86.04% of the positive Friedmann terms, and positive curvature offsets 51.13% of their sum. [Internal transfer into a radiation fluid](RADIATION_EXIT.md)
