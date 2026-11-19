# What the first elastic collision does to the particle spectrum

For the whole emitted population on the quiet wall, the initial second energy moment is **1.10429** in units of the mean energy squared. A classical thermal gas of massless particles with the same mean energy has a second moment of **1.33333**. In the stipulated contact interaction, collisions initially broaden the spectrum toward that value on all four walls. The amount of broadening also depends on whether collisions can occur during the available time. [First-collision protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7016401650) · [All values and numerical controls](SHAPE.json) · [Calculation](shape.py)

## One collision in the expanding gas

Let $M_j=\int k^j dN/y^2$ denote one-sign comoving moments, with $k$ measured in units of the reference energy $H_0$. Both signs have the same distribution. We use the isotropic, homogeneous, dilute Boltzmann closure and independent incoming particles. Only opposite signs scatter, with $d\sigma/d\Omega_*=\kappa\hat s/(4\pi)$ in their center-of-mass frame, where $\hat s$ is the physical two-particle invariant mass squared. The outgoing participants conserve energy and momentum in the local freely falling frame.

An event samples the two incident comoving momenta $x,y$ independently with probability measure $k\,dN/(y^2 M_1)$. For incident opening cosine $u$, the relative flux times the cross section is proportional to $x y(1-u)^2$. The event-weighted probability density of $u$ is $3(1-u)^2/8$, with mean $-1/2$. The outgoing energies expressed as comoving momenta are

$$k_\pm=\frac{x+y\pm z\sqrt{x^2+y^2+2xyu}}{2},$$

where $z$ is uniform on $[-1,1]$. This variable is the direction cosine in the collision center-of-mass frame relative to its boost. Expansion redshifts all four physical momenta by the same local scale factor at the collision event.

Let $A=M_2/M_1$, $B=M_3/M_1$, $C=M_4/M_1$, and $\bar k=M_1/M_0$. Integrating $(k_+^m+k_-^m-x^m-y^m)/2$ for $m=2,3$ gives the change per participating particle, averaged over collisions.

$$G_2=\frac{5A^2-4B}{12},\qquad G_3=\frac{3AB-2C}{4}.$$

Number and energy give $G_0=G_1=0$. The initial slope of the dimensionless moment $R_m=\langle k^m\rangle/\bar k^m$ per mean particle collision depth is $G_m/\bar k^m$. Thus the initial second moment rises precisely when $B/A^2-1<1/4$. This is the relative variance of the energies of particles *entering collisions*. The incident distribution of a classical massless thermal gas gives equality. A monoenergetic distribution at $k=1$ gives $G_2=1/12$ and $G_3=1/4$. The replay tests both limits.

| Prescribed wall, full spectrum | Initial $R_2$ | Initial $R_3$ | Initial slope of $R_2$ per mean collision | Initial slope of $R_3$ per mean collision |
|---|---:|---:|---:|---:|
| Quiet | 1.104287 | 1.321595 | +0.067573 | +0.247649 |
| Plus | 1.104311 | 1.321527 | +0.067617 | +0.247851 |
| Minus | 1.104028 | 1.321338 | +0.067420 | +0.246748 |
| Far | 1.104223 | 1.321363 | +0.067591 | +0.247668 |
| Classical thermal reference | 4/3 | 20/9 | 0 | 0 |

On the quiet wall through shell time 2.16, the mean collision depth is $25.4908\Lambda$, with $\Lambda=y^2\kappa H_0^4$. The first-order changes are $\Delta R_2=1.72248\Lambda$ and $\Delta R_3=6.31278\Lambda$, each with corrections of order $\Lambda^2$ for this fixed history. These derivatives apply at small depth. The later distribution and the collision rate change together once rescattering matters. The conditional fraction of incoming collisions between two particles from the lower source-frequency band is 91.89%. Mixed low-high collisions account for 7.94%, and high-high collisions for 0.172%. Band labels refer to the original *pair* frequencies. Scattering redistributes individual momenta. [Computed optical depths](../transport/CONTACT.json)

The angular control rebuilds the one-particle distribution from the original positive pair weights and four direction nodes for the single-daughter projection. It constructs a positive three-node Gaussian rule from those particle momenta. The rule preserves all collision-biased raw energy moments through degree five. Independent integrals over the incident and outgoing angles evaluate $k_+^m+k_-^m$ directly. For $m$ through three, the sum is a polynomial in both incoming energies after averaging over the outgoing direction. The replay checks the slopes, the separate low-low, low-high and high-high contributions, and the vanishing changes in number and energy.

## A particle count for chemical equilibrium

As a secondary interpretation of the reconstructed spectrum, consider an isolated, charge-symmetric gas of one massless Dirac flavor. Each sign has $g=2$ spin states. Its emitted proper number and energy densities are $H_0^3y^2 M_0/a^3$ and $H_0^4 y^2 M_1/a^4$. With pair-changing reactions in chemical equilibrium and zero charge asymmetry, the Fermi–Dirac comparison has zero chemical potential. Its per-sign equilibrium number and energy densities are

$$n_{\mu=0}=\frac{3g\zeta(3)}{4\pi^2}T^3,\qquad
\rho_{\mu=0}=\frac{7 g\pi^2}{240}T^4.$$

Equating its energy density to the emitted energy fixes $T$. The common redshift cancels in the resulting number ratio. On the quiet wall,

$$\frac{n_{\mu=0}(\rho_{\rm emitted})}{n_{\rm emitted}} =
 \frac{3g\zeta(3)}{4\pi^2 M_0}
 \left(\frac{240 M_1}{7 g\pi^2}\right)^{3/4}\frac{1}{\sqrt{|y|}}
 =\frac{2.58595}{\sqrt{|y|}}.$$

At the illustrative $y=0.10$ the zero-chemical-potential reference needs **8.18 times as many particles of each sign** as the modeled emission supplies. The four walls give ratios of 8.14–8.19. Elastic scattering preserves both particle counts. A number-preserving Fermi–Dirac distribution with the same total energy could instead have nonzero fugacity. At this same coupling the two conserved quiet-wall moments fix fugacity $0.000173936$ and $\mu/T=-8.65682$. Its $R_2=1.333326$ and $R_3=2.222192$ agree closely with the dilute classical values in the table.

These are equilibrium *comparisons* using only the particles emitted within the modeled cutoff. Additional species, an ambient plasma, energy exchange, and number-changing processes can alter the accounting. The leading $y^2$ emission approximation restricts interpretation of the coupling scaling to its perturbative domain. Whether collisions reach a thermal distribution at fixed number requires a kinetic evolution.

## Reproduction

From the study root, regenerate the daughter data before running the shape calculation.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB shape/shape.py --repo ../.. --output /tmp/hdblast-shape.json
~~~
