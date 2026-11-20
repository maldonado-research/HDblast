# An incoming Maxwell pulse on a spherical wall

A pulse can transfer most of its energy to a scalar on a spherical wall. It can also set the scalar moving in the opposite direction from the planar prediction. In the model studied here the signed impulse depends on the wall radius, the travel distance from the emitter, and the timing of the waves. A large scalar kinetic energy alone does not establish that the field will climb a potential at positive field values.

## The two bulk regions and the wall

Take two copies of the five-dimensional Minkowski annulus $r_s\le r\le R$ with metric $-dt^2+dr^2+r^2d\Omega_3^2$. Join their outer boundaries to the same fixed round wall. On each side take three independent Maxwell fields $A^{(I)}=A(t,r)\sigma_I$, where $\sigma_I$ are orthonormal invariant one-forms on the unit three-sphere with $d\sigma_I=\pm2\,\star_3\sigma_I$. Equal amplitudes in the triplet make the tangential stress isotropic. Their magnetic eigenvalue is four. In this probe problem the metric and radius are held fixed; the action has a canonical wall scalar, no scalar potential, and gauge coefficient $K(\psi)>0$ on the wall.

Put $q=\sqrt3\, A/R$ and $w=r/R$. Dividing by the proper wall three-volume $2\pi^2R^3$ gives the two bulk energy densities and the surface energy
$$E_{\rm bulk}^{(i)}=\tfrac12\int_{r_s}^R w\left(q_t^2+q_r^2+\frac{4q^2}{r^2}\right)dr,\qquad
E_{\rm wall}=\tfrac12\dot\psi^2+\tfrac12 K\left(\dot q_b^2+\frac{4q_b^2}{R^2}\right).$$
Since the wall is compact, multiplying by its three-volume gives a finite total energy. The exact radial equation and the two dynamical surface conditions are
$$r q_{tt}=(rq_r)_r-\frac{4q}{r},\qquad
(K\dot q_b)^{\cdot}=-\sum_i q_{i,r}(R)-K\frac{4q_b}{R^2},\qquad
\ddot\psi=\tfrac12K_{,\psi}\left(\dot q_b^2-\frac{4q_b^2}{R^2}\right).$$
The field takes the same value on both wall faces. The radial derivative on each face points toward increasing areal radius. The angular magnetic term affects both the radial propagation and the force on the scalar. The wall has zero velocity, so its mechanical support performs no work.

The initial fields, velocities and scalar are all zero. At each emitting tube impose the incoming characteristic $I_i=(q_{i,t}-q_{i,r})/2$. The returned characteristic is $O_i=q_{i,t}-I_i$. The true net flux into the annulus is $w_s(I_i^2-O_i^2)$, where $w_s=r_s/R$. This is a local radiating boundary condition on a finite sphere. Both angular scattering within the annulus and energy leaving it through the ports are resolved. A source at a different radius can give a different field at the wall despite supplying the same energy.

## Protocol, signed impulse and retained energy

Set the initial kinetic length $K_0=1$, propagation speed one and incident energy per unit proper wall volume one, counting both ports. Lengths and times are reported in these units. Thus a wall of physical radius 2 with $R/K_0=32$ in the table has physical $K_0=0.0625$. The dimensionless coupling $gK_0\sqrt{E_{\rm in}}$ is $\sqrt8$ unless varied explicitly.

We continue the exponential $K/K_0=\exp(g\psi)$ to finite upper and lower values. Write $s=|g\psi|$, $u=\ln20-0.2$ and $x=(s-u)/0.4$. Define $h=s$ below $u$, $h=\ln20$ above $u+0.4$, and $h=u+0.4(x-2.5x^4+3x^5-x^6)$ between them. Then $\ln(K/K_0)={\rm sgn}(\psi)h$. The coefficient is positive and nondecreasing in $\psi$ throughout. The lower cap $K/K_0=1/20$ matters for negative-velocity trajectories. Changing either cap defines a different model.

The source uses the compact positive rising pulse of the planar benchmark, normalized at the emitting ports. Its shape is proportional to $\sqrt{x_t}(1-x_t)$ with $x_t=[1+\exp(-4t/d)]^{-1}$, times smooth cutoffs. It turns on at $t=-3d$ and switches off near $t=d\ln(19)/4$. Here $d$ is the pulse-width parameter. The complete support lasts $3.79861d$. The default run uses $d=1$, travel distance $R-r_s=1$, equal in-phase incident waves and readout $t=20.79861$. The bulk field must propagate across this distance before it drives the scalar.

Starting at rest, the scalar velocity is the difference of two nonnegative integrated forces,
$$\dot\psi=J_E-J_B,\qquad J_E=\tfrac12\int K_{,\psi}\dot q_b^2\,dt,\qquad
J_B=\tfrac12\int K_{,\psi}\frac{4q_b^2}{R^2}\,dt.$$
The table gives the final scalar energy as a fraction of the unit incident dose. In each row the scalar has crossed a cap and is moving further into the constant-$K$ region. Its subsequent acceleration is therefore exactly zero even if some radiation is still near the wall. The finite readout determines its asymptotic kinetic energy under these assumptions.

| $R/K_0$ | Electric impulse $J_E$ | Magnetic impulse $J_B$ | Scalar velocity | Retained scalar energy |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 1.43537 | 1.28102 | +0.15435 | 1.1911% |
| 4 | 1.42204 | 2.12226 | −0.70022 | 24.5151% |
| 8 | 1.30305 | 0.55080 | +0.75225 | 28.2941% |
| 16 | 1.33600 | 0.10897 | +1.22704 | 75.2811% |
| 32 | 1.34306 | 0.02617 | +1.31689 | 86.7106% |
| 128 | 1.34515 | 0.00162 | +1.34353 | 90.2537% |
| 4096 | 1.34529 | 0.00000158 | +1.34528 | 90.4895% |

The flat characteristic calculation gives 90.4898%. At $R=32$ only 5.5816% of the incident dose has returned to the ports by the reported readout. At that time 7.7078% is still in bulk or surface Maxwell fields. The returned fraction and retained scalar fraction can therefore sum to less than one at a finite time.

## Dependence on the source

* **Pulse width.** At $R=8$, halving $d$ gives a positive velocity $+0.75944$ and 28.8375% scalar energy. Doubling $d$ gives $J_E=1.64502$, $J_B=2.63693$, velocity $-0.99192$ and 49.1949%. Equal input energy can drive opposite field directions.
* **Emitter radius.** At $R=4$, changing the travel distance from $0.5$ to $1$ to $1.5$ gives velocities $-1.15185$, $-0.70022$ and $+0.14824$. At $R=32$ doubling the distance from one to two changes the retained fraction from 86.7106% to 86.3625%.
* **Polarity.** Adding a smooth opposite-polarity lobe of duration 40 after the primary pulse makes the signed incident characteristic integral on each port zero. The full packet is renormalized to unit energy. It gives $-0.60624$ at $R=4$ with 18.3761% retained, and $+1.26067$ at $R=32$ with 79.4638% retained.
* **Relative phase.** A perfectly odd pulse drives a Dirichlet field at the shared wall. Its scalar response vanishes. Half the energy in a single-sided input is in this odd sector, which bounds one-sided scalar capture by 50% in the symmetric geometry. For one-sided incidence at $R=32$ the measured scalar fractions are 38.4403% for $g=\sqrt8$ and 43.3553% for $g=4$. Constant $K$ or no incident pulse gives zero scalar energy.

## Independent checks

The weak radial discretization uses degree-four Legendre–Gauss–Lobatto elements. The surface gauge kinetic term adds exactly $K$ to the shared node's bulk mass, and its magnetic term adds $4K/R^2$ to the stiffness. The time integrator advances total nodal momentum and reconstructs velocity using the instantaneous mass. It separately accumulates the two incoming energies, the two returned energies and the two scalar impulses. On the archived nonlinear trajectories the maximum absolute energy residual is below $10^{-8}$ of the input unit. For five representative runs, changing from twelve to eighteen elements per side and simultaneously tightening the integration tolerance changes the scalar energy by at most $1.5\times10^{-7}$. The tested large-radius energy differs from the separately integrated planar limit by $2.7\times10^{-6}$.

With $K$ fixed and time dependence $\exp(-i\omega t)$ the bulk radial solution is an exact linear combination of $J_2(\omega r)$ and $Y_2(\omega r)$. The frequency response and the energy transmitted by finite pulses agree with independent numerical evolution. [Bessel tests and broadband comparison](BESSEL.md) · [Pulse solver](fixed.py) · [Numerical record and time traces](FIXED.json)

These examples give an isotropic, finite-volume mechanism on a prescribed spherical wall. A scalar receiving a negative kick moves away from a potential rise placed at positive $\psi$. The calculations here omit that potential and gravitational evolution. A cosmological trajectory requires a coupled field and geometry solution with the appropriate signed loading.
