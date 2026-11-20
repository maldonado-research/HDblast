# What can the pressure of an incoming wave do?

In the earlier experiment the boundary followed a specified Gaussian trajectory. Here the wave and the boundary evolve together. This answers a more limited physical question. Can a passive mechanical wall acquire energy from a packet whose initial projection on the massless mode is zero? The additional spring, inertial mass and damping specify an illustrative material. A separate calculation below checks what gravity and the original scalar boundary condition would require.

## Wave and boundary

Use the central flat wall with unit warp factor at rest and the same even transverse field as in the [scattering experiment](../scattering/SCATTERING.md). The domain is one retained half-line of conformal depth $z\ge b(t)$. The action per tangential mode is $\tfrac12\int dt\,dz\,w(z)[q_t^2-q_z^2-k^2q^2]$. The two choices of gauge kinetic term give $w=e^A(1-\phi^2)^2$ or $w=e^{3A}$. The packet occupies $6<z<10$ at launch, $k=1$, and both weighted initial moments vanish. In length and time units of the central-wall width, use mechanical mass per normalized tangential volume $M=1$, natural frequency $\omega=1$ and damping rate $\Gamma=0.8$. Energies are in units of $Mc^2$ for that same normalization.

At a freely reflecting moving boundary the natural field condition is $q_z=-\dot b\,q_t$. Conservation of field energy on its time-dependent domain then gives

$$
P=\tfrac12 w(b)[(1-\dot b^2)q_t(t,b)^2-k^2q(t,b)^2],\qquad
\dot E_{\rm field}=\dot b\,P.
$$

The force on this lower boundary is $-P$. Advance the mechanical equation and accumulated damping work,

$$
M(\ddot b+\Gamma\dot b+\omega^2 b)=-P,\qquad
D(t)=M\Gamma\int_0^t\dot b^2\,dt'.
$$

The continuum conserves $E_{\rm field}+M(\dot b^2+\omega^2b^2)/2+D$. The wall starts at rest. A fixed-wall control can feel pressure while receiving zero mechanical work. A mobile wall with $\Gamma=0$ moves without energy entering the damping channel. The absorber in this experiment is a phenomenological mechanical degree of freedom. Identifying a physical radiation bath would require a material and a thermalization calculation.

## Measure the light that stays near the wall

Write $I=\int_b^\infty w q_t\,dz$, $J=\int_b^\infty w q\,dz$ and $N(b)=\int_b^\infty w\,dz$. The field equation, moving Neumann condition and Leibniz rule give $\dot I=-k^2J$ and $\dot J=I-\dot b\,w(b)q(t,b)$. Consequently

$$
\frac{d}{dt}\left[e^{-ikt}(I+ikJ)\right]=-ik e^{-ikt}\dot b\,w(b)q(t,b).
$$

Computing both sides by different routes measures the instantaneous projection energy $E_0=(I^2+k^2J^2)/(2N(b))$. At rest it is the energy of the normalizable massless mode. When the boundary continues moving the projection can change. Measuring $E_0$ directly is important because a total-energy residual can exceed this small quantity.

Let $\epsilon$ denote incident energy in the chosen units. For a balanced packet of fixed shape, its wave amplitude is proportional to $\sqrt{\epsilon}$. Pressure drives a displacement of order $\epsilon$. Hence the first nonzero damping energy is of order $\epsilon^2$ and the source in the rotating-moment identity is of order $\epsilon^{3/2}$. When its leading integral is nonzero, the projected energy is of order $\epsilon^3$. The coefficients depend on the packet, the measurement time and the material. Holding the wall fixed to obtain its leading forcing provides an independent check on both coefficients.

## Numerical checks

Use the map $z=x+b(t)h(x)$, with $h=(1-x/L_m)^4$ for $0\le x<L_m$ and zero above it. Evolve $q$ and the canonical momentum density $w(z)(\partial z/\partial x)q_t$ with centered fluxes on $0\le x\le60$. A second flux term accounts for transport by the grid. The boundary flux vanishes in these coordinates. The input is rebalanced and normalized on each mesh before evolution. The output includes independent spatial and temporal integrals for the massless projection, a frozen-boundary forcing prediction, an alternate mesh map, a later readout, energy closure and the far-boundary field. The principal readout at $t=36$ and the later one at $t=44$ precede even the earliest possible arrival of the launch packet at $z=60$.

At incident energy 0.04 and time 36, extrapolating the two finest meshes ($\Delta z=0.0125,0.00625$) gives the following fractions of the incoming energy. The change from the preceding extrapolation is at most 0.50% of the projected energy across the two materials.

| Gauge kinetic term | Damping work | Projected massless light |
| --- | ---: | ---: |
| $(1-\phi^2)^2$ | 0.496% | 4.41 parts per million |
| $e^{2A}$ | 0.850% | 1.48 parts per million |

The material with greater damping work in these runs retains less bound light. The finest direct projection and source integral agree within 0.0062% for either material. The largest instantaneous relative energy imbalance on those two finest example runs was $1.04\times10^{-6}$. At weak incident energy 0.0025 the independently predicted leading dissipated and projected energies agree with the coupled result to within 0.37% on the 0.025 mesh. These figures describe this packet family and the specified boundary law.

[Detailed runs and convergence](RESPONSE.json) · [Evolution code](response.py) · [Gravitational and scalar junction check](COMPATIBILITY.md)
