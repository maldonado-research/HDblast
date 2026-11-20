# Exact reflection and transmission with spherical Maxwell waves

A stationary Maxwell control tests the curvature term and both boundary conditions without integrating the time-domain equations. The scalar coupling is off and the surface kinetic coefficient is fixed at $K=1$. Each copy of the five-dimensional Minkowski annulus has a source port at $r=R-1$ and shares the same wall at $r=R$. Frequencies and radii below use the same units as the fixed-wall calculation.

For time dependence $e^{-i\omega t}$, the radial amplitude of the lowest transverse vector harmonic on the three-sphere satisfies

$$q'' + r^{-1}q' + (\omega^2-4/r^2)q=0.$$

Its independent solutions are $J_2(\omega r)$ and $Y_2(\omega r)$. We solve for their coefficients using the wall data. For an even incoming channel choose $q(R)=1$ and $2q'(R)=K(\omega^2-4/R^2)q(R)$. For an odd channel choose $q(R)=0$ and $q'(R)=1$. The radial coordinate increases from the port toward the wall on both copies.

At a port set $v=-i\omega q$. In the local characteristic convention of the time-domain drive, the incident amplitude is $I=(v-q')/2$ and the outgoing amplitude is $O=(v+q')/2$. The real wall data enforce $|O/I|=1$ separately for the even and odd eigenchannels. A single incoming wave from one port is a coherent superposition of both channels. If $s_+$ and $s_-$ are their complex reflection amplitudes, the fractions of power returned to the same port and delivered to the opposite port are $|s_++s_-|^2/4$ and $|s_+-s_-|^2/4$. These add to one because the two ports have identical weights. These fractions refer to time-periodic waves at the specified finite ports.

| Wall radius | Frequency | Fraction transmitted to the other port | Wall speed per unit even incident amplitude |
| ---: | ---: | ---: | ---: |
| 2 | 0.45 | 0.001932 | 0.18054 |
| 4 | 0.45 | 0.56956 | 1.72827 |
| 16 | 0.45 | 0.95498 | 1.96180 |
| 2 | 4 | 0.18769 | 0.66017 |
| 4 | 4 | 0.18570 | 0.76566 |
| 16 | 4 | 0.19618 | 0.86149 |

At frequency $0.45$ the transmitted fraction ranges from about $0.19\%$ to $95.5\%$ across these three geometries with the same port-to-wall distance. At frequency 4 all three transmit roughly one fifth. This control has a fixed wall and no transfer to the scalar. It illustrates why a pulse can show strong radius and frequency dependence before scalar dynamics enter the model.

An independent frequency-domain solve uses the mass, gradient, and angular matrices of the time-domain discretization, with incident characteristic drives and outgoing port conditions. Twelve combinations of radius and frequency test both symmetry channels. On the fine grid the maximum absolute discrepancy in the complex reflection coefficient or in the complex wall speed gain is below $10^{-8}$. [Exact and discrete calculations](bessel.py) · [Numerical results](FIXED.json)

For the unit-energy compact pulse injected through the left port, weight the exact single-frequency transmission by the Fourier power of its incident characteristic and integrate over frequency. This predicts right-port energy fractions of **0.15617124** at $R=2$ and **0.50666446** at $R=4$. A time-domain calculation gives 0.15617125 and 0.50666447. At readout the total energy left in the two annuli and wall is less than $5\times10^{-9}$ of the incident dose. That leftover bounds the amount of additional energy that can still reach the right port. Doubling the sampled temporal bandwidth and refining frequency resolution changes the spectral prediction by less than $4\times10^{-10}$. These pulse fractions average over the incident spectrum and can differ appreciably from the transmission at one frequency.
