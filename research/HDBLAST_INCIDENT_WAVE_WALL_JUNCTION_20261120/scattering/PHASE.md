# The return time of a reflected packet

A static wall reflects all the incoming continuum flux in this single, lossless channel. The two electromagnetic couplings can still be distinguished by the shape and timing of the return. We compare the fixed central cut of the same tuned geometry, with no brane kinetic term, no incoming massless mode and the same detector convention. The earlier weight is $G=(1-\phi^2)^2$. The weight selected by the fixed-profile null-sector condition is $G=e^{2A}$ with $A(0)=0$. Masses and wave numbers are given in inverse wall widths. The time unit is $L_0/c$. [Preceding spectral calculation](../prior/retarded.py)

Write $z=\int_0^y e^{-A}dy$, $j=\sqrt{G e^A}$, $u=j\psi$, and $D=u_z-(j_z/j)u$ for a transverse mode. The brane requires $D(0)=0$. Use a time dependence $e^{-i\omega t}$ and an outward reference phase $e^{imr}$ with $m>0$ and $\omega=\sqrt{k^2+m^2}$. Here $r=e^{(5/9)y-A_\infty}/(5/9)$ and $A_\infty=\lim_{y\to\infty}[A(y)+(5/9)y]$. The asymptotic offset between $r$ and $z$ on this geometry is common to both weights.

If $f$ is the outgoing radial solution, the reflected field with a unit-amplitude incident wave is $f^*+\mathcal R f$. Write $d_f=f_z-(j_z/j)f$ at the wall. Then $\mathcal R=-d_f^*/d_f$, so $|\mathcal R|=1$. Real coefficients and the reflecting boundary conserve flux. The two weights change the complex phase. The quantity compared between the two systems is $\mathcal R_{\rm selected}/\mathcal R_{\rm earlier}$ at the same $m$.

In the large-$r$ tail the outgoing solution is proportional to $\sqrt{\pi mr/2}\,H^{(1)}_\nu(mr)$, with $\nu=2$ for the selected weight and $\nu=4.6$ for the earlier weight. We remove the order-dependent phase $-\pi\nu/2-\pi/4$ to refer both outgoing amplitudes to $e^{imr}$. Differentiating the relative phase with respect to $\omega$ gives the relative centroid delay of a sufficiently narrow flux-normalized packet. The factor $d m/d\omega=\omega/m$ is included at fixed $k=1$. Negative values mean that the energy centroid arrives earlier in the selected model. The causal front propagates on the common light cone.

| $m L_0$ | Relative phase (rad, principal branch) | Relative delay ($L_0/c$) |
| ---: | ---: | ---: |
| 0.25 | 1.64477 | −7.16313 |
| 0.5 | 1.08810 | −5.89119 |
| 1 | −0.51965 | −5.13734 |
| 2 | 2.95930 | −1.84071 |
| 4 | 1.41519 | −0.37405 |

There is also a finite-width test. The incident flux spectrum per unit $dm$ is proportional to $\exp[-1/(1-x^2)]$ for $|x|<1$ and vanishes outside that support, with $x=(m-m_c)/(0.15 m_c)$. The returning packet has the same spectral flux at infinity. Its relative centroid displacement is the flux-weighted average of the relative phase derivative. This follows from the first time moment of the Fourier transform of a flux-normalized packet, or by integrating the time moment of the asymptotic energy flux. At $m_c=0.5,1,2$ we obtain **−5.89709, −5.12441 and −1.86551** in units of $L_0/c$. At very different bandwidths the mean need not equal the central group delay or the peak delay.

One integration propagates the Neumann solution outward and extracts its Wronskian with the outgoing wave. A second propagates the outgoing wave inward from a deeper surface, at tighter integration tolerance and halved phase-derivative step. Their relative group delays differ by less than $9\times10^{-8}$ in the stated units over the five benchmark masses. The packet integration uses two Gaussian quadrature orders. For the packet centered at $m=1$, the [replay](../reproduce.py) also transforms the reflected amplitudes into the time domain on a uniform frequency grid and computes the energy centroid of each return. Subtracting the two centroids gives $-5.1244134$ in the same units. This timing comparison concerns asymptotic incoming continuum packets on a fixed central slice. [Calculation](phase.py) · [Full numerical results](SCATTERING.json)
