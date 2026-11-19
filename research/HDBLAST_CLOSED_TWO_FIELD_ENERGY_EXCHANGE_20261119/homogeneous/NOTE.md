# A homogeneous daughter can exchange energy with the modulus

Take the four-dimensional Einstein-frame scalar kinetic term to be $-[(\partial\Theta)^2+Z(\Theta)(\partial\chi)^2]/2$, in reduced Planck units, with the archived fitted parent potential and a massless shift-symmetric daughter. A minimally coupled Jordan-frame daughter gives $Z=1/f$. A constant rescaling of its field sets $Z=f_*/f$. The results below use that normalization. They hold in the homogeneous sector on a round closed three-sphere.

With proper-time derivatives and physical radius $a$, the shift symmetry gives
$$
P=a^3 Z\dot\chi,\qquad \rho_\chi=p_\chi=\frac{P^2}{2a^6Z},\qquad
\dot\rho_\chi+6H\rho_\chi=Q_\chi=-\rho_\chi\frac{d\ln Z}{dt}.
$$
The parent obeys $\ddot\Theta+3H\dot\Theta+V_{,\Theta}=\rho_\chi(\ln Z)_{,\Theta}$, so $\dot\rho_\Theta+3H\dot\Theta^2=-Q_\chi$. Thus between any two times $a^6\rho_\chi$ changes by the exact factor $Z_i/Z_f=f_f/f_i$. Returning to the same parent field restores this quantity even if the Universe has expanded. The sign of $Q_\chi$ reverses when the parent reverses direction at fixed nonzero coupling. A charge-free daughter remains empty. A nonzero homogeneous charge describes coherent stiff matter with $w_\chi=1$. Radiation, entropy production and a temperature require further physics.

## Backreacted launches

We use $x=H_0t$, $a=H_0a_{\rm phys}$, $v=d(\Theta/M_{\rm Pl})/dx$, $h=H/H_0$, $u=V/V_0$, $r=\rho_\chi/V_0$ and $V_0=3M_{\rm Pl}^2H_0^2$. The code evolves $h'=a^{-2}-v^2/2-3r$ and enforces $h_i=0$, $a_i^{-2}=u_i+r_i$, $v_i=0$ initially. A bounce then has $a''_i/a_i=u_i-2r_i>0$. The Friedmann constraint $h^2+a^{-2}=u+v^2/6+r$ is monitored throughout. The whole trajectory uses the tabulated nonlinear fit.

The chosen launch displaces the parent toward the specified minimum by the archived amplitude $H_0/(2\pi M_{\rm Pl})=0.00140270$ at detuning $10^{-3}$. Each branch has $H_0/M_{\rm Pl}=0.00881344$. Comparisons between $Z=f_*/f$ and $Z=1$ fix initial density, parent position, parent velocity, radius and expansion rate. The free charge normalization changes between these comparisons. Controls with zero initial density use their own constrained bounce radius. In the table $J=\int_0^{x_{\rm cross}}Q_\chi\,dx/(H_0 V_0)$. All columns refer to the first passage through the potential minimum, for $r_i/u_i=0.01$.

| Minimum $\phi_m$ | Daughter metric | $x_{\rm cross}$ | $a^6r/(a_i^6r_i)$ | $J$ | $r(x_{\rm cross})$ | $r'(x_{\rm cross})$ |
|---|---|---:|---:|---:|---:|---:|
| +0.5 | $f_*/f$ | 0.520171 | 1.229546 | +0.00143779 | 0.00657043 | +0.00229670 |
| +0.5 | constant | 0.516528 | 1 | 0 | 0.00539769 | −0.00529462 |
| −0.5 | $f_*/f$ | 0.535344 | 0.899098 | −0.00061340 | 0.00430981 | −0.00954897 |
| −0.5 | constant | 0.541080 | 1 | 0 | 0.00471465 | −0.00787102 |

On the positive-minimum branch, the coupling deposits energy throughout the first descent. At the crossing, its rate $0.00885231$ exceeds the simultaneous dilution $6hr=0.00655561$, giving a brief increase in physical density during expansion. Starting near $r_i=0.01$, the daughter still arrives at the minimum with a lower density than at launch. The negative-minimum branch transfers energy in the reverse direction through its first descent. Lowering the seed to $10^{-5}u_i$ reduces the respective first-passage integrals to $+1.43572\times10^{-6}$ and $-5.99961\times10^{-7}$. For fixed field endpoints the weighted-density factor is independent of charge, while the exchanged amount tends to zero with the seed. At the first expansion maxima of the 1% coupled runs, $x=1.25098,1.49633$ respectively, the daughter fractions are 0.556% and 0.332%.

With the parent initially *exactly* at the crest, a 1% coherent daughter seed gives $v'_i=3r_i(\ln Z)_{,\Theta}=+0.0136269$ in both fits. It pushes the parent to $\Theta-\Theta_*=+0.01$ at $x=0.404875$ and $0.463435$ respectively. This initial force points toward the specified minimum only for $\phi_m=-0.5$. Reversing the sign of the conserved charge leaves the force unchanged. The charged constant-metric and uncharged varying-metric controls keep the parent fixed at the crest through their $x=0.8$ endpoint.

## Replay and checks

From the deposit root, run `python homogeneous/exchange.py --output /tmp/homogeneous.json`. An optional `--gravity` sets the path to the pinned predecessor's gravity directory. NumPy and SciPy are required. The 1% displaced runs were repeated with an independent state representation and integrator, evolving $a,a',\chi'$ with the spatial Einstein equation and the daughter wave equation. The largest difference in the reported times, fields, velocities, radii, expansion rates, densities and integrated sources is $5.04\times10^{-10}$. Over 251 sampled times per displaced run the Friedmann and parent-balance residuals are below $8.1\times10^{-11}$. The weighted daughter-balance residual is below $1.4\times10^{-14}$. The direct representation conserves dimensionless charge to $4.6\times10^{-14}$. These classical tests cover the specified initial states and their first expansion, without an assumption about occupation of short-wavelength daughter modes.
