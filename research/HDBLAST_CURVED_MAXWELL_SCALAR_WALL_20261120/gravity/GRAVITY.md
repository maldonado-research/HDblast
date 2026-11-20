# A closed gravitating wall driven by a finite Maxwell pulse

An observer fixed at areal radius $r_s=1.65$ launches an outgoing electromagnetic pulse in each of two identical five-dimensional spacetime regions. The waves reach a moving spherical wall. We solve the Einstein and Maxwell equations between the emitter and the wall, including the reflected waves, curvature, gravitational mass and scalar backreaction. A gauge-invariant interaction on the wall determines the energy transfer.

Two cases distinguish loading the scalar from reversing contraction. In the first, the wall starts at rest and 76% of the electromagnetic energy incident at the wall does net work on its scalar in the simulated interval. The scalar retains kinetic energy and the wall contracts. In the second, a wall **already expanding at launch** first contracts and then expands again. The supplied field-dependent coupling and positive potential are essential to the second example. The initial spacetime, pulse and wall law remain input choices.

## Bulk and wall equations

The AdS length and tuned wall tension are one and $\kappa_5^2=6$. Each retained interior has metric

$$ds^2=-N^2 f\,dt^2+f^{-1}dr^2+r^2d\Omega_3^2,\qquad f=1+r^2-m/r^2.$$

The time coordinate equals wall proper time. Three Abelian potentials point along the three orthonormal left-invariant one-forms on the unit $S^3$, each with common amplitude $A/\sqrt3$. Their total stress is isotropic. The bulk Lagrangian per unit three-sphere is $-r(DA)^2/2-2A^2/r$ with the two-dimensional volume element. Let $P=A_t/(Nf)$ at fixed radius, $D=A_r$ and $c=Nf$. The wave and radial Einstein equations become

$$D_t=(cP)_r,\qquad P_t=(rcD)_r/r-4NA/r^2,$$
$$m_r=2rf(P^2+D^2)+8A^2/r,\qquad (\log N)_r=2(P^2+D^2)/r,\qquad
 m_t=4rNf^2PD.$$

For example, $m$ decreases at a fixed radius for a purely outward pulse, where $D=-P$. The last equation advances the mass enclosed by the emitter. We independently check it inside the annulus and at the wall.

Write $a(t)$ for the wall radius, $v=\dot a$, $H=v/a$ and $X=\sqrt{f_b+v^2}$. The wall normalization fixes $N_b=X/f_b$. Both normals point toward larger radius within their respective retained regions. Let $w=\dot A_b=XP+vD$ and $z=nA=vP+XD$. The outward-travelling characteristic arriving at the wall is $U=P-D$, which gives $z=w-(X-v)U$. The boundary condition for the returning characteristic is $P+D=[2w-(X-v)U]/(X+v)$. On the wall the invariants are $E^2=w^2/a^2$ and $B^2=4A_b^2/a^4$.

For surface Lagrangian $-(\partial\psi)^2/2-V(\psi)-K(\psi)\sum F_i^2/4$, with $\pi=\dot\psi$, the sheet and scalar equations are

$$\frac{d(Ka w)}{dt}+4KA_b/a=-2az,\qquad
\dot\pi+3H\pi+V'=K'(E^2-B^2)/2.$$

The surface Maxwell energy $e_g=K(E^2+B^2)/2$ contributes $e_g/3$ to pressure. The scalar contributes $e_s=\pi^2/2+V$ and $p_s=\pi^2/2-V$. With $e=e_g+e_s$ and $p=e_g/3+p_s$, the spatial and temporal Israel equations are

$$X/a=1+e,\qquad
\frac{\ddot a}{a}=(1+e)(1-2e-3p)-1-m_b/a^4-2T_{nn},\qquad
 T_{nn}=\frac{w^2+z^2}{2a^2}-\frac{2A_b^2}{a^4}.$$

The two-sided influx is $Q=-2wz/a^2$. The areal mass evaluated at the moving wall also satisfies $\dot m_b=4a^3(Xwz/a^2-vT_{nn})$. The negative magnetic term in $T_{nn}$ follows from the Maxwell stress. It can give a negative normal pressure while Maxwell energy remains nonnegative.

## Source and controls

Initially $a=2$ and all Maxwell and scalar fields vanish. The relation $m_i=4+16H_i^2$ enforces the spatial junction. The stationary tube supplies $U_s=S\sin^4(\pi \tau_s/0.083)$ for $0<\tau_s<0.083$ and zero outside. The clock on that tube advances by $d\tau_s=N_s\sqrt{f_s}\,dt$. Its inward-travelling characteristic leaves the annulus through the tube. Thus the interior receives any returning energy and the enclosed mass changes accordingly. The physical orthonormal outgoing characteristic amplitude is $\sqrt{f_s}U_s$. The initial packet is zero and the source has finite duration and energy. Its internal power supply is unspecified.

The coupling $K=0.055\exp[\log(C)\tanh(\psi/0.017)]$ is positive and bounded. For the initially resting case use $S=4.5$, $C=20$ and $V=0$. At $t=0.5$ the wall registers incident proper energy 0.98572 and returning energy 0.20214. The net work on the scalar is 0.75296 in units of tension times unit three-sphere volume. The scalar density is 0.14162 while the surface electromagnetic density is 0.00577. The final contraction speed is 0.71029. The surviving gauge energy, pressure work during contraction and scalar work enter the wall energy ledger. A constant $K$ gives zero scalar work from the same initially empty state and the same prescribed source. The input amplitude 2.3 with varying $K$ yields a smaller work fraction of approximately 32%.

The second case uses $H_i=0.12$, $S=4.8$ and $C=40$. The potential is zero for $\psi\le0.035$. It rises to 0.10 over $0.035<\psi<0.065$ according to $10y^3-15y^4+6y^5$ times 0.10, with $y=(\psi-0.035)/0.030$. It stays constant above the ramp. The field acquires its kinetic energy from the electromagnetic interaction. The resulting potential energy produces negative wall pressure. On the 320-cell mesh the radius turns from expansion to contraction near 2.049 and has a subsequent minimum near **1.8543**, safely outside the source. At wall time 2.2 the radius is 1.8840 and the velocity is positive. At the second turn the total energy density is about 0.1141, the pressure about −0.090 and the acceleration divided by the radius about +0.1114. By the end 0.85605 units of electric work and −0.02323 units of magnetic work give 0.83281 net scalar work. Dividing by the local wall-incident energy 0.96225 gives **86.5%**. The stationary source registers 1.22027 units of emitted energy. Source and wall energies refer to different observers. Magnetic curvature opposes scalar loading in the integrated budget for this pulse. Coherent scalar work and the electromagnetic sheet are the only wall matter included here. This calculation does not assign a particle temperature to either.

## Numerical checks

We evolve the two characteristic fields $P-D$ and $P+D$ on $x=(r-r_s)/(a-r_s)$. We integrate $A$ inward from its wall value. The two radial gravitational constraints use cumulative Simpson quadrature with the inner mass and outer lapse as boundary data. The unsquared spatial Israel condition provides $X=a(1+e)$ for the characteristic boundary values. The geometric $\sqrt{f_b+v^2}$ from the bulk supplies an independent check and sets the lapse. The wall acceleration evolves through the temporal Israel equation. Third-order upwind-biased differences in the interior and second-order closures near the boundaries are paired with fourth-order Runge–Kutta time steps. Each step is constrained by the local characteristic cell-crossing time. We stop each history at its listed final time or before reaching the fixed source, at approximately $a=r_s+0.025$.

The radiation absorbed by a surface with a prescribed efficiency was replaced here by a sheet boundary condition. For each run we integrate the incident and returning powers $a(w-z)^2/2$ and $a(w+z)^2/2$, their difference, scalar work, and the pressure work. We check $a^3e+\int 3a^2p\,da=\int a^3 Q\,dt$ and its scalar counterpart. At three internal radii moving with the grid, we check the local temporal Einstein mass equation and the time derivative of the reconstructed field. The empty-bulk control reproduces $a^2=4-t^2$ for an initially stationary radiation-free wall.

In the positive-potential case the second minimal radius on 160, 240 and 320 cells is approximately 1.8568, 1.8548 and 1.8543 respectively. A further 640-cell check gives a minimum radius of 1.8540 at time 1.6656. The maximum spatial-junction discrepancy falls to $3.7\times10^{-6}$ and the interior mass discrepancy to $8.0\times10^{-5}$ on that grid. The 320-cell peak spatial-junction discrepancy is below $3.0\times10^{-5}$ and the temporal mass discrepancy at the three internal radii is below $6.4\times10^{-4}$. The source mass stays above 0.089 and the computed annulus remains untrapped. A marginal case with the same $H_i$ and incident amplitude, $C=20$ and plateau height 0.09 illustrates why that convergence check matters. A 48-cell mesh falsely predicts a bounce. Finer grids reach the inner cutoff while the wall is still contracting. Beyond that cutoff the model would require additional information about the emitting interior. All numerical experiments here assume exact spatial homogeneity of the wall and isotropy of the summed triplet. They address a class of Einstein–Maxwell–wall solutions and leave the origin and observational consequences of the higher-dimensional pulse open.

[Coupled evolution and reproduction](gravity.py) · [Numerical results](GRAVITY.json)
