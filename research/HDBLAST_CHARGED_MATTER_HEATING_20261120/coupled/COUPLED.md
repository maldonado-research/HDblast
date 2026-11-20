# Can coherent wall radiation heat charged matter before turnaround?

The calculation adds a small population of equal-mass particles to the five-dimensional finite-source experiment. The electromagnetic current, gravitational wall junction, scalar, metric, both bulk Maxwell waves, and charged-fluid state evolve together. The scalar decay and the separate radiation fluid stay zero. Opposite charges can gain ordered momentum even in a collisionless fluid. We measure how much energy instead appears as random motion in the local rest frames of the charged fluids. A specified scattering law supplies the irreversibility in this experiment.

## Initial data and field convention

The round wall begins at radius 2 and zero expansion rate. The proper volume factor used throughout is $a^3$ with the constant unit three-sphere volume divided out. The initially cold population has wall-frame rest-mass density $D_0=10^{-4}$ in units of the tuned tension. It is included in the gravity calculation before it conducts. The initial bulk mass on **each** side is fixed by the exact Friedmann constraint, $\mu_0=4-16D_0(2+D_0)=3.99679984$. The zero-seed check restores the original value of 4. Each source, a proper depth 0.002 from the wall, launches the same finite pulse as in the earlier model, with dose 0.123 unless specified. The surface coefficient grows from 0.001 to 13,200.04. The earlier study gives the bulk equations and scalar potential. [Original solver and equations](https://github.com/maldonado-research/HDblast/tree/aee9b3ed9c1f25257df2fcecbd8e47ceb2c7656a/research/HDBLAST_CURVED_MAXWELL_SCALAR_WALL_20261120/loading)

There are three independent Abelian fields. Each has physical electric amplitude $X/\sqrt{3}$ along its own invariant orthonormal axis of the three-sphere. The physical magnetic field is parallel to it. Each field couples to its own equal pair of positive and negative streams. Each of these six streams has one sixth of the total lab-frame number density. Their velocities have equal magnitude, opposite signs within each pair, and point along the corresponding field axis. This arrangement produces zero net spatial momentum and isotropic total spatial stress. It selects a special matter content. Define $Q=|q|/(\sqrt{3}m)$ in the dimensionless field units and $D=D_0(2/a)^3$. Each individual carrier has mass $m$. The total electromagnetic power transferred to the gas per unit wall volume is $JX$, where $J=DQv$ and $v$ is the signed velocity of the positive stream. The collinear bulk flow has zero magnetic Lorentz force. The pressure tensor within each stream is assumed isotropic in its own rest frame.

Two preparations are compared. For pre-existing mobile charges their current responds from the beginning. For a delayed conversion the same conserved rest mass starts in cold electromagnetically neutral precursors. At a prescribed wall proper time they convert into balanced pairs at zero drift and zero initial temperature. In this ideal mass-matched limit no energy is added or subtracted at conversion. The field equations and the gas stress are continuous there. No binding energy, ionization dynamics or reaction rate is supplied. The reference time 10.5017 was chosen from the first zero of the magnetic field after the scalar descends in the **neutral-seed** control. This selects a phase with strong electric field, which can do work on charges. Earlier and later switch times test this choice. The proper-time switch does not arise from a calculated chemical process.

## Relativistic charged gas and drag

Write $u=\gamma v$, $\gamma=\sqrt{1+u^2}$, $\Theta=p_0/(n_0m)$ for the proper pressure divided by the proper rest-mass density of one species, and $h$ for specific enthalpy including rest energy. We use the Taub–Mathews ideal-gas closure
$$
 h(\Theta)=\frac52\Theta+\sqrt{1+\frac94\Theta^2}.
$$
It recovers the cold and ultrarelativistic ideal-gas limits. For $\Theta\geq0$ the adiabatic squared sound speed is $\Theta h'/[h(h'-1)]\leq1/3$. The mixture's energy and isotropic spatial pressure in the wall frame are
$$
 \rho_g=D(h\gamma-\Theta/\gamma),\qquad
 p_g=D\left(\frac{h u^2}{3\gamma}+\frac{\Theta}{\gamma}\right).
$$
Each stream's spatial momentum per unit wall-frame rest density is $hu$. The lab relaxation force is $-\nu hu$; its opposite-sign partner experiences the opposite force. The summed collisional energy exchange with external matter vanishes. A positive drag rate converts ordered motion into internal energy through
$$
 \frac{d(hu)}{d\tau}+Hhu=QX-\nu hu,\qquad
 (h'-1)\dot\Theta+\frac{\Theta u}{\gamma^2}\dot u
 = -3H\Theta+\nu h u^2 .
$$
For conducting runs these coupled equations are solved simultaneously at every Runge–Kutta stage. Particles do no electric work until activation. The temperature parameter of this closure is expressed relative to the chosen rest mass. It gives no independent temperature in kelvin. The closure assumes each stream equilibrates internally, approximates the massive equation of state, and puts relative momentum exchange into a prescribed wall-frame drag rate. It does not derive a collision kernel, a photon distribution, chemistry or a temperature shared by all species. Radiation from accelerating charges and pair creation/annihilation would need separate terms at these high energies. [Taub–Mathews closure, Mignone and McKinney, Eq. 13](https://arxiv.org/html/0704.1679#S2.E13)

The current $-J$ is added to $d(KX)/d\tau$ and $\rho_g,p_g$ enter the junction and the radial geometry. These choices give two independent energy identities,
$$
 \frac{d(a^3\rho_g)}{d\tau}+3Ha^3p_g=a^3JX,\qquad
 \frac{d(a^4 e_{\rm EM})}{d\tau}=a^4\left[2(I^2-O^2)+\frac{\dot K}{2}(B^2-X^2)-JX\right].
$$
Here $I,O$ are the incoming and outgoing Maxwell characteristics at one of the two identical wall faces. The prefactors in the two energy integrals differ, so their raw totals have different normalizations. A third check uses
$$
 W=(a/2)^2\gamma^{2/3}\Theta\exp[\mathrm{asinh}(3\Theta/2)],\qquad
 \dot W=\frac23 (a/2)^2\gamma^{2/3}\exp[\mathrm{asinh}(3\Theta/2)]\nu h u^2\geq0 .
$$
This is an entropy variable obtained by integrating the per-particle proper first law; it stays regular for the initial cold state. Its direct state value is compared with its positive time integral during the simulation. In the collisionless limit, $\Theta$ remains zero and the charged gas obeys the analytical impulse relation $u=Q[C-a_{\rm on}C_{\rm on}/a]$, where $\dot C=X-HC$. The simulation checks this relation and the gravity and electromagnetic energy constraints independently.

## What counts as random energy

The main diagnostic is $D(h-1-\Theta)$, the proper thermal energy above rest mass per particle, multiplied by the number of particles per unit **wall** volume. This excludes the energy of bulk counterstreaming. The full lab energy associated with thermal internal energy and pressure is $\rho_g-D\gamma$, which also includes the boost of hot material and is at least as large. We report both as fractions of the total added wall energy $e$, including all particle rest energy, scalar, and electromagnetic energy. The background tension is excluded. A rise in a gas energy density alone can also arise from reversible drift. This distinction matters in the collisionless controls and near a plasma oscillation.

## A successful finite-interval conversion and its controls

In the selected experiment the wall bounces at radius **1.88104**, the scalar completes its descent at wall time **4.51322**, and the neutral carriers convert at time **10.5017**. The chosen projected charge-to-mass ratio is $Q=10^8$ and the drag rate is $\nu=6$. Both are free dimensionless parameters. The first post-exit magnetic zero on the independently evolved neutral precursor is at time **10.5017103** on the 48-element early mesh. It was used to choose a phase and was not fitted to a measured cosmic signal. The parameter and timing surveys reported here are exploratory.

In the 64-element primary run, intrinsic random particle energy reaches **91.6638%** of all added wall energy, and the lab thermal excess at its peak is **91.6643%**. At the peak of the lab thermal fraction the drift speed is only **0.00256** times light speed. The closure gives $\Theta=1851$ there and $p_g/\rho_g=0.33333332$. The intrinsic fraction exceeds 90% at wall time **11.03911**, radius **9.65411**. It stays above this level during a single expanding interval until the first maximum radius **9.76159**, time **12.42350**, corresponding to **1.113%** growth (0.0110714 e-folds). The run subsequently contracts. This short expanding interval and the assumptions about matter preparation rule out reading the example as the thermal history of the observed universe.

At the conversion event $a^4e_{\rm EM}=43.360626$. By the peak of the lab thermal fraction the net later bulk flux adds **−0.000786**, the scalar adds **zero**, work on the gas removes **43.357336**, and **0.002503** remains in the coherent field. These numbers close to within $2\times10^{-10}$ in the stated comoving units. About 99.992% of the field's comoving energy at activation has gone into the work integral. The rest-frame random-energy diagnostic and the independent gas energy/pressure work ledger distinguish that electric work from reversible charged drift. The numerical value of $\Theta$ tests an ultrarelativistic closure; an actual plasma at such energy could radiate, create pairs or require a different collision law. A normalized microphysical model would have to calculate those processes.

For the following sweep the initial mass, pulse, $Q=10^8$, $\nu=6$ and 24-element early grid are held fixed. The table gives the maximum intrinsic random-energy fraction during post-exit expansion.

| Wall clock when the neutral precursors become mobile | Peak intrinsic share of added wall energy |
|---:|---:|
| 4.7 | 0.86% |
| 7.0 | 46.47% |
| 9.5 | 87.90% |
| 10.5017 | 91.65% |
| 11.0 | 90.53% |
| 11.5 | 87.39% |

**How much energy is available at each phase?** After the scalar descends, $K$ is constant. Define $M=a^2B$, $Y=a^2 X$, and $F=\int 2a^4(I^2-O^2)d\tau$ between activation and the comparison time. Faraday's law gives $\dot M=2aX$. With sign conventions fixed above the exact ledger is
$$
 \int a^4JX\,d\tau=\frac{K}{2}\left(Y_{\rm on}^2-Y_*^2+M_{\rm on}^2-M_*^2\right)+F.
$$
If the remaining electric field and net bulk flux are small and high conductivity approximately freezes $M$, the predicted fraction of the initial comoving wave energy converted to electrical work is $X_{\rm on}^2/(X_{\rm on}^2+B_{\rm on}^2)$. In table order its values are **0.00905, 0.49623, 0.95547, 1.00000, 0.98927, 0.95773** on the 48-element neutral-seed trace. The 24-element coupled runs give **0.00946, 0.50706, 0.95896, 0.99993, 0.98782, 0.95609** for the comoving work $\int a^4JX\,d\tau$ divided by the sheet energy $a_{\rm on}^4 e_{\rm EM,on}$, evaluated at the peak expanding post-exit lab thermal fraction in each run. They agree to within 1.09 percentage points. This is a conditional phase prediction for the chosen highly responsive conductor. It may fail at weaker coupling or with a significant continuing input flux. The reference fields are sampled at the specified wall clock (linear interpolation between 10.5 and 10.6 for 10.5017). Run `python -I -B phase.py` to generate all six predictions and compare them with the coupled result. [Phase calculation](phase.py)

An electric field does work directly on moving charges. Magnetic energy can only reach them through subsequent field evolution. Activating early also changes the later electric field through the current, so phase comparisons require the coupled Maxwell calculation. Switching off collisions at time 10.5017 gives finite charged drift and exactly zero internal energy. Reducing $Q$ by a factor of ten at that same time gives 63.99% intrinsic random energy. Keeping $Q=10^8$ and changing the drag rate to 2 or 20 gives maxima of 90.87% and 91.59%, respectively. Equal total energy in waves can lead to substantially different heating histories.

With charges mobile from the start ($Q=5{,}000$, $\nu=10$), the particles reach **6.56%** intrinsic random energy during expansion. The wall bounces but its scalar does not complete the descent by radius 16. The otherwise identical collisionless gas stays internally cold. Raising the source dose from 0.123 to 0.130 restores a complete scalar exit and subsequent contraction. The overall hot fraction peaks at 6.62% during expansion, while the peak **after the exit** is 2.21%. This is a different preparation and part of its heating happens before the late electromagnetic field becomes available.

## Reproduction and numerical accuracy

The C++17 program retains the parent radial constraints, proper-time source, returning wave, wall Israel equations, upwind cubic elements and Runge–Kutta integration. It explicitly aligns the time step with delayed conversion and locally limits the fractional change in the gas variables while the accelerated particles are still cold. The default source-normalization quadrature uses 32,768 Simpson panels. Its factor is 0.4998200679213546. The 64-element primary uses six late elements; the 48-to-64 change in the peak intrinsic fraction is $8.4\times10^{-6}$, and the inferred duration shifts by $7.3\times10^{-7}$ e-folds. At a fixed 24 early elements, increasing the late resolution to ten elements and decreasing the time-step factor from 0.08 to 0.05 changes the peak intrinsic fraction by $1.3\times10^{-6}$. Compilation with separate optimization and floating-point contraction disabled reproduces the primary peak and duration to the reported precision.

During the primary run the minimum bulk mass is 0.096536, the minimum lapse is 0.997781 and the radial squared gradient stays positive. Peak Friedmann, gas comoving energy, total wall energy and Maxwell comoving energy residuals are approximately $2.4\times10^{-10}$, $1.1\times10^{-7}$, $1.1\times10^{-7}$ and $3.5\times10^{-10}$ in model units. The normalized error in the independent integral of the entropy variable is $1.1\times10^{-6}$. The tube mass discrepancy is below $3.8\times10^{-6}$. The collisionless impulse is verified to $7.0\times10^{-11}$ relative to $1+|u|$ in the delayed control. The matter remains in a regular positive-mass region throughout the reported primary trajectory.

Run from this directory with Python 3 and a C++17 compiler:
~~~sh
python -I -B charged.py --output /tmp/charged-result.json --jobs 2
~~~

The default parameter file lists all 22 controls and resolution checks. Local temperatures and share maxima use the proper wall clock. Absolute temperatures, microscopic scattering rates, photon production and comparison with cosmological data would require further inputs. [Parameters](cases.json) · [Coupled solver](charged_forward.cpp) · [Run driver](charged.py) · [Machine-readable results](CHARGED.json)
