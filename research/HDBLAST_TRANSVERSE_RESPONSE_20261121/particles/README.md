# Can a supplied transverse wave turn the actual particles?

A local nonlinear test follows charged particles with their finite birth momenta through electric and magnetic fields varying across the expanding wall. The stored axial electric field, parallel magnetic guide, mean expansion and accepted births come from the earlier wall trajectory. New nonzero spatial field modes respond to the moving particles. We compare the result with a second calculation in which the imposed wave propagates in vacuum while test particles respond to it. The spatial Maxwell model here has periodic boundaries and no extra incident bulk drive. Binary scattering is switched off in this experiment. [Full results and controls](RESULTS.json)

At the 1% electrical-work marker a standing magnetic wave with initial peak field 0.04 and initial physical wavenumber $50/L$ produces a transverse RMS momentum of about **2,284 $mc$** in the dedicated case and **2,186 $mc$** in the shared case on the refined runs. Despite the large kick, their transverse cumulative momentum distributions still miss the old thermal targets by approximately **27.5 and 27.3 percentage points**. At some momentum threshold roughly 27 more or fewer particles out of 100 lie below the threshold than thermal equilibrium would predict. The axial electric field initially has magnitude about 0.207. The seed is supplied as an initial condition. The maximum spatial RMS field induced by turning on the current, measured relative to the same vacuum evolution, is about 0.00010 and 0.00016 in this test. [Required transverse impulse](../history/RESULTS.md) · [Seeded linear response](../modes/MODES.json)

| Shared case, standing transverse magnetic seed | Transverse RMS momentum at 1% work ($mc$) | Maximum difference from the thermal transverse CDF |
|---|---:|---:|
| Zero seed, finite angles at birth | 0.272 | Nearly 100 points |
| Peak 0.005, self-consistent wave | 275 | 78.2 points |
| Peak 0.04, self-consistent wave | 2,180 | 27.5 points |
| Peak 0.04, prescribed vacuum wave | 2,177 | 27.6 points |

The table uses the base spatial and birth mesh. Separate runs vary physical wavenumber (5, 50, 250 and 1,000 in units of $1/L$), seed size, polarization, time step, reference wall step, number of harmonics and spatial and birth-angle quadrature. At wavenumber 50, doubling birth and time resolution moves the 0.04 shared CDF discrepancy from 0.275 to 0.273. Increasing retained harmonics from three to seven then gives 0.273. These tests describe prescribed families of perturbations. They do not optimize over possible fields or imply that a particular seed existed.

## What evolves

Use conformal time $d\eta=dt/a$, phase $\phi=k_c x$, comoving momentum $u=ap/m$, conformal field amplitudes $\mathcal E=a^2 E$ and $\mathcal B=a^2 B$, and $v=u/\sqrt{a^2+|u|^2}$. For a positive representative,

$$
\phi'=k_c v_x,\qquad u'=Q[\mathcal E+v\times\mathcal B].
$$

The negative charge at the same spatial phase has momentum $(u_x,-u_y,-u_z)$. Both trajectories then obey the opposite-charge Lorentz forces in the same fields $E=(0,E_y,E_z)$, $B=(0,B_y,B_z)$. This symmetry enforces zero charge density, zero longitudinal-to-wave current $J_x$, and permits both electromagnetic polarizations. The modes satisfy
$$
\mathcal E_y'=-\partial_x\mathcal B_z-a^3J_y/K,\quad
\mathcal E_z'=\partial_x\mathcal B_y-a^3J_z/K,\quad
\mathcal B_y'=\partial_x\mathcal E_z,\quad \mathcal B_z'=-\partial_x\mathcal E_y .
$$

Use $K=2.54$. For a real cosine coefficient the comoving current is $2Q\sum w v\cos(n\phi)$, with the analogous sine projection. The weights sum to the modeled $D$. This is the same effective normalization as the symmetric wall current. Only nonzero harmonics evolve in response to this current. An external mean drive maintains the prescribed homogeneous history; its work is accounted for separately. Density variations arise from particle orbits. The source injects spatially uniform cohorts at quantiles of the recorded accepted mass. The birth polar directions sample the archived angle-dependent acceptance and the azimuths form uniform quadrature rings. Offsetting each ring's spatial sampling makes the unperturbed current cancel in each retained nonzero Fourier mode.

## Energy and an exact particle check

At each step the Maxwell modes propagate freely for a half step. After a half spatial drift and magnetic rotation, the electric impulse and its Fourier current are solved together using the discrete velocity $(u_{+}+u_{-})/(\sqrt{a^2+u_+^2}+\sqrt{a^2+u_-^2})$. This preserves exchange between wave energy and particle energy during the electric update. A second magnetic rotation, drift and free wave step completes the split. Magnetic rotation preserves the particle momentum magnitude. The energy ledger includes injection, work by the prescribed mean field, energy exchanged with the nonzero modes, and the change associated with the mass term as $a$ varies. The largest relative particle ledger residual across the archived cases is below $4\times10^{-13}$. The spatial wave ledger is also checked. [Numerical evolution](patch.cpp) · [Histories and sampling](run.py)

A separate test sets $a=1$, switches off the current response and launches the exact traveling vacuum wave $E_z=A\cos\phi$, $B_y=-E_z$, with $\phi=k(x-\eta)$. A charge started at rest at phase $\phi_0$ satisfies $p_z=(QA/k)(\sin\phi_0-\sin\phi)$, $p_x=p_z^2/2$ and $\gamma-p_x=1$. Its elapsed time follows from $k\Delta\eta=\int_\phi^{\phi_0}(1+p_z^2/2)\,d\phi$. An independent implicit inversion of this integral reproduces the simulated mean energy, transverse RMS momentum and axial pressure share to within six parts per million at the finest step. [Plane-wave check](check.py) · [Measured convergence](CHECK.json)
