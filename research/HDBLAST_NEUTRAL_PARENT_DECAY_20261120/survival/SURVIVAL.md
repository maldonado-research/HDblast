# How long after opening the decay threshold do the charged daughters appear?

Opening a decay channel does not immediately convert the parent population. On the published cold-precursor wall, a particular smooth mass law opens at proper time 10.341. At the dimensionless rate normalization $\kappa=50$, half of the population converts by 10.964. At $\kappa=1000$, half converts by 10.475. Both calculations use exactly the same function of the local scalar field. A separate finite-time quantum calculation checks when it is reasonable to use a local decay rate.

## Two-body phase space and survival

Assume a neutral real scalar parent $\chi$ at rest with a Yukawa interaction $-g\chi\bar f f$, producing one Dirac particle and its antiparticle of fixed mass. There is no color multiplicity. For three degenerate daughter flavors, independent final states add at this order and the summed width uses $g^2=\sum_i y_i^2$. Time is in units $L/c$ and the daughter mass parameter is $\mu=mcL/\hbar$. Throughout this note $M$ denotes the parent mass in the same dimensionless units. Choose a smooth positive mass function,

$$
M(\psi)=2\mu R(\psi),\qquad R=0.95+0.10\,S[(\psi-\psi_*)/w+1/2],
$$

where $S=0$ for argument $z\le0$, $S=1$ for $z\ge1$, and $S=[1+\exp(1/z-1/(1-z))]^{-1}$ between the endpoints. The scalar crossing is $\psi_*=0.77$ and $w=0.03$. These are chosen near the heating window of the previous charged-fluid study. The parent is initially stable against the stipulated vacuum two-body channel, including during the early electromagnetic pulse.

The local flat-space tree width in wall time is

$$
\Gamma=\frac{g^2 M}{8\pi}\left(1-4\mu^2/M^2\right)_+^{3/2}
=\kappa R(1-R^{-2})_+^{3/2},\qquad \kappa=\frac{g^2}{4\pi}\mu .
$$

The power $3/2$ includes the scalar-to-Dirac spin matrix element and two-body phase space. This width assumes low occupation, a sufficiently dense final spectrum, a locally adiabatic background and negligible electromagnetic effects during formation. For independent parents maintained at rest, $F(t)=1-\exp[-\int_0^t\Gamma dt']$ is the exact solution **of that rate equation**. It sums depletion. It has yet to incorporate microscopic memory or Pauli blocking.

We use the $D_i=10^{-4}$, electrically inert precursor trajectory, including its rest-mass effect on gravity. The table imports that trajectory without a change to its geometry. On this wall $R$ first reaches 1 at 10.34058, reaches its upper plateau at 10.82063, and expansion ends at 12.42343. The integral of $R(1-R^{-2})_+^{3/2}$ by turnaround is **0.0573111**.

| $\kappa$ | Fraction converted at 10.5017 | Fraction at 10.8293 | Fraction by turnaround | Time of 50% conversion |
|---:|---:|---:|---:|---:|
| 10 | 0.0106 | 0.0939 | 0.436 | Not reached |
| 50 | 0.0520 | 0.389 | 0.943 | 10.964 |
| 100 | 0.101 | 0.627 | 0.997 | 10.731 |
| 1000 | 0.656 | 0.99995 | $>0.99999$ | 10.475 |

Releasing 90% by turnaround requires $\kappa=40.18$. Releasing 90% by 10.8293 requires 233.5. Releasing it by the earlier constitutive model's clock time of 10.5017 requires 2,158. Changing $\psi_*$ from 0.74 to 0.79 changes the threshold time from 9.433 to 10.984 and raises the turnaround requirement from 27.19 to 60.62. Changing $w$ from 0.015 to 0.06 at fixed $\psi_*=0.77$ preserves the threshold time and raises the turnaround requirement from 38.59 to 43.81.

At birth each daughter has $p/\mu=\sqrt{R^2-1}\le0.3202$, so its kinetic energy is at most $0.05\mu$. The source of its rest and kinetic energy is the parent, including work by the prescribed scalar-dependent mass. At fixed initial comoving parent count one,

$$
(1-F_f)R_f+\int R\,dF=R_i+\int(1-F)\,dR ,
$$

in energy units of $2\mu$. The numerical output checks this identity and redshifts each birth momentum by $a_{\rm birth}/a_f$ to compute the remaining kinetic energy without electric acceleration or collisions. At $\kappa=50$ the field supplies 0.09750 in these units, 0.04465 becomes daughter kinetic energy at birth and 0.04365 remains at turnaround. Large late charged-gas heating, if realized, needs subsequent energy transfer.

## Matching the requested force strength

The chosen force parameter of the main charged calculation is $Q=q/(\sqrt3\mu)=10^8$, where $q$ is the norm of a charge to a single gauge component in the published field normalization. The Yukawa coupling $g$ is independent of this gauge charge. At $q=1$, matching $Q$ gives $\mu=5.77\times10^{-9}$. Even allowing $g^2/(4\pi)=1$, the rate equation gives only $3.31\times10^{-10}$ conversion by turnaround. Choosing $q=10$ yields $3.31\times10^{-9}$. Applying these figures as literal quantum predictions would extrapolate the local approximation far outside its range. The accumulated zero-momentum detuning phase after the gate opens is only $0.19519\mu$. For these ultralight parameters it is far below one radian across the entire available interval, and the closed-space spectrum matters.

For $g^2/(4\pi)\le1$, 90% by turnaround instead needs $\mu\ge40.18$ and hence $q\ge6.96\times10^9$ at the same $Q$. Electromagnetic distortion is another consistency check. Small field action over a Compton time requires $|qX|/(\sqrt3\mu^2)=Q|X|/\mu\ll1$. With late $|X|$ of order $9\times10^{-4}$ and the same $Q$, even that weak-field condition requires $\mu\gg9\times10^4$. The condition is necessary for the vacuum-width interpretation in this regime and is insufficient to establish an adiabatic threshold or thermalization. Identifying a microscopic species would also require a gauge normalization, a parent abundance and a consistent gravitational scale.

## Finite-time continuum check

For a distinct weak-decay test, inject one parent with no daughters at $t_b=9$ and project onto daughter momentum states at turnaround. We prescribe an external smooth test envelope $s(t)=S[(t-t_b)/d]S[(t_f-t)/d]$, $d=0.45$, which multiplies the Yukawa interaction only in this test. The finite endpoints make the measurable protocol explicit. Daughters are treated as freely redshifting local adiabatic plane waves during the coupling. Local spinor spin sums, one parent at rest, and physical normalization of its spatial wavefunction give the following leading-order probability,

$$
\frac{P^{(2)}}{g^2}=\frac1{2\pi^2}\int_0^\infty dk\ k^2
\left|\int_{t_b}^{t_f} dt\,
\frac{s(t)p(t)}{a(t)^{3/2}E(t)\sqrt{M(t)}}
\exp\left[i\int_{t_b}^{t}(2E-M)\,dt'\right]\right|^2 ,
\quad p=k/a,\qquad E=\sqrt{\mu^2+p^2}.
$$

The absolute square retains interference between emission times. In a time-independent background the long-time replacement $|\int dt\,e^{i(2E-M)t}|^2/T\to2\pi\delta(2E-M)$ gives the preceding tree width, including its coefficient. For the finite protocol we compare against $P_\mathrm{local}^{(2)}=\int s^2\Gamma\,dt$, before exponentiating. At $g=0.01$ all probabilities in the table are below $10^{-4}$, keeping parent depletion negligible.

| Daughter mass $\mu$ | $P^{(2)}/g^2$, varying parent and scale factor | Coherent / local probability | $P^{(2)}/g^2$, fixed $a$ and subthreshold $R=0.95$ |
|---:|---:|---:|---:|
| 8 | 0.06776 | 2.157 | 0.02761 |
| 32 | 0.13496 | 1.074 | 0.00530 |
| 128 | 0.50496 | 1.005 | 0.0000993 |

The static subthreshold control has zero infinite-time vacuum width. Its finite response depends on work and spectral broadening from switching the external interaction. At $\mu=32$ halving the switching time raises that response from 0.00530 to 0.01375 in units of $g^2$. These transient amplitudes cannot establish an on-shell decay below threshold. At $\mu=128$ the varying-background finite-window probability agrees with the first-order local prediction within 0.5%. The code checks time and momentum refinement, higher momentum cutoff, and a fixed above-threshold mass. This continuum test neglects Dirac mode mixing from the expanding metric and background gauge field during emission. It treats the correspondence to the local width in the heavy, weak-field regime and does not solve the quantum heating problem.

Run `python -I -B survival/local_decay.py --output /tmp/SURVIVAL.json` from the capsule root. Supply `--source PATH` to locate another copy of the pinned input. NumPy and SciPy are required. [Rate and energy calculation](local_decay.py) · [Finite-time amplitudes](finite_time.py) · [Numerical output](SURVIVAL.json)
