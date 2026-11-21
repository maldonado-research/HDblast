# What happens when the electric current first reverses?

When the fluid's charge current passes through zero, how much momentum is still carried by particles that have never collided? In **both** fixed-pole wall histories, there are too few collision incidences by this point to produce the pressure assigned by the fluid. Even giving every collision to the particles that would otherwise carry the largest axial stress fails this test. [Machine-readable results](REVERSAL.json)

| At the first current reversal | Dedicated vector pole | Shared vector, inclusive channel |
|---|---:|---:|
| Wall time in units of $L$, native run | 6.57480 | 6.57264 |
| Fraction of whole-run absolute electrical work completed | 53.21% | 39.25% |
| Collision incidences per particle predicted by the thermal rate | 0.425 | 0.221 |
| Smallest fraction of particles that must have collided to fit the fluid axial stress | 60.2% | 67.1% |
| Minimum stress carried by the unhit population, divided by the entire fluid axial stress | 1.99 | 4.86 |

Consider 100 particles already born by the reversal. In the dedicated case the thermal collision calculation supplies about 42.5 participations. Repeated hits would reduce the number of *different* particles touched. Yet it takes at least 60.2 different particles to have been deflected to achieve the reference pressure. In the shared case the same comparison is about 22.1 participations versus 67.1 particles that must have been touched. Rounding to one decimal place is for readability.

We define the event as the **first change from negative to positive in the positive-charge stream's mean axial velocity after particles are present**. This is a sign change of the charge-weighted current because the opposite-sign stream is its reflection. The instantaneous electric field and the cumulative electric work have different zeroes. When the axial mean crosses zero, individual particles can still be moving rapidly along the axis. An isotropic rest-frame fluid assigns $P_{zz}/(n_{\rm lab}m)=\Theta$ at this event. The stress check permits the full range of flow and temperature in the two native endpoints surrounding the zero.

## The count and the stress

The event budget to the marker is evaluated independently of the thermal interpolation table. At each integration time the single-particle thermal event coefficient $C_{\rm event}(\Theta,u)$ comes from direct invariant-mass quadrature. If $D$ is cumulative comoving rest mass summed over the six sign/axis populations, the mean incidence count for *all particles present at the marker* is

$$
\bar N_{\rm hit}(t_*)=
\frac{1}{D(t_*)}\int^{t_*} D(t)\,
\frac{b_{\rm ch} S_0 D(t)}
 {6(mL)^3a^3\sqrt{1+u^2}}\, C_{\rm event}(\Theta,u)\,dt .
$$

The factor $b_{\rm ch}$ is one for the dedicated pole and one-third for the shared inclusive pole. The six populations are equivalent under their axis and charge symmetries. Each physical encounter counts once per participant. This convention permits direct comparison with the fraction of particles that could ever have been hit. [Thermal event integral](../history/RESULTS.md)

For an unhit birth, the recorded electric impulse and expansion enclose its axial momentum $p_z/m$ at reversal for every possible birth direction and time in its native birth cell. If $z_0$ is the smallest allowed $|p_z|/m$ in that cell, the contribution to the axial stress is at least $z_0^2/\sqrt{1+z_0^2+q_0^2}$ with $q_0=\sqrt{1.05^2-1}$. This uses the maximum possible transverse momentum of a freely moving birth. Count at most one new distinct participant per incidence, sort the cohorts by this cost and keep only the cheapest compulsory unhit fraction. Every other particle may contribute zero stress. This gives the lower bounds in the table. The sorted cumulative costs can also be inverted to find the smallest number that must have been deflected to meet the reference stress. [Characteristic and sorting derivation](STRESS.md)

The shared case permits an independent check without integrating to the zero. Its first reversal occurs **before** half of the total absolute electrical work. Grant the system all the future incidences from that later half-work stage and divide by the earlier, smaller particle population. Even this deliberately excessive allowance is below 0.484 incidences per particle. The minimum unhit stress still exceeds **2.35 times** the full fluid stress at reversal. [Earlier encounter budgets](../history/HISTORY.json)

## Numerical and physical scope

The archive includes native and halved wall time steps, two time and invariant-mass quadratures for the event count, the precise native bracket of each current zero and a comparison to the archived half-work count. The stress calculation widens its native endpoint ranges conservatively. A further sensitivity run enlarges each birth momentum interval by $\pm500m$ and increases the total collision budget by 1%. The halved wall step gives remaining stress factors of 1.990 and 4.861. With both additional padding choices the factors still exceed 1.93 and 4.78, respectively. The two temporal and invariant-mass quadratures agree within 0.0013% on the reversal counts and reproduce the archived half-work counts within 0.00071%.

This is a consistency test of the specified history, thermal rate, axial pressure, birth spectrum and locally parallel fields together. It identifies a time when each fluid description requires more redirection than its own collision count allows. A kinetic solution can develop a different velocity distribution, collision rate and electric current. Subsequent expansion, heat, and the cosmological origin of the source require separate calculations.

Reproduce with `python -I -B ballistic/reversal.py --traces /path/to/native --history /path/to/HISTORY.json --output /tmp/REVERSAL.json --jobs 2`.
