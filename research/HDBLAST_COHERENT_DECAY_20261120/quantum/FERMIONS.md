# Finite-time response of fermions to the coherent parent

A time-local decay width assumes that a nearly monochromatic parent acts long enough for daughter energies to resolve the decay threshold. We tested this assumption by evolving Dirac modes on the archived wall. We use a single test species, zero electromagnetic charge, an undepleted parent and fixed geometry. This isolates the Yukawa interaction. The full charged model also contains a gauge background and a current. The physically matched mass modulation may be strong enough to invalidate a weak-coupling expansion. [Parent supply study](https://github.com/maldonado-research/HDblast/pull/191) · [Weak response](FERMIONS.json) · [High-mass response](NARROW.json)

## What is counted

For a canonically normalized physical parent $\Phi=\sqrt{\sigma}L\chi$ with Yukawa interaction $y\Phi\bar\Psi\Psi$, write $\lambda=y\sqrt{S_0}$, $S_0=\sigma L^4$ and $\mu=mL$. After rescaling the Dirac spinor by $a^{3/2}$, the proper wall-time equation for one particle–antiparticle slot is

$$
 i\frac{d u_n}{d\tau}=\left[\frac{n+3/2}{a}\sigma_x+
       (\mu+\lambda w\chi)\sigma_z\right]u_n .
$$

A round unit three-sphere has $2(n+1)(n+2)$ slots per massive Dirac shell, including both intrinsic eigenvalue signs. A homogeneous scalar preserves each slot. The following sums count one pair per occupied slot. No sign, flavor or antiparticle factor needs to be appended for this one Dirac species. [Dirac eigenspinors on spheres](https://arxiv.org/abs/gr-qc/9505009)

All runs start in the instantaneous negative-energy state. Let $u_0^-(\tau),u_0^+(\tau)$ be the negative and positive initial states both evolved with $\lambda=0$, and let $u_\lambda^-$ evolve with the same metric and the specified coupling. The interaction-induced amplitude in that transported basis obeys

$$
 \left.\partial_\lambda\langle u_0^+(T)|u_\lambda^-(T)\rangle\right|_0=
   -i\int w\chi\,\langle u_0^+(\tau)|\sigma_z|u_0^-(\tau)\rangle\,d\tau.
$$

Our weak pair response $C$ sums the squared absolute values of these amplitudes over the finite band with their multiplicities. This choice exactly cancels free propagation in the amplitude and has a nonnegative coefficient. The files also give symmetric finite-$\lambda$ projections at $\lambda=\pm1$, final probabilities in the common instantaneous energy basis, and the gravity-only count. The latter observables contain interference. For example, on the $\mu=1$ interval between natural zeros of $\chi$, the 512-shell instantaneous count is 3.59340 with gravity alone, 3.52682 with $\lambda=+1$ and 3.66220 with $\lambda=-1$. A detector convention matters when attributing an excess to the source.

For a parent mass $M$ and $\beta^2=[1-4\mu^2/M^2]_+$, the stationary one-flavor width is $y^2 M\beta^3/(8\pi)$. Multiplying by the instantaneous action count $2\pi^2S_0a^3\rho_\chi/M$ yields a local prediction per $\lambda^2$. Damping the coherent field by $\Gamma\dot\chi$ releases energy proportional to $\dot\chi^2$. Dividing that energy by $M$ per pair yields a second local comparison. Both are integrated over the same interval and weighted by $w^2$,

$$
 C_{\rm action}=\frac{\pi}{4}\int a^3 w^2\rho_\chi\beta^3d\tau,\qquad
 C_{\rm damping}=\frac{\pi}{4}\int a^3 w^2\dot\chi^2\beta^3d\tau.
$$

In a static continuum, a weak source proportional to $\cos(M\tau)$ reproduces these normalizations after cycle averaging. A finite static-sphere calculation with 2,304 shells and a Hann envelope gives $C/C_{\rm action}=1.00031$ from the exact Fourier transform. The time-evolved spinors reproduce that finite transform to 0.031%. [Mode solver and normalization test](fermions.py)

## Which comparisons support a local rate?

For the archived $\mu=1,16$ stable-parent histories, the short protocol runs from 0.5 wall-time units before the threshold to 0.5 after it. The long protocol starts at the same time and stops at the wall's first maximum. Smooth windows multiply the interaction by the product of compact ramps $S[(\tau-\tau_i)/d]S[(\tau_f-\tau)/d]$ with widths of 0.125 in the short test and 0.45 or 0.9 in the long test. The final protocol keeps $w=1$ between the first and last post-opening zero of the actual parent field before maximum size. Because $\chi$ vanishes at both endpoints, its Yukawa mass shift also vanishes there.

| Protocol | $C/C_{\rm damping}$ at $\mu=1$ | $C/C_{\rm damping}$ at $\mu=16$ |
| --- | ---: | ---: |
| Short threshold crossing | 341.62 | 4.706 |
| Through turnaround, 0.45 ramps | 5.512 | 1.040 |
| Through turnaround, 0.9 ramps | 3.481 | 1.030 |
| Between natural post-opening zeros | 2.593 | 1.080 |

The large short-window value needs a control. Moving an otherwise identical window entirely below threshold gives $C=0.006916$ at $\mu=1$, comparable to $0.007194$ in the crossing window. Both local on-shell integrals vanish in the closed control. Finite exposure, switching and the evolving geometry can all contribute spectral weight. Interpreting the factor of 342 as an intrinsic decay-rate enhancement would misidentify the observable.

Each displayed low- and intermediate-mass result includes the first 512 shells. For the smooth long windows, doubling the band changes $C$ by at most 0.013% for the tested 0.45 ramps. Between natural zeros the finite-band changes are larger, 0.73% and 1.59% respectively. Halving the integration step and comparing exact sign-symmetric evolution against its derivative test the numerical solution separately. The reconstructed parent oscillator matches its archive with maximum phase-space relative errors of $1.8\times10^{-8}$ at $\mu=1$ and $2.0\times10^{-6}$ at $\mu=16$. The calculation reconstructs the oscillations using the measured geometry and mass schedule, including where the archived heavy field advances by more than half a cycle between samples.

## The high-mass candidate

The narrow-source $\mu=400$ test uses an independently sampled continuation of the parent-only equations through a short post-threshold interval. Initially the source and its gauge/geometry history match the new joined preparation, whose final gauge increment is 2.5. In this continuation the parent survives and the daughter current vanishes. The source fields, radius and proper velocities are sampled every 0.0002 wall-time units. We also ran a stable continuation of the previously published preparation with final gauge increment 26. The two histories permit a direct comparison of the preparation dependence. [Generation of the histories](narrow_trace.py) · [High-mass mode solver](narrow.py)

The matched stable continuation opens at $\tau_*=6.24180457$, the same time as the joined primary. The short window spans $\tau_*\pm0.035$ with ramps of 0.0125. The longer window spans $\tau_*-0.12$ through $\tau_*+0.50$ with ramps of 0.045 or 0.085.

| High-mass history | Short crossing $C/C_{\rm damping}$ | Longer interval $C/C_{\rm damping}$ | Longer interval with gentler ramps |
| --- | ---: | ---: | ---: |
| Matched preparation, final $K=2.54$ | 2.2702 | 1.0101 | 0.9982 |
| Earlier preparation, final $K=26.04$ | 2.1161 | 1.0079 | — |

The matched short-window count is $2.670\times10^{-6}$ pairs per $\lambda^2$ in 12,288 shells. The longer window gives $2.355\times10^{-5}$ in 6,144 shells. Doubling the respective bands changes them by 0.00991% and 0.000365%. Halving the step changes them by 0.0948% and 0.0696%. Exact unitary evolution at $|\lambda|=1$ and its derivative agree to the precision displayed. The maximum spinor norm error is below $5.3\times10^{-14}$.

An equal-length subthreshold window again has finite response and zero local decay rate. In the longer high-mass protocol, the final mass excess $M-2\mu$ approaches 40 and the post-threshold duration is 0.50, or about 20 radians of threshold resolving phase at the final excess. Suitably long weak off-shell tests can approach the local law. This agreement does not validate the microscopic charged, finite-coupling evolution of the joined wall. Electric production, fermion occupation and depletion, and the feedback of all populated modes on the fields and geometry would be needed to predict an absolute energy transfer. The finite-band sums reported here apply to the chosen time intervals and initial quantum states.

Run the fermion calculations from the study root with Python, NumPy and SciPy.

~~~sh
python -I -B quantum/fermions.py --output /tmp/FERMIONS.json
python -I -B quantum/narrow.py --output /tmp/NARROW.json
~~~

Both input histories are included. The following commands regenerate them from the pinned files in this study. Generating the narrow histories requires a C++17 compiler. The output sampling interval is 0.0002 in proper wall time. [Archived parent study](https://github.com/maldonado-research/HDblast/tree/dfdbad8b4859d331fab1dcfa4a09051d2113d906/research/HDBLAST_DRIVEN_PARENT_SUPPLY_20261120/coupled)

~~~sh
python -I -B quantum/extract.py prior/PARENT.json.gz --output /tmp/HISTORY.json
python -I -B quantum/narrow_trace.py prior --output /tmp/NARROW_HISTORY.json
~~~
