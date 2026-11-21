# Which changes in the particle distribution require an exchange of energy?

Imagine giving every particle a perfect steering wheel. Its direction can change by any angle, at any time, while it keeps its energy in the chosen frame. How close could this bring the particle counts to those expected in a hot gas? We tested the actual distributions at the two-polarization study's early work marker with both sideways particle currents included. Each charge is compared with its corresponding drifting thermal stream. All counts use one common wall-time slice. [Input study](https://github.com/maldonado-research/HDblast/pull/198)

## What the particle counts permit

| Fraction of particles, in percent | Separate interactions, fitted pulse | Shared interaction, fitted pulse | Separate, small pulse | Shared, small pulse |
|---|---:|---:|---:|---:|
| Initial discrepancy across the original 20 energy-direction bins | 47.16 | 55.70 | 61.91 | 61.74 |
| Smallest possible 20-bin discrepancy if energy stays fixed in the shell frame | 0 | 0 | 0 | 0 |
| Discrepancy that survives in 19 shell-frame energy bins | 18.97 | 19.66 | 19.24 | 19.70 |
| Smallest possible 20-bin discrepancy if energy stays fixed in the target stream frame | 14.97 | 21.55 | 16.11 | 16.77 |
| Discrepancy that survives in 19 target-stream-frame energy bins | 26.43 | 25.94 | 20.87 | 19.27 |

Entries use 16,384 equal-weight birth pairs per run. The small pulse has two transverse components with combined launch amplitude 4% of the reference and added emitted energy 0.16%. The fitted source cost about 1,356 and 492 times the original emitted energy. The discrepancy between two histograms is half the sum of the absolute count differences. Thus the 19 shell-frame energy bins alone distinguish the fitted populations from the target even if directions can be distributed to reproduce all 20 original counts. A mathematical ability to reproduce those 20 counts requires neither a physical scattering law nor momentum conservation for the scatterers. The table quantifies the most any external redirection could do under individual energy conservation. [All finite-bin counts and optimal bounds](ELASTIC.json)

Let the scattering center move at **any fixed speed parallel to the streams**, including a different speed for each charge. For each sign at least one energy threshold still gives a count discrepancy. Adding the two sign-resolved discrepancies with their original equal weights gives lower bounds of **8.75% and 9.18%** for the pressure-fitted pulse, or **6.75% and 7.05%** for the small pulse in the two histories. These are rounded down across the original and refined samples. The bounds cover every axial boost continuously. The threshold may depend on speed and sign. This is an obstruction for the empirical snapshots under an immediate redistribution in one fixed frame per charge. Binary collisions, scattering centers moving in several different frames, simultaneous electrical work, and scattering during birth can change individual energies in every single fixed frame. [Threshold witnesses for every speed interval](CONTINUUM.json)

## Thermal reference on the wall's time slice

Use units with daughter rest mass and light speed equal to one. Fold each species' axial momentum onto its positive target-stream direction. Write $r_t$ for the magnitude of the reference drift rapidity, $\beta=\tanh r_t$, $g=\gamma'$ for the energy in the target-stream frame and $p'=\sqrt{g^2-1}$. With $x=(g-1)/\theta$ and polar cosine $\mu$ in that frame, the target counting measure is

$$
dP = \frac{e^{-x}g p'}{2 e^{1/\theta}K_2(1/\theta)}
       (1+\beta p'\mu/g)\,dx\,d\mu .
$$

The angular tilt follows from the invariant momentum measure $d^3p/\gamma$. It converts the sample in stream coordinates to counts on the common wall-time slice. Reassigning particles to directions preserves these initial sampling weights. At another observation time the positions, the relevant hypersurface, or the thermal target may have changed.

## Exactly optimizing the 20 old counts

For a center at axial rapidity $r_s$ relative to the wall, put $d=r_t-r_s$. The energy preserved by one elastic turn is

$$
g_s=\cosh(d)\,g+\sinh(d)\,p'\mu .
$$

For the five tabulated frames $0\le r_s\le r_t$. Each of the 20 old bins is a connected rectangle in $(g,\mu)$, so its image under this equation is an interval. The two cosine endpoints give the extrema at fixed $g$. The minimum over $g$ can also occur at $g=1/\sqrt{1-\mu^2\tanh^2 d}$ when $\mu<0$. We test this point and the endpoints. Limiting faces of adjacent cells share their boundaries. The sampled energy cutoffs used here have zero empirical mass on a shared reachability endpoint. A particle can feed any bin whose image contains its energy. For a fixed energy in the scattering-center frame, the boosts and inverse boosts of these momenta lie on the same constant-energy sphere, and a rotation of the momentum can connect them.

If $S$ is a subset of bins, let $L(S)$ be the fraction of particles whose allowed destinations all lie in $S$, and $T(S)$ the desired thermal fraction in $S$. The sharp minimum discrepancy is $\max_S[L(S)-T(S)]$. This follows by capping the flow into each target bin at its target mass. The deficit of the maximum accepted flow is the minimum unmatched mass, and a minimum cut gives the subset expression. Any leftover mass can then go to a permitted bin. We evaluate all $2^{20}$ subsets and independently construct an allocation using a linear program. Both give the same answer. Turning individual particles according to randomized angular distributions realizes the result in expectation. The finite list of indivisible tags would have additional count granularity of one tag. In the target-stream frame all four directions of one energy band are available to every particle in that band. The optimal residual then equals exactly the five-bin energy-marginal discrepancy. The additional energy-spectrum measurement uses the fixed edges in the output JSON for $(g_s-1)/(\theta\cosh d)$.

## A bound covering all axial speeds

Let $b=\tanh(r_t-r_s)$ and use the signed axial momentum $p'_\parallel$ in the target-stream frame. Every subluminal axial frame lies in $-1<b<1$. For each particle the quantity

$$
Y_b=\frac{g_s}{\theta\cosh d}=\frac{g+b p'_\parallel}{\theta}
$$

is affine in $b$. On a whole speed interval, its minimum and maximum occur at the two endpoints. Counting samples whose maximum energy is below a threshold bounds the observed CDF from below everywhere on the interval. Counting those whose minimum energy is below the threshold bounds it above. The program uses 128 adjacent intervals spanning the closed interval from −1 to 1, and 96 positive energy thresholds. The endpoint speeds are included by continuity.

The comparison thermal probability is enclosed independently. In target-rest coordinates let $q=p'/\theta$ and $\epsilon=1/\theta$. The normalized radial density is proportional to $q^2 e^{-\sqrt{q^2+\epsilon^2}}$. Relative to the massless radial density $q^2e^{-q}/2$, their normalized ratio lies between $e^{-\epsilon}$ and $e^\epsilon$. Their total variation is bounded by $(e^\epsilon-1)/2$. The exact conditional angular density is $(1+\beta v'\mu)/2$ with $v'=q/\sqrt{q^2+\epsilon^2}$. Using $(1+\beta\mu)/2$ adds total variation at most $\epsilon/8$, since $1-v'\le\epsilon/q$ and $\mathbb{E}[1/q]=1/2$ under the massless radial law. The triangle inequality combines these bounds. The normalized energy differs from $q(1+b\mu)$ by $\sqrt{q^2+\epsilon^2}-q\in[0,\epsilon]$.

We partition $\mu$ into 1,024 intervals and integrate their angular weights exactly. The extrema of $1+b\mu$ over each speed-angle rectangle occur at its corners. The massless radial CDF is $1-e^{-q}(1+q+q^2/2)$. Its monotonicity, the shift in the energy threshold by at most $\epsilon$, and the total-variation allowance enclose the massive thermal CDF. We then maximize the guaranteed empirical-minus-thermal or thermal-minus-empirical CDF gap over the 96 thresholds in each interval and minimize over intervals. This method keeps the finite width of every interval and the nonzero particle mass. It supplies a conservative bound for the tested empirical population. The calculation neither fits a scattering angle nor assumes a particular microscopic collision rate.

The original 4,096-pair and refined 16,384-pair measurements repeat the count and frame comparisons. Their resolution agreement is an empirical check on sampling the prescribed birth distribution. The sample-specific certificates do not supply a rigorous sampling error for its continuum limit.

## Reproduction

The two scripts consume the charge-tagged NPZ files produced by the kinetic component from the archived two-polarization solver. NumPy and SciPy are required. The outputs contain all binned probabilities and the sign of each threshold witness.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B bound/elastic.py --snapshots inputs --output output/ELASTIC.json
OPENBLAS_NUM_THREADS=1 python -I -B bound/continuum.py --snapshots inputs --output output/CONTINUUM.json
~~~
