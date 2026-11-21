# Can different birth directions and times use the same fermion state?

They can. A sphere of allowed birth momenta can cross a given momentum label on both its front and back surfaces. An electric reversal can produce further crossings. A state occupied at the first crossing can prevent a later birth. On the heavier of the two recorded trajectories, accounting for repeated crossings increases the estimated blocked fraction from **0.00290%** to **about 0.0058%** of the particles requested in the measured interval. The calculation assumes collisionless propagation between births, a sharp two-body energy, a locally continuous density of states and an isotropic source in the wall frame. The histories come from the previously published fluid model. [Earlier angle-integrated calculation](https://github.com/maldonado-research/HDblast/pull/192) · [Results and numerical controls](HISTORY.json)

## Which directions belong to a pair?

Write $k=a\mu\sqrt{R^2-1}$, $\mu=mL$, and choose the axis along the positive-charge field convention, with $\dot b=\mu Q aX$. The kinetic momenta of daughters born to a parent at rest point in directions $\boldsymbol n$ and $-\boldsymbol n$. The corresponding canonical labels in the local Weyl approximation are
$$
\boldsymbol K_+=k\boldsymbol n-b\hat z,\qquad
\boldsymbol K_-=-k\boldsymbol n+b\hat z=-\boldsymbol K_+ .
$$
Thus the opposite charges and their reversed outward normals give both partners the same signed normal sweep speed $v=\dot k-\dot b u$, $u=\boldsymbol n\cdot\hat z$. A change of axis reverses the sign assigned to the dipole. The acceptance of a physical pair tests the mirrored labels together. Treating the two speeds as $\dot k\pm\dot b u$ at the *same* physical polar angle would impose a different pair geometry.

The model has three equal charged flavors with both charges and two spin states. On the unit three-sphere the leading phase volume per spin and species is $2\pi^2/(2\pi)^3$ per $d^3K$. The available daughter rate per unit $du/2$ is therefore $12k^2|v|$ if both sides of the moving surface can supply empty states. The bare total requested daughter rate is $\dot N=2\pi^2 S_0\dot D_c/\mu$. Counting pairs divides both rates by two. For a vacuum scalar Yukawa decay at rest, the spin amplitude is proportional to $\boldsymbol\sigma\cdot\boldsymbol p$ up to a unitary convention on antiparticle spin. At nonzero momentum its two singular values are equal. The reduced spin state has equal weights and the direction rate is isotropic in that limit. The external field can change the matrix element, dispersion and allowed spins. Here both spins are granted.

Set $T=\dot N/(12k^2)$. The bare demand per crossed state is $h=T/|v|$. A sharp isotropic single-passage source loses a fraction $\tfrac12\int_{-1}^1[1-|v|/T]_+\,du$ even if every crossing finds empty states. For $|\dot b|>|\dot k|+T$ the missed fraction is exactly $T/(2|\dot b|)$. In the same regime the first moment of the accepted flux, divided by the entire requested rate, is $-(\dot k/\dot b)T/(2|\dot b|)$. Some directions can saturate even while the sum over all available directions exceeds the total demand. The dipole injects signed momentum into the two counterstreaming charges. Their total momentum cancels, and the second moments of identical axisymmetric distributions on three orthogonal axes sum to an isotropic tensor. This symmetry leaves the time-dependent distribution and its relaxation open. [Local direction-resolved calculation](../science/DIRECTIONAL.json)

## One ceiling covers every visit to a label

At fixed canonical cylindrical coordinates $(r,z)$, each allowed visit solves
$$
r^2+[z+b(t)]^2=k(t)^2,\qquad u_j=[z+b(t_j)]/k(t_j).
$$
The upper occupancy of one fermion state is one over all of its visits. Denote its total bare load by $H=\sum_j h_j$. Pushing the isotropic birth measure into canonical momentum space gives the density $H$ per state. If sources can always take the free part, the maximum accepted amount is $\min(1,H)$. The continuum loss per unit state volume decomposes as
$$
[H-1]_+=\sum_j[h_j-1]_+
       +\left[\sum_j\min(1,h_j)-1\right]_+ .
$$
The second term counts additional loss after granting each visit its best individual filling. The integration starts with empty states, including labels that might have been occupied before the measurement window. It gives an optimistic estimate conditional on the no-clearance prescription.

For the same states, two one-way rate models illustrate the difference between a capacity and a probability. An irreversible Poisson source whose spin-paired sector is either empty or occupied has $df/dH=1-f$, hence $f=1-e^{-H}$. A molecular-chaos closure that factors two equal final-state vacancies has $df/dH=(1-f)^2$, hence $f=H/(1+H)$. These functions describe different assumptions about correlations. Quantum phase interference, inverse transitions and cross-label collisions can change both.

A separate sensitivity calculation starts a paired label empty, accepts up to its vacant fraction at each visit, then allows occupation to decay by $e^{-18\Delta\tau}$ until the next visit, with zero in-scattering. Clearing both partners together is deliberately favorable. The number 18 is the fluid drag coefficient used in the predecessor. The actual rate at which collisions clear an individual state requires a microscopic collision kernel. The hard collisionless result applies to particles retaining their canonical labels throughout the measured interval. [Microscopic transport assumptions](https://github.com/maldonado-research/HDblast/tree/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120/kinetics)

## Delivery on the two measured backgrounds

All entries are percentages of the **requested daughter births inside the selected window**, from 0.001% to 99.999% of cumulative production. The source comes from the original, unblocked fluid history. The table shows results at halved evolution step, rounded to the precision supported by the time and angle quadratures.

| $mL$ and action choice | Single passage, angular capacity | Added loss from revisits | Total capacity loss | Hard cap with label clearing at 18 | One-way correlated loss | One-way independent-vacancy loss |
|---|---:|---:|---:|---:|---:|---:|
| $4\times10^{10}$, $S_0=3.609\times10^{48}$ | 0.8204% | 1.073% | 1.89% | 1.67% | about 6.4% | about 9.5% |
| $4\times10^{12}$, $S_0=1.042\times10^{54}$ | 0.002899% | 0.00291% | approximately 0.0058% | approximately 0.0058% | about 0.051% | about 0.092% |

The five quadrature and integration settings give total hard-cap losses of 1.89345% to 1.89384% in the first example and 0.005808% to 0.005821% in the second using the decomposed integral. The equivalent direct integral gives 1.89275% to 1.89373% and 0.005805% to 0.005839%. These are sensitivity ranges from deterministic numerical checks. Statistical or rigorous numerical confidence intervals would need a separate construction. Small missing birth fractions alone do not measure changes in gas temperature or wall expansion.

## Reconstructing the label occupation

The code uses native time steps containing $k,b,\dot b,D_c$ from the same archived solver that supplies the local angular test. It linearly interpolates $k$ and $D_c$ and uses a cubic Hermite interpolation of $b$ from its measured endpoints and derivatives. It subdivides rapidly changing electric slopes, including electric reversals, then solves the intersection equation exactly on each piecewise-linear $(b,k)$ segment. Scaling $b$ and $k$ by $\mu$ conditions the arithmetic. A tree of conservative interval bounds rejects segments unable to intersect the sampled label. A translated sphere can contribute two roots on one segment, which are counted separately. Grazing double roots have zero measure in momentum space. The sampled birth itself is retained exactly, and its companion root on that segment is obtained by factoring out the known root.

Time quantiles use a mixture of bare production and the analytically integrated single-passage deficit, with inverse weighting to restore the intended birth distribution. At each time the Gaussian angular integration is split at the grazing direction $u=\dot k/\dot b$ and at distances proportional to $T/|\dot b|$ on either side. A geometric progression of angular panels resolves the narrow loss and its tails. Evaluating the quadratic coefficients relative to the sampled birth radius and center limits cancellation near grazing. The scalar identity $\int dV\,[H-1]_+=\int dN\,[1-1/H]_+$ checks the sum of the separately measured single-passage and revisit terms. For the latter integral, distributing each label's extra loss among its visits in proportion to $\min(1,h_j)^2$ improves sampling of secondary crossings without changing the integral. The numerical controls vary the time quantile count from 8,192 to 32,768, the time offset, the angular Gaussian rule from four to eight points per panel, the electric-field subdivision and the native evolution step. Three manufactured histories test a single monotone passage, two radial passages whose integrated loss has a closed form, and two roots during one translation.

Small changes in the underlying trajectory can move a grazing crossing of another part of the history past an angular quadrature point. For this reason, portable comparisons of the integrated secondary-crossing terms use their resolved numerical precision. The continuum integral is conditional even when the discrete quadratic roots are found accurately.

The trajectory and canonical-label count use the leading continuum approximation to the charge-coupled Dirac spectrum on the round wall. The exact finite-volume state constraints, magnetic branches and the finite energy width of a real interaction can matter in thin angular bands. Neither the hard cap nor the rate examples calculate the true pair-production amplitude in the changing gauge field. The hot fluid's drag and possible scattering during production are independent assumptions. [Implementation](labels.cpp) · [Replay entry point](history.py)

With `PRIOR` set to the pinned #192 study directory, run from this capsule root.

```sh
python -I -B theory/history.py --prior "$PRIOR" --output /tmp/history.json
```
