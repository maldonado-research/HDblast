# How much does the probe mass change the angular splitting?

A massive charged scalar can be excited by the expansion as well as the electromagnetic field. Repeat the finite-band crossing test for a massive conformally coupled complex scalar of charge $q=1$ and for a neutral twin of the same mass. The conformal frequency becomes

$$\omega^2=(n+1)^2-2kqA/\sqrt3+q^2A^2/3+\mu^2a^2.$$

Here $\mu$ is the mass in inverse bulk curvature lengths. Both species start in the instantaneous positive-frequency vacuum on the initial stationary slice, where $A=A'=0$ and $a'=0$. They use the same geometry and readout. A massive field has a changing conformal frequency even at $q=0$, and this initial-state prescription matters. At the actual crossing $A=0$, the three $n=2$ labels have the same instantaneous conformal energy $\sqrt{9+\mu^2 a^2}$. The neutral twin preserves that angular degeneracy in its occupation. The charged field need not.

| $\mu$ | Neutral occupation in each $n=2$ state | Largest to smallest charged occupation in $n=2$ | Total neutral excitation energy | Total charged excitation energy |
|---:|---:|---:|---:|---:|
| 0 | below $10^{-28}$ | 53,306 | below $10^{-22}$ | 0.05109 |
| $1/(\sqrt3\,10^8)$ | below $10^{-28}$ | 53,306 | below $10^{-22}$ | 0.05109 |
| 0.1 | 0.00001268 | 8.983 | 0.04415 | 0.09509 |
| 1 | 0.00085199 | 1.01491 | 1.30395 | 1.36042 |

The energies include particles and antiparticles with $0\le n\le32$ on the full three-sphere in conformal units. They use the instantaneous basis at wall time 10.48052 in the original trajectory. They exclude shorter wavelengths and any interpretation as total heat density. The second mass matches the illustrative charge-to-mass ratio of the coupled stream calculation at $q=1$. It shifts no tested occupation by more than $1.85\times10^{-18}$. This compares mass scales only. The quantum test starts in the vacuum and uses the unmodified electromagnetic background.

Increasing the mass allows the expansion to populate all directions, strongly reducing this angular diagnostic. At $\mu=1$ the largest and smallest occupations differ by only about 1.49%. Similarity within one angular shell does not establish a thermal distribution. Subtracting the neutral energy from the charged energy measures the change between two evolutions and includes interference between gravitational and electromagnetic excitation. It cannot be assigned as a separately counted population created by light.

From the study directory run

~~~sh
python -I -B quantum/mass.py --trace background/trace.csv.gz --output /tmp/MASS.json
~~~

[Full massless test and definitions](PAIRS.md) · [Program](mass.py) · [Results](MASS.json)
