# What changes when the pressure and the collision model describe the same particles?

The previous fluid calculation used a convenient pressure law, the Taub–Mathews equation of state (TM), and computed collisions assuming each particle stream had a relativistic Maxwell–Jüttner energy distribution. At very high temperatures these two descriptions give almost identical energy and pressure. Equal low-order moments alone still leave room for important differences in scattering. We tested two ways to make the descriptions internally consistent. Each run starts with empty particle sectors and evolves the pulse, wall, electric field, parent source, particles' fluid moments and gravitational equations together. [Predecessor](https://doi.org/10.5281/zenodo.24408924) · [Non-unique momentum transfer at matched moments](../moments/RESULTS.md)

| Description inside each stream | Hot radius factor, dedicated elastic channel | Hot radius factor, one pole shared among three flavors |
|---|---:|---:|
| Relativistic thermal distribution, exact Synge thermodynamics | 2.15583 | 2.06749 |
| Isotropic particles of one rest-frame energy, exact TM thermodynamics | 2.10717 | 1.99871 |

The table uses the archived spatial grid, halved wall time step, and a refined collision table for the single-energy shell. The peak intrinsic energy shares were 97.62% and 97.42% for the thermal examples and 97.51% and 97.25% for the shell examples (dedicated, then shared). Replacing the previous approximate thermal pressure law with the exact one changed growth by less than one part in a million. The larger change came from changing the spread of speeds used to compute the collisions. Both distributions are assigned locally at every time step. These runs test the effect of that assignment; the particle simulations instead let collisions reshape a population. [Full wall runs](FLOW.json) · [Distribution-resolved particle replay](../driven/RESULTS.md)

On the original mesh the shared single-energy run misses exact radius doubling by 0.065%. Increasing the number of early bulk radial elements from 24 to 48 gives factors of 2.0381 (shell) and 2.1085 (thermal). Increasing the late grid on its own from six to twelve elements changes either result by less than three parts in one hundred million in absolute growth. These checks apply to the specific finite pulse and the transfer between early and late grids. With 72, 96 and 144 early elements the shell results are 2.04218, 2.04270 and 2.04288. The corresponding thermal results are 2.11278, 2.11333 and 2.11352. Both assigned-distribution examples remain above doubling on the refined grids. A continuum error bound has yet to be derived. [Radial-grid tests](SPATIAL.json) · [Additional early resolutions](CONVERGENCE.json)

## Turning the recorded fluid state back into a particle distribution

Use particle mass $m=1$, rest-frame pressure parameter $\Theta=P/n$ and mean specific rest-frame energy $e$. For any isotropic distribution $e=\langle E\rangle$ and $e-3\Theta=\langle 1/E\rangle$. Cauchy–Schwarz gives $e(e-3\Theta)\ge1$, with equality precisely when all energies are the same. TM attains equality for $e=1.5\Theta+\sqrt{1+2.25\Theta^2}$. A natural distribution realizing this stress is uniform directions and momentum magnitude $k=\sqrt{3\Theta e}$. In a thermal gas $h=e+\Theta=K_3(1/\Theta)/K_2(1/\Theta)$. These are properties of relativistic kinetic theory. [Thermodynamic sources](../SOURCES.md)

Let $s_d=S/D$ be the tracked signed momentum per unit rest and $U/D$ the tracked excess above a cold beam with the same momentum. The Synge implementation recovers the temperature from

$$
U/D=\frac{h^2-1}{\sqrt{h^2+s_d^2}+\sqrt{1+s_d^2}}-\frac{\Theta h}{\sqrt{h^2+s_d^2}} .
$$

We used this difference formula to retain small internal energies in fast streams. Given $h$, the four-velocity is $u=s_d/h$ and $\gamma=\sqrt{1+u^2}$. We also converted the archived thermal force to momentum drag using the Synge enthalpy consistently. The independent primitive-recovery check spans cold, hot and boosted limits. [Implementation](synge.hpp) · [Checks](CHECK.json)

## Force for an isotropic momentum shell in each rest frame

Put $r=\operatorname{asinh}|u|$, $v=\tanh r$, and $t=\operatorname{asinh}k$. In the positive stream's rest frame, one incident negative particle of negative-stream rest-frame direction cosine $z$ has energy $e_2'=e\cosh(2r)-kz\sinh(2r)=\cosh t_2$. At fixed $z$, isotropy of the first particle makes $q=s-4$ uniform between $4\sinh^2[(t-t_2)/2]$ and $4\sinh^2[(t+t_2)/2]$. Write $B(q)=F(s)\sigma_T(s)$, with $F=\sqrt{s(s-4)}/2$, and $\overline B(z)$ for this uniform average. Elastic even-angle scattering removes half the incident pair momentum difference on average. Interchanging the beams changes the sign of that difference. The dimensionless wall-time relaxation coefficient is therefore

$$
C=\frac{1}{\gamma h e^2}\int_{-1}^1\frac{dz}{2}
      (e-kz/v)\overline B(z).
$$

At zero drift the removable singularity gives $C=[\overline B+2B(4k^2)]/(3 e h)$ with $q$ uniform on $[0,4k^2]$. This follows by differentiating the two-boost collision integral at $v=0$ and integrating the term in $dB/dq$ by parts. The cold beam yields $B(4u^2)/\gamma$. The solver uses the same fixed pole as the previous study. We evaluated the inner integral with a tabulated primitive and the outer integral by Gauss quadrature in $t_2$, resolving crossings of the resonance. A separate calculation directly integrates two independent emission directions and their relative azimuth. The largest sampled difference from the 72-point angular reference is $1.6\times10^{-5}$ relative. [Rate calculation](shell.hpp) · [Independent angular integration](check.py)

Reproduce from the study root with

```sh
python -I -B closure/flow.py --prior prior --output /tmp/FLOW.json
python -I -B closure/flow.py --spatial --prior prior --output /tmp/SPATIAL.json
python -I -B closure/flow.py --convergence --prior prior --output /tmp/CONVERGENCE.json
python -I -B closure/check.py --output /tmp/CHECK.json
```
