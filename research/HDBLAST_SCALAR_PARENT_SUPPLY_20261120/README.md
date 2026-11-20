# Can vacuum production supply the neutral parents assumed by HDblast?

HDblast's neutral-parent decay study starts with a cold spin-zero population on a closed wall. Its mass changes as the wall scalar evolves, eventually allowing it to decay to charged daughters. This calculation asks whether that same mass history, acting on a quantum field initially in vacuum, could generate the assumed number of parents by the opening of the decay channel. It borrows the TOE project's method of comparing a varying mass with a constant-mass control on the same geometry. The quantum equation and state counting here are those of a real scalar. [HDblast decay study](../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/README.md) · [TOE mass-pulse calculation](https://github.com/maldonado-research/Unified-Theory-of-Everything/blob/d9fa25bbca056135c7e808ca006ac2807a9190d5/checkpoint_2026_09_06/dynamic_eft_review.md)

**For one free conformally coupled real scalar and a daughter mass $mL=1$, fewer than 525 parent quanta are possible in expectation on the whole sphere at the opening**, in the specified instantaneous particle basis. The original fluid's initial rest-density and charge normalization imply over **1.50 billion** parents for a late effective gauge coupling below one. Even the conservative upper bound falls short by a factor greater than **2.86 million**. This compares inventories at the point where decay becomes possible. The precursor geometry already includes the original cold matter. A replacement history beginning without those parents would require a new gravitational solution. [All-mode certificate](all_shell_bound.py) · [Fluid normalization](../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/coupled/DECAY.md)

## Field, geometry and counts

Use natural units, curvature radius $L$, dimensionless wall proper time $t$ and dimensionless radius $a$. Take a neutral real scalar with curvature coupling $\xi=1/6$, vacuum initial data at $t=0$, and no interactions or decay during propagation. Write $\mu=mL$ for the **daughter** mass and $ML=2\mu R(\psi)$ for the parent. The principal HDblast mass ratio rises smoothly from $0.95$ to $1.05$, reaching $R=1$ when $\psi=0.77$. The trial uses $\mu=1$ to compare directly with the earlier TOE wall calculation. HDblast leaves this scale free. Both comparisons use the neutral-precursor trajectory, a fixed reference solution with the original seed. The interpolation of $a$ and $\psi$ uses the same Hermite values and slopes as the published survival calculation. The constant-mass control freezes $ML$ at its initial value of 1.9 for the whole evolution. [Mass function, input and timing](../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/survival/local_decay.py)

For scalar harmonics of degree $n\ge0$ on the unit round three-sphere, set $j=n+1$. There are $j^2$ real harmonic oscillators per shell. Conformal time obeys $d\eta=dt/a$. The rescaled oscillator and its vacuum initial conditions satisfy

$$
u_j''+\omega_j^2 u_j=0,\qquad \omega_j^2=j^2+c^2,\qquad c(t)=a(t)M(t)L,
\qquad u_j(0)=1/\sqrt{2\omega_j(0)},\quad u_j'(0)=-i\omega_j(0)u_j(0).
$$

Define the instantaneous expected occupation $N_j=|\omega_j u_j-i u'_j|^2/(2\omega_j)$. The expected number over the sphere is $\sum_j j^2N_j$. The real scalar carries no independent antiparticle copy. At proper time $t=10.34058195$ the first channel opens on the specified trajectory. Results at that instant are

| Parent mass history | Homogeneous mode ($n=0$) | All modes through $n=32$ | Through $n=64$ | Through $n=128$ |
|---|---:|---:|---:|---:|
| HDblast $2R(\psi)$ | 0.00008641 | 1.22177 | 1.63683 | 1.68820 |
| Constant $ML=1.9$ | 0.00100134 | 0.03639 | 0.03823 | 0.03857 |

The changing mass raises the resolved total while lowering the occupation of the homogeneous mode. The latter is the spatial mode closest to the assumed parent at rest. An independent adaptive evolution in the instantaneous basis reproduced all 24 reported counts within $2.54\times10^{-7}$. Each entry is a finite-time expectation. Counts at later times in the receipt continue the free field for comparison and omit decays, inverse reactions and backreaction. [Propagator](parent_modes.py) · [Numerical results](MODES.json)

## Bounding the contribution from every mode

For a vacuum-initialized oscillator the Bogoliubov coefficients obey $|\alpha|^2-|\beta|^2=1$ and $|d|\beta|/d\eta|\le |\omega'|\,|\alpha|/(2\omega)$. Consequently $r_j=\operatorname{asinh}|\beta_j|$ is at most one half the total variation of $\log\omega_j$. On the archived interpolant $c=aML$ first decreases to a minimum and then increases up to the opening. The checker verifies the sign of the derivatives of the Hermite polynomials and accounts for the smooth mass change. The gate is constant during contraction. It verifies $c_i=3.8$, $c_b>3.5$ at the minimum, and $c_f<19.2$ at opening. Hence

$$
r_j\le \tfrac14\log\frac{(j^2+3.8^2)(j^2+19.2^2)}{(j^2+3.5^2)^2}.
$$

Summing $j^2\sinh^2 r_j$ with outward-rounded rational bounds through $j=1000$ gives less than 516.908. For larger $j$, the inequality $\log(1+x)\le x$ bounds $r_j\le D/j^2$ with $D=(3.8^2+19.2^2-2\cdot3.5^2)/4$. Using $\sinh x\le x/(1-x)$ for $0\le x<1$ and $\sum_{j>J}j^{-2}\le1/J$ bounds the entire remaining number by 8.038. The resulting 525 ceiling needs no numerical mode integration or assumption about cancellation of phases. Extending the sign and magnitude checks on the interpolant through turnaround gives a corresponding bound of 699 through that epoch for a freely propagated field. These bounds concern occupations on slices. An interacting calculation could deplete and repopulate the field between slices. It applies to the prescribed Hermite curve through the pinned data. [Bound checker](all_shell_bound.py) · [Exact-arithmetic bounds](BOUND.json)

## Comparing with the parent supply

The fluid starts with parent rest density $D_0=10^{-4}$ in wall-tension units, $a_i=2$, $M_i L=1.9$ and wall action $S_0=\sigma L^4$ in natural units. Its initial comoving number is $N_{P,i}=2\pi^2a_i^3D_0 S_0/(M_iL)=0.0083112458\,S_0$. The charge requirement concerns its eventual daughters. Keeping $Q=|q|/(\sqrt3\,\mu)=10^8$ at $\mu=1$ and requiring $\alpha=q^2/(4\pi S_0K)<1$ at the late $K=13200.04$ gives $S_0>1.80857\times10^{11}$ and $N_{P,i}>1.50314\times10^9$. This comparison uses the surface-dominated local gauge normalization. Both tests are conditional on the specified particle content and mass scale. A pre-existing parent population, different coupling to curvature, direct interactions, different masses or a different expansion history can alter the production test. The calculation does not establish an absolute particle yield for HDblast or a temperature for its plasma. [HDblast microscopic assumptions](../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/consistency/CONSISTENCY.md)

## Reproduction

From this directory, use the archived input distributed in the preceding study. The mode solver uses NumPy and SciPy, checks the oscillator Wronskian and halves its time step. The bound checker uses the Python standard library. A second solver uses adaptive DOP853 transport of the Bogoliubov coefficients and independently constructs the reference history in the Bernstein basis. All three check the input digest.

```sh
python parent_modes.py ../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/prior/CHARGED.json
python all_shell_bound.py ../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/prior/CHARGED.json
OPENBLAS_NUM_THREADS=1 python control.py ../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/prior/CHARGED.json
```

The numerical receipt records both mass histories and three readout times. Halving the maximum time step from 0.00025 to 0.000125 changed the recorded counts by less than $7.6\times10^{-7}$. Source trajectory and mass-law data © 2025–2026 Ricardo Maldonado, distributed under CC BY 4.0. This calculation reads those data and computes the additional mode histories and bounds. [Data license](../../LICENSE) · [Archived HDblast study](https://doi.org/10.5281/zenodo.24406309)
