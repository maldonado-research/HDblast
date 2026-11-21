# Re-evolving the wall with restricted particle births

A daughter pair that cannot form leaves energy in its neutral parent. Keeping this energy available changes later decay opportunities, and the delayed charges can draw a different amount of energy from the electric field. Here each rate law is inserted into the time evolution of the parent, gas, electric field and expanding wall. Both tested masses retain a 90%-hot phase through a radius increase of about **2.25**. This conclusion applies to the modeled energy flow. The collision rate, actual occupations and angular distribution of previously created particles still need a kinetic evolution. [Live runs](FEEDBACK.json) · [Earlier coupled model](https://github.com/maldonado-research/HDblast/pull/192)

## Four local source laws

The parent is resolved while its minimum is driven, and its subsequent free oscillations are cycle averaged. For radius $a$, daughter mass $\mu=mL$, mass ratio $R=M/(2m)$ and $z=\sqrt{R^2-1}$, the birth shell in comoving momentum space has radius $k=a\mu z$. Let $P=\Gamma\rho_\chi$ be the effective parent power before applying an angular constraint, and $D$ the comoving daughter rest energy. The unconstrained rest source is $\dot D_{\rm try}=a^3 P/R$. With $u$ the positive daughter's cosine relative to the chosen fixed electric axis, define the normal speeds in units of $a\mu$ by
$$
c=Hz+R\dot R/z,\qquad d=QX,\qquad
s=\frac{\pi^2 S_0 P}{6\mu^4 R z^2},\qquad w(u)=|c-du|.
$$
Here $w$ measures swept state capacity in a given direction and $s$ is the requested speed for equal emission per solid angle. The formula applies after the kinematic threshold opens. The opposite charge at the reflected angle has the same normal speed. Write $\langle F\rangle=\frac12\int_{-1}^1F(u)\,du$. We compare four acceptance factors $A(u)$.

| Run | Acceptance | Assumption |
|---|---|---|
| Unlimited | $1$ | Original fluid law |
| Angle average | $\min(1,\langle w\rangle/s)$ | Pooling directions gives a generous total cap, with zero birth dipole |
| Directional cap | $\min(1,w/s)$ | Uniform trial emission independently fills each direction up to its fresh-state capacity |
| Joint vacancy | $w/(w+s)$ | Initially empty states meet the shell once, with equal and independent probabilities that each member of a pair finds a vacancy |

The last law follows from $df/dh=(1-f)^2$, $f(0)=0$, over an attempted filling $h=s/w$. The solution $f=h/(1+h)$ accepts the fraction $f/h=1/(1+h)$. It is a local comparison prescription. Correlations between the two daughters, a residual occupation when the same state returns, collisions and finite energy width can all change it. [Implemented angular integrals](angular.hpp) · [Directional counting conventions](../science/RESULTS.md)

## Energy and directed momentum

Each trial transfer first acquires a factor $A_0=\langle A\rangle$. The parent action $I=a^3\rho_\chi/M$ decays as $\dot I=-A_0\Gamma I$ after the drive. Only $A_0a^3P$ enters the comoving daughter energy, and $\dot D=A_0\dot D_{\rm try}$. With $A_1=\langle u A(u)\rangle$, newborns inject directed comoving momentum at the rate $J_{\rm birth}=\dot D_{\rm try}z A_1$ in the sign-oriented two-stream variable. Opposite charges have opposite physical momenta. In the code the gas energy above $\sqrt{D^2+S^2}$ receives the term $-v_c J_{\rm birth}$, where $S$ is directed momentum and $v_c=S/\sqrt{D^2+S^2}$. This accounts for the momentum added to that streaming energy at fixed total birth energy. The gas rest-frame mixing ledger also includes the corresponding $-\gamma v J_{\rm birth}$ term. Charge abundance and drift then determine the electric current and its work. [Live solver and ledgers](feedback.py)

Here “refused” means integrated trial rest energy less accepted rest energy, divided by all trial rest energy over the run. Energy retained by the parent can be tried again later, so this quantity differs from a reduction in the final particle count. The inputs on the two rows of each mass are the same as the frozen-history examples. The half-step reductions for the two directional laws and the native controls give the following results.

| $mL$ | Law | Refused trial rest, native / half step | Hot growth factor, native / half step |
|---|---|---:|---:|
| $4\times10^{10}$ | Unlimited | 0% | 2.25282 |
| $4\times10^{10}$ | Angle average | 0.10543% | 2.25282 |
| $4\times10^{10}$ | Directional cap | 0.83039% / 0.83028% | 2.25282 / 2.25231 |
| $4\times10^{10}$ | Joint vacancy | 5.05617% / 5.05629% | 2.25277 / 2.25226 |
| $4\times10^{12}$ | Unlimited | 0% | 2.25286 |
| $4\times10^{12}$ | Angle average | 0% | 2.25286 |
| $4\times10^{12}$ | Directional cap | 0.002923% / 0.002865% | 2.25286 / 2.25241 |
| $4\times10^{12}$ | Joint vacancy | 0.04335% / 0.04340% | 2.25286 / 2.25241 |

At the higher mass, joint vacancy leaves about 97.19% of the wall's added energy in intrinsic random particle motion when the hot radius has doubled. Its final comoving rest energy is 0.00212786. At the smaller mass, even the 5.06% refused fraction leaves a final rest near 0.00212786 because later attempts supply the energy. The finite proper-time growth ends at a maximum radius near 34.6, followed by contraction. The higher-mass signed birth impulse is very small; for the joint-vacancy case its net sign changes with step refinement, so the sign of that integrated impulse is unresolved. The postprocessed calculation resolves short electric reversals inside native steps and supplies the more accurate frozen-history angular peak. The live rates use native Runge–Kutta stages. [Direction- and time-resolved peak](../science/RESULTS.md)

Over all 12 runs the maximum absolute errors in the gas, oscillator and total-surface energy ledgers are $1.34\times10^{-10}$, $2.87\times10^{-9}$ and $1.64\times10^{-8}$ in code units. The Friedmann constraint residual stays below $4.90\times10^{-9}$. The unconstrained run reproduces the earlier hot-growth benchmark. These checks test conservation and reproducibility for a prescribed local model. The source remains an ansatz for an occupied, interacting quantum field; the gas pressure and effective drag still come from the fluid closure. Replacing the prescribed drag with a fixed elastic scattering model throughout the higher-mass evolution gives a 2.156-fold hot radius increase and a peak intrinsic random-energy fraction of 97.62%. [Explicit scattering result](../boost/RESULTS.md)

The astronomical implications require a separate heat and entropy budget. [Cosmological comparison](../cosmology/BUDGET.md)

From the study root run `python -I -B feedback/feedback.py --base prior --output /tmp/FEEDBACK.json --jobs 2`.
