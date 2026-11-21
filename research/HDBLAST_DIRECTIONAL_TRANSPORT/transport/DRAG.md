# How much hot drag can the converting scalar supply?

The heaviest HDblast decay benchmark needs over **16.45%** of the reference elastic unitarity ceiling to deliver its prescribed hot drag while passing its published, angle-summed fermion supply test. An isolated spin-zero resonance coupled equally to three unpolarized Dirac flavors supplies at most **2.78%** in an elastic channel that preserves flavor. Including every allowed outgoing flavor gives a ceiling of **8.33%** for the charge-odd momentum relaxation described below. These inequalities hold at every incoming energy in the single-resonance model. They test whether the interaction already used to convert the parent can also supply the fluid's drag. [Archived state-supply limits](../../HDBLAST_COHERENT_DECAY_20261120/kinetics/SWEEP.md) · [Reference collision integral](../../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/consistency/CONSISTENCY.md)

## Spin, flavor and partial-wave unitarity

Consider a local, field-free binary collision of distinguishable particles $i,\bar i$ in their center-of-momentum frame, with momentum magnitude $k$. Assume negligible final-state occupation for the hot scattering, equal masses, three equal diagonal Yukawa couplings and one spin-zero neutral resonance. An isolated scalar pole forms one linear combination of the available spin and flavor channels. In a flux-normalized $J=0$ partial wave, write its scattering matrix as

$$
S-I=2iA(s)\,|v\rangle\langle v|,\qquad |A(s)|\le1,\qquad \|v\|=1 .
$$

For visible flavor $i$, let $B_i=\|P_i v\|^2$, summing its helicities. For this model $B_i=1/3$. A scalar at rest has a rotationally invariant decay distribution after summing the daughter spins. Hence the transport weight $1-\cos\theta_*$ gives $\sigma_T=\sigma$. Tracing the squared $J=0$ amplitude over both initial spins of each particle, with weight $1/4$, gives

$$
\frac{\sigma_{T,i\bar i\to j\bar j}}{4\pi/k^2}
  =\frac{2J+1}{(2s_i+1)(2s_{\bar i}+1)}\,|A|^2 B_i B_j
  \le\frac{B_i B_j}{4}.
$$

Only $i\bar i$ forms the diagonal resonance for an incident flavor $i$. Three outgoing choices share its width. For same-flavor elastic scattering the coefficient is at most $1/36$. Summing the three outgoing flavors gives at most $1/12$. The term $4\pi/k^2$ is the normalization used by the HDblast hot-gas screen. The parity-even fermion pair that couples to the scalar is $^{3}P_0$. Positive pair parity requires odd orbital angular momentum and total angular momentum zero requires $L=S=1$. Its magnetic degeneracy is $2J+1=1$. [PDG cross-section formula, Eq. 51.1](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-cross-section-formulae.pdf) · [PDG angular-momentum assignments, Table 15.2](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-quark-model.pdf)

As an independent normalization check, evaluate the archived spin-averaged tree resonance at $s=M^2$, with equal daughter masses $m$ and $\Gamma_i=y_i^2 M(1-4m^2/M^2)^{3/2}/(8\pi)$. Its cross section divided by $16\pi/(s-4m^2)$ gives $B_iB_j/4$ for both equal and unequal couplings. Quantitative saturation of the upper bound across a thermal distribution would require the largest possible partial-wave amplitude at every sampled energy. The proof uses the upper bound throughout. [Executable check](drag.py) · [HDblast spin-averaged cross section](../../HDBLAST_NEUTRAL_PARENT_DECAY_20261120/consistency/CONSISTENCY.md)

## Comparing with the fluid

The reference HDblast drag kernel uses counterstreaming dilute ultrarelativistic Boltzmann distributions of particles and antiparticles, small opposite drifts, common temperature $T$, and the density $n_t$ of opposite charges of the same flavor. Its integrand for $\nu_{\rm drag}$ weights $\sigma_T$ by the nonnegative factor $v_{\rm M}(p_{1z}-p_{2z})^2$. Thus a pointwise fraction of $4\pi/k^2$ gives the same upper fraction of the thermal drag ceiling. The inclusive statement concerns the charge-odd directed momentum, in the symmetric three-flavor setup. Each outgoing isotropic equal-mass pair has zero *mean* particle-minus-antiparticle four-momentum at fixed incoming total momentum. This holds in the fluid frame by linearity of the boost. Thus re-emission into a different flavor also erases at most the incident directed difference, just as isotropic elastic randomization does. The target density per incident flavor is still $n_t$; summing over the three incident flavors scales their total initial directed momentum and its relaxation equally. This moment calculation requires the stated symmetries. A general reacting mixture needs the full collision operator.

Let $S_{\rm drag}$ be the smallest action at which the reference ceiling can deliver $\nu L=18$ on the published hot slices. The required fraction of that ceiling at a given action is $S_{\rm drag}/S_0$. Let $S_{\rm sweep}$ be the maximum action passing the archived test of the angle-summed, instantaneous, sharp supply of initially vacant states in the central 99.998% of births. At $mL=4\times10^{12}$, the favorable endpoint $S_0=S_{\rm sweep}=1.26591\times10^{54}$ still requires a fraction 0.16455. The inclusive resonance would require $S_0\ge 12 S_{\rm drag}=2.49974\times10^{54}$; the elastic channel would require $S_0\ge36 S_{\rm drag}=7.49923\times10^{54}$. The minimum shortfall factors are **1.975** and **5.924**, respectively. The published example at $S_0=1.04156\times10^{54}$ requires 20% of the reference ceiling. [Both mass cases and action limits](DRAG.json)

Consequently one spin-zero pole cannot account for both prescribed rates within this particular sharp particle-source and local thermal-scattering interpretation. The limits on the pole follow from unitarity and the model's quantum numbers. The trajectory-dependent floors and ceilings come from numerical histories. Additional collisional channels, a different spin distribution, plasma forces, finite formation width and departures from the reference gas can all change this conclusion. Off-resonance crossed exchange and gauge processes are outside the single-pole proof. The fluid heating history has not been recomputed with these rates.

From the study root run the following command to reproduce the numerical comparison with the archived half-step sweep.

~~~sh
python -I -B -O transport/drag.py ../HDBLAST_COHERENT_DECAY_20261120/kinetics/SWEEP.json
~~~
