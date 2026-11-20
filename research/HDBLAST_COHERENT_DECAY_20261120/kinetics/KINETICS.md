# Can the charged particles occupy the available states?

Only one magnetic branch of the instantaneous Dirac spectrum can receive sharp, energy-conserving births throughout the measured 1%–99% production interval of the successful primary run. In the part delivering the 75th through 90th percentiles of the particles, the fluid requires about 34,000 times the available states if the wall action is just large enough for one elastic channel to provide the assumed hot drag. Occupied birth states would have to be emptied repeatedly. The fluid evolution prescribes the daughter source and its drag and therefore cannot test this requirement on its own. [Numerical results](KINETICS.json)

## Fermion-state count on the round wall

Assume three Dirac flavors, each charged under one of the three Abelian fields. Each flavor has particles and antiparticles, daughter mass $\mu=mL=400$, and classical $Q=q/(\sqrt3\mu)=2\times10^8$. Use natural units with $S_0=\sigma L^4$. The homogeneous field has potential $A=a^2 B/2$ in a unit invariant direction. In a temporal gauge, the instantaneous spatial two-spinor operator on a unit three-sphere is

$$
\mathcal D_b=2\boldsymbol\sigma\cdot\boldsymbol J+3/2-b\sigma_3,
\qquad b=\mu Q A.
$$

The $3/2$ is the unit-sphere spin connection. Reversing the sign of the charge preserves the spectrum. Put $\ell=2J+1$. The two unpaired extremal states have eigenvalues $\ell+1/2\pm |b|$. The remaining two-by-two blocks have eigenvalues

$$
\lambda_{\ell j,\pm}=\tfrac12\pm\sqrt{(\ell-|b|)^2+4j|b|},
\quad 1\le j\le\ell-1.
$$

Each eigenvalue in this list carries multiplicity $\ell$ from the other commuting $SU(2)$ action. For instance, the block on $|J,m,\uparrow\rangle,|J,m+1,\downarrow\rangle$ has diagonal elements $2m+3/2-b$ and $-2(m+1)+3/2+b$ and off-diagonal element $2\sqrt{(J-m)(J+m+1)}$. Restoring the four-component Dirac Hamiltonian gives one positive-energy state per intrinsic eigenspinor for each charge, with energy $E/m=\sqrt{1+(\lambda/a\mu)^2}$. At $b=0$ this gives $2(n+1)(n+2)$ positive-energy states per sign and flavor at $|\lambda|=n+3/2$, including spin. The code checks direct spin-block spectra against the integer counting formula for 25 field and cutoff combinations. [Round-space Dirac spectrum](https://doi.org/10.1007/BF01199016) · [Calculation](kinetics.py)

Consider a local decay at rest with measured $R=M/(2m)$. Even counting all states with $E\le M/2$ only admits $|\lambda|\le k_b=a\mu\sqrt{R^2-1}$. The first paired branch has absolute eigenvalue at least $\sqrt{4|b|}-1/2$. On the measured background this bound exceeds the largest birth cutoff in each reported window by a factor of at least 56.5. Also $|b|$ exceeds $k_b+1$ throughout. Thus only the extremal branch $\ell+1/2-|b|$ contributes. In a homogeneous time-dependent Abelian background an extremal state keeps its representation labels. Electric acceleration changes its physical momentum through $b(t)$. The count therefore includes the new birth states swept out by the electric field, the changing radius, and the moving mass threshold. A microscopic scalar decay could impose selection rules that further reduce the rate into these states.

For each production window, let $b_-,b_+$ bound $|b|$ and $k_+$ bound $k_b$. Counting every integer from $\lceil b_- -1/2-k_+\rceil$ through $\lfloor b_+ -1/2+k_+\rfloor$ gives a generous candidate set. Its total size is $N_{\rm cap}=6\sum\ell$. The factor six counts three flavors and both charges. The gas rest energy per unit comoving coordinate volume is $D_c$ in the coupled output; the particle number added over the whole three-sphere is $\Delta N=2\pi^2 S_0\Delta D_c/\mu$. The archive uses $\dot D_c=a^3\Gamma\dot\chi^2/R$. The comparison uses its integrated number directly, so rapid variation in the instantaneous rate is retained. [Coupled model](../coupled/joined_forward.cpp)

| Share of final number produced | Proper wall time interval | Collisionless upper $S_0$ | Minimum label clearance per $L$ at the hot drag floor |
|---|---|---:|---:|
| 1%–10% | 6.254–6.334 | $1.07\times10^{25}$ | 0 |
| 10%–25% | 6.334–6.487 | $1.57\times10^{24}$ | 0 |
| 25%–50% | 6.487–6.827 | $1.11\times10^{23}$ | 16.7 |
| 50%–75% | 6.827–7.410 | $2.50\times10^{21}$ | 505 |
| 75%–90% | 7.410–8.180 | $2.18\times10^{19}$ | 44,100 |
| 90%–99% | 8.180–10.11 | $4.06\times10^{19}$ | 9,415 |

Every accepted native step intersecting a bin contributes both endpoint values, including the exterior bracketing endpoints. The recorded range of $|b|$ is expanded on both ends by 1% of its span plus two maximum step changes. The upper birth cutoff is increased by 1%. These margins and a second grid test numerical sensitivity. The spectral counting statement is exact for the stated input intervals. The time-dependent enclosures are numerical estimates. The middle 98% is the measured domain of this diagnostic.

## Comparing particle supply with elastic drag

The independent collision test considers two oppositely drifting, dilute ultrarelativistic Maxwell–Boltzmann gas populations of a single flavor at a time. There is one opposite-charge target population. The temperature proxy is $T/m=\Theta$. During hot expansion only, select $\Theta\ge100$ and drift $|v|\le0.1$. The $s$-wave elastic unitarity bound $\sigma_T\le4\pi/k_{\rm cm}^2$ gives a thermal **momentum-drag** ceiling $L\nu\le\pi n_t L/T^2$. The proper target density in curvature units is $n_t L^3=S_0 D_c/(6\gamma a^3\mu)$, with $\gamma=(1-v^2)^{-1/2}$. This uses the relativistic dilute-gas kernel and assumes local hard scattering can be approximated by its field-free form. It describes transport among hot particles. The low-energy production states require the magnetic spectrum derived above. [Derivation of the drag ceiling](https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_NEUTRAL_PARENT_DECAY_20261120/consistency/CONSISTENCY.md)

The primary trajectory spends $32.88 L$ of wall time in this hot low-drift regime. Supplying $L\nu=18$ throughout requires $S_0\ge7.38\times10^{23}$. The most restrictive point has $T/m=14,252$, drift $-0.00190$, and $D_c/a^3=6.05\times10^{-7}$. Supplying just the integrated drag optical depth over the same interval requires $S_0\ge2.54\times10^{23}$. These are necessary conditions for this one-channel thermal benchmark. The elementary collision frequency for a particular newly created particle can differ substantially from the average hot momentum-relaxation rate.

If every occupied birth label can clear at most $L\nu_c$ times per unit dimensionless wall time, Pauli exclusion gives $\Delta N\le N_{\rm cap}(1+L\nu_c\Delta\tau)$. This even allows every candidate label to participate during the entire bin. With $L\nu_c=18$, the 75%–90% bin allows at most $S_0=3.23\times10^{20}$. The pointwise hot-drag requirement exceeds that upper bound by a factor of **2,284**. Even matching the integrated hot drag alone exceeds it by a factor of **786**. At the pointwise floor the same bin requires average label clearance of at least **$4.41\times10^4/L$**. The 64-cell control gives **$4.51\times10^4/L$**, a pointwise floor $7.56\times10^{23}$, and integrated floor $2.54\times10^{23}$. Both use the same physical interactions and six late cells.

This rules out identifying the single constant $18/L$ with both the hot drag and a ceiling on birth-state turnover in this one-channel, sharp-birth model. The classical source does not specify such a turnover ceiling. Fast scattering of newly born particles, a finite formation bandwidth, a different set of final states, and medium corrections to the fermion mass or dispersion each change the test. The magnetic branch restriction itself is a reason to recalculate the vacuum, three-dimensional decay width in the actual background. Observational conclusions require separate evidence.

From the capsule root, reproduce the primary and alternate measurements with

```sh
python -I -B kinetics/kinetics.py --study coupled/JOINED.json --output /tmp/KINETICS.json
```

The runner adds read-only diagnostics to a temporary copy of the supplied coupled solver, checks the solver hash against the archived data, and verifies that the particle number and hot-expansion growth reproduce. [Readout](readout.inc) · [Measured windows and counts](KINETICS.json)
