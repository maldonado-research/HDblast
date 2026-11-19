# Energy flux and directional pressure of the outgoing pairs

For a single traveling parent at comoving momentum $q=8$ along the $z$ axis, the computed daughters deposit 38.45% of their energy in the longitudinal pressure moment and 30.77% in each transverse pressure moment on the quiet shell with the wide switch. For a parent at rest all three fractions equal one third. Averaging over an isotropic distribution of incident directions also gives three equal fractions. These are moments of outgoing particle occupations at the common readout $H_0\tau=2.10$. [Beam results](BEAM.json)

## Angular distribution

Use the conformal momenta of a massless Dirac particle and antiparticle with total energy $\omega$ and net momentum of magnitude $q$. The spin-summed square of the scalar Yukawa matrix element on the two-body final-state phase space is proportional to $\omega^2-q^2$. At fixed $(\omega,q)$ it is independent of the direction in the pair center-of-momentum frame. Thus $\mu=\cos\theta_*$ has uniform measure $d\mu/2$. The boost speed from that frame to the shell frame is $\beta=q/\omega$; the geometry and interaction switch enter through the distribution in $\omega$.

For the positive-$z$ pair the two conformal daughter energies are $e_{1,2}=\omega(1\pm\beta\mu)/2$, and their signed longitudinal momenta are $p_{1,2,z}=\omega(\beta\pm\mu)/2$. Locally measured energies and momenta are these quantities divided by $a_f$. A free particle contributes energy $e$ and longitudinal momentum flux $p_z^2/e$ in conformal units. Adding the two particles and performing the angular average gives
$$
 L(\beta)=\frac{\langle p_{1z}^2/e_1+p_{2z}^2/e_2\rangle}{\omega}
  =\frac12\int_{-1}^1\frac{(\mu+\beta)^2}{1+\beta\mu}\,d\mu
  =2-\frac1{\beta^2}+\frac{(1-\beta^2)^2}{\beta^3}\operatorname{atanh}\beta .
$$
The endpoint limits are $L(0)=1/3$ and $L(1)=1$. Near zero, $L=1/3+8\beta^2/15+8\beta^4/105+\cdots$. The coefficient of $\beta^{2n}$ is $8/[(2n-1)(2n+1)(2n+3)]$ for $n\ge1$, so the longitudinal fraction increases with the boost. Each transverse fraction is $(1-L)/2$, because the massless particles satisfy $\sum_i p_i^2/e=e$ and the system is azimuthally symmetric. The threshold $\omega=q$ has zero pair-production measure in this scalar channel.

## Which direction carries momentum?

Spatial translation symmetry fixes the net comoving momentum of each pair. Annihilating a parent in the traveling $+\mathbf q$ mode deposits $+\mathbf q$. Excess stimulated creation of another parent in that same mode deposits $-\mathbf q$. The sign reversal exchanges the two polar directions and leaves the diagonal pressure fractions unchanged. Let $dN_-$ and $dN_+$ be the positive contributions of annihilation and of excess stimulated creation in the one-parent minus switched-vacuum experiment. With the $y^2$ factor understood, the energy, momentum, and longitudinal pressure moment at the cut are
$$
 {\cal E}_f=\frac{1}{a_f}\int\omega(dN_-+dN_+),\qquad
 {\cal P}_{z,f}=\frac{q}{a_f}(N_- - N_+),\qquad
 {\cal S}_{zz,f}=\frac1{a_f}\int\omega L(q/\omega)(dN_-+dN_+).
$$
The net daughter comoving momentum equals the comoving momentum depleted from the occupied parent mode. A spatially uniform switch can supply energy while supplying zero net spatial impulse. For both transverse components together, ${\cal S}_{zz,f}+2{\cal S}_{xx,f}={\cal E}_f$. For a box-normalized tag, divide these energy moments by the physical box volume to obtain contributions to the occupation-diagonal stress tensor. An instantaneous local quantum stress on an anisotropic state can also involve particle–antiparticle coherences outside this occupation measurement. In an isotropic orientation average the net daughter momentum vanishes. The resulting outgoing massless energy has isotropic kinetic pressure $\rho/3$.

## Measurements

Each row uses the wide switch $1.15\le s\le2.05$ on the fine quiet history. The energy ${\cal E}_f$ and signed momentum ${\cal P}_{z,f}$ are reported per $y^2 H_0$; the counts are per $y^2$.

| $q$ | Annihilation count | Stimulated count | ${\cal E}_f$ | ${\cal P}_{z,f}$ | ${\cal S}_{zz,f}/{\cal E}_f$ | Each transverse fraction |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.339193 | 0.154948 | 10.528530 | 0 | 0.333333 | 0.333333 |
| 4 | 0.281920 | 0.103063 | 8.368647 | 0.578977 | 0.347997 | 0.326002 |
| 8 | 0.195070 | 0.027871 | 5.082772 | 1.082474 | 0.384531 | 0.307735 |

For $q=8$ on this wall, restricting to total comoving daughter energy $\omega\le60$ changes the longitudinal fraction to 0.385107. The narrow switch $1.25\le s\le1.95$ gives 0.386941 for the full frequency calculation and 0.387388 with the same cutoff. The finite-window boost mixture therefore matters when interpreting a directional measurement.

A local boost tests how much of the single beam's pressure anisotropy is explained by its net energy flux. Write $F={\cal P}_{z,f}/{\cal E}_f$, $Z={\cal S}_{zz,f}/{\cal E}_f$, and $X=(1-Z)/2$. To set the flux of the local kinetic tensor to zero, boost by $v=2F/(1+Z+D)$ along $z$, where $D=\sqrt{(1+Z)^2-4F^2}$. In this energy frame $\rho_{\rm L}/{\cal E}_f=(1-Z+D)/2$, $p_{\parallel,{\rm L}}/{\cal E}_f=(-1+Z+D)/2$, and $p_{\perp,{\rm L}}/{\cal E}_f=X$, with a common physical-volume normalization understood throughout. A boost of an isotropic radiation fluid would satisfy $F^2=(1+X)(Z-X)$. On the fine quiet wide-window $q=8$ history the right side exceeds $F^2$ by 0.055074. The flux-removing speed is 0.157643, and the resulting pressure fractions are $p_{\parallel,{\rm L}}/\rho_{\rm L}=0.363150$ and $p_{\perp,{\rm L}}/\rho_{\rm L}=0.318425$. Thus the one-direction excess has residual shear in its energy frame. Mixing all initial directions isotropically gives a different distribution.

The code evolves the three parent momenta in a shared conformal clock. It evaluates the two signs of the Fourier transform on frequency grids beginning exactly at each light-cone threshold. On both wall meshes it samples 2,048 points in each time window and integrates frequencies through 512 with step 0.125. Repeating the fine quiet and nearby-positive histories with 4,096 samples, frequency cutoff 1,024 and step 0.0625 changes the largest longitudinal fraction by less than $4.3\times10^{-10}$ in absolute value. This tests numerical resolution over the sampled range. A 256-node independent two-daughter angular integration verifies the analytic longitudinal and transverse fractions at nine boosts within $4.5\times10^{-15}$ on this implementation. The stable small-$\beta$ series avoids cancellation. The original wall archive and crossing normalization are authenticated before evolution.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B beam_stress.py /path/to/HDblast > beam-replay.json
~~~
