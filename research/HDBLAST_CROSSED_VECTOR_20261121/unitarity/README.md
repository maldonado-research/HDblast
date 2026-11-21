# A probability check for crossed vector exchange

The published HDblast collision kernel keeps a neutral vector resonance formed and destroyed by a charged particle–antiparticle pair (the $s$ channel). The simplest local field theory with the vertex $gV_\mu\bar\chi\gamma^\mu\chi$ would also let the same pair exchange the vector (the $t$ channel). Here we test a specific extension that uses the original timelike pole and an undressed spacelike propagator $1/(t-M^2)$. This section concerns the **dedicated single-flavor benchmark** with $M/m=8300$ and $g^2=12\pi(0.30)$. It asks where that prescription would require sizeable corrections to satisfy probability conservation. The original pulse histories still use the isolated pole. The conserved equal-mass vector current removes the propagator’s longitudinal numerator. This test omits the distinct long-range electromagnetic exchange, for which screening or an infrared prescription would be needed. [Published interaction assumptions](https://github.com/maldonado-research/HDblast/blob/9467ffbc67e9540da6762fecd1c44e5e22c80522/research/HDBLAST_DIRECTIONAL_BIRTH_20261121/scattering/RESONANCE.md)

## One helicity state suffices for a necessary test

Take a distinguishable particle and antiparticle with the same physical helicity before and after scattering. The differences of helicities in both pairs vanish. Under a common rotation of the final momenta about the $y$ axis, the $J=0$ projection therefore uses $d^0_{00}=1$. Put $z=\cos\theta$ and write $s,t$ in mass-squared units. In the spinor phase convention used in the code, contraction of the two vector currents gives the numerators

$$
 N_t=2s-6m^2+2m^2 z,\qquad N_s=-4m^2 z.
$$

To see the angular dependence, choose incoming Pauli spinors $(1,0)$ and rotate the outgoing spinors together to $(\cos(\theta/2),\sin(\theta/2))$. The antiparticle carries the conjugate spin label. Inserting the positive- and negative-energy Dirac spinors normalized to $u^\dagger u=v^\dagger v=2E$ gives $N_t=4E^2\cos^2(\theta/2)+4p^2[2-\cos^2(\theta/2)]$. The $s$ term integrates to zero regardless of its energy-dependent width. With no identical-state factor for the unlike charges, the magnitude of the real diagonal tree partial wave is

$$
 |a_0|=\frac{\beta g^2}{32\pi}\int_{-1}^{1}
 \frac{2s-6m^2+2m^2 z}{M^2+(s-4m^2)(1-z)/2}\,dz,\qquad
 \beta=\sqrt{1-4m^2/s}.
$$

The code evaluates the integral analytically, checks it against an independent angular quadrature and checks the numerators against complex four-component spinors. At $m=0$ it reduces to $(g^2/8\pi)\log(1+s/M^2)$. [Helicity construction](https://doi.org/10.1016/0003-4916(59)90051-X) · [Helicity projection and unitarity, Eqs. 2.1–2.3](https://arxiv.org/abs/1510.02110) · [Finite-mass normalization, Eq. 2](https://arxiv.org/abs/1805.07306)

The same check also applies to the **identical equal-helicity channel**. For two like charges take the initial Pauli spinors $(1,0)$ and $(0,1)$ along the opposite momenta and rotate both outgoing spinors together. Direct spinor contractions give the same $N_t(z)$ and $N_u(z)=-N_t(-z)$. Antisymmetrizing the two fermions gives the sum of the $t$-exchange function at $z$ and $-z$. Its $J=0$ integral doubles. In the full-angle convention the incoming and outgoing identical two-particle states each carry $1/\sqrt{2}$. The normalized diagonal $|a_0|$ consequently equals that of the unlike channel. This normalization also follows by including the $1/2!$ two-body phase space in the optical theorem. Both signs therefore have a failing diagonal elastic wave above the same threshold in this dedicated benchmark. [Equal-helicity fermion states, footnote 3](https://arxiv.org/abs/2108.00013) · [Identical-state normalization, Eqs. 2.54–2.58](https://arxiv.org/abs/2409.16790)

## When would this amplitude need to change?

With $S_{ii}=1+2i a_{ii}$, probability conservation requires $|S_{ii}|\le 1$ even if many other final states are open. Hence any physical diagonal partial wave $a$ lies in $|a-i/2|\le1/2$. The real part of a *physical* amplitude cannot have magnitude greater than $1/2$. The dedicated tree term first exceeds that value at **$\sqrt{s}/M=1.42749146$**, corresponding to $s/m^2=140379347.51$. Its values at $\sqrt{s}/M=1,2,3$ are 0.312, 0.724 and 1.036 in magnitude. This one diagonal element already detects the breakdown. Below the crossing, other helicity or angular-momentum channels and higher-order terms could also matter. [Partial-wave bounds in spinning theories](https://arxiv.org/abs/2108.00013)

For a real tree value $A$ the shortest distance to the allowed disk is $\sqrt{A^2+1/4}-1/2$. This is a lower bound on the size of a **complex amplitude correction**. It gives no unique corrected cross section. At $|A|=1/2$ the minimum correction is 41.4% of $|A|$. For each mass value in an archived encounter the calculation also asks whether every physically allowed diagonal amplitude differs from the proposed value by more than 10%, 50% or 75%. Any purely real tree amplitude generally needs an imaginary correction at higher order. The size of the required correction is what makes the diagnostic informative here.

In the finest late sample, **29.16%** of the previously accepted elastic encounters exceeded the real-part crossing and carried **38.52%** of the incoming conformal energy summed over encounters. For **18.74%** of all encounters a physically admissible diagonal amplitude would require more than a 50% change in the complex plane. A particle can collide repeatedly, so encounter-energy fractions count some particles several times. Population-level transport and the shell-frame stress derivative are measured separately. These findings limit the use of the added tree-level diagrams. Different couplings, spacelike form factors and a complete unitary interaction can lead to different transport. [All cases and controls](PARTIAL.json)

Run the script on the event tapes reconstructed from the collision-photon study.

~~~sh
python -I -B unitarity/partial_wave.py --events EVENTS --output PARTIAL.json
~~~
