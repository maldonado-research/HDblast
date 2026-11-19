# A microscopic finite-time decay control

In Minkowski spacetime, take a real scalar of mass $m$ at rest and a massless Dirac field in its vacuum, coupled by $y f(t)\phi\bar\psi\psi$, where $f(t)=\exp[-t^2/(2T^2)]$. Subtract all asymptotic number and energy changes obtained from the initial vacuum under the same switching protocol. The remaining observables are excesses per additional prepared rest particle to leading order $y^2$. The sign convention for external work is positive when the external switch supplies energy to the fields.

Spatial homogeneity forces a pair connected to the rest mode to have momenta $\mathbf{p},-\mathbf{p}$ and total energy $\omega=2|\mathbf{p}|$. The spin sum is $4y^2p_1\cdot p_2=2y^2\omega^2$. The two-body massless phase space is constant and the external scalar normalization is $1/(2m)$. Thus the pair spectral weight away from resonance acquires the factor $(\omega/m)^2$, with overall scale fixed by the rest width $\Gamma_0=y^2m/(8\pi)$. The annihilation branch uses detuning $\omega-m$. Creating a rest scalar together with the pair uses $\omega+m$. The first branch has initial occupation factor one. In the second branch, stimulation raises the initially occupied rest mode's creation factor from one to two. Vacuum subtraction leaves one copy. Contributions from creation in previously empty scalar modes cancel.

Write $D$ and $R$ for the excess annihilation and creation branch yields. Since $\widetilde f(\nu)=\sqrt{2\pi}T e^{-T^2\nu^2/2}$, their leading spectra are

$$
 d(D,R)=\frac{\Gamma_0}{2\pi}(\omega/m)^2\big|\widetilde f(\omega\mp m)\big|^2\,d\omega.
$$

Set $u=mT$, $x=\omega/m$, $L=\Gamma_0\int f^2dt=\Gamma_0\sqrt\pi T$. Dividing by $L$ gives $(u/\sqrt\pi)x^2e^{-u^2(x\mp1)^2}dx$ for $x\ge0$. Define $A=(D+R)/L$ for excess pair count, $B=(D-R)/L$ for signed parent depletion, and $C=\int\omega\,d(D+R)/(mL)$ for daughter energy. With $g=e^{-u^2}/(u\sqrt\pi)$ the half-line moments yield

$$
 A=1+\frac{1}{2u^2},\quad
 B=A\,\mathrm{erf}(u)+g,\quad
 C=\left(1+\frac{3}{2u^2}\right)\mathrm{erf}(u)
       +\left(1+\frac{1}{u^2}\right)g.
$$

At the in/out endpoints the interaction vanishes, so the asymptotic free-energy difference equals work from the switch. Its excess is $W_{\rm ext}/(mL)=C-B=[\mathrm{erf}(u)+g]/u^2>0$.

| $mT$ | Excess pairs / $L$ | Signed depletion / $L$ | Daughter energy / $mL$ | External work / $mL$ | Stimulated creation / $L$ |
|---:|---:|---:|---:|---:|---:|
| 1 | 1.5 | 1.47160494 | 2.52185948 | 1.05025454 | $1.4197531\times10^{-2}$ |
| 2 | 1.125 | 1.12490429 | 1.37502655 | 0.25012225 | $4.7852758\times10^{-5}$ |
| 4 | 1.03125 | 1.031249999974 | 1.093750000002 | 0.062500000028 | $1.3126935\times10^{-11}$ |
| 8 | 1.0078125 | 1.0078125 | 1.0234375 | 0.015625 | $6.5975715\times10^{-34}$ |

At large $mT$ the normalized three matter observables approach one. The external work per $mL$ falls as $1/(mT)^2$ even when stimulated creation is exponentially suppressed. The off-resonance energies of the annihilation products also contribute to switching work. A controlled golden-rule limit requires both many oscillation times and weak enough coupling that the order-$y^2$ transition yields remain small. For a very short pulse one must also check that $LA\ll1$. These in/out calculations assume a prescribed coupling in flat spacetime. The pulse duration alone supplies no physical decay or thermalization timescale on the cosmological shell.

The standard-library script evaluates the elementary moments and independently integrates both spectra and their energy moments at two resolutions. The maximum relative difference over the 24 quantities at 8192 panels is $2.27\times10^{-10}$. The tiny creation branch is evaluated directly with a rescaled positive integrand and through a complementary error function. For each case it obeys $R/L\le e^{-u^2}/(4\sqrt\pi\,u^5)$.

```sh
python -I -B gaussian.py > reproduced.json
```
