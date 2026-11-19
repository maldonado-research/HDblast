# Which particles can the oscillating shell modulus produce?

The archived matter metric fixes the small-amplitude coupling of the modulus to the trace of the matter stress tensor. A single real, massless, minimally coupled scalar has a nonzero two-body decay width. The coefficient is $1/(32\pi)$ times the dimensional estimate in folder 146. A massless Dirac field, a conformally coupled massless scalar, and a classical four-dimensional massless gauge field have zero linear tree-level width through this universal coupling. A choice of daughter field therefore changes the physical prediction. [Archived estimates](../archive/BLAST_REHEATING_RESULTS.json) · [Computed widths](WIDTHS.json) · [Replay](widths.py)

## Normalization of the source

On the two archived tanh-branch minima the four-dimensional Einstein and matter metrics obey $g_E=f g_J$ in Planck units. Thus $g_J=A^2 g_E$ with $A=f^{-1/2}$. The archive tabulates the dimensionless canonical coordinate $\Theta$ using the **reduced** $M_{\rm Pl}=2.435\times10^{18}\,{\rm GeV}$. Define the canonically normalized fluctuation $\sigma=M_{\rm Pl}(\Theta-\Theta_m)$. Then
$$
\frac{A(\Theta)}{A_m}=1+\frac{\alpha\sigma}{M_{\rm Pl}}+O(\sigma^2/M_{\rm Pl}^2),\qquad
\alpha=\left.\frac{d\ln A}{d\Theta}\right|_m=-\frac12\left.\frac{d\ln f}{d\Theta}\right|_m .
$$
With $\delta S_m=\tfrac12\int\sqrt{-g_J}\,T_J^{\mu\nu}\delta g^J_{\mu\nu}$, the linear interaction in locally canonical Einstein-frame matter variables is $(\alpha/M_{\rm Pl})\sigma T^\mu{}_\mu$. The overall sign follows the stress-tensor convention and drops out of a single partial width. All daughter masses below are the physical Einstein-frame masses at this minimum, $m_i=A_m m_{i,J}$. This constant rescaling matters if a mass is supplied in matter-frame units. A convention using $d\ln A^2/d\Theta$ gives $2\alpha$ and must retain the $1/2$ in the metric variation. A direct mass-independent Yukawa term $y\sigma\bar\psi\psi$ is an additional interaction.

## Decay amplitudes and phase space

Let $M=m_\sigma$ and $x=(m_i/M)^2$. The width vanishes for $x\ge1/4$. For one real scalar $s$ with matter-frame curvature coupling $\xi R_J s^2/2$ in the conventional action, the improved stress tensor contains the term $\xi(\eta_{\mu\nu}\Box-\partial_\mu\partial_\nu)s^2$. Combining its trace with the scalar kinetic and mass terms gives the on-shell two-particle amplitude, up to a common sign,
$$
{\cal M}(\sigma\to ss)=\frac{\alpha}{M_{\rm Pl}}\left[(1-6\xi)M^2+2m_s^2\right].
$$
The two-scalar matrix element of $s^2$ contains both field contractions. For a second check, set the local Einstein scale factor to one and vary only the matter factor $A$. The Fourier mode $v_k=A s_k$ satisfies
$$
v_k''+\left[k^2+A^2m_{s,J}^2+(6\xi-1)A''/A\right]v_k=0 .
$$
Expanding $A/A_m$ and using $\sigma''=-M^2\sigma$ in the harmonic rest frame gives exactly the mass-squared insertion $\alpha[(1-6\xi)M^2+2m_s^2]\sigma/M_{\rm Pl}$. This checks the relative sign between the improvement and the mass term.

A Dirac field rescales as $\psi_E=A^{3/2}\psi_J$, leaving the mass term $-A m_{f,J}\bar\psi_E\psi_E$. Its effective Yukawa magnitude is $|\alpha|m_f/M_{\rm Pl}$. The spin-summed squared matrix element for a distinguishable fermion and antifermion is $2(\alpha m_f/M_{\rm Pl})^2(M^2-4m_f^2)$. Four-dimensional gauge kinetic terms contain the factor $\sqrt{-g_J}\,g_J^{\mu\rho}g_J^{\nu\sigma}$, whose Weyl factors cancel. The classical massless gauge trace thus gives zero at this order. A quantum trace anomaly, running gauge coupling or explicit modulus-dependent gauge kinetic function needs separate microscopic input.

The integrated Lorentz-invariant two-body phase space is $\beta/(8\pi)$, $\beta=\sqrt{1-4x}$. Multiply by $1/(2M)$ and the squared amplitude, with a further $1/2!$ for identical real scalars. One colorless Dirac flavor carries its full spin sum and requires no identical-particle factor. The resulting widths in the parent rest frame are
$$
\Gamma_0=\frac{\alpha^2 M^3}{32\pi M_{\rm Pl}^2},\qquad
\frac{\Gamma_s}{\Gamma_0}=(1-6\xi+2x)^2\sqrt{1-4x},\qquad
\frac{\Gamma_f}{\Gamma_0}=4x(1-4x)^{3/2}.
$$
Independent open species add their widths, with a complex scalar containing two real degrees of freedom and a color multiplicity multiplying a Dirac result. An energy fluid with $p=\rho/3$ has zero **mean** trace. Decay probes the vacuum-to-two-particle matrix element at invariant mass $M$, which is nonzero for a minimal real scalar.

| One species | Width relative to $\Gamma_0$ at zero daughter mass | Maximum relative width over open masses | $x$ at maximum |
|---|---:|---:|---:|
| Minimal real scalar, $\xi=0$ | 1 | 1.115419 | 0.1 |
| Conformal curvature coupling, $\xi=1/6$ | 0 | 0.071554 | 0.2 |
| Colorless Dirac fermion | 0 | 0.185903 | 0.1 |

The massive Dirac rate vanishes both at zero mass and at pair threshold. Maximizing its cubic phase-space factor times $x$ gives $1-10x=0$. For the scalar the interior nonzero stationary point is $x=(1+6\xi)/10$ when it lies in the open interval. A global maximum also compares the endpoints. The scalar improvement coupling is a microphysical choice.

## Calibrating both archived minima

The primary calibration uses the archive's tabulated $\alpha$ and $M/H_{\rm top}$ for the modified potentials with interior minima. The original ledger used $\Gamma_{\rm estimate}=\alpha^2M^3/M_{\rm Pl}^2$ and left the prefactor unspecified. At the same field parameters the single minimal-scalar value is $\Gamma_0/\Gamma_{\rm estimate}=0.00994718$. A model in which many independent real minimal massless fields are open would need at least 101 of them to reach this original estimate.

| $\phi_m$ | Archived $\alpha$ | Archived $M/H_{\rm top}$ | $\Gamma_0/H_{\rm top}$ divided by $(H_{\rm top}/M_{\rm Pl})^2$ | $\Gamma_0/H_{\rm top}$ at $t=10^{-3}$ |
|---|---:|---:|---:|---:|
| $+0.5$ | 0.2999263 | 7.683604 | 0.4059048 | $3.15294\times10^{-5}$ |
| $-0.5$ | 0.1260048 | 16.542539 | 0.7149581 | $5.55356\times10^{-5}$ |

Here $t$ controls detuning and $H_{\rm top}/M_{\rm Pl}=0.27870549\sqrt{t}$. Thus the ratios in the last column scale linearly with $t$ for massless minimal daughters. At $t=10^{-10}$ they are $3.15294\times10^{-12}$ and $5.55356\times10^{-12}$. The machine-readable table gives all five original scale choices, in both relative units and GeV. These are perturbative flat-vacuum rest-frame widths around a harmonic minimum. Applying a single width to a cosmological condensate requires oscillations whose amplitude and expansion rate vary slowly on a period. The calculation supplies neither the matter content nor a thermalization rate.

For an independent background check, write $\phi=-\tanh y$, $W=1-\phi+\phi^3/3$. A primitive of $W$ is $F(y)=y+(2/3)\ln\cosh y+(1/6)\tanh^2y$. The warp integral gives
$$
f(y)=2\int_0^\infty du\,e^{-2[F(y+u)-F(y)]/3},\quad D=1-Wf/3,\quad Z_y=2WD/f+6D^2/f^2 .
$$
Since $d\Theta/dy=+\sqrt{Z_y}$, the logarithmic coupling is $\alpha=D/(f\sqrt{Z_y})$. The Einstein potential at the minimum obeys $M^2/H_{\rm top}^2=3 f(0)^2 V_{,\phi\phi}(1-\phi_m^2)^2/(t f_m^2 Z_y)$ with $V/t=(1-\phi/\phi_m)^2(1+b\phi)$. Integrating directly at the nominal $\phi_m$ gives $\Gamma_0$ coefficients 0.4056137 and 0.7157419. Their relative offsets from the historical values are $-0.072\%$ and $+0.110\%$, consistent with using the exact minimum in place of the archived nearest mesh point and numerical second derivative. Both comparisons are retained. [Archived geometry definitions](../archive/blast_common.py)

## Replay

From the study root, run
~~~sh
python -IB widths/widths.py --ledger archive/BLAST_REHEATING_RESULTS.json --output /tmp/WIDTHS.json
~~~
The replay integrates the two-body radial delta function with its Jacobian and averages an explicit four-component Dirac trace over 40 outgoing directions. The largest absolute discrepancy from the dimensionless analytic width ratios is below $10^{-14}$. It independently optimizes the three massive channels, checks the threshold and conformal zeros, reconstructs both minima by warp quadrature and compares every scale to the archived ledger.
