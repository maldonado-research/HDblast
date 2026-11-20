# How a wall source couples to the photon continuum

The scalar-weighted Maxwell field on each static tuned leaf has a massless photon and a continuous set of heavier modes that extends down to zero mass. A wall source can excite these additional modes. Their coupling at very small mass falls as the **6.2 power of the mass** for the chosen five-dimensional interaction. The spectral weight determines how much extra energy the same wall current delivers to the continuum. A short and a slow current produce very different mixtures.

## Outgoing boundary response

We use one retained half-leaf, in units of the wall length $L_0$, and set $g_5^2=1$. Four-dimensional masses and frequencies are measured with the warp factor normalized to one at the shell. The bulk scalar is $\phi=-\tanh y$, the shell is at $b=-\operatorname{arctanh}\phi_b$, and the specific added Maxwell weight is $w=G=(1-\phi^2)^2$. The tuned warp is $A(y)=-[P(y)-P(b)]/3$ with $P(y)=y+(2/3)\log\cosh y+\tanh^2 y/6$. Define $p=w e^{2A}$, $dz=e^{-A}dy$ and $j=(we^A)^{1/2}$. For a transverse vector of four-dimensional mass $m$, its radial profile $\psi$ solves $-(p\psi_y)_y=m^2w\psi$. With $u=j\psi$, the Schrödinger potential is $j_{zz}/j$. Primes below refer to $z$.

Choose time dependence $e^{-i\omega t}$. Take $u$ to have unit outgoing plane-wave amplitude toward the far throat and define $d=u'-(j'/j)u$. Its phase is arbitrary. The one-copy response of the boundary value to an applied current is

$$
B(-[m+i0]^2)=-\frac{u_b}{j_b^2 d_b},\qquad
\rho(m)=\frac{2m}{\pi}\operatorname{Im}B
       =\frac{2m^2}{\pi j_b^2|d_b|^2}\quad(m>0).
$$

The last equality follows from the conserved outgoing flux $\operatorname{Im}(u^*d)=m$. It keeps the small positive spectral density numerically resolvable even when the Green function is predominantly real. $\rho$ is per unit mass, with units fixed by the wall length convention. The complete static spectral measure also contains a discrete massless pole of residue $1/N$, where $N=\int_b^\infty w\,dy=(1+\phi_b)^2(2-\phi_b)/3$. Doubling the action for two equivalent sides halves both residues for the same total source. All continuum-to-pole ratios here are invariant under that common rescaling. [Original five-dimensional reduction](https://github.com/maldonado-research/HDblast/blob/1475602d67399699bf4461d72c7a8c1742c8ab08/research/HDBLAST_GAUGE_KINETIC_FIRST_DESCENT_20261119/origin/ORIGIN.md)

## Why slow drives see very little leakage

At large positive $y$ the warp is $A\sim a_0-ky$, where $k=5/9$ and $a_0=[P(b)+(2/3)\log 2-1/6]/3$. The convenient far radial coordinate is $r=e^{ky-a_0}/k$, for which $dr/dz\to1$. Then $j\sim J_\infty r^{-s}$ with

$$
s=(1+4/k)/2=4.1,\qquad J_\infty=4e^{-2a_0/k}k^{-s}.
$$

A unit-amplitude outgoing shape (flux $m$) for the leading inverse-square potential is $\sqrt{\pi mr/2}\,H^{(1)}_{s+1/2}(mr)$ up to phase. At $1\ll r\ll1/m$ its magnitude approaches $C_\nu (mr)^{-s}$ with $C_\nu=2^{\nu-1/2}\Gamma(\nu)/\sqrt\pi$, $\nu=s+1/2$. Its decaying zero-energy profile on the interior is proportional to $j$. The integrated radial equation gives $p_b\psi_y(b)=m^2\int_b^\infty w\psi\,dy$, with the flux derivative at infinity vanishing for any fixed $m>0$ because $j$ decays. Matching the leading zero profile in this integral gives $|j_b d_b|\sim (C_\nu N/J_\infty)m^{2-s}$ and hence

$$
\rho(m)\sim c m^{2s-2},\qquad
c=\frac{2J_\infty^2}{\pi C_\nu^2 N^2}.
$$

| Wall scalar $\phi_b$ | $N$ | Threshold coefficient $c$ | $\rho(0.25)$ | $\rho(1)$ | $\rho(4)$ |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.666667 | 0.0834056 | 0.00001528 | 0.061746 | 0.593853 |
| +0.5 | 1.125000 | 0.0786815 | 0.00001486 | 0.094817 | 1.037990 |
| −0.5 | 0.208333 | 0.1642741 | 0.00002954 | 0.093093 | 1.038357 |

## Energy from a prescribed transverse current

Apply a weak, conserved spatially transverse current at the static wall, with Gaussian proper-time profile $J(t)\propto\exp[-t^2/(2\sigma^2)]$, fixed spatial momentum $k_\parallel>0$, and zero incoming classical field. The energy deposited in a normalized oscillator is proportional to its squared boundary amplitude times $|\widetilde J(\sqrt{k_\parallel^2+m^2})|^2$. The common factors cancel against the normalized massless oscillator. This gives the continuum-to-pole energy ratio and continuum share of the **deposited** energy,

$$
R(\sigma)=N\int_0^\infty\rho(m)e^{-\sigma^2m^2}\,dm,\qquad f_{\rm cont}=\frac{R}{1+R}.
$$

At long duration, $R\sim Nc\,\Gamma((2s-1)/2)/(2\sigma^{2s-1})$, with $2s-1=7.2$. The energy comes from the external current. This experiment describes a classical response on a static background. It does not calculate photon pair creation from vacuum, the movement of a self-gravitating shell, or the conversion of cosmological matter energy into radiation.

| $\sigma/L_0$ | Continuum share at $\phi_b=0$ | At $+0.5$ | At $-0.5$ |
|---:|---:|---:|---:|
| 0.25 | 44.457% | 70.094% | 29.526% |
| 0.5 | 16.523% | 36.521% | 8.932% |
| 1 | 2.0773% | 5.6859% | 0.9856% |
| 2 | 0.05374% | 0.12570% | 0.02693% |
| 4 | 0.000461% | 0.000822% | 0.000265% |

## A second scalar weight changes the threshold power

As a separate calculation take $w=e^{2A}$, with the shell normalization $A(b)=0$. This changes the five-dimensional action and differs from the original $w=\operatorname{sech}^4y$. A constant rescaling of the action on a single leaf leaves continuum-to-pole ratios unchanged. Here $j=e^{3A/2}$, so $s=3/2$, $J_\infty=k^{-3/2}$ and $C_{2}=2^{3/2}/\sqrt\pi$. The same outgoing matching gives

$$
\rho(m)\sim \frac{m}{4k^3N^2},\qquad N=\int_b^\infty e^{2A}dy.
$$

At mass 0.0003 the observed ratios to this prediction are 0.99999793, 0.99999833 and 0.99999770 for $\phi_b=0,+0.5,-0.5$. The respective analytic coefficients are 1.3590324, 0.8978140 and 1.6833279. The continuum begins at zero for both weights but has substantially different strength at low mass. For the second weight its inverse-mass-squared moment diverges logarithmically. Its positive spectral representation implies $B(Q)-1/(NQ^2)\sim [1/(4k^3N^2)]\log(1/(QL_0))$ at small spacelike $Q$.

## Numerical and spectral checks

The program propagates $u_y=(\log j)_y u+e^{-A}d$ and $d_y=-(\log j)_y d-e^{-A}m^2u$ inward. It starts with outgoing Hankel data at $y=8$, keeps the $dr/dz$ Jacobian and fixes the flux exactly there. It compares matching surfaces at 8 and 10 and a second calculation propagating a real, regular Neumann solution outward to $y=9$ and taking its Wronskian with the outgoing pair. On the reported masses 0.005–8 the maximum spectral density discrepancy in either check is below $9\times10^{-7}$ relative. The density remains positive over the sample. The analytic sign for all positive masses follows from the flux.

A different check integrates the measured timelike spectrum at spacelike momentum $Q$. It compares $\int_0^\infty\rho(m)/(Q^2+m^2)\,dm$ directly with the independently computed spacelike boundary response after removing its massless pole. At $Q=0$ the independent value is $\int_b^\infty [\int_y^\infty w(t)dt]^2/p(y)\,dy\,/N^2$. We use 24- and 16-point Gaussian rules in each indicated mass interval through $M=32$. Above mass 8 the quadrature uses a matching surface at $y=6$. For the dispersion tail we use the large-$m$ expansion $\rho=2[1+(H'-H^2)_b/(2m^2)+O(m^{-4})]/(\pi j_b^2)$ with $H=j'/j$. It follows by inserting the outgoing logarithmic derivative in its Riccati equation. Changing the tail match from 32 to 16 is a numerical sensitivity check. For the original weight the timelike dispersion and elliptic continuum agree within $1.4\times10^{-8}$ relative at $Q=0,0.3,1,3$ on all three leaves. For the Gaussian drive no such asymptotic tail is needed. Positivity and the $Q=0$ moment bound its omitted contribution to $R$ by $N M^2\exp(-\sigma^2M^2)\int\rho(m)/m^2\,dm$ for $M\sigma\geq1$. [Numerical results](RETARDED.json) · [Replay code](retarded.py) · [Independent spacelike response](../euclidean/EUCLIDEAN.md)
