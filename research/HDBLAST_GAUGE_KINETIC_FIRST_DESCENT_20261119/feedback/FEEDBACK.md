# When the photons act back on the rolling field

The fixed-background calculation can predict more energy in light than an unchanged parent could supply. We test that approximation by letting the quantum oscillators take energy from the parent and including their energy in the expansion rate. The calculation defines a specific **finite Hamiltonian model**. Its upper frequency cutoff changes the physics while the coupling is moving.

## Initial data and variables

We take one U(1) field, both polarizations, every complete round three-sphere harmonic from $L=2$ to a fixed integer $K$, and the conformal in-vacuum on the constant-coupling plateau. The gauge action is $-I(\Theta)^2F_{\mu\nu}F^{\mu\nu}/4$. We vary the end-to-end change $\delta=\log(I_f/I_i)$ with the prescribed smooth step over 20% to 80% of the field displacement between launch and first minimum. If the field reverses inside that region the same field-dependent action operates on its way back. The scale factor and field are classical. Units are $x=H_0t$, $a=H_0a_{\rm physical}$, $h=H/H_0$, $v=d\Theta/dx$, $\epsilon=(H_0/M_{\rm Pl})^2$, and densities are measured in $V_0=3M_{\rm Pl}^2H_0^2$. Each oscillator has instantaneous occupation $n$ and anomalous covariance $m=\langle bb\rangle$. A prime here denotes $d/dx$. Set $g=d\log I/d\Theta$, $s=gv$, and $w_L=L(L^2-1)$. The photon equations are

$$ n'_L=2s\,\mathrm{Re}\,m_L,\qquad m'_L=-2i(L/a)m_L+s(2n_L+1),\qquad
|m_L|^2=n_L(n_L+1). $$

In this finite system we subtract the half quantum of each oscillator for all $x$. Its energy is $L/(2a)$, independent of $I$. This prescription gives energy and stress zero for the conformal vacuum with constant $I$. The reference zero-point energy redshifts as free radiation. We omit its gravitational source. With $B=\epsilon/(3\pi^2 a^4)$ the remaining equations read

$$ r_\gamma=B\sum_{L=2}^K w_L n_L,\qquad
S=2Bg\sum_{L=2}^K w_L\mathrm{Re}\,m_L,\qquad Q=vS, $$
$$ v'=-3hv-3u_{,\Theta}-3S,\qquad h'=a^{-2}-v^2/2-2r_\gamma,\qquad a'=ah. $$

We evolve $v'$ using $S$ directly, including when $v=0$. The two energy balances are $r'_\gamma+4hr_\gamma=Q$ and $(u+v^2/6)'=-hv^2-Q$. The code integrates the source without dividing by the parent velocity, checks the Friedmann constraint and oscillator purity, and compares a separate set of position–momentum covariance equations at both weak and strong coupling. For the decoupled-background control the radiation equations still run while the corresponding loss and gravity terms are turned off in the parent system. That control serves only as a comparison.

## How much survives after the switch?

The following are photon percentages of the total scalar-plus-photon energy at the **first field minimum**, after the coupling has become constant. In every row the photon and parent energy affect the evolving geometry. The weak change is $I_f/I_i=e^{0.1}\simeq1.105$. The exploratory strong change is $e^3\simeq20.1$. The cutoff applies to harmonic label $L$, with proper energy $H_0L/a$. Changing $K$ changes which oscillators can take energy from the field.

| Highest harmonic $K$ | Weak, positive minimum | Weak, negative minimum | Strong, positive minimum | Strong, negative minimum |
| ---: | ---: | ---: | ---: | ---: |
| 48 | 0.582% | 0.703% | 63.98% | 52.29% |
| 64 | 0.692% | 0.956% | 60.97% | 48.55% |
| 96 | 0.724% | 1.068% | 44.60% | 30.83% |
| 128 | 0.720% | 1.081% | 36.60% | 25.34% |
| 192 | 0.685% | 1.020% | 21.56% | 20.40% |

For $K=48$, $\delta=3$ the field reverses inside the switch. Its total variation in $\log I$ by first minimum is 3.530 and 3.965 on the two branches. A photon-dominated outgoing state is possible in those two finite systems. The same contrast with a larger band gives a different outcome. Common local terms in an absolute stress do not cancel between different *backreacting* trajectories. A cutoff-independent claim would require matching or renormalizing the parent and gravity terms while the field moves and then demonstrating a stable late energy. The closed oscillators have no collisions. Radiation domination at one event establishes neither a temperature nor a surviving hot expanding era. [Subsequent interval with constant coupling](COAST.md)

## A rigorous energy bound for a specified range of wavelengths

From $|m|^2=n(n+1)$ and $|n'|\le2 |s|\sqrt{n(n+1)}$, vacuum initial data imply $n_L\le\sinh^2 V$, where $V=\int|d\log I|$ is the actual path length, including any reversals. The Friedmann equation gives the total density $h^2+a^{-2}$. On the first expanding segment, $a\ge a_i$, so the fraction in these modes has the bound

$$ F_{L\le K}\le \frac{\epsilon\,K(K-1)(K+1)(K+2)}{12\pi^2 a^2(1+a^2h^2)}\sinh^2 V
\le\frac{\epsilon\,K(K-1)(K+1)(K+2)}{12\pi^2 a_i^2}\sinh^2 V. $$

This uses the finite sum $\sum_{L=2}^K L(L^2-1)=K(K-1)(K+1)(K+2)/4$. It applies even if other energies also contribute to homogeneous Einstein gravity, provided the included photons obey these oscillator equations with the specified initial vacuum. An interaction which scatters energy into these modes calls for another bound. On either archived launch and for $V=0.1$, the modes whose **initial** proper frequencies are at most $0.25M_{\rm Pl}$ have $K=28$ and can supply at most **0.433%** of total energy throughout this expansion segment. Allowing initial frequencies through $0.5M_{\rm Pl}$ gives $K=56$ and at most **6.71%**. Even to avoid excluding one-half requires reaching at least $K=93$, an initial photon frequency of about $0.820M_{\rm Pl}$. This is a necessary frequency reach for this variation and normalization. It does not imply that the four-dimensional theory is valid to those energies or that this much energy is actually produced. The same proof can be reapplied with a different $H_0/M_{\rm Pl}$ or a different total variation. [Bound calculation](budget.py)

## Why increasing the cutoff changes a coupled answer

For high frequencies while $I$ is moving, the leading instantaneous polarization of a vacuum mode is $n_L\simeq (\partial_\eta\log I)^2/(4L^2)$ in conformal time $d\eta=dx/a$. Thus $\sum w_Ln_L$ contains a transient term proportional to $K^2$, even though the smooth, completed interaction on a fixed trajectory gives a fast-decaying outgoing spectrum. At mid-switch in the weak $K=192$ runs, the energy in the upper fifth of the retained labels is 1.035 and 1.039 times this leading estimate. It is misleading to count such an instantaneous contribution as particles that survive after the coupling settles. The calculation illustrates the effect of a specified regulator. The cutoff comparison provides no error bar for a renormalized continuum answer. [Solver and independent checks](feedback.py) · [Numerical output](FEEDBACK.json)
