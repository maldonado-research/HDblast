# Independent transfer calculation for the alternative bath

A Gaussian bath matched to the slowly varying conformal mode can add energy to the physical field on the recorded rolling shell. If all of its signed energy loss is assigned to an initially empty massless sector, the inferred density becomes negative. This occurs even when the massive field starts in this bath's own vacuum. The positive particle-count rule gives a different radiation ledger. Here we evaluate both rules on identical fixed histories using finite oscillator transfers.

## Frozen channel

Use the fine quiet and nearby positive histories in the source-work release. Preparation occurs at $T=0.70$ and $s=H_0\tau=0.76902765$. Normalize $a$ by the original quiet crossing, retain $0\le k\le8$, and use the published $y=1$ rate $\gamma=m^2/(8\pi E)$ with $E=\sqrt{(k/a)^2+m^2}$. At the midpoint of each time step use the physical phase-space pair $x=\sqrt{a^3E}\,\chi$, $p=\Pi/\sqrt{a^3E}$. Define $n=(\langle x^2\rangle+\langle p^2\rangle-1)/2$, $c=(\langle x^2\rangle-\langle p^2\rangle)/2$, and $d=\langle xp+px\rangle/2$. The shear $j=D/E$ for the conformal target follows from $u=a\chi$ and $u'/u=-i\Omega-\Omega'/(2\Omega)$, with $\Omega=aE$. With $h=\dot a/a$ and signed $\mu=100(\phi-0.5)$,

$$
D=\tfrac12(3h+\dot E/E)=h+\frac{h\mu^2+\mu\dot\mu}{2E^2}.
$$

At a fixed midpoint the target covariance has $n_t=j^2/4$ and $z_t=-j^2/4-ij/2$ where $z=c+id$. The original physical bath corresponds to $j=0$. For the sheared target the number that counts departures from its vacuum is $n_B=n+j d+(j^2/2)(n+1/2+c)$. The signed removal of physical Hamiltonian energy per mode is $\gamma E(n-n_t)$, while the positive count prescription assigns $\gamma E n_B$. The two agree for the original physical bath. All reported densities use $d\nu=k^2dk/(2\pi^2)$ and a massless redshift law.

Treating $a,E,\gamma,j$ as constant on a step of length $\delta$, the linear channel gives
$$
n_1=n_t+(n_0-n_t)e^{-\gamma\delta},\quad
z_1=z_\infty+(z_0-z_\infty)e^{-(\gamma+2iE)\delta},\quad
z_\infty=\frac{\gamma z_t}{\gamma+2iE}.
$$
The code integrates these exponentials analytically over each step for both powers. It increments $a^4R$ by $a\int d\nu\,E\int\gamma(n-n_t)ds$ or $a\int d\nu\,E\int\gamma n_B ds$. It converts the final moments back to $(\chi,\Pi)$ before changing the midpoint basis. This construction avoids numerical time integration of the covariance differential equations. Midpoint freezing is second-order for a smooth history and introduces a finite time-discretization error. The proper-time interpolants for $\ln a$ and $\phi$ use cubic Hermite data from the recorded $h$ and field velocity. Time steps are eight times denser within $0.035$ on either side of the mass crossing to resolve the low-momentum modes.

## Results

$R_H$ is the density inferred from signed Hamiltonian loss. $R_{\rm count}$ weights the specified positive number by the instantaneous mode energy. The first line uses the physical vacuum at preparation and as the bath target. The next two lines vary the preparation for the same sheared bath. All use the original $y=1$ feedback-on wall histories prescribed throughout the calculation.

| Initial vacuum and bath | Quiet $R_H$ at $s=2.10$ | Quiet $R_{\rm count}$ | Nearby positive $R_H$ | Nearby positive $R_{\rm count}$ | Sampled minimum quiet $R_H$ |
|---|---:|---:|---:|---:|---:|
| Physical, physical | 37.15761 | 37.15761 | 37.09392 | 37.09392 | 0 |
| Physical, sheared | 37.03162 | 37.02787 | 36.96815 | 36.96336 | −0.06450 |
| Sheared, sheared | 36.84129 | 36.83954 | 36.77811 | 36.77531 | −0.08615 |

The minimum in the sheared-preparation quiet experiment occurs near $s=1.5394$ among sampled steps. At preparation its signed and counted powers both vanish exactly. Later negative density therefore persists when the two initial states match. At the final cut the nearby positive contrast in this case is −0.17151% for signed Hamiltonian loss and −0.17434% for the count rule, each compared with quiet. The negative intermediate density is a consistency problem for interpreting signed loss as energy deposited in passive, initially empty daughters. An external drive or an interacting reservoir would need its own energy accounting.

The finest calculation uses 48 Gauss–Legendre momentum nodes, a maximum ordinary time step of $10^{-4}$ and a maximum near-crossing step of $1.25\times10^{-5}$. Each final radiation density changes by less than $6.1\times10^{-5}$ from the preceding 24-node calculation at half the time resolution. The independent adaptive covariance calculation agrees to within $2.1\times10^{-5}$ in each of the displayed endpoint densities. The comparison uses the same recorded background and bath assumptions, so it tests the evolution and deposition algorithms.

## Replay

With NumPy and SciPy installed, point the program at the repository root or the release's historical-input tree containing the two cited source files. It checks the SHA-256 signatures of the wall archive and the quiet crossing summary before reconstructing the proper clock. The report verifier uses Python's standard library.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B bath_transfer.py /path/to/repository > transfer-replay.json
python -I -B verify_transfer.py transfer-replay.json
~~~

[Transfer solver](bath_transfer.py) · [Computed results](TRANSFER.json)
