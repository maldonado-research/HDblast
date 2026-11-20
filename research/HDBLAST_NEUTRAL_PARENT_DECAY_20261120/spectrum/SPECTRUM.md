# Which fermion pairs fit on a closed wall?

A neutral scalar at rest couples through $y\chi\bar\Psi\Psi$ to a Dirac field of mass $\mu$. We test the elementary decay that underlies the fluid's adopted reaction rate. At the first late magnetic zero in the neutral precursor control, a particle matched to $Q=10^8$ and $q=1$ is too light for an energy-conserving decay into two free fermions on the sphere. At much larger charge the applied electric field can make a charged pair directly in one tested quantum mode. That process has no Yukawa vertex and would have to be included separately in a microscopic theory.

## Measure the phase space at the proposed time

We use the published neutral-seed wall with the same initial seed energy $D_i=10^{-4}$. We interpolate the tabulated radius, scalar, expansion rate and electric and magnetic fields with shape-preserving cubics in wall time. This supplies a reference history before any reaction or particle current is added. The proposed mass law is $M=2\mu R(\psi)$ with $R_i=0.95$, $R_f=1.05$, central $\psi=0.77$ and smooth width $0.03$. Energies and masses below are proper energies in units of $\hbar c/L$. Times are in units of $L/c$. The scalar-dependent parent changes from subthreshold to superthreshold according to the flat-space criterion at wall time 10.34058. At the following magnetic zero, $t=10.5017102$, $a=9.55331$, $R=1.0322505$ and $X=-0.00088811$. The reference wall starts contracting again at $t=12.42343$.

The intrinsic two-component Dirac eigenvalues on the unit round three-sphere are $\pm(n+3/2)$, each sign with multiplicity $(n+1)(n+2)$ for $n=0,1,\ldots$. In a four-component massive Dirac field each intrinsic eigenspinor supplies one positive and one negative energy solution, with $E_n=\sqrt{\mu^2+[(n+3/2)/a]^2}$. Thus one positive-energy shell has multiplicity $g_n=2(n+1)(n+2)$. A homogeneous scalar parent is a singlet of the spatial symmetry. Its Yukawa overlap is diagonal in the spatial eigenspinor basis; the antiparticle carries the conjugate label. The sum counts $g_n$ independent pair slots at that energy. The number of pair slots over the first $b$ shells is $2b(b+1)(b+2)/3$ per Dirac field, before any filling or rescattering. These facts describe charged fermions. Scalar daughters require their own spectrum and curvature coupling. [Dirac eigenspinors of spheres, Camporesi and Higuchi](https://arxiv.org/html/gr-qc/9505009v1)

For any shell to fit within the parent's energy budget in the free frozen spectrum,

$$\mu\sqrt{R^2-1}\ \geq\ \frac{3}{2a}.$$

Actual stationary resonance also requires $2E_n=M$. At the reference magnetic zero the first inequality gives **$\mu\geq0.61331$**. The earlier material parameter $Q=q/(\sqrt3\,\mu)=10^8$ then requires **$|q|\geq1.0623\times10^8$**. With $q=1$ the lightest pair costs 0.31403 while the parent has only $1.19194\times10^{-8}$. The electromagnetic potential vanishes on this slice even though the electric field persists.

The photon coupling measured by short-distance, slowly varying charged probes in the wall-dominated regime would be $\alpha=q^2/(4\pi S_0 K)$, where $S_0=\sigma L^4/(\hbar c)$. Requiring $\alpha<1$ at the kinematic mass floor gives $S_0>6.80\times10^{10}$ at the late $K=13200.04$, or $S_0>8.98\times10^{17}$ if imposed also at the initial $K=0.001$. These are conditional bounds on a common microscopic interpretation. The initial-energy inventory alone gives $N_\chi=\pi^2 a_i^3 D_i S_0/(\mu R_i)$ neutral parents. At the floor, one whole initial parent requires $S_0\geq73.8$. Many parents, emission into vacant final states, and a kinetic treatment are additional assumptions. [Earlier unit normalization and charge tests](https://github.com/maldonado-research/HDblast/pull/187)

## How the gauge field shifts the threshold

Write $s=qA/\sqrt3$, with $A=a^2 B/2$. In a left-invariant spin frame the intrinsic operator on the unit sphere is $D_s=2\boldsymbol\sigma\cdot\mathbf J+3/2-s\sigma_3$, up to orientation conventions. For one left spin $j=0,1/2,\ldots$ there is an independent right multiplicity $2j+1$. The two extreme states have eigenvalues $2j+3/2\pm s$. The internal two-state blocks at total spin projection $h=-j+1/2,\ldots,j-1/2$ have eigenvalues

$$d_\pm=\tfrac12\pm\sqrt{(2j+1)^2+s^2-4hs}.$$

Indeed, the diagonal entries of such a block are $1/2\pm(2h-s)$ and the off-diagonal entry is $2\sqrt{(j+1/2)^2-h^2}$. The internal square root is at least $\sqrt{8j}\geq2$ for $j\geq1/2$, and there is no internal block for $j=0$. The exact minimum magnitude is therefore
$$d_{\min}(s)=\min_{n\geq0}\left|n+\tfrac32-|s|\right|.$$
Direct Hermitian diagonalization independently checks degeneracies, both signs of charge, internal branches and edge crossings.

Using the bounded interpolated reference history throughout its open, expanding interval, $q=1$ has $d_{\min}/a\geq0.15270$. Even the largest allowed daughter momentum for that matched mass is only $1.8485\times10^{-9}$. The frozen-geometry channel stays closed. For a different $q$ matched to the same $Q$, the magnetic field can move an edge through zero. A necessary floor at a general time is $1.5/[Q|A|+a\sqrt{R^2-1}]$ for $\mu$. At late times it can be much smaller than the free floor at $A=0$. It tests energetic access and gives neither a rate nor an upper bound on the number of pairs in a changing field.

## Discrete transitions and time resolution

For a static sphere at the reference magnetic zero, resolve one initial homogeneous parent and an empty Dirac sector. A normalized spatial parent contributes $1/\sqrt{2MV}$, $V=2\pi^2a^3$. In each Dirac block $H=\mu\sigma_z+p_n\sigma_x$ with $p_n=(n+3/2)/a$. Scalar mass modulation couples through $\sigma_z$ with squared positive–negative overlap $p_n^2/E_n^2$. A Gaussian temporal window $w(t)=\exp[-t^2/(2T^2)]$ thus gives the first-order transition probability divided by its exposure $\int w^2dt$,

$$\frac{\Gamma_T}{y^2}=\frac{\pi}{MV}\sum_n g_n\frac{p_n^2}{E_n^2}\,\delta_T(2E_n-M),\qquad \delta_T(z)=\frac{T}{\sqrt{\pi}}e^{-T^2z^2}.$$

This is a controlled time-frequency experiment on a fictitious frozen background. A finite window can inject energy. Its off-shell weight is not a spontaneous decay rate. The flat comparison uses $g_n/V\to p^2dp/\pi^2$, giving $\Gamma_\infty/y^2=M(1-4\mu^2/M^2)^{3/2}/(8\pi)$ in the continuum and long-time limit. At fixed finite volume the long-time spectrum consists of individual lines.

| $\mu$ | $T$ | Finite sphere / flat with same width | Broad flat / sharp flat rate | Distinct lines within $1/T$ |
|---:|---:|---:|---:|---:|
| $5.77\times10^{-9}$ | 2 | 0.9123 | $2.88\times10^{16}$ | 1 |
| $5.77\times10^{-9}$ | 16 | $1.84\times10^{-9}$ | $4.51\times10^{14}$ | 0 |
| 1 | 2 | 0.9936 | 6.835 | 7 |
| 10 | 2 | 0.9996 | 1.118 | 21 |
| 100 | 2 | 0.999996 | 1.0011 | 19 |
| 100 | 128 | 2.936 | 1.0000003 | 1 |

At $\mu=10^4,T=1$ with gauge and geometry held fixed the continuum rate is recovered to $4.4\times10^{-7}$ relative. Resolving the threshold requires $2\mu(R-1)T\gg1$. Sampling many nearby lines requires $T(2\beta/a)\ll1$ with $\beta=\sqrt{1-R^{-2}}$, away from the first shell. At this epoch there are only 1.92172 proper time units before the turnaround. Even assigning that entire interval to one observation, the first criterion reaches unity only at $\mu=8.07$. The Gaussian examples are frequency-resolution controls and may extend beyond the physical interval. Slow variation of the parent mass, metric and electromagnetic Hamiltonian must be justified separately.

## Direct electric transition in one resolved fermion mode

Large charges change the problem again. In the extremal $n=0$ spin state the spatial eigenvector is independent of $A$, and the moving FRW metric can be accounted for by rescaling the fermion field by $a^{3/2}$. In proper time an exact two-state test Hamiltonian is

$$i\frac{d}{dt}\begin{pmatrix}u\\v\end{pmatrix}=
\begin{pmatrix}p(t)&\mu\\\mu&-p(t)\end{pmatrix}
\begin{pmatrix}u\\v\end{pmatrix},\qquad p(t)=\frac{3/2-Q\mu A(t)}{a(t)}.$$

It has no parent and no Yukawa decay. Start in the normalized instantaneous negative-energy state at $t_c-Y/\sqrt v$, where $p(t_c)=0$ and $v=|\dot p(t_c)|$. Read out the positive-energy projection at $t_c+Y/\sqrt v$. This prescribes a vacuum for this one mode at the start of a selected passage. It does not track a vacuum from before the initial pulse. A linear, infinite sweep gives $P_{\rm LZ}=\exp(-\pi\mu^2/v)$. [Zener's two-level crossing](https://doi.org/10.1098/rspa.1932.0165)

| $\mu$ at fixed $Q=10^8$ | Linear sweep | Direct history, $Y=16$ |
|---:|---:|---:|
| 0.61331 | 0.999978306 | 0.999978304 |
| 1,000 | 0.9652464 | 0.9652430 |
| 100,000 | 0.0290950 | 0.0291035 |

At the kinematic mass floor the edge crosses $2.88\times10^{-6}$ proper time units before the magnetic zero. Its half-interval is 0.06855 for $Y=16$. With $Y=8$ the result is 0.999978329. Reversing the sampled real Hamiltonian gives the same probability to $6\times10^{-13}$; holding the momentum constant gives a transition probability below $5\times10^{-34}$. The derivative inferred from the measured $A$ agrees with the electric field to about $6\times10^{-5}$ relative. Independently reconstructing $A$ with a Hermite spline using the measured electric field for its endpoint slopes changes the occupation by $1.3\times10^{-9}$ or less. All sub-sample crossing times refer to the interpolated history.

Adjacent extremal angular branches pass through zero about $1.92\times10^{-6}$ time units apart locally, or only $0.000449$ of the scale $1/\sqrt v$. The primary single-mode occupation is conditional on a prescribed empty positive-energy state and external fields. Many other modes, pre-existing occupation, particle currents, backreaction and Pauli blocking matter for a total yield. A large one-mode transition indicates that a static final-state vacuum for a Yukawa width is a consequential physical assumption. It does not measure a cosmological production rate or a temperature.

## Reproduce

~~~sh
python -I -B spectrum/spectrum.py --history prior/CHARGED.json --output /tmp/SPECTRUM.json
~~~

[Program](spectrum.py) · [Computed shells, bounds and controls](SPECTRUM.json)
