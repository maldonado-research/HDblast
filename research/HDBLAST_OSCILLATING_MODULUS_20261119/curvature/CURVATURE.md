# Does radiation overtake the modulus before turnaround?

The two four-dimensional closed launches at detuning $t=10^{-3}$ reach their first maximum before the modulus completes the prescribed harmonic preparation. With the calibrated width active from launch, radiation supplies only 0.00256% and 0.00581% of the energy there. The modified-detuning branches with minima at $+0.5$ and $-0.5$ therefore both fail to reach radiation equality on their first expansion. The post-oscillation envelope gives a separate quantitative condition for initial states that do reach the dust regime.

## Matter preparation and allowed curvature

Write $x=H_{\rm top}t_{\rm proper}$, $h=H/H_{\rm top}$, $\gamma=\Gamma/H_{\rm top}$, and $S,R=\rho_m,\rho_r$ in units of $3M_{\rm Pl}^2H_{\rm top}^2$. The one-real-scalar width is $\Gamma=\alpha^2m^3/(32\pi M_{\rm Pl}^2)$. Define $M=S_i a_i^3$ and $r_i=R_i/S_i$. The system evolves a pressureless oscillation envelope, massless radiation and spatial curvature $K$,

$$
S_x+3hS=-\gamma S,\qquad R_x+4hR=\gamma S,\qquad
h^2=S+R-K/a^2,\qquad a_{xx}/a=-(S+2R)/2.
$$

Radiation equality means $R=S$ with $h>0$. A simultaneous equality and first maximum is the marginal case. In the early pure-dust limit the unique critical curvature is

$$
K_{\rm crit}=1.12766656\,(\gamma M)^{2/3}.
$$

Positive curvature below this threshold permits equality during expansion. Higher curvature produces a maximum first. Since the acceleration is strictly negative throughout the envelope, even a universe reaching equality can subsequently recollapse. Spatially flat expansion from early pure matter reaches equality at $\Gamma\Delta t=1.05380240$, when $H/\Gamma=0.58009667$.

The table gives the critical initial curvature fraction $K/[a_i^2(S_i+R_i)]$. Each scan fixes $a_i,S_i,R_i$ and recomputes the initial expansion rate from Friedmann as $K$ varies. The first three rows specify $a_i=1$, zero daughters, and matter density corresponding to **flat** $H_i=m/50$. These are separately prepared envelope controls. The final row uses each measured flat harmonic-gated trajectory after capture and two qualified cycles.

| Preparation | Detuning | $\phi_m=+0.5$ | $\phi_m=-0.5$ |
|---|---:|---:|---:|
| $H_i=m/50$ reference | $10^{-3}$ | 0.00392314 | 0.00343167 |
| $H_i=m/50$ reference | $10^{-10}$ | $8.45138\times10^{-8}$ | $7.39277\times10^{-8}$ |
| $H_i=m/50$ reference | $10^{-20}$ | $1.82079\times10^{-14}$ | $1.59272\times10^{-14}$ |
| Measured flat harmonic anchor | $10^{-3}$ | 0.0204800 | 0.0295948 |

The measured flat anchors have $a_i=25.1591,26.8207$. Their physical $K_{\rm crit}$ in the archive normalization is 0.00215444 and 0.00364103. Positivity of $h_i^2$ already requires $K<0.10520,0.12303$ respectively. The original $K=1$ trajectories thus require their own field evolution from the bounce. On those trajectories the harmonic-gated first maxima occur at $x=1.260889,1.616293$, with $h_x=-0.904,-1.151$ and zero gated radiation. The always-on first-maximum fractions reported above use the same initial geometry. A bound covers their entire measured expansion. Since $Ra^4$ increases and total density decreases until turnaround, $f(x)\le f_{\rm turn}(a_{\rm turn}/a_{\rm bounce})^4$. The respective upper bounds are 0.00546% and 0.0167%. The smaller-detuning rows in the table vary the physical width on a specified matter preparation. Tracking a new curved hilltop launch at either smaller detuning would require its own field evolution.

For all eight finite-anchor scans, an integration at $0.9K_{\rm crit}$ reaches equality during expansion. At $1.1K_{\rm crit}$ it reaches turnaround first. In the flat continuation of the measured anchors the scale factor grows by factors 70.25 and 48.63 before equality. These results concern energy fractions in the four-dimensional effective cosmology. Thermalization and coupling to observed matter require additional interactions.

## Event calculation and independent checks

For $K>0$ define $b=Ka/M$, $g=\gamma M/K^{3/2}$, $d\theta=\sqrt K\,dx/a$, $P=Sa^3/M$ and $Q=Ra^4K/M^2$. Primes here denote $d/d\theta$. The code evolves

$$
b''=P/2-b,\qquad P'=-gbP,\qquad Q'=gb^2P,\qquad
(\Gamma\Delta t)'=gb .
$$

The conserved Friedmann constraint is $b'^2+b^2=Pb+Q$. Initial data are $b_i=Ka_i/M$, $P_i=1$, $Q_i=r_i b_i$, and $b'_i=\sqrt{b_i(1+r_i-b_i)}$. The regular limiting matter solution starts at $b=b'=Q=0$. The code locates the first downward crossing of $b'=0$ and the first upward crossing of $Q/b-P=0$. Proper-time acceleration at the envelope maximum is proportional to $-P/(2b^2)-Q/b^3<0$. Its sign distinguishes turnaround from a bounce even though the Friedmann constraint alone gives $h=0$ at both.

The limiting critical value is $g_*=0.83508166118$. At the simultaneous event $b=0.81006079$ and $\Gamma\Delta t=0.90379316$. The second integrator eliminates velocity using the Friedmann constraint, sets $f=R/(S+R)$, $u=\Gamma\Delta t$ and takes $\psi=\operatorname{arccot}(aH/\sqrt{k_{\rm phys}})$ from its initial angle to $\pi/2$ at turnaround. In the normalized variables $\cot\psi=b'/b$, $b=e^{-u}\sin^2\psi/(1-f)$, and

$$
\frac{df}{d\psi}=\frac{2(1-f)}{1+f}(g b-f\cot\psi),\qquad
\frac{du}{d\psi}=\frac{2g b}{1+f}.
$$

Both coordinates locate $f=1/2$ at the same $g_*$ within $10^{-12}$. A different Runge–Kutta method and looser error control shift $g_*$ by less than $2\times10^{-12}$. The maximum absolute conformal Friedmann residual on the three reported control integrations is below $3\times10^{-14}$. The initial-size sequence $b_i=0.01,10^{-4},10^{-6}$ with zero initial radiation gives $g_*=0.8346273,0.8350812,0.83508166$.

## Why the early-matter boundary is unique

Write $U=\Gamma/\sqrt{(\rho_m+\rho_r)/(3M_{\rm Pl}^2)}$ and $\cos\psi=H/\sqrt{(\rho_m+\rho_r)/(3M_{\rm Pl}^2)}$. From initially pure matter, $f$ increases up to turnaround. Indeed $\dot f=(1-f)(\Gamma-Hf)$ starts positive, and at any candidate stationary point its derivative is $-(1-f)\dot H f>0$ because $\dot H=\ddot a/a-H^2<0$. Taking $f$ as the independent variable gives

$$
\frac{d\psi}{df}=\frac{(1+f)\sin\psi}{2(1-f)(U-f\cos\psi)},\qquad
\frac{dU}{df}=\frac{U(3+f)\cos\psi}{2(1-f)(U-f\cos\psi)}.
$$

On expansion both off-diagonal derivatives in coordinates $(\psi,-U)$ are positive. The small-$U$ solution has $\psi\sim(U/g)^{1/3}$ and $f=2U/5+g^{-2/3}U^{5/3}/7+O(U^2)$. The latter expansion follows by substituting into $df/d\ln U=2(1-f)(U/\cos\psi-f)/(3+f)$. At fixed small $f$, the lower-$g$ solution starts at larger $\psi$ and smaller $U$. The cooperative differential equations preserve this strict ordering. The lower-$g$ solution must reach $\psi=\pi/2$ at a smaller $f$. Continuous dependence on $g$ and the bounds below then give a unique critical value for the early pure-matter family.

For an analytic bracket, use $u=\Gamma\Delta t$ and the exact integral $R/S=e^u\int_0^u [a(z)/a(u)]e^{-z}dz$. Monotone expansion gives $R/S\le e^u-1$. Concavity from the limiting origin gives $a(z)/a(u)\ge z/u$ and hence $R/S\ge(e^u-1-u)/u$. Equality during expansion requires $u\ge\ln2$ and occurs by the positive root $u_*=1.25643121$ of $e^{u_*}=1+2u_*$ if expansion lasts that long. Before turnaround, $P\le1$ bounds $b(\theta)\le(1-\cos\theta)/2$ on $[0,\pi]$ by the forced-oscillator integral. Integrating by parts in $b'$ gives $b'(\pi)\le0$ for decreasing $P$, so $u_{\rm turn}\le\pi g/2$. Thus $g<2\ln2/\pi$ suffices for turnaround first. Conversely $b'$ stays positive through $\pi/2$ by its cosine-kernel representation. While $u\le u_*$ the sine-kernel integral gives $b\ge e^{-u_*}(1-\cos\theta)/2$. It follows that $g>2u_*e^{u_*}/(\pi/2-1)=15.465$ suffices for equality before turnaround. These broad bounds and the ordering proof are independent of the numerical root.

## Reproduction

Run from the study directory after reproducing the microscopic widths and measured scalar backgrounds.

```sh
python -IB curvature/curvature.py --background background/BACKGROUND.json --widths widths/WIDTHS.json --output /tmp/CURVATURE.json
```

The result records the independent turning-point ratios on both sides of the critical value, constraints, finite anchors, controls and microscopic rate for each detuning and minimum.
