# A compact incident pulse on the original rolling B1 shell

The original B1 branch already crosses its particle-production threshold. This experiment asks whether adding a constraint-matched incoming scalar disturbance changes the predicted number of particles, and whether a local crossing rule predicts the sign obtained from propagating the free quantum modes. It uses a fixed mass law and a common early vacuum. The measurements are numerical model outputs.

## Same-band comparison

In the common comoving band $0\leq\kappa=k/(a_{*,0}H_0)\leq8$, $a_{*,0}$ is the unperturbed scale factor at the crossing and $s=H_0\tau$. Each entry for a number is per comoving volume in units of $a_{*,0}^3 H_0^3$. The direct modes have minimal curvature coupling and a time-dependent mass. Both pulse signs and the zero control start from the same rolling shell. The sign labels describe the front that first reaches the shell. The full initial scalar data also contain a compensating pulse farther from the shell.

| Fine-mesh quantity | Zero pulse | Positive front | Negative front |
|---|---:|---:|---:|
| First crossing time, $T_*$ | 1.230400 | 1.229991 | 1.230806 |
| Crossing speed, $d\phi_b/ds$ | 1.283190 | 1.278746 | 1.289862 |
| Original event rule, band number target | 3.609286 | 3.752859 | 3.471790 |
| Leading event rule with derivative factor set to one, same histories | 3.683268 | 3.671131 | 3.700226 |
| Direct free modes at $s=2.10$, band number | 3.514078 | 3.506934 | 3.548930 |
| Direct free modes, change relative to zero | — | −0.203% | +0.992% |
| Decaying B1 shell gas, local radiation variable $R$ at $s=2.10$ | 0.00005053 | 0.00005238 | 0.00004884 |

The derivative-corrected rule gives the opposite sign of the number change for both matched histories. The positive front also fails B1's own local validity test, which requires $G\ge G_{\min}$. Here $G=100$ and $G_{\min}$ is approximately 43.7, 216, and 33.4 for the zero, positive and negative branches, respectively. The primary packet spans about 0.063 in shell proper time on arrival, compared with $1/\sqrt q\simeq0.088$ for the original crossing. A one-sided quartic fit over the previous 0.2 in $s$ therefore smooths over a feature shorter than the characteristic particle formation time. The positive branch has $v_\mathrm{fit}=1.3270$ while the crossing velocity is 1.2787. The corresponding values for the negative branch are 1.2404 and 1.2899.

An additional six coupled evolutions disabled only the derivative factor in the gas prescription, starting with the same pre-crossing data. The fine-mesh radiation at $s=2.15$ becomes $5.58198,5.54734,5.63418$ in units of $10^{-5}$, ordered zero, positive and negative. Its relative changes switch to −0.621% and +0.935% from +3.67% and −3.35% in the corrected gas. This comparison lets the shell dynamics and matter backreaction respond to the change.

These branches support a conclusion about the limitations of an event approximation for these fixed, time-dependent histories. The late free-mode results depend on the full matched scalar data. The fitted correction also uses an expansion in inverse $q$, whose applicability varies across the branches. The same simulations retain the historical phenomenological $\chi$-gas energy transfer and decay. The direct modes are computed on these saved histories, without feeding their stress tensor back into the spacetime. Their occupation estimates the free, undecayed production. Adding it to the gas density would double count the same modeled degree of freedom. No inference about a hot Big Bang, measurable cosmological spectrum, or an ultraviolet completion follows from this calculation.

## Constraint-matched initial data

The common background evolves with the original B1 dynamics until $T_0=1.11$. The scalar seed has $dc=0.01$, the production threshold is $\phi_*=0.5$, and $G=100,y=1,b=6.49\times10^{-6},\mu_0=0$. From the checkpoint, three separate Cauchy problems evolve under the original RK4, wall junctions, quadratic mass law and gas prescription. Both meshes use $L=16$, $\kappa_\mathrm{damp}=10$, fourth-order spatial differences and CFL 0.5. Fine/coarse spacings near the wall are 0.0005/0.001. The unperturbed fine replay gives the original archived $T_*=1.2303997675$ and $s_*=1.4861216180$.

Define $f(x)=\exp[1-1/(1-x^2)]$ for $|x|<1$ and zero elsewhere. Set

$$
F(Z)=\epsilon f((Z+0.12)/0.02)+c f((Z+0.25)/0.105),\qquad
\delta\phi=F,\quad\delta\phi_T=-F'.
$$

The sign is for a right-moving characteristic in the local conformal coordinates. Both scalar data and the geometry equal the rolling background in the collar adjacent to the wall. For $\epsilon=+0.001,-0.001$, the fine mesh gives $c=-0.00122564,+0.00502900$. The negative branch requires a larger compensator. This asymmetry matters when comparing late responses. Compensator width and location were chosen by searching for a matching root for both signs before inspecting spectra. They define one particular family of globally matched initial data.

For background variables $Q=A_Z,u=A_T,p=\phi_T,r=\phi_Z,B_0=B$, hold $A,u$ fixed, set $B=B_0+g$, and write $K=F'(r-p+F')$, $D=Q-u^2/Q$, $S=(1+u/Q)K/3$. Subtract the original Hamiltonian and momentum constraints and solve for $z=\exp(-2g)-1$ and the change in $B_T$. With $E=\exp(2B_0)$, $U_0=U(\phi)$ and $\Delta U=U(\phi+F)-U_0$, the resulting linear radial equation is

$$
z' = \frac{2}{D}\left[(EU_0/3-S)z-S-E\Delta U/3\right],\qquad
\delta B_T=(-u g'-K/3)/Q .
$$

Integrate from $Z=-0.10$ with $z=0$. Adjust $c$ until $z(-0.355)=0$, then continue by the identity outside the two ends. $D$ and $Q$ stay bounded away from zero. Both endpoints match within $5\times10^{-16}$ in $z$. For the positive branch the maximum change in the original finite-difference Hamiltonian/momentum residuals, relative to the larger local sum of absolute terms, is $1.43\times10^{-5}/7.84\times10^{-5}$ on the coarse grid and $2.37\times10^{-6}/1.45\times10^{-5}$ on the fine grid. The negative branch shows similar refinement. These are two-ended matched data on a numerical background with finite residuals.

In this chart radial characteristics obey $dZ/dT=\pm1$. The front can first reach the wall at $T=1.210$ and the remote compensating data at $T=1.255$. All three B1 event times lie strictly between them. Through $T=1.20$ the coarse grid has a pre-arrival difference in field speed of at most $1.31\times10^{-7}$ and in field value below $8.8\times10^{-11}$. On the fine grid the speed difference is below $4\times10^{-12}$. The numerical leakage decreases on refinement. Since radial matching integrates outward from the wall edge, the metric inside $Z>-0.145$ is independent of the distant compensation coefficient. These facts isolate the change at the first event to the incoming front. The late spectra and radiation include both the front and the compensator.

## Free quantum calculation and physical scale

Propagate the minimally coupled real scalar $X_k=a^{3/2}\chi_k$, using standard unit-Wronskian normalization,

$$
X_k''+\left[\kappa^2/a^2+G^2(\phi_b-\phi_*)^2-\tfrac94h^2-\tfrac32 h'\right]X_k=0,\qquad h=H/H_0 .
$$

Primes here denote $d/ds$. Interpolate the shell fields in proper time and integrate $h$ for the common scale factor with $a_{*,0}=1$. A positive-frequency second-order WKB state starts at $s=0.80$ on the shared prehistory. A matching adiabatic projection at $s=2.10$ gives $|\beta_\kappa|^2$. The comparable event-rule occupation is $(1+D/q)\exp[-\pi(\kappa/a_*)^2/q]$, where each branch has its own crossing scale $a_*$ in the common normalization. Integrate both with $\kappa^2d\kappa/(2\pi^2)$. The factor $a_*^3$ is thus included when changing variables to physical momenta at the event. The curvature coupling differs from a conformally coupled, massless control.

The band was selected to remain below $M_5/H_0=b^{-1/3}=53.61074$ in the earlier-start test as well. For every retained trajectory over $s=0.50$ to 2.15 the physical single-particle energy $\sqrt{(\kappa/a)^2+G^2(\phi_b-\phi_*)^2}$ is below $52.48H_0$ throughout this band. The bound uses all extrema of the interpolated mass and scale factor, and conservative interval extrema. The small margin to the chosen effective gravity scale leaves sensitivity to ultraviolet model corrections. Out-of-band integrals up to 60 in the receipts are numerical diagnostics. The in-band physical quasiparticle energy densities at $s=2.10$ in units of $H_0^4$ are 95.267, 95.081 and 96.232. These estimates omit vacuum polarization and the interacting energy-transfer problem.

The ordering survives halving the mesh spacing and momentum quadrature spacing, tightening the time integrator, moving the initial vacuum to $s=0.50$, changing the endpoint basis from order two to order zero and measuring at $s=1.95,2.02,2.10,2.15$. All comparisons use the same basis and momentum interval across branches. Zeroth-order WKB here retains the derivative of its amplitude. Late absolute particle number is finite-time and basis dependent. The largest fine-mesh unit-Wronskian residual in the baseline integrations at $s=2.10$ is below $1.5\times10^{-8}$. The full numerical outputs include the relative endpoint adiabaticity, which is below 0.011 at the default initial time and 0.007 at the primary final time across sampled momenta. Two mesh resolutions establish numerical agreement for this comparison and do not establish an asymptotic convergence order.

## Reproduction

The compact files in `runs/` contain the measured proper-time prehistories, shell time series, original event reports and matched-data receipts for both resolutions. The prehistory files in this bundle retain the three history columns required by the quantum calculation. A full Cauchy checkpoint for evolving the bulk can be regenerated in a new directory. Python 3.12, NumPy and SciPy with `cumulative_simpson` suffice.

```sh
python summarize_rolling.py runs replay-summary.json
OPENBLAS_NUM_THREADS=1 python rolling_modes.py runs replay-plus.json --grid fine --branch=0.001
python run_rolling.py /path/to/original-project scratch --grid fine --branch capture
python run_rolling.py /path/to/original-project scratch --grid fine --branch quiet
python run_rolling.py /path/to/original-project scratch --grid fine --branch=0.001
python run_rolling.py /path/to/original-project scratch --grid fine --branch=-0.001
```

Repeat with `--grid coarse` for the other resolution. Add `--leading-only` to the three restart commands for the six coupled derivative-factor controls. Their shell histories appear under `runs/{grid}/leading-{branch}`. The file `leading_coupled.json` lists matched proper-time samples. The driver authenticates the four original B1 source files by SHA-256. [Original source tree](https://github.com/maldonado-research/HDblast/tree/fe35172f7c94dbd616c677045d59c6616c194d9b/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source) Provide a tree with their original `research/HDBLAST_CHECKPOINT_20260930` and `research/HDBLAST_CHECKPOINT_20260927` paths. `rolling_results.json` contains both grids, numerical constraints and arrival controls. The `modes/` receipts record both grids and all finite-state sensitivity checks. The strict calculation uses $\Delta\kappa=0.125$ and relative tolerance $2\times10^{-10}$. For an independent derivation of the same canonical scalar equation, the original project's source path is `research/HDBLAST_CHECKPOINT_20260927/preheating/verify/rerun/ph_lib.py`. All field amplitudes, compensator supports, vacuum times and bands are model choices. The registration and historical unperturbed sources predate this exploratory matched-pulse study.
