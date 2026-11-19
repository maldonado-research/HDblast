# Evolving two possible decay prescriptions on the same four wall histories

At the common cut $s=H_0\tau=2.10$, the signs of all three pulse responses survive a change of oscillator used to define a decay. During the preceding mass crossing, the two energy definitions disagree even on the sign of the accumulated energy. A reservoir that starts empty and receives exactly the physical Hamiltonian energy lost to the alternative bath has a negative balance for part of its evolution.

The calculation uses the published fine-mesh reciprocal $y=1,b=1$ histories as prescribed functions, adopts $a=1$ at the original quiet reference crossing, and keeps $0<\kappa<8$. All three experiments begin at the native preparation time $s_i=0.769027655$ with empty reservoirs. Changing the bath or the initial state does not recompute the recorded geometry. The model contains a single minimally coupled scalar with the archived dimensionless signed mass $\mu=100(\phi_b-0.5)$. Frequencies, densities and powers below use $H_0$ as the unit.

## What is evolved

For each comoving momentum let $E^2=(\kappa/a)^2+\mu^2$, $v=a^3$, $h=\dot a/a$, and $\gamma=\mu^2/(8\pi E)$, with dots in $s$. In canonical variables $(\chi,\Pi=v\dot\chi)$ evolve their real symmetrized covariance $(X,C,P)$. Damping adds $-\gamma(\Sigma-\Sigma_b)$ to the ordinary Hamiltonian covariance equation. The instantaneous physical Hamiltonian ground state H has $X_0=1/(2vE), C_0=0, P_0=vE/2$.

The alternative B is the local WKB tangent of $u=a\chi$ at $W=\Omega=aE$ with $u'/u=-i\Omega-\Omega'/(2\Omega)$ in conformal time. In physical time its real shear is

$$
 D=h+\frac{h\mu^2+\mu\dot\mu}{2E^2},\qquad
 X_B=X_0,\quad C_B=-\frac D{2E},\quad P_B=P_0+\frac{vD^2}{2E}.
$$

This choice specifies a pure reference state at each instant. The minimally coupled conformal equation itself contains the curvature term $-a''/a$. The present comparison uses $W=\Omega$ as a freely specified jump prescription.

Let $U=(P-P_0)/(2v)+vE^2(X-X_0)/2=E n_H$. Two powers per mode enter the comparison,

$$
 L_H=\gamma\left(U-\frac{D^2}{4E}\right),\qquad
 L_B=\gamma E n_B=\gamma\left(U+D C+\frac{vD^2X}{2}\right).
$$

The first is signed energy removed by the dissipator from the physical Hamiltonian. The second counts positive B oscillator quanta with assigned energy $E$. Their difference $L_d=L_B-L_H$ is the power required from an additional drive if a radiation reservoir is fed by the quantum count. For the physical bath, set $D=0$ throughout. The code integrates $L_d$ directly from its expanded expression. A change of basis applied only in reporting would leave the evolution intact; these runs actually change its damping target.

The stored rates incorporate the measure $d\kappa\,\kappa^2/(2\pi^2)$. For each power, the code computes both $J=\int L\,ds$ and $Z=a^4\rho=\int aL\,ds$. The latter assumes radiation redshift. The field obeys the independent check $\Delta\!\int U = S+A-J_H$, where $S=\int v\mu\dot\mu(X-X_0)$ is the explicit mass work and $A$ integrates the expansion term $-3h$ times the ground-subtracted comoving pressure. Mode measure and proper time are implicit in both work integrals. Fitting a thermal distribution or evolving a microscopic daughter sector lies beyond this calculation.

## Fixed-wall results

Entries give quiet endpoint densities in $H_0^4$ and pulse endpoint density changes relative to quiet in percent. Each pair lists signed physical loss followed by the positive count. The near-positive, near-negative, and far-positive runs use the archived injection signs and locations.

| Preparation / bath | Quiet densities | Near positive (%) | Near negative (%) | Far positive (%) |
|---|---:|---:|---:|---:|
| H / H | 37.157633 / 37.157633 | −0.171407 / −0.171407 | +1.016698 / +1.016698 | +0.099361 / +0.099361 |
| H / B | 37.031636 / 37.027887 | −0.171382 / −0.174208 | +1.024878 / +1.034943 | +0.102033 / +0.101511 |
| B / B | 36.841314 / 36.839557 | −0.171511 / −0.174340 | +1.024278 / +1.034352 | +0.101649 / +0.101128 |

On the quiet H/B history the initial comoving signed and counted powers are −0.13354848 and +0.13354848. The signed power first crosses upward at $s=0.78005343$. The deeper minimum of the accumulated signed balance is $Z_H=-0.07362851$ at $s=1.54064317$. At the cut the integrals $(J_H,J_B,J_d)=(73.2785861,73.2970067,+0.0184206)$ close. The radiation-weighted drive density is −0.00374874. Different time weights account for the sign difference between the integrated drive and its final radiation-weighted density.

Preparing B in its own vacuum costs 0.07506558 in initial comoving field excitation relative to H, or 0.49951245 in density. The initial dissipative powers vanish. Later the quiet signed balance reaches $Z_H=-0.09833568$ at $s=1.54081047$. Thus initial state mismatch alone cannot account for the negative excursion. The original H/H control and the positive-count reservoirs have minimum accumulated balance zero. Both B protocols have positive late signed balances. If one interprets $Z_H$ as energy in a preloaded free-radiation reservoir on the identical fixed histories, keeping it nonnegative throughout would require initial quiet densities of at least 0.92155 (H/B) or 1.23079 (B/B). This conditional bound is $-\min Z_H/a_i^4$.

At a simple mass zero the nonadiabatic term in the target injection scales as $\gamma D^2/(4E)\sim\dot\mu^2\mu^4/[128\pi(p^2+\mu^2)^3]$, with $p=\kappa/a$. A low-momentum mode samples a peak scaling as $p^{-2}$ over a time of order $p/|\dot\mu|$; multiplying its time integral by the three-dimensional momentum measure gives order $p\,dp$. This local estimate addresses that specific term. A finite band and a positive quadrature grid are essential qualifications on the numerical results.

## Reproduction and checks

Use a GitHub HDblast checkout containing the 2026-11-18 source-work wall archive and the pulse-response quiet rolling summary. The script checks both file SHA-256 digests before evolving anything. Only NumPy and SciPy are needed. Run from this directory with the checkout path in `REPO`.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -I -B evolve.py "$REPO" > EVOLUTION24.json
python -I -B evolve.py "$REPO" --order 48 --tol 2e-10 --curves QUIET48.csv > EVOLUTION48.json
python -I -B evolve.py "$REPO" --order 96 --tol 2e-10 > /tmp/hdblast-evolution96.json
python -I -B compare.py EVOLUTION24.json EVOLUTION48.json /tmp/hdblast-evolution96.json > CHECK.json
```

The wall interpolation uses cubic Hermite polynomials for log scale factor and mass with their recorded first derivatives. DOP853 uses maximum proper-time step 0.001, relative tolerance $10^{-9}$ for 24 nodes and $2\times10^{-10}$ for 48 and 96 nodes, and absolute tolerance one thousandth of each relative tolerance. Comparing all 12 evolutions, the maximum absolute shift in final signed, counted or drive density is $3.73\times10^{-8}$ from 24 to 48 nodes and $5.65\times10^{-8}$ from 48 to 96. The original archived radiation, linearly interpolated to the common cut, agrees with the H/H replay within $4.01\times10^{-5}$. Regenerating the curves evaluates all 3 quiet solutions at 401 equally spaced times. Its columns contain $a^4\rho$; extrema reported in the JSON use continuous-time power-zero events.
