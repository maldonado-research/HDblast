# Matching initial pressure and scalar force limits pair-count reversal

Fixing the incoming parent occupation and free-field energy permits a phase correlation to reverse three excess pair-count contrasts in the moving-parent shell model. Here we additionally match the initial isotropic pressure and the scalar force to their values in the incoherent preparation. On the finest momentum grid, the distant positive pulse retains at least **0.3828% of its original positive contrast**, even at the maximum Gaussian covariance allowed by the occupation. The two near pulses can each reverse under these same initial constraints when optimized individually. [Pre-calculation plan](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-6996927941) · [Parent calculation](../../HDBLAST_MOVING_PARENT_20261119/phase_feasibility/FEASIBILITY.md)

## Sources at preparation

For the isotropic population $n_q=1$ on $0<q<8$, write $c_q=\langle a_{\boldsymbol q}a_{-\boldsymbol q}\rangle$ in the instantaneous oscillator basis with real positive starting field amplitude. The difference relative to $c_q=0$ in the initial field variance of a mode is proportional to $\mathrm{Re}\,c_q/E_q$, where $E_q^2=m_i^2+q^2/a_i^2$. The change in the kinetic variance is proportional to $-E_q\mathrm{Re}\,c_q$ with the same positive prefactor. They cancel in the free Hamiltonian. The mass depends on the background scalar, so matching its force requires matching the field variance. For the minimally coupled scalar, the pressure change is proportional to $-[m_i^2+2q^2/(3a_i^2)]\mathrm{Re}\,c_q/E_q$. Thus the two additional source conditions are

$$
\int_0^8 q^2 \,dq\; \frac{\mathrm{Re}\,c_q}{E_q}=0,
\qquad
\int_0^8 q^2 \,dq\; \frac{q^2\mathrm{Re}\,c_q}{E_q}=0.
$$

Imaginary $c_q$ satisfies both conditions mode by mode and can change the initial field–momentum correlation. Pressure and force are matched at the starting clock. Later stresses evolve on the prescribed wall and mass histories.

## How the bound works

A row $i$ denotes a pulse-minus-quiet excess count per mean parent in the comoving band, with each wall's own switched vacuum subtracted. With $d_i$ its signed incoherent contrast, the phase-dependent result is $d_i(1-m_i)$ where $m_i=\mathrm{Re}\sum_q V_{iq}c_q$. The response coefficients $V_{iq}$ and Gauss weights are defined in the parent calculation. Simultaneous reversal requires $\min_i m_i>1$.

At $n=1$, the classical separable covariance disk has radius 1. The occupation bound $|c_q|\le\sqrt{n(n+1)}=\sqrt2$ allows general two-mode Gaussian states, including mixed states. Anomalous covariances inside the disk need not correspond to pure states. Both source constraints are homogeneous, so the maximum common margin scales exactly with the disk radius.

For the numerical unit-disk problem let $G$ be a row basis for the momentum-weighted source constraints $G\,\mathrm{Re}\,c=0$. For any simplex weights $p_i$ and any real moment multipliers $\lambda$, a computable upper bound and a feasible-state lower bound are

$$
\min_i\mathrm{Re}\sum_q V_{iq}z_q
\;\le\;\max_{\substack{|c_q|\le1\\G\,\mathrm{Re}\,c=0}}\;\min_i m_i(c)
\;\le\;\sum_q\left|(pV)_q+(\lambda G)_q\right|.
$$

The upper bound follows by weighting the rows, adding the vanishing moment term, and maximizing separately in each complex disk. Minimizing the upper expression gives the finite-dimensional dual. We find candidate states and multipliers with an inscribed 256-sided polygon per node. We project the trial real covariance back into the equality nullspace, scale the state to respect the circular disk, and directly evaluate both bounds in the *original circular problem*. With $\mathrm{Re}\,c_q=0$ each remaining variable lies in a real interval and the corresponding dual expression is $\sum_q|\mathrm{Im}(pV)_q|$. The JSON files record both witnesses for every case, including the original physical moment residuals.

## Three pulses and 24 selections

The following margins use the fine wall, wide switch, full pair-energy spectrum and occupation bound $\sqrt2$.

| Momentum nodes | Feasible common margin | Upper bound on common margin |
|---:|---:|---:|
| 24 | 0.99742746 | 0.99745892 |
| 48 | 0.99631479 | 0.99634626 |
| 96 | 0.99614396 | 0.99617164 |

The upper certificate puts all its weight on the distant positive pulse. Its original contrast is about $+2.06\times10^{-4}$ in units of the squared coupling, so the certified residual is at least $+7.8\times10^{-7}$. At 96 nodes the near positive and near negative pulses attain margins of at least 1.12199 and 1.06289 respectively, each with its own source-matched state. The corresponding fine three-arm optimum for separable states lies in $[0.70438015,0.70439972]$. An independent canonical-transfer and direct Fourier calculation places the distant-pulse quantum margin at 0.99616484. [Independent transfer and momentum refinement](../source_check/INDEPENDENT.md) · [96-node certificate](HIGH_Q.json)

| 96-node constraints and responses | Separable bounds | General Gaussian bounds |
|---|---:|---:|
| Fine moments, all 24 selections | 0.70258347–0.70260321 | 0.99360307–0.99363099 |
| Fine and coarse moments, all 24 selections | 0.64188321–0.64189892 | 0.90775994–0.90778216 |
| Imaginary only, three fine wide full pulses | 0.64407257 | 0.91085616 |
| Imaginary only, all 24 selections | 0.64159799 | 0.90735658 |

The 24 selections combine three pulses, two wall grids, two time windows and full or screened pair energies. Matching only the two fine-grid moments supplies an upper bound for every state that matches both grids. The separate four-condition result uses the small difference between the archived initial scale factors and masses. It describes the additional requirement that one covariance match two slightly different numerical initializations.

For that calculation write $u=\sqrt{q^2+(a_fm_f)^2}$, $v=\sqrt{q^2+(a_cm_c)^2}$, $s=u+v$, and $d=u^2-v^2\ne0$. The four moment functions span $u,1/u,v,1/v$. They can be evaluated stably in the equivalent basis $s,1/s,1/(uvs),1/(uvs^3)$. This follows from $u-v=d/s$, $1/v-1/u=d/(uvs)$ and $1/u+1/v-4/s=d^2/(uvs^3)$. We scale each weighted basis row and orthogonalize with QR before solving. The largest normalized physical source residual of the finest both-mesh witness is below $6\times10^{-16}$.

These are finite-grid certificates for a band-limited population on the archived model trajectories. They concern the matched-vacuum excess count. The numerical refinements leave the physical threshold close, and do not give a continuum error bound. Determining a source-matched state's effect on later shell evolution requires evolving the driving scalar and wall together with the parent state.

## Reproduce

From the parent directory of this appendix, run the following with NumPy and SciPy and set `repo_root` to the repository checkout containing the previous moving-parent study.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B source_match/source_match.py "$repo_root" > /tmp/hdblast-source24.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B source_match/source_match.py "$repo_root" --nodes 48 --samples 2048 --padding 32 > /tmp/hdblast-source48.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B source_match/source_match.py "$repo_root" --nodes 96 --samples 2048 --padding 32 > /tmp/hdblast-source96.json
```

The baseline uses 1024 switch samples and Fourier padding 16. The other two momentum grids use 2048 switch samples and padding 32. The driver authenticates its response helpers, which authenticate the historical wall inputs. [Baseline](BASELINE.json) · [48 nodes](REFINED.json) · [96 nodes](HIGH_Q.json) · [Driver](source_match.py)
