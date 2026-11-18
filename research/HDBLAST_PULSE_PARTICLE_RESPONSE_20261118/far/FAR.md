# A compensated pulse whose second scalar disturbance cannot reach the readout

The rolling B1 wall crosses $\phi_*=0.5$ after receiving a positive incident pulse. We compare it to the quiet crossing with the second scalar disturbance placed too far away to arrive during any measured part of the outgoing evolution. The directly evolved low-momentum spectrum increases by 0.05115% on the finer grid. The model's native event prescription predicts a 3.978% increase in the same band for the same background pair.

## What reaches the wall

On the common numerical slice $T=1.11$, the primary scalar bump has amplitude +0.001 and support $[-0.14,-0.10]$. A second bump centered at −0.85 with half-width 0.25 restores the exterior constraint data. The two continuous constraints determine the metric response after fixing the scalar data. The signed second amplitude is solved numerically. On the finer grid it is +0.03151452. The two zero-denominator checks give $\min|A_Z|=1.1840$ and $\min|A_Z-A_T^2/A_Z|=0.6568$. Its left metric mismatch is $4.4\times10^{-15}$ after shooting. The changes in the normalized discrete Hamiltonian and momentum constraint residuals have maxima $2.37\times10^{-6}$ and $1.45\times10^{-5}$, respectively. Both improve on refinement.

The null characteristics obey $|dZ/dT|=1$ in this conformal bulk chart. The primary front can reach the wall at $T=1.21$. The positive fine-grid wall crosses at $T=1.22999135$. The earliest arrival of the remote control is $T=1.71$. The outgoing cut $s=H_0\tau=2.15$ occurs at $T=1.610320$, leaving 0.09968 coordinate-time units. The scalar pulse has a metric perturbation required by the constraints. Integrating the first-order equation for $y=\exp(-2\delta B)-1$ with zero data at the right endpoint fixes its solution at any $Z$ from data to its right. The metric and velocity perturbation between the two bumps is therefore independent of the remote coefficient. The entire perturbation produced by changing that remote control lies outside the past of the readout, while the primary pulse and its associated geometry are in its past.

As a direct numerical check, shifting the control center to −0.90 and increasing its half-width to 0.30 preserves its right edge and changes its required amplitude to +0.0359075. Through the full run the two wall histories differ by at most $7.3\times10^{-14}$ in $\phi$ and $1.5\times10^{-11}$ in $H/H_0$.

## Two predictions on the same band

All profiles use the common comoving scale factor normalization in which the quiet crossing has $a_*=1$. The rolling mode solver uses minimal coupling, $G=100$, the same second-order adiabatic basis, and $0\le\kappa\le8$ as in the other rolling runs. The table gives the percentage change of the comoving band occupation from quiet at $s=2.15$. The event column applies the native curvature-corrected filling prescription on that same momentum band.

| Local amplitude | Rolling spatial grid | Direct free modes | Native event |
|---|---|---:|---:|
| +0.001 | Coarse | +0.049249% | +3.98163% |
| +0.001 | Fine | +0.051147% | +3.97787% |
| +0.0005 | Coarse | +0.025792% | +1.92716% |
| +0.0005 | Fine | +0.026713% | +1.97316% |

At fine resolution the full-pulse mode response ranges from +0.05113% to +0.05116% over four outgoing cuts from $s=1.95$ to $s=2.15$. The shifted control yields +0.05114720% at the final cut. Its unshifted counterpart gives +0.05114719%. The full-pulse native response is about 78 times the directly evolved response. These are relative changes of prescribed band occupations on native coupled histories. Absolute excitation, energy deposition by the incoming pulse and feedback of the computed modes require additional calculations. The native shell gas uses its own event rule in these runs.

## The missing mirror

Reversing the primary to −0.001 gives an uncompensated endpoint $y_{\rm left}=+0.1793870$ on the fine grid. We sampled 151 coefficients across $[-0.15,0.15]$ and refined sampled interior minima for four compact remote changes with the same center and width. The smallest resolved endpoint was 0.179330 for an incident scalar, 0.179203 for a scalar velocity, 0.179337 for a metric velocity, and 0.179323 for a metric displacement. The displacement search retained 64 admissible samples after screening for characteristic degeneracy. None of these families supplied the compact mirror in the tested interval. Consequently this remote-control result is one-sided in the primary amplitude. The signed runs with two nearer bumps answer a separate question about their combined disturbance.

## Reproduction

From this folder run either check below. The quick check uses the Python standard library. The full check reads the compressed numerical wall histories, recalculates the time margins, reruns all five mode integrations with NumPy and SciPy, and rebuilds the report. Both use the sibling classical directory for the common quiet histories and prehistory.

    python -I -B replay_far.py
    python -I -B replay_far.py --full

To regenerate a native branch, supply the HDblast repository root and a new output location to run_rolling.py. Run capture on each requested grid first. The following example regenerates the fine nominal, half-amplitude and relocated arms.

    python -B run_rolling.py REPOSITORY WORK/runs --grid fine --branch capture
    python -B run_rolling.py REPOSITORY WORK/runs --grid fine --branch quiet
    python -B run_rolling.py REPOSITORY WORK/runs --grid fine --branch 0.001
    python -B run_rolling.py REPOSITORY WORK/runs --grid fine --branch 0.0005
    python -B run_rolling.py REPOSITORY WORK/runs --grid fine --branch 0.001 --control=-.9 --spread=.30 --label shifted
    python -B probe_far.py WORK --grid fine --centers=-.85,.25 --amps=-.001 --kind=incident --profile

Repeat the final command with velocity, metric and metric_field for the three other profile diagnostics. FAR_RESULTS.json gives all four cut times, comoving mode occupations and event changes. The modes and runs directories hold the individual receipts and sampled trajectories.

The fine-grid full positive response also survives shifting the initial vacuum to $s=0.50$, using the leading adiabatic analyzer, or halving the momentum spacing and tightening the solver tolerance. At $s=2.10$ these alternatives give +0.05106%, +0.05150% and +0.05111%, respectively. [Sensitivity calculation](sensitivity.py) · [All four endpoint checks](SENSITIVITY.json)
