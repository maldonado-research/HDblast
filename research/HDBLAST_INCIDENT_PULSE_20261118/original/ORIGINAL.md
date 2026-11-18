# Cross-check with the published B1 solver

We ran the historical B1 fourth-order finite-difference and Runge–Kutta solver with its registered massive-particle gas and scalar junctions. The wrapper authenticates the original four source modules, supplies the matched initial deviations after reference initialization, and leaves the original evolution, matter creation event handling and wall laws in place. The pulse and quiet control use the identity conformal chart, zero initial displacement, a uniform grid to $Z=-0.50$ and zero constraint damping. Both outer layers remain at the original quiet data. The continuum light-crossing time from that boundary to the shell is 0.50, beyond the comparison endpoint 0.18. The independent shrinking-domain calculation eliminates artificial outer-boundary inputs entirely in its retained numerical domain. [Runner](compare_solver.py) · [Source hashes](../SOURCES.json)

| Historical solver spacing | Peak scalar change from quiet | Peak recorded time | Maximum absolute difference from finest shrinking-domain trace / RK4 peak |
|---:|---:|---:|---:|
| 0.001 | 0.0010504104 | 0.12050 | 1.287% |
| 0.0005 | 0.0010492809 | 0.12025 | 0.247% |
| 0.00025 | 0.0010492710 | 0.12025 | 0.233% |

The last column compares the whole trace at the historical output times and uses linear interpolation of the reference trace where necessary. The finest peak heights agree much more closely than the maximum difference over the time series. Historical output is sampled at intervals of 0.00025 or 0.0005. All six runs reached $T=0.18$ with zero particle creation events. Every retained sample has $R=X=0$ for radiation and massive particles. These calculations use the same laws and initial pulse. Their agreement tests the evolution and field readout and does not provide an independent observation. [Comparison data and replay](COMPARISON.json)

The archived source tree and the new study are both included in the Zenodo ZIP. From this directory, with NumPy and SciPy installed, choose a fresh output directory and run

```sh
python -I -B compare_solver.py /path/to/extracted/archive /tmp/hdblast-rk4 --spacing 0.00025
```

The tabulated historical outputs in this directory are sufficient for the standard-library comparison. [Raw runs](runs) · [Waveform comparison script](compare_histories.py)
