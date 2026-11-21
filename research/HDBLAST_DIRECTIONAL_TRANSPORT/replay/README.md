# Isotropic fermion injection on the moving birth shell

A sharp shell creates particles with comoving momentum magnitude $k$. Under the local collisionless assumptions of the HDblast sweep, a patch with $u=\cos\theta$ can receive new states at a rate proportional to $|\dot k+\dot b u|$. An isotropic source requests the same rate from every solid angle. At $mL=4\times10^{12}$ the angularly summed capacity covers the request on every tested step, but individual directions leave **28.9926 parts per million of total production unfilled** on the half-step grid. The largest instantaneous unfilled share is **25.5057%**. Both numbers use the archived action and histories. [Upstream inputs and conventions](https://github.com/maldonado-research/HDblast/tree/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120/kinetics)

## Two orderings of the angular average

For an archived step define $h=2\pi^2 S_0\dot D_c/(12\mu k^2)$, $C=|\dot k|/h$ and $D=|\dot b|/h$. Let $\langle\cdot\rangle$ denote the uniform average over $-1\le u\le1$. Directional depletion under isotropic injection leaves a fractional shortage $A(C,D)=\langle[1-|C+D u|]_+\rangle$. Summing capacity first gives the archive’s $[1-\langle|C+D u|\rangle]_+$. Convexity of the positive part puts the angularly resolved shortage at least as high as the summed shortage for the same histories.

For $D>0$ the angular integral is the area of the triangle $[1-|v|]_+$ inside $[C-D,C+D]$, divided by $2D$. The fraction of solid angle seeing any shortage is the overlap length of that interval with $[-1,1]$, divided by $2D$. At $D=0$ the missing fraction is $[1-C]_+$. The minimum over direction of the normalized capacity is $[C-D]_+$. For $D\ge C$ a tangent direction has zero capacity. The integrated shortfall counts the unmet part of the request across a band of nonzero solid angle.

The replay uses the archive’s accepted updates, birth-shell speeds, step-centered $k$, cumulative rest increments, clipping to the two production windows and signed field interpolation. The temporal integral splits at the field zero and at absolute fields $D=|C-1|$, $C$ and $C+1$. Above $C+1$ the shortage is exactly $1/(2D)$, giving a logarithmic integral in time. Below that final edge an eight-point Gauss rule integrates each smooth piece, with further subdivision whenever the nonzero endpoint ratio exceeds two. The code tests constant fields, reversals, a one-sided ramp, a broad reversal and the interior peak against analytic controls. A direct replay of the archived sweep verifies the summed results independently of the angular readout. [Angular and temporal integration](directional.cpp) · [Trajectory replay](replay.py)

For $C\le1/2$ a maximum of $A$ occurs at zero field. For $C>1/2$ the maximum occurs at $D=\sqrt{(1+C)^2-2}$ when the step crosses that value. To locate a peak, the replay checks the endpoints and both signed stationary fields. At the higher mass, the half-step peak occurs at $\tau=6.6183483613$ and $6.6183484828$, on either side of the reversal at $6.6183484220$. The radial expansion alone suffices at that exact reversal, where $C=1.21540$. Moving away from zero creates a band of deficient directions as the electric sweep cancels the radial sweep. At the lower mass the maximum includes the reversal at $\tau=6.6183518247$.

## Production fractions at the archived actions

The lower-mass action is $3.6085895820\times10^{48}$. The higher-mass action is $1.0415592373\times10^{54}$, five times its reference drag floor. Every cumulative percentage in the table divides by the archived reference final production. The peak divides by the instantaneous requested rate. Normalizations conditional on the tested births are available in the machine-readable output. [Complete results](results.json)

| $mL$ | Central production window | Clock | Summed unmet, % of total | Isotropic unmet, % of total | Peak instantaneous isotropic unmet, % |
|---|---:|---|---:|---:|---:|
| $4\times10^{10}$ | 98% | native | 0.1051637 | 0.8199783 | 99.64920 |
| $4\times10^{10}$ | 98% | half | 0.1051631 | 0.8200053 | 99.64919 |
| $4\times10^{10}$ | 99.998% | native | 0.1051637 | 0.8203913 | 99.64920 |
| $4\times10^{10}$ | 99.998% | half | 0.1051631 | 0.8204183 | 99.64919 |
| $4\times10^{12}$ | 98% | native | 0 | 0.002898008 | 25.50596 |
| $4\times10^{12}$ | 98% | half | 0 | 0.002898070 | 25.50571 |
| $4\times10^{12}$ | 99.998% | native | 0 | 0.002899200 | 25.50596 |
| $4\times10^{12}$ | 99.998% | half | 0 | 0.002899262 | 25.50571 |

The higher-mass half-step calculation places 58.0392 parts per million of total requested births in directions with some shortage in the widest window. Of that population, 28.9926 parts per million of total production exceed the local capacity. The nominal tails contain 20 parts per million of the archived reference yield. The higher-mass half-step endpoint exceeds that reference by 0.1291 parts per million. Counting the untested tails at the two limiting shortage fractions, zero and one, brackets the full-history shortfall between 28.9926 and 49.1217 parts per million of the actual half-step yield. The calculation assumes favorable availability on both sides of the surface with a sharp momentum and no collisions during injection. Previously occupied states, finite spectral width, anisotropy and feedback on the fluid require additional dynamics.

## Action dependence on the archived higher-mass history

The local request $h$ is proportional to $S_0$ with the trajectory fixed. The table uses the half-step grid and central 99.998% of production. Multiples of the archived drag floor $S_d=2.0831184746\times10^{53}$ are conditional rate benchmarks. All actions, both masses, cumulative fractions and peak times are included in the output. [Action scans](results.json)

| Action | Maximum summed demand / capacity | Isotropic unmet, ppm of total | Peak instantaneous isotropic unmet, % |
|---|---:|---:|---:|
| $0.1$ times archived action | 0.08228 | 2.8967 | 3.8122 |
| $S_d$ | 0.16455 | 5.7936 | 7.1371 |
| Archived action, $5S_d$ | 0.82277 | 28.9926 | 25.5057 |
| Summed action ceiling, $1.2659138001\times10^{54}$ | 1.00000 | 35.2540 | 29.2893 |
| $12S_d$ | 1.97465 | 69.9712 | 49.3741 |
| $36S_d$ | 5.92396 | 209.646 | 83.1194 |
| $10$ times archived action | 8.22773 | 288.020 | 87.8460 |

## Replay

Use Python 3.9+ and a C++17 compiler available as `g++`. The `--archive` argument points to the upstream capsule in the indicated revision. The script verifies the hashes of its seven inputs before running either history. Temporary files hold roughly 300 MiB of accepted-step records. The optional `--slices DIRECTORY` retains them for independent checks.

```sh
git clone https://github.com/maldonado-research/HDblast.git
git -C HDblast checkout --detach dc99ff5c5894abcff3e11e02e4b97e19da31de00
python -I -B replay.py \
  --archive HDblast/research/HDBLAST_COHERENT_DECAY_20261120 \
  --output /tmp/directional-results.json
```

Slice records use native-endian binary64 values in the order given by `COLUMNS` in the replay. The inclusive endpoints record shell radii and signed center speeds. The original step supplies the unsigned radial speed and centered radius for each window. A filename ending in `-central98.bin` identifies the 98% window. The other binary records cover the central 99.998%.
