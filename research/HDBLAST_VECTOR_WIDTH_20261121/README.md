# A narrower vector interaction can still reduce the beam pressure

At shell time 7.05, a width parameter of **0.043** gave an initial pressure-change rate of **−0.4342 ± 0.0012** when the exchange diagrams were included. The isolated resonance at that same width gave **−0.01039 ± 0.00068**. Only **0.000854%** of the additional decrease came from collision energies above the tested partial-wave bound. For comparison, the original isolated resonance at width 0.30 gave −0.3595 ± 0.0057 on this same population. This identifies a narrower coupling worth evolving through the pulse. [Complete numerical results](RESULTS.json) · [Replay code](narrow.py)

The initial distribution contains 16,384 positive and 16,384 negative particles of one species. Its normalized directional pressure is $Q=(P_{zz}-\mathrm{tr}P/3)/\rho=0.43039$. Each entry below is the instantaneous collision contribution to $dQ/d(t/L)$. More negative values give faster initial loss of pressure along the primary field direction. The mass ratio is $M/m=8300$ throughout. Changing $w$ changes both the coupling $g^2=12\pi w$ and the timelike resonance width in the propagator. [Published stress generator and definitions](https://github.com/maldonado-research/HDblast/blob/378a5223fa764b29a991b8df6ac825a051f86277/research/HDBLAST_CROSSED_VECTOR_20261121/crossed-stress/README.md)

| Width $w$ | Threshold $\sqrt{s}/M$ | Isolated resonance, total | With exchange, total | Added decrease from above threshold |
|---:|---:|---:|---:|---:|
| 0.065 | 12.9504 | −0.02126 ± 0.00093 | −0.98981 ± 0.00221 | 13.72% |
| 0.043 | 48.2190 | −0.01039 ± 0.00068 | −0.43420 ± 0.00116 | 0.000854% |
| 0.035 | 116.9642 | −0.00751 ± 0.00062 | −0.28829 ± 0.00089 | 0% |

| Width | Isolated below | Isolated above | With exchange below | With exchange above |
|---:|---:|---:|---:|---:|
| 0.065 | −0.02037 | −0.000896 | −0.85604 | −0.13377 |
| 0.043 | −0.010394 | −1.6706 × 10⁻⁸ | −0.434201 | −3.6343 × 10⁻⁶ |
| 0.035 | −0.007509 | 0 | −0.288286 | 0 |

The threshold solves $(w/0.30)|\mathrm{Re}\,a_0(s;0.30)|=1/2$ for the specified equal-helicity diagonal partial wave. At width 0.043 every recorded encounter in the earlier finest time-integrated sample passed this necessary check. The *possible* pairs at the final time still include 825 unlike pairs and 731 unordered same-sign pairs above the threshold. We counted and integrated every one of these rare pairs. Increasing the tail integration from 18 polar nodes and 32 azimuths to 24 nodes and 48 azimuths changed its full pressure contribution by under $2\times10^{-10}$ relative. At width 0.035 even the conservative ceiling $\sqrt{s}/M\le2\max\gamma/(M/m)=71.8311$ is below the threshold, so the above-threshold contribution vanishes on this finite snapshot. [Partial-wave derivation](https://github.com/maldonado-research/HDblast/blob/378a5223fa764b29a991b8df6ac825a051f86277/research/HDBLAST_CROSSED_VECTOR_20261121/unitarity/README.md)

The ordered above-threshold counts at width 0.065 were 7,824,425 opposite-charge, 7,677,916 positive–positive and 7,927,644 negative–negative, out of $16{,}384^2$ entries for each charge product. At width 0.043 the corresponding counts were 825, 760 and 702. The two same-sign products each receive the usual one-half incoming-pair factor. The exact below- and above-threshold counts are distinct from the sampled stress integrals at width 0.065.

The uncertainties are integration standard errors from eight matched randomized batches, each with 32,768 incoming pairs per charge combination. Exhaustive rare-tail contributions are added separately. The original particle trajectories were generated at width 0.30 with the isolated pole. Evolving the narrower interaction would change the particles, their current and their subsequent collisions. Satisfying this one necessary real-part test does not establish unitarity of every helicity eigenchannel, perturbative convergence, or equilibration. Other helicity amplitudes and angular momenta, higher-order corrections and a self-consistent time evolution each remain to be tested.

## Reproduce

The calculation uses the files in commit `378a5223fa764b29a991b8df6ac825a051f86277`, Python 3.12, NumPy 2.3.5 and SciPy 1.17.0. In a checkout containing the parent studies, run these commands from the new study directory that contains `narrow.py`. They regenerate the parent inputs and run the measurement. With existing event tapes, only the last command is needed. The JSON records each batch, both regions for both amplitudes, the exhaustive tail counts, angular refinements and the SHA-256 digest of the input momentum matrix.

~~~sh
research=..
work=/tmp/hd-width-replay
mkdir -p "$work"
python -I -B "$research/HDBLAST_AXIAL_CURRENT_20261121/fields/extract.py" \
  --study "$research/HDBLAST_POLARIZATION_20261121" --cases dedicated_half --end 7.0505 \
  --output-dir "$work" --output "$work/fields.json" --jobs 2
OPENBLAS_NUM_THREADS=1 python -I -B "$research/HDBLAST_COLLISION_PHOTONS_20261121/events/replay.py" \
  --axial "$research/HDBLAST_AXIAL_CURRENT_20261121" --polarization "$research/HDBLAST_POLARIZATION_20261121" \
  --scattering "$research/HDBLAST_SCATTERING_20261121" --history "$work/dedicated_half.npz" \
  --artifacts "$work/events" --output "$work/events.json" --jobs 2
OPENBLAS_NUM_THREADS=1 python -I -B "$research/HDBLAST_VECTOR_WIDTH_20261121/narrow.py" --study "$research/HDBLAST_CROSSED_VECTOR_20261121" \
  --snapshot "$work/events/t7.05000-p14-n8192-r91283-particle-pole.npz" \
  --history "$work/dedicated_half.npz" --output "$work/widths.json"
~~~

[Reproduction script](narrow.py) · [Complete results and sampling checks](RESULTS.json)
