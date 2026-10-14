# Geometric shell-contact schematic

`contact-support.png` and `.pdf` illustrate the support bounds for the
constructive finite-angular theorem, using x=pi/2-chi, a=tau*/4 and
d0=tau*/3. The particular drawing uses tau*=pi/8; it is not a numerical
certificate of the background's negative-s extension or an evolved waveform.
The theorem may require a smaller scale, without changing the contact logic.

Blue shading is possible causal support, not a claim of a nonzero field at
every shaded point. Smooth data are flat at their support edge; the proven
boundary trace is nonzero somewhere in each sufficiently short interval after
arrival, and need not be nonzero at the exact arrival instant.

The dashed curve is T=T0=b tan(tau*/2), expressed as
 tau/tau*=1/2+asin(sin(tau*/2) sin(x))/tau*.
The whole moving-ball ambient slab through T0 lies within the controlled
hemisphere rectangle. Arrival at tau=a occurs strictly before T0 at the shell.

`make_contact_figure.py` regenerates the figures and the support-bound CSV.
It needs NumPy and Matplotlib and is separate from the exact symbolic replay.
No scientific conclusion depends on raster appearance.

For a cloud workspace with a read-only default cache location, use task-specific
cache paths (the home directory need not be changed):

```sh
MPLCONFIGDIR=/tmp/hdblast-contact-mpl XDG_CACHE_HOME=/tmp/hdblast-contact-cache python make_contact_figure.py
```
