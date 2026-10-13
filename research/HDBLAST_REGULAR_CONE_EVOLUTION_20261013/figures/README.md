# Finite-propagation schematic

`support-geometry.pdf`/`.png` show a geometric support bound, not numerical
wave amplitudes or a computed background. Units are the freely chosen
local scale t0, and the illustrative shell ratio is r_b/t0=4.5. The proof
separately requires a sufficiently small t0 for the analytic neighborhood.
CSV columns specify the plotted bounds. Generate optional replacement
figures in a scratch directory with NumPy and Matplotlib installed:

    MPLCONFIGDIR=/tmp/hdblast-mpl python make_support_figure.py --output /tmp/hdblast-cone-figure

Figure generation is not part of the pinned algebra replay. PDF metadata
may change on regeneration. Do not overwrite the frozen research capsule.
