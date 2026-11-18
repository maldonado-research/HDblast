# The area of two finite light cones on the HDblast shell

Two patches of a receiving observer’s sky each cover one tenth of the full sky. Trace their light rays into the past by the same observer-normalized affine length $L$ and measure the areas they cover. The registered shell metric predicts that **one patch has a larger area than its de Sitter value and the other has a smaller area** for every sufficiently small, nonzero canonical perturbation. Reversing the perturbation reverses the two signs. This is a GitHub supplement to the archived joint-curvature study. [Underlying model and data](https://doi.org/10.5281/zenodo.24352632) · [Optical derivation](DERIVATION.md)

For a sky patch of solid angle $\Omega$, define $D_{\rm cap}=\sqrt{\mathrm{area}/\Omega}$. In the background $D_{\rm cap}=L$. Normalize rays to unit frequency in the receiving observer’s rest frame. The observer's dimensionless elapsed clock is $q_0=49\ell/32$, with $\ell=10^{-6}$. The chosen affine length ends at $q_1=99\ell/64$ in the background. In the perturbed geometry the endpoint is selected by $L$ itself. The positive cap faces the first of the two registered quadrupole axes. Its negative partner faces the second.

| First variation at zero canonical amplitude | Outward numerical enclosure for the positive direction |
| --- | --- |
| Central pencil, in the model normalization $\mu\rho_b^4/(\ell L^3)$ | $(0.0811747,\;0.0826894)$ |
| Finite cap, in the same normalization | $(0.0584458,\;0.0595364)$ |
| Fractional distance $L^{-1}\,dD_{\rm cap}/d\epsilon$ per unit canonical amplitude | $(3.7210,\;4.1805)\times10^{-25}$ |

A third check combines the optical endpoint and the weighted integral before applying the input-cell enclosures. It tightens the central-pencil interval to **$(0.0817743,\;0.0820899)$** in the model normalization. [Combined-kernel check](kernel_check.py) · [Exact result](KERNEL.json)

The last row gives the scale of the prediction. At small amplitude $\epsilon$, the fractional distance change is this coefficient times $\epsilon$, plus $o(\epsilon)$. The fractional area change has twice the distance coefficient. The calculation establishes the sign for sufficiently small amplitudes but gives no numerical amplitude range or instrumental sensitivity. The much larger numbers in the first two rows reflect an explicit mathematical rescaling. [Exact coefficients](BEAM.json) · [Independent sign calculation](OPTICS.json)

**Why the sign matters.** Four-dimensional Einstein gravity with positive Newton constant and matter satisfying the pointwise null energy condition bounds $D_{\rm cap}$ from above by $L$ for an initially regular cone before its first conjugate points. The cap above $L$ violates this bound for the specified shell propagation and calibration. Four-dimensional models with different gravitational equations, matter violating the condition, or a different metric for light remain possible. Turning the construction into an astronomical test would require an independent affine-length calibration and information about physical source areas or equivalent distance measurements. Equal redshift does not in general set equal affine length. [Focusing equation and comparison assumptions](DERIVATION.md)

## Reproduce

Use Python 3.12 and the standard library. Keep the sibling directory from the archived joint-curvature study at its repository location. Run the following commands from this directory.

~~~sh
python -I -B optics.py ../HDBLAST_JOINT_CURVATURE_20261118/METRIC.json > /tmp/hdblast-optics.json
cmp OPTICS.json /tmp/hdblast-optics.json
python -I -B beam_distance.py > /tmp/hdblast-beam.json
cmp BEAM.json /tmp/hdblast-beam.json
python -I -B kernel_check.py > /tmp/hdblast-optical-kernel.json
cmp KERNEL.json /tmp/hdblast-optical-kernel.json
~~~

The first calculator uses uniform curvature enclosures over the archived time window. The second follows the individual waveform cells, propagates their interval bounds through the light-cone integral and converts the cap result into a fractional distance. The third calculator integrates one positive kernel against the same waveform cells. The published fractions give the exact endpoints. The decimal display is rounded. [Uniform replay](optics.py) · [Cellwise replay](beam_distance.py) · [Combined-kernel replay](kernel_check.py)
