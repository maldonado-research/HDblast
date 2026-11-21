# Can an anisotropic birth law preserve each stream’s zero directed momentum?

At the archived $mL=4\times10^{12}$ action, the source requests **1.680 times** the maximum angular birth rate compatible with zero directed momentum in each charge stream at two instants near an electric reversal. Summing over all angles at either instant gives demand equal to 0.790 times capacity. This constrains even a freely adjustable angular birth law on the prescribed sharp, collisionless history. [Accepted-step witness](MOMENT.json)

Write $c=|\dot k|$, $d=|\dot b|$, $h=P/(12k^2)$ and $u=\cos\theta$, reflecting the angular coordinate as needed. In units of $c$ the available directional flux is $g(u)=|1+ru|$, with $r=d/c$. For an admitted density $0\le j\le g$, the stream’s source moment is proportional to $\langle u j\rangle$ and its total to $\langle j\rangle$, where the average uses $du/2$. Choose $z\in(0,1)$ such that $\int_{-1}^z u g(u)\,du=0$. Among all zero-moment laws the maximum mean is $F_0/c=\tfrac12\int_{-1}^z g(u)\,du$. The bound follows by writing $\langle j\rangle=\langle(1-u/z)j\rangle\le\langle(1-u/z)_+g\rangle$. Filling all available angles up to $z$ attains it. An azimuthal average preserves the axial constraint and total rate. The attaining axisymmetric law also has vanishing transverse moments. [Angular variational argument](../ANGULAR.md)

The replay finds all five signed field reversals in the widest production window and selects the crossing step with the largest $h/c$. Both signed values $\dot b/c=\pm4/3$ occur inside that step. At $r=4/3$ the cusp is at $u=-3/4$. Integrating the linear and quadratic pieces of $g$ and $ug$ gives

$$
0=\frac{z^2+1}{2}+\frac r3(z^3-1)-\frac1{3r^2},\qquad
\frac{F_0}{c}=\frac12\left[z-1+\frac r2(z^2+1)+\frac1r\right].
$$

Rational bisection gives $z=0.4360871053855383$ and $F_0/c=0.489767540520615$. Separate two-point Gaussian quadratures of the unsimplified clipped functions reproduce the capacity and give a mean axial moment of $-1.4\times10^{-17}$. On the half-step grid the witness occurs at proper times 6.6183483643 and 6.6183484798. With $h/c=0.8227726384$, the angularly summed and zero-moment capacities give the two demand ratios above. The native grid gives 1.67995 for the zero-moment ratio.

At fixed history, $h$ scales with $S_0$. This pair of instants implies $S_0\le6.20003\times10^{53}$ on the half-step grid. At that upper limit the reference drag already requires **33.5985%** of its thermal single-channel ceiling. The equal-flavor scalar pole’s inclusive ceiling lies lower by a factor of **4.03182**, and its same-flavor elastic ceiling by **12.0955**. The computed scalar fluctuation energy from $m/8$ through $64m$ is **1.50454%** of the classical outgoing parent energy at that action. Both required fractions increase as $S_0$ decreases. This action bound uses one ratio of field to radial speed. Additional times, prior occupancy and other physical constraints can further restrict the source. [Transport assumptions](../transport/DRAG.md) · [Scalar pulse normalization](../pulse/PULSE_MODES.md)

From this directory, save the accepted-step traces from the pinned archive and run the moment check with Python’s standard library.

```sh
python -I -B replay.py --archive "$capsule" --slices "$slices" --output /tmp/directional-results.json
python -I -B moment.py --slices "$slices" --output /tmp/directional-moment.json
```

Here `capsule` is the archived HDblast coherent-decay study directory and `slices` is a writable output directory. The witness receipt records each trace’s SHA-256 hash, the crossing step, both interpolation times, the exact rational bracket of the angular cut and the native-grid control. [Reproduction code](moment.py) · [Results](MOMENT.json)
