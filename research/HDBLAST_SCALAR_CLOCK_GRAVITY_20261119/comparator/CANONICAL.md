# Which four-dimensional sources can reproduce the expanding shell?

On the simulated interval, the scale factor and expansion rate together exclude a broad class of flat four-dimensional sources. They also allow an explicit four-dimensional canonical scalar that reproduces the same homogeneous expansion. This exploratory comparison fits the measured simulation histories separately. It uses the same physical proper time as the shell calculation.

## A bound using only two endpoints

Assume spatially flat, homogeneous and isotropic four-dimensional Einstein gravity, constant reduced Planck mass $M_4$, expansion $H>0$ and a total isotropic source. Here $a$ denotes the ordinary three-dimensional scale factor, all constant vacuum terms are included in the total stress, and a dot denotes a physical proper-time derivative. Combining the Friedmann and acceleration equations gives

$$
\frac{d(Ha^3)}{dt}=a^3\left(3H^2+\dot H\right)=\frac{a^3}{2M_4^2}(\rho-p).
$$

The dominant energy condition for a perfect fluid includes $p\leq\rho$. Any such source present throughout the interval must satisfy $(Ha^3)_2/(Ha^3)_1\geq 1$. This includes every expanding mixture of positive-density radiation, dust and canonical scalar matter with nonnegative total scalar potential. The inequality uses only endpoint geometry and applies even when the components exchange energy. The result follows directly from the two gravitational equations. The corresponding energy-condition bounds have also been studied in observational cosmology. [Integrated bounds and their interpretation](https://arxiv.org/html/0802.0706#S2)

| Retained quiet wall, $s=1.60$ to $2.10$ | $(Ha^3)_2/(Ha^3)_1$ | Effective $w=p/\rho$ averaged over $d\log a$ |
|---|---:|---:|
| Fine, $\lambda=0$ | 0.540246 | 3.811225 |
| Fine, $\lambda=0.01$, static calibration, cutoff 120 | 0.540245 | 3.811235 |
| Coarse, $\lambda=0.01$, same calibration and cutoff | 0.540270 | 3.811692 |

For the fine positive-step run the native rates $H/H_0$ at those times are 0.506829 and 0.176688, and $\Delta\log a=0.146017$. The average is $-1-2\log(H_2/H_1)/(3\Delta\log a)$ from four-dimensional energy conservation. On the nearby positive, nearby negative and far pulsed histories, the corresponding fine ratios with $\lambda=0.01$ are 0.540118, 0.540867 and 0.540287. The table interpolates $H$ and $\log a$ independently to common endpoint clocks. As a derivative-free check that also dispenses with interpolation, take the first *recorded* sample at or after 1.60 and the last at or before 2.10. The ratio is between 0.54029 and 0.54121 across all 16 grid, arm and coupling combinations. The large violation of the lower bound is already present at zero Yukawa coupling. Its meaning concerns a four-dimensional interpretation of the total expansion. It does not diagnose a five-dimensional energy-condition violation or measure the pressure of the emitted pairs. [Numerical histories and bounds](CANONICAL.json)

The exclusion extends to open four-dimensional FLRW spatial sections. With curvature convention $H^2+k/a^2=\rho/(3M_4^2)$, the derivative of $Ha^3$ gains $-2ka$, which preserves its lower bound for $k\leq0$. An independent check follows from matter conservation, $d(a^6\rho)/dt=3Ha^6(\rho-p)$. For $r=a_2/a_1>1$, the two endpoints therefore require $k/(a_1^2H_0^2)\geq(h_1^2-r^6h_2^2)/(r^4-1)$ if the dominant energy condition holds throughout. The fine quiet positive-step endpoints give a strictly positive lower bound of 0.22929. A closed model would also have to satisfy the energy condition at intermediate times and would have a different intrinsic spatial curvature from the simulated flat shell. The scalar reconstruction below uses the flat geometry.

## A canonical scalar construction

The effective four-dimensional total stress obeys the weaker null energy condition on a decreasing $H(t)$ history. Write $s=H_0t$, $h=H/H_0$, $u=\Phi/M_4$ and $U=V/(M_4^2 H_0^2)$. Choose a canonically normalized minimally coupled homogeneous scalar as the only source and take the increasing-field branch. The inverse equations are

$$
u'(s)=\sqrt{-2h'(s)},\qquad U(u(s))=3h(s)^2+h'(s).
$$

Here a prime in these two expressions means a derivative with respect to $s$. Given $h'<0$, the field is strictly monotonic. The parametric function $U(u)$ is then single-valued on the swept interval. Its energy density is positive, $\rho/(M_4^2H_0^2)=3h^2$. When $U$ becomes negative, the canonical field can have $p>\rho$ and violate the dominant energy condition while retaining positive kinetic energy. The squared propagation speed of a canonical scalar is one. The ratio $p/\rho$ inferred from the background has a different physical meaning. Potential reconstruction from expansion data has a long history. Starobinsky also considers separately conserved dust and observational inverse problems. [Reconstruction method](https://arxiv.org/abs/astro-ph/9810431)

We interpolated each native $h(s)$ with a cubic spline and located the extrema of $h'$ on $1.35\le s\le2.10$ from both the knots and the interior roots of $h''$. Across all 16 histories the interpolated derivatives lie between −1.243 and −0.1725. A second shape-preserving cubic interpolation has negative derivative throughout the same interval. These are checks on interpolants of the numerical histories. In the fine quiet positive-step run $U$ crosses zero at $s=1.500441$ and reaches a minimum near −0.598615. All 16 cubic fits have one zero on the specified interval. The minimum and zero are derivative-based estimates. The two-endpoint exclusion above requires neither a derivative estimate nor a sign for an interpolated potential.

We then integrated the four-dimensional Klein–Gordon equation and the Raychaudhuri equation forward for the two fine quiet histories. We used only the fitted $U(u)$, the initial scalar velocity, initial expansion and initial proper time in this forward evolution. In dimensionless variables, $v=u'$, $v'=-3hv-dU/du$, $h'=-v^2/2$, and $(\log a)'=h$. We use $u$ as the independent integration variable so each Runge–Kutta step respects the knots of $U$. The Friedmann constraint supplies a separate residual. In both fine quiet runs the maximum absolute difference in $h$ at the fitted scalar knots is below $2\times10^{-13}$, and the Friedmann residual is below $2\times10^{-13}$ at both integration resolutions. These checks validate the reconstruction of the fitted history; its physical precision is limited by the original five-dimensional simulation.

Each potential is determined from the whole expansion history it reproduces. A predictive four-dimensional competitor would need a specified common potential, a description of how incoming disturbances are prepared, and predictions for additional observables. The homogeneous metric alone permits the canonical construction above. The spot checks use the zero-coupling history and the physical positive step of the finite-cutoff response calculation. [Effective four-dimensional descriptions and primary literature](../literature/REFERENCES.md)

To reproduce after unpacking the bulk histories, run the following from this study directory.

~~~sh
python -I -B bulk/pack_runs.py unpack --runs /tmp/hdblast-native-runs
OPENBLAS_NUM_THREADS=1 python -I -B comparator/canonical.py --runs /tmp/hdblast-native-runs --output /tmp/canonical.json
~~~
