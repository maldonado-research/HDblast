# When the equatorial sphere stops crossing all at once

In the radial Einstein–scalar simulation the equatorial two-sphere reaches the potential minimum at one proper time. Its crossed volume grows initially as $C_{\rm eq}\sqrt{\delta x}$, with $C_{\rm eq}=2.083$–$2.087$ for the two negative-amplitude runs. To test the sensitivity of this result to angular coherence, prescribe a spread $\tau$ in the onset times along that two-sphere. This calculation varies the local geometry of a scalar level set. The angular timing pattern and $\tau$ are inputs. [Radial first contacts](../patches/PATCHES.md)

Write $\xi$ for signed radial angular distance from the equator and $\mathbf n$ for a unit vector labeling a point on it. Retain the velocity $v>0$, transverse curvature $c>0$ and proper measure from the radial solution. A prescribed shift $\tau g(\mathbf n)$ in the local onset time gives the leading local profile

$$
q-q_m=v[\delta x-\tau g(\mathbf n)]-c\xi^2/2,\qquad
f\simeq C_{\rm eq}\sqrt{\tau}\,A_g(s),\qquad
A_g(s)=\left\langle\sqrt{[s-g(\mathbf n)]_+}\right\rangle_{S^2},
\qquad s=\delta x/\tau.
$$

Here $\delta x$ is measured from the earliest contact and $0\leq g\leq1$. The square root gives the half-thickness of the crossed interval in the direction transverse to the sphere. The angular average supplies the fraction and distribution of directions already in motion beyond the minimum.

## A closed-form example

An axisymmetric quadrupole plus a constant gives $g=1-n_z^2=\sin^2\theta$. Two antipodal points reach the minimum first. With $u=|\cos\theta|$, the proper area average on the round contact sphere gives $A(s)=\int_0^1\sqrt{[s-1+u^2]_+}\,du$. Its closed form is

$$
A(s)=\tfrac12
\begin{cases}
\sqrt s-(1-s)\operatorname{atanh}\sqrt s, &0<s<1,\\
1,&s=1,\\
\sqrt s+(s-1)\operatorname{atanh}(1/\sqrt s),&s>1 .
\end{cases}
$$

For $s\ll1$, $A(s)=s^{3/2}(1+s/5+O(s^2))/3$. The initially crossed volume is proportional to $(\delta x)^{3/2}/\tau$. For $s\gg1$, $A(s)=\sqrt s[1-1/(3s)+O(s^{-2})]$. A window for the square-root law therefore requires the time spread to be much shorter than the elapsed time, while the elapsed time is still short enough for the radial contact expansion to apply. The order of the two limits matters. At $s=1$ the last-to-cross set is an entire circle and $A'(s)$ has a logarithmic singularity. This feature depends on the axisymmetric timing pattern.

| Elapsed time divided by the prescribed spread | Crossed volume relative to the simultaneous-sphere prediction |
|---:|---:|
| 0.1 | 3.403% |
| 1 | 50.000% |
| 10 | 96.597% |
| 100 | 99.666% |

These percentages compare two predictions at the same elapsed time. They are ratios to $C_{\rm eq}\sqrt{\delta x}$. The value of $\tau$ for nonspherical HDblast data would require independent initial data and a new Einstein–scalar evolution.

## A triaxial check

A general real quadrupole can be diagonalized. After shifting its minimum and normalizing its full range, consider $g_\eta=\eta n_y^2+n_z^2$, $0<\eta\leq1$. The two earliest points are $\pm\hat x$. In local tangent coordinates their angular delays are quadratic with eigenvalue ratio $\eta$. Integrating the square-root transverse thickness over the two growing angular ellipses gives $A_\eta(s)\sim s^{3/2}/(3\sqrt\eta)$ for fixed $\eta$ and $s\ll\eta$. At late scaled times the angular mean $\langle g_\eta\rangle=(1+\eta)/3$ gives $A_\eta(s)=\sqrt s[1-(1+\eta)/(6s)+O(s^{-2})]$.

For $\eta=1/2$, the numerical crossed volume reaches 4.866%, 67.532% and 97.459% of the simultaneous prediction at $s=0.1,1,10$ respectively. Both $\eta=1$ and $\eta=1/2$ give log slopes within $1.2\times10^{-6}$ of $3/2$ between $s=10^{-6}$ and $10^{-5}$. For each azimuth the calculation can use the closed form above with the local angular amplitude, or directly integrate the original square-root integrand. At 13 scaled times the two methods agree to $9.0\times10^{-14}$ relative or better. Neither the timing spread nor the actual angular pattern follows from the radial data. [Numerical checks](UNFOLDING.json) · [Replay](unfolding.py)
