# What happens at the limiting photon load?

The closed-universe constraint has a double root at the limiting radiation load. A double root determines the initial radius. The evolution equations determine whether the universe initially expands or contracts. On both rolling backgrounds in this study, it contracts to fourth order in time at the limiting load. A slightly lighter load allows a short expansion whose duration obeys a quarter-power law.

## Expansion at the double root

Let $x=H_0t$, let $U_i$ be the initial dimensionless potential energy, and put $F_i=-3U'_i$. The initial scalar speed and expansion rate vanish. The prepared radiation density is $r=d/a^4$. Introduce the relative load $f=4U_i d$ and $\delta=\sqrt{1-f}$. The larger radius solution of the constraint has

$$ a_i^2=\frac{1+\delta}{2U_i},\qquad \alpha=\left.\frac{a''}{a}\right|_i=U_i-r_i=\frac{2U_i\delta}{1+\delta}. $$

The acceleration equation is $a''/a=U-v^2/3-r$, where $v=q'$ and $v'(0)=F_i$. At $f=1$, $\alpha=0$. The first time derivative of the acceleration also vanishes. At that instant $r''=0$, and the second time derivative of the acceleration is $U'_i F_i-2F_i^2/3=-F_i^2$. Time-reversal symmetry then gives

$$ \frac{a(x)}{a_i}=1-\frac{F_i^2}{24}x^4+O(x^6). $$

Both archived releases have nonzero $F_i$. The result is a strict local maximum. A stationary field with $F_i=0$ would instead permit an exactly static equilibrium at this load.

## Turnaround just below the limiting load

Write $N=\log(a/a_i)$. Near the release and for small positive $\alpha$,

$$ N(x)=\tfrac12\alpha x^2-\tfrac1{24}F_i^2x^4+O(\alpha x^4+x^6). $$

Solving $N'(x_{\max})=0$ for its first positive root gives

$$ x_{\max}\sim \frac{\sqrt{6\alpha}}{|F_i|}\propto(1-f)^{1/4},\qquad N_{\max}\sim\frac{3\alpha^2}{2F_i^2}\propto(1-f). $$

Thus the gain in radius shrinks linearly in the distance from the limiting load while the duration of expansion shrinks more slowly. The prefactors depend on the archived potential and release point. These asymptotic formulas apply close to the limiting load.

## Numerical check

The solver evolves the field and the exact equation $h'=a_i^{-2}e^{-2N}-2r_ie^{-4N}-v^2/2$. Time, $N$, field displacement and velocity are rescaled to avoid subtracting nearly equal radii. At exact criticality it integrates to four positive times and compares $-N$ to $F_i^2 x^4/24$. At $x=0.0003$ the measured-to-predicted ratios are 1.0000025 on the positive branch and 1.0000017 on the negative branch. The subsequent table compares the numerically measured turnaround time and log-radius gain with the two subcritical predictions. A ratio approaching one tests the prefactor as well as the power. The calculation accepts the gap $1-f$ directly to retain precision.

| $1-f$ | Time ratio, $\phi_m=+0.5$ | Time ratio, $\phi_m=-0.5$ | Gain ratio, $\phi_m=+0.5$ | Gain ratio, $\phi_m=-0.5$ |
| ---: | ---: | ---: | ---: | ---: |
| $10^{-8}$ | 0.823492 | 0.776952 | 0.759147 | 0.695876 |
| $10^{-12}$ | 0.997113 | 0.995791 | 0.996179 | 0.994450 |
| $10^{-16}$ | 0.999971 | 0.999957 | 0.999962 | 0.999944 |

The comparison at gap $10^{-8}$ also uses the independent unscaled solver. These statements hold for the homogeneous closed Einstein equations with conserved radiation, one canonical scalar, and the specified turnaround initial data. This calculation assumes the radiation has already been prepared. See [the critical calculation](critical.py) and [the radiation evolution](photons.py).
