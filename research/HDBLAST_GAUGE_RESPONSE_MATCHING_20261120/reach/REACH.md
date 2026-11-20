# How large a coupling change fits inside the calibrated model?

One way to produce more light is to change its coupling to the rolling field more strongly. The same change also reduces the kinetic coefficient left over after we calibrate the response of the short waves. This gives a quantitative restriction on how much the coupling can change along a specified path.

## A sharp bound along one descent

Write $g(q)=d\log I/dq$, $\delta=\log(I_{\rm out}/I_{\rm in})$, and $Z=1-cg^2/a^2$, with $c>0$ fixed by the added band. Let $z_*$ be the minimum permitted kinetic coefficient, with $0<z_*<1$. On a monotone crossing between the endpoints of the switch, $a(q)$ is the radius when the field passes $q$. The local condition $Z\ge z_*$ implies

$$
|\delta|=\left|\int g(q)\,dq\right|
\leq \int |g(q)|\,|dq|
\leq \sqrt{\frac{1-z_*}{c}}\int a(q)\,|dq|.
$$

For a specified positive continuous radius history, this upper bound is sharp. A slope of constant sign proportional to $a(q)$ saturates the pointwise inequality. Turning it smoothly to zero in successively thinner endpoint layers approaches the bound while matching constant couplings outside. Requiring smooth endpoints and a strict reserve generally leaves a supremum. The construction treats the radius history as given. Changing the coupling changes the history once the light and field exchange energy.

With $z_*=0.05$, the original peaked profile uses 38.45% and 38.46% of the possible contrast along the two **archived, prescribed** paths. For example, at $K=192$, even the best shape along those paths allows only $|\delta|\leq0.9333$ and $0.7313$. The strong example requires $|\delta|=3$.

## A necessary test even if the trajectory changes

On the first expanding segment the largest radius up to switch completion is $a_{\rm end}$. Every field value must be crossed at least once, including if the field temporarily reverses. The coupling is a single-valued function of the field. At some crossing of each value, $|g(q)|\le a_{\rm end}\sqrt{(1-z_*)/c}$. A necessary radius is therefore

$$
a_{\rm end}\ \geq\ a_{\rm req}
=\frac{|\delta|}{|q_{\rm off}-q_{\rm on}|}\sqrt{\frac{c}{1-z_*}} .
$$

Total conservation implies that energy density cannot increase during this expanding segment. Indeed $\rho+p=(2Z+aZ_a/3)v^2/6+4r_\gamma/3\ge0$ while the reserve holds. In the dimensionless time $x=H_0t$, let $u_i$ be the total density at the stationary launch with $a_i=1/\sqrt{u_i}$. The closed Friedmann equation gives $(da/dx)^2\le u_i a^2-1$. Integrating from launch yields $a(x)\le a_i\cosh(\sqrt{u_i}x)$. Thus

$$
x_{\rm end}\ \geq\ \frac{\operatorname{arcosh}(\max(1,a_{\rm req}/a_i))}{\sqrt{u_i}} .
$$

At $K=192$ the strong excursion requires radii of at least 3.487 and 4.536. Its earliest possible completions under these assumptions are $x=1.921$ and $2.193$. For scale, the two original prescribed trajectories had already reached their first minima at $x=0.517$ and $0.541$. This test leaves open a sufficiently delayed crossing on a different trajectory. Because $c$ grows quadratically with $K$, this specific positive-reserve prescription cannot keep a fixed nonzero contrast on a fixed finite-radius, finite-width crossing while sending $K$ to infinity.

## Trying a broader coupling in the coupled equations

For this exploratory check we used $u=(q-q_{\rm on})/(q_{\rm off}-q_{\rm on})$ and chose $d\log I/du$ proportional to $a_{\rm ref}(u)T(u)$. The normalization fixes $\int_0^1 (d\log I/du)\,du=\delta$. The smooth taper $T$ rises from zero over the first fifth of the switch, stays one in the middle, and falls to zero over the final fifth. The radius comes from the archived uncoupled path and is interpolated once. The resulting field-only function stays fixed while the coupled model evolves its own radius. The aim of this choice was to spread the required slope. We did not select the shape by optimizing the photon share.

On the reference histories it uses 80.05% and 80.06% of the upper bound. In the fully coupled runs with $\delta=3$, we obtained the following photon percentages of total scalar-plus-photon energy at the first potential minimum.

| Highest harmonic | Positive minimum | Lowest $Z$ | Negative minimum | Lowest $Z$ |
| ---: | ---: | ---: | ---: | ---: |
| 48 | 81.40% | 1.000 | 65.01% | 1.000 |
| 56 | 74.09% | 0.632 | 61.56% | 0.406 |
| 60 | 82.15% | 0.429 | 69.81% | 0.072 |
| 64 | 76.32% | 0.209 | stopped | 0.050 |
| 72 | stopped | 0.050 | stopped | 0.050 |

At $K=60$ both outgoing states stay above half until the first maximum radius or the field's first re-entry into the switching interval. Over these respective intervals the radius grows by 2.59% and 10.53%. The second run comes close to the chosen kinetic reserve. Lowering the reserve would change the stopping criterion. In the original peaked profile the same strong contrast had already reached the reserve on both branches by $K=54$. Changing the shape can extend the finite range where a radiation majority occurs. The five-dimensional model must still supply the interaction and its shape before the percentage becomes a prediction of that model.

The checks include integration in the field coordinate for the prescribed geometry, independent variables for the coupled quantum and gravitational trajectories, and separate propagation after the coupling stops changing. [Calculation](reach.py) · [Results](REACH.json) · [Calibrated equations](../matched/MATCHED.md)
