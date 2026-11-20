# Conservative release into the same three Maxwell fields

During the bounce calculation, an incident electromagnetic pulse sets a wall scalar in motion. Here the scalar subsequently returns energy to the same three gauge fields. Their stress has the exact isotropic equation of state $p_A=e_A/3$. With the initial pulse density 0.123 and the profile below, they account for over 90% of the wall's material energy while its radius grows from 3.4290 to 3.9083. The model then recollapses. The independent radiation fluid and damping coefficient both vanish throughout.

## What drives the release?

Use the timelike gravitational wall equations in the [loading solver](loading.py). On each of the two matched faces the incoming and outgoing characteristic amplitudes are $I$ and $O$. The gauge triplet has wall-frame electric amplitude $X=I+O$, magnetic amplitude $B=2A/a^2$, and surface density $e_A=K(X^2+B^2)/2$. The initial kinetic coefficient 0.001 already rises to 0.04 as the early pulse loads the scalar. The retained scalar climbs the same positive plateau of height 0.111 as in the preceding branch.

Now tie two changes in the action to the field value. From $\psi=0.36$ to 0.38 the old potential falls by $0.111S$ while the kinetic coefficient increases by $\Delta K S$. Both use the same smooth unit step $S((\psi-0.36)/0.02)$. The coupling has no explicit time dependence. There is a small tangential magnetic field on the expanding sphere at the start of the drop. Raising the gauge coefficient increases the energy carried by this field. Including the electric term, the proper power transferred from scalar to Maxwell per unit wall volume is
$$T=\tfrac12\dot K(B^2-X^2).$$
A magnetic field therefore gives positive transfer while $K$ increases. The scalar, gauge and bulk terms obey
$$\dot e_\psi+3H(e_\psi+p_\psi)=-T,\qquad
\dot e_A+4He_A=T+2(I^2-O^2).$$
The last term includes both faces. No additional friction or fluid source appears. Three equally excited Abelian fields give isotropic stress on the round three-sphere. This is a coherent electromagnetic state; a temperature or a thermal spectrum requires further physics.

## Main trajectory and energy budget

With $\Delta K=2000$, the wall bounces at $a=1.87793$ and reaches the potential descent at $a=2.72544$. Its scalar velocity there is $+0.04815$. The gauge seed has $B=0.0111097$, $X=-0.00016416$ and $e_A=2.46905\times10^{-6}$ in the wall-density units. The incident pulse has already ended.

| Event | Proper time | Radius | Fraction in wall Maxwell fields |
| --- | ---: | ---: | ---: |
| Start of potential descent | 3.05527 | 2.72544 | 0.00220% |
| End of potential descent | 3.30732 | 2.96608 | 87.0709% |
| First post-descent crossing of 90% | 3.92464 | 3.42903 | 90.0000% |
| Peak during expansion | 5.68256 | 3.90830 | 92.1052% |
| Expansion ends | 5.69580 | 3.90832 | 92.1051% |

Once the descent finishes, $K$ is constant. Scalar kinetic energy dilutes as $a^{-6}$, whereas the Maxwell density dilutes as $a^{-4}$ apart from bulk leakage. This accounts for the later increase from 87.07% to over 90%. During the expanding interval above 90%, $\ln(a_{\rm end}/a_{\rm start})=0.13083043$. The minimum bulk mass parameter on the full trajectory is approximately 0.0680. Continuing through the maximum gives $a=1.6$ at proper time 8.92716.

Integrating from the beginning to the end of the potential descent gives a comoving scalar-to-gauge transfer $\int a^3T\,d\tau=2.34527523$, outward leakage $\int 2a^3O^2 d\tau=0.0000006071$, and gauge pressure work $\int Ha^3e_A d\tau=0.05879399$. The gauge energy $a^3 e_A$ increases by 2.28648063. The scalar energy $a^3e_\psi$ decreases by 1.93111085 and its pressure work is $-0.41416438$. In these comoving units the common unit three-sphere volume is factored out. By the expansion maximum the total post-descent-start leakage has reached 0.00405302. The early pulse is supplied through the timelike wall characteristics; independently sourcing those characteristics from finite-distance bulk data is a separate gravitational initial-boundary-value problem.

## How much change in the coupling is necessary?

The gauge equation implies a useful exact identity,
$$\frac{d}{d\tau}\left(\frac{a^4e_A}{K}\right)=
\frac{a^4}{K}\left[2(I^2-O^2)-\dot K X^2\right].$$
After the external pulse $I=0$, and $a^4e_A/K$ is nonincreasing whenever $\dot K\ge0$. For any target fraction $f_*$ reached on such a segment from $i$ to $f$,
$$\frac{K_f}{K_i}\ \ge\ \frac{a_f^4 f_* e_f}{a_i^4 e_{A,i}}.$$
The main run satisfies these assumptions from the start of the potential descent to the first crossing of 90%. The required ratio for these endpoints is at least **49,778.57**; the realized ratio is **50,001**. Equivalently, the required seed is at least $2.45806\times10^{-6}$ and the actual seed is $2.46905\times10^{-6}$. The final coefficient $2000.04$ is **2,000,040 times its initial value**. A smaller seed, different cosmological trajectory or continued incident flux changes the bound. In particular, with residual illumination the maximum possible increase of $a^4e_A/K$ on a segment with increasing $K$ is $\int 2a^4I^2/K\,d\tau$.

## Controls and limits

| Incident dose | $\Delta K$ | Largest gauge fraction during expansion after the start of descent |
| ---: | ---: | ---: |
| 0.123 | 0 | 0.0026% |
| 0.123 | 1,500 | 78.4390% |
| 0.123 | 1,850 | 89.0672% |
| 0.123 | 1,900 | 90.2136% |
| 0.123 | 2,000 | 92.1052% |
| 0.123 | 2,150 | 93.4041% |
| 0.123 | 2,400 | 84.8366% |
| 0.123 | 8,000 | 1.1974% |
| 0.12275 | 2,000 | 80.7447% |
| 0.12325 | 2,000 | 84.7112% |

A larger change can interrupt the descent. On the ramp, at an instant of zero scalar velocity, the sign of the acceleration is the sign of $0.111-\Delta K(B^2-X^2)/2$. At $\Delta K=2400$ the field temporarily reverses and later finishes the descent. At 8000 it reverses and never reaches $\psi=0.38$ before the calculation reaches its inner radius. The raw characteristic dose also matters: with $\Delta K=2000$, changing that dose by about 0.2% in either direction removes the 90% expanding epoch in the tested runs. A universal conversion efficiency cannot be inferred from this scan.

High-order quadrature on each accepted evolution step checks the scalar and gauge comoving budgets separately, including leakage, as well as the identity for $a^4e_A/K$. Tightening the primary ODE tolerance from $2\times10^{-12}$ to $2\times10^{-13}$ shifts the duration by less than $8\times10^{-10}$ e-folds. This is a conditional solution for a specified effective action with a large change in its gauge kinetic coefficient. The origin of that function, thermalization, and a connection to observations remain open. [Release calculation](return_field.py) · [Numerical results](LOADING.json)
