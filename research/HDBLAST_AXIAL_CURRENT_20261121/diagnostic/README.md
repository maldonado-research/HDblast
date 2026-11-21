# How soon can a different particle current reverse the axial electric field?

The negative electric field initially accelerates opposite charges in opposite directions. At the first work measurement the fluid already moves close to light speed. Even a redistribution that sends every positive and negative charge at the speed limit in the most favorable directions can produce only a small extra current that weakens the negative field. A reduction of that current can instead strengthen the negative field. These limits use the recorded radius, incoming wave and birth history, and hold for any homogeneous current satisfying the speed limit.

| History | Upper limit on field weakening at time 6.30944 | Upper limit at time 6.45 | First field zero with the original fluid current | Earliest field zero in the extremal current calculation |
|---|---:|---:|---:|---:|
| Separate interactions | 0.080% | 6.07% | 6.47163 | 6.46984 |
| Shared interaction | 0.027% | 2.14% | 6.47044 | 6.46984 |

The percentages are relative to the magnitude of the original field at that time. Both original fields and both extremal upper fields are negative at 6.45. In the separate-interaction history the sign of the reference field changes shortly afterward. A relative percentage near a reference zero is poorly conditioned, so the archived table also gives the absolute field changes. A conservative interpolation bound gives a negative axial field through time 6.46982 in both cases. These statements concern the measured inputs with the wall's incident characteristic fixed. Photons or changes in geometry and production could alter those inputs. The [finite-collar calculation](../return/README.md) examines the effect of changing the incident characteristic when the wall sends waves back into the bulk. [All bounds and both numerical integrators](CURRENT.json)

## Derivation

Let $F=a^2 E$ and $G=a^2 B$ be the conformal axial electric and magnetic fields, with wall proper time $t$. Compare two histories beginning with the same fields at first birth. Write $v_r=u/\sqrt{1+u^2}$ for the reference signed electric-current velocity and $w$ for the new charge-weighted axial velocity, with $|w|\leq1$. Each charge sign is included. The difference satisfies

$$
\frac{d}{dt}\binom{\delta F}{\delta G}
=
\begin{pmatrix}
-2/K &-2/a\\
2/a &0
\end{pmatrix}
\binom{\delta F}{\delta G}
+
\binom{Q D(v_r-w)/(K a)}{0}.
$$

Here $K=2.54$, $Q=10^7$, and $D$ is the cumulative comoving rest in the normalized triad. The common radial incident characteristic cancels. For a unit electric impulse at $s$, write $(f,g)$ for the homogeneous solution with $(f(s),g(s))=(1,0)$. The derivative of $f^2+g^2$ is $-4 f^2/K$, so $|f|\leq1$ and $|g(t)|\leq 2(t-s)/a_{\min}$. Variation of constants gives

$$
 f(t)\geq e^{-2(t-s)/K}-2[(t-s)/a_{\min}]^2 > 0
$$

through time 6.67 on both histories. Consequently $g\ge0$ and $f(t)\le e^{-2(t-s)/K}$. Positivity makes $w=-1$ and $+1$ the extremizers for the axial electric component at every time in this interval. The equality cases trace limiting light-speed currents. Massive particles may approach these limits only subject to their dynamics and finite energies.

For a conservative upper integral on each measured time cell, the script maximizes $D$, $1/a$ and $1+v_r$ separately using the endpoint values of the interpolated data. For the bound in the opposite direction it uses $1-v_r$. It propagates the exponentially damped envelope analytically and checks whole cells for the earliest possible positive field. A midpoint calculation propagates both signs of the extremal current, and an independently implemented adaptive differential-equation solver checks both field components. These are bounds on filtered current responses for specified input histories. They do not enclose numerical errors of the underlying five-dimensional evolution. [Code](current_bounds.py)

## Which currents account for the late difference?

The field solver allows this question to be answered directly, using the recorded current in each time step. We propagated the final electric measurement backward through the exact linear field updates, then multiplied its response to a current impulse by the actual difference between the particle and fluid currents in each step. The signed contributions add to the measured field difference. This is a retrospective energy-independent check of the field update and identifies when the difference arose. It does not predict what the current would be for a different field. [Adjoint calculation](memory.py) · [All four paired numerical refinements](MEMORY.json)

In the finest separate-interaction calculation with collisions, at time 6.57491 the dynamic axial field was 0.06580 lower than the field driven by the original fluid current. At 6.67 the sign of this difference had changed. The contributions to the 6.67 field from current differences before 6.57481 summed to **−0.05803**. The later currents added **+0.08879**, giving **+0.03076** overall. Both histories receive the same incoming radial wave. The later current changes are large enough to reverse the accumulated difference. Each sum uses physical electric field units at the time of measurement.

~~~sh
python -I -B diagnostic/current_bounds.py --inputs /path/to/field-tables --output /tmp/CURRENT.json
python -I -B diagnostic/memory.py --inputs /path/to/field-tables \
    --snapshots /path/to/live/snapshots --output /tmp/MEMORY.json
~~~
