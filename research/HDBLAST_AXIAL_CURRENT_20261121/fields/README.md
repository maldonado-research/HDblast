# Which field and current drive the shell?

The spherically averaged solver follows the outgoing and incoming electromagnetic waves between the outer source and the wall. We recorded the inward characteristic at the wall itself and the charged fluid's current, then used those as independent inputs to a second solution of the wall Maxwell equations. On both benchmark histories it recovers the axial electric and magnetic fields through two electric-field reversals and one current reversal. [Recorded checks and energy balances](FIELDS.json)

## Where the energy comes from during the first reversal

The following quantities are conformal energies per project-normalized volume from the first particle birth to wall time 6.6805. Incoming and outgoing energies refer to passage through the inner wall of the radial collar. Emission by the outer source has a separate normalization and acts at a different time. The values shown are for the parent fluid and its axial field. Particle-current experiments use this fluid case as a reference.

| Quantity | Separate interactions | Shared interaction |
|---|---:|---:|
| Stored axial electromagnetic energy at first birth | 613.39 | 613.39 |
| Additional axial energy arriving through the wall after first birth | 0.0000747 | 0.0000773 |
| Energy leaving through the wall | 137.87 | 142.59 |
| Gross positive work by the field on the fluid | 760.91 | 821.22 |
| Gross energy returned by the fluid to the field | 306.84 | 403.96 |
| Stored axial electromagnetic energy at the endpoint | 21.46 | 53.54 |

On the separate-interaction history the axial electric field crosses zero at time 6.471635. Until the current reverses at 6.574812 the particles do **285.26** units of work on the field. The stored field energy grows by **266.30** and **18.95** escape back into the radial collar. A second electric-field zero occurs at 6.660851. The shared-interaction history gives the corresponding times 6.470436, 6.572655 and 6.653651. This separates energy stored from the earlier illumination, energy returned by particles, and new external illumination during the reversal. A calculation that replaces the fluid by individual particles must recompute the current, including its sign. The table gives dimensionless model energies. It is not an estimate of the universe's observed energy density.

## Equations and normalization

Let $I(t)$ be the incoming axial characteristic, $X$ the axial electric field, $B$ the axial magnetic field, $a$ the wall radius and $K$ the gauge kinetic coefficient. The reflected characteristic is $O=X-I$. Between first birth and wall time 6.6805 both numerical histories have $K=2.54$ at every accepted step. Set $F=a^2X$, $G=a^2B$ and $d\eta=dt/a$. The homogeneous angular mode has curl magnitude 2, with signs fixed by the original wall solver. Its equation is

$$
\frac{dF}{d\eta}=-2G-\frac{2aF}{K}+\frac{4a^3 I-J}{K},
\qquad \frac{dG}{d\eta}=2F.
$$

The physical current of the archived fluid is $j=QD\,u/[a^3\sqrt{1+u^2}]$, where $Q=10^7$, $D$ counts total accepted comoving daughter rest and $u$ is the fluid's axial spatial four-velocity. The recorded source term is $J=a^3j$. The particle solver on one normalized representative beam uses $J=QD\langle s v\rangle$, where the two charges have signs $s=\pm 1$. This is the normalization of the symmetric reduced wall equations. In the kinetic cross section the target number density for a particular opposite charge and flavor includes the factor $1/6$ from the three flavors and two signs. The angular mode energy, particle work, and wall source are reported in the parent study's common normalization. If identical rotated copies are assigned equal shares of the total energy, each copy receives one third of the aggregate energy entries. The first-hit probability for that flavor still uses its own target density. The accompanying three-flavor particle calculation extends this normalization to unequal instantaneous flavor populations by using each flavor’s own tagged charges as the source of its own field. [Independent populations](../independent/README.md)

The Maxwell energy identity in this interval is

$$
U(t)-U(t_b)=2\int_{t_b}^t a^4(I^2-O^2)\,dt'
                  -\int_{t_b}^t a X J\,dt',
\qquad U=\tfrac12 K(F^2+G^2).
$$

The final term is energy received by the particles when it is positive under the integral. A negative integrand returns energy to the field. Before first birth, changes in $K$ can also supply energy. The postbirth calculation starts with the original stored field and integrates each sign of work directly. For the two additional polarizations the same equations apply to their own incoming characteristics and field strengths. Their reference radial solutions have zero added transverse particle current. The earlier and delayed profiles are normalized separately to the emitted energy of the original reference pulse. Their actual coefficients are applied by the kinetic solver.

## Numerical tests

The program uses the authenticated two-polarization study's C++ builder, preserving its solver and equations, and adds a recorder at each accepted step of the postbirth solution. It retains the existing radius, source energy, fluid moments, and angular birth law as well as all three inward characteristics. It also records the fluid current directly and checks its reconstruction from $D,a,u$. The first 5% of the electrical-work history agrees with the compact reference samples. The original early and delayed mode normalization factors and field values at all three earlier markers pass an independent comparison.

The wall replay uses an implicit midpoint step in proper time for $F,G$, with the tabulated inputs evaluated at midpoints. The maximum discrepancy in axial $F$ divided by the peak axial field norm is $4.5\times10^{-9}$ for the separate-interaction history and $6.9\times10^{-9}$ for the shared one. The narrow late transverse pulse gives the larger relative peak discrepancy, below $2.7\times10^{-6}$. Integrating the incident power, reflected power and particle work directly over the native time knots balances the axial field energy to better than $3\times10^{-9}$ of the peak in each history. The independently evolved radial energy ledger also agrees. These figures check the specified numerical trajectories and their field convention. The returning radiation from a new particle current can change the incoming data in a dynamical bulk. [Radial driver and replay](extract.py)

~~~sh
python -I -B fields/extract.py --study /path/to/HDBLAST_POLARIZATION_20261121 \
    --output-dir /tmp/axial-inputs --output /tmp/FIELDS.json --jobs 2
~~~

The command writes transient history arrays for each interaction case. The arrays include the physical axial fields and inward characteristic `X, B, incoming`, conformal current `J`, normalized transverse components `E1, B1, I1, E2, B2, I2`, and cumulative fraction of the reference electrical work. The time interval starts with the last recorded zero daughter population immediately before birth. The trajectory, temperature and angle-dependent birth data have the names used in the parent study. [Parent fields and source normalization](https://github.com/maldonado-research/HDblast/tree/173a1c13f9b4722b91fd2346d837d129c3836d4f/research/HDBLAST_POLARIZATION_20261121)
