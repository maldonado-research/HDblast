# Can the electromagnetic field give its energy to particles without collisions?

Yes, in this calculation of the final expansion. Charged daughters born at different times each accelerate in the electric field and contribute to its current. Before turnaround, positive electrical work transfers an amount equal to **152.4%** of the initial electromagnetic energy to the particles. Their current returns **61.6%** to the field during other parts of the motion. The net transfer is **90.85%** of the initial field energy. These percentages concern the same energy moving back and forth.

At turnaround the daughter particles carry **83.30%** of the model wall's added energy and electromagnetic fields carry **8.39%**. With neutral daughters the corresponding percentages are **0.0181%** and **91.67%**. The charged particles remain strongly directional. The result distinguishes energy gained from the field, energy spread among particles, and thermodynamic equilibrium.

## Inputs and scope

The test imports the 64-element, $\kappa=100$ primary trajectory from the continuous decay study. We start at $\tau=10.34$, immediately before the scalar-controlled decay threshold is crossed at $10.34110$. The archived initial electromagnetic potential, displacement and gauge kinetic coefficient are $A=a^2B/2=0.00137388$, $P=a^2X=-0.0810095$ and $K=13{,}200.04$. We follow the archived radius $a$, neutral mass ratio $R=M/(2m)$, converted fraction $f$, and scalar energy through its expansion maximum at $12.42354$. At that time 99.675% of the parents have converted. The scalar density included in the quoted wall fractions is the imported potential plus scalar kinetic energy. The background tension is excluded.

**The field and charged particles evolve together**. The geometry, scalar and decay history are prescribed. We turn off exchange with the surrounding bulk in this test. In the prior coupled-fluid trajectory the recorded bulk transfer over the same interval was $-0.00438$ in comoving units, or $-0.0101\%$ of the initial field energy. Its scalar-to-electromagnetic transfer was zero. A bulk responding to the collisionless currents could give a different flux. The test cannot establish how the modified stress changes the bounce or turnaround.

There are three equal families with orthogonal electromagnetic axes. Each neutral decay produces charges of opposite signs with opposite birth momenta of magnitude $p_b/m=\sqrt{R_b^2-1}$ and a uniform birth direction in the parent rest frame. Families equally share the stipulated total daughter rest density. No second factor for partners or families is applied. Set $c=1$, measure momenta in units of the daughter mass and write $D=8D_i/R_i$, with $D_i=10^{-4}$, $R_i=0.95$ and initial radius 2. The positive member of a family has $Q=10^8$ in the projected field normalization. Its negative counterpart has the same energy and contributes the same sign to the electrical current.

## Evolution and angular integrals

For a packet born at $a_b,A_b$ and angle cosine $z$, its subsequent longitudinal and transverse momenta are

$$
p_\parallel=x+sz,\qquad p_\perp^2=s^2(1-z^2),\qquad
x=Q(A-A_b)/a,\qquad s=a_b\sqrt{R_b^2-1}/a.
$$

The curved-wall magnetic field rotates the transverse momentum about the same axis. This leaves these observables fixed between electric impulses apart from the expansion factors. The homogeneous round-three-sphere connection also leaves this axisymmetric momentum evolution intact. The signed velocity is $v_\parallel=(x+sz)/\gamma$, where $\gamma=\sqrt{1+x^2+s^2+2xsz}$. With $u_\pm=\sqrt{1+(x\pm s)^2}$ and $U=u_++u_-$, the exact angular means are

$$
G=\langle\gamma\rangle=\frac23\left(U-\frac{u_+u_-}{U}\right),
\qquad\langle v_\parallel\rangle=\partial_xG,\qquad
\langle1/\gamma\rangle=2/U .
$$

These formulas stay finite at $x=0$ and $s=0$. Integrating over the daughter births with measure $df$ gives $J=QD a^{-3}\int \langle v_\parallel\rangle\,df$. In wall time the closed Maxwell equations become

$$
\dot A=P/a,\qquad
\dot P=-4A/a-a^2J/K .
$$

The comoving electromagnetic energy is $W=K(P^2+4A^2)/2$. It obeys $\dot W=-a^4JX$. As an independent check on the particles and the expansion, define $U_{\rm m}=Da[(1-f)R+\int G\,df]$. Direct differentiation gives

$$
d(W+U_{\rm m})=
 D\left[(1-f)R+\int\langle1/\gamma\rangle\,df\right]da+
 Da(1-f)\,dR .
$$

At birth each pair draws exactly its assigned parent energy. This balance includes work associated with expansion and variation of the surviving parent mass. Field work, daughter energies and the two terms on the right are integrated independently.

## Streams and resolution

Near the end the instantaneous electromagnetic share changes substantially. Over the final quarter unit of wall time its sampled range is 0.093% to 8.39%, with **3.26%** as the time average. In the representative family, 45.45% of the positive carriers move in the positive axial direction at turnaround. The mean signed velocity is $-0.091c$, while the mean absolute longitudinal velocity is $0.99996c$. The 10th, 50th and 90th percentiles of particle Lorentz factor are about 848, 3,915 and 11,148. Approximately 59% of the daughter energy exceeds the cold-stream energy at the same mean momentum.

An independent quadrature constructs the distribution's energy–momentum tensor. In the energy rest frame of this one charged species, over 99.99999% of its summed pressure is parallel to its field axis. An equilibrium gas in its rest frame has one third along each axis. Summing three rotated families yields isotropic **total** pressure while retaining the directional structure of each family. Gas temperature, collisional entropy production and a photon spectrum are outside this collisionless example.

The solver uses equal-width intervals in birth time, weighted by the actual increments in $f$, and a centered drift–kick–drift field evolution in conformal time. The baseline has 4,096 intervals and 131,072 field steps. Doubling the number of birth intervals at fixed time step changes the final electromagnetic share by $2.70\times10^{-5}$ and its final-quarter mean by $8.96\times10^{-6}$ in absolute share. Halving the field step instead changes the final share by $9.82\times10^{-7}$. The peak combined energy-balance error is $6.1\times10^{-6}$ of initial field energy and falls to $1.4\times10^{-6}$ at the smaller step. Linear interpolation of the archived mass and conversion data produces a final electromagnetic share of 8.385%, close to the shape-preserving interpolation. The charge-free solution matches the analytic conformal oscillator with absolute field-variable error below $1.3\times10^{-11}$. Direct 8- and 16-node angular integration checks the analytic current and spectrum.

Finite birth bins approximate a continuous distribution. The percentile and instantaneous current values are more sensitive to the relative phases of streams than the time-averaged field energy. A reversible collisionless evolution can dephase and later rephase. This experiment follows only the stated late window and gives no long-time irreversible damping rate.

Run `python -I -B kinetic/vlasov.py --source coupled/DECAY.json --output /tmp/KINETIC.json` from the capsule root. Python with NumPy and SciPy, and a C++17 compiler named `g++`, are required. [Runner](vlasov.py) · [Current and field evolution](vlasov.cpp) · [Controls and spectra](KINETIC.json)
