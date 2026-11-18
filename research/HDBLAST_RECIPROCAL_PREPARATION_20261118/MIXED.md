# What a mixed displacement and velocity preparation must cancel

The higher-dimensional packet has exactly zero initial scalar readout and a positive late mean in the normalized model units. A four-dimensional observer measuring the same quadrupole can use these facts to bound how the displacement and velocity of the comparison fields are prepared. This comparison uses a release at the registered clock origin and free propagation thereafter. The actual higher-dimensional packet has nonzero initial data away from the shell. [Earlier shell and modal construction](https://github.com/maldonado-research/HDblast/blob/4028ccb085adf56d7a74c737ab8dee11e85d5ea4/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/TOWER.md)

Let $c$ be the real constant detector vector, $b$ the initial displacement and $v$ the initial velocity in elapsed time $x=q/\ell$. Let $C_j,S_j$ be the physical free scalar solutions with respective initial pairs $(1,0)$ and $(0,1)$ for dimensionless squared mass $P_j\geq0$. The readout is

$$
F(x)=\sum_j c_j[b_j C_j(x)+v_j S_j(x)].
$$

This includes countably many fields if $c,b,v\in\ell^2$. Nonnegative instantaneous oscillator energy gives $|S_j(x)|\le x$ on this expanding clock. Energy divided by the instantaneous squared frequency decreases and gives $|C_j(x)|\le1$. The associated sums therefore converge uniformly. Write $L$ for the rational lower endpoint of the target late physical mean on $W=[5/4,11/8]$, with $L>4.01324558$. Let $\varepsilon$ bound each of the absolute errors on the selected measurements. [Window certificate](https://github.com/maldonado-research/HDblast/blob/4028ccb085adf56d7a74c737ab8dee11e85d5ea4/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SHAPE.json)

## Displacement and detector angle

If $v=0$, the later mean requires $M=\sum_j|c_j b_j|\geq L-\varepsilon$. At release $|\langle c,b\rangle|=|F(0)|\leq\varepsilon$. Hence for nonzero vectors and $L>\varepsilon$,

$$
|\cos\angle(c,b)|\leq \frac{\varepsilon}{L-\varepsilon},
\qquad
\min\left(\sum_{c_jb_j>0} c_jb_j,\; -\sum_{c_jb_j<0} c_jb_j\right)
\geq\frac{L-2\varepsilon}{2}.
$$

Indeed $M\leq\|c\|\|b\|$ by Cauchy–Schwarz, and the smaller of the two signed totals is $(M-|F(0)|)/2$. At $\varepsilon=1/2000$ the cosine magnitude is below $0.000125$ and each signed total exceeds $2.0061$. The rational cosine bound places the angle within $0.00714$ degree of perpendicular. At exact agreement $c$ is perpendicular to $b$ when the latter is nonzero, including for a mixed preparation. If the initial summed positive oscillator energy is finite, $F'(0)=\langle c,v\rangle$ by convergence in the velocity norm. Exact initial silence then also makes a nonzero $v$ perpendicular to $c$. Errors on amplitude alone provide no corresponding derivative tolerance.

## A bound including both parts of the preparation

Let $T$ be the normalized early triangular exposure $T[F]=\int_0^1 2(1-x)F(x)\,dx$ and $W$ the normalized late window mean. The per-mode velocity result gives $T[S_j]>0$ and $|W[S_j]|\le9T[S_j]$. Define $N_d$ as the smaller of the positive and negative totals of $c_jb_j$ in absolute value. Define $N_v$ from the signed values $c_jv_jT[S_j]$ by the same rule. Both quantities are expressed in the registered amplitude units. Their decomposition refers to the chosen independent canonical species. A change of basis in a degenerate mass sector can change its individual signed contributions, but the inequality holds for every such decomposition.

The displacement contribution is bounded in both exposures by $M=|F(0)|+2N_d\le\varepsilon+2N_d$. Since the early target vanishes, the absolute early exposure of the velocity contribution is at most $\varepsilon+M$. Its late absolute exposure is at most $9(\varepsilon+M+2N_v)$ by summing the per-mode bound. Matching the late mean thus requires

$$
L\leq20\varepsilon+20N_d+18N_v .
$$

In particular, if each part has a common sign across species, then $N_d=N_v=0$ even when the *two parts* have opposite signs. An arbitrary mixture of these releases needs at least $L/20>0.20066$ error in one of the initial value, the early triangular exposure and the late mean. This lower bound is one twentieth of the corresponding calibrated target amplitude. Absolute tolerances scale under a common gain rescaling. If errors on these three measurements are all at most $1/2000$, an unrestricted preparation must satisfy $20N_d+18N_v>4.0032$. This statement uses one point readout at the known release and two finite exposures.

The exact four-dimensional replica can meet this requirement with correlated, signed initial data. Its ambient dimension cannot be inferred from the sign decomposition. Applying the bound to an experiment requires physical justification for the allowed source preparation and the same angular readout, clock and calibration.

[Exact arithmetic checker](angle.py) · [Rational bounds](ANGLE.json)
