# Which source acts first when the interaction turns on?

The final energy in emitted particles gives limited information about the source acting on the shell as the interaction turns on. Here a short-time expansion yields a quantitative test of the source inserted into the three junction conditions. It applies at a fixed, finite pair-energy cutoff and with the parent and bath initially independent. All changes are coefficients of $y^2$ relative to the uncoupled evolution on the same trajectory.

## A short-time ratio independent of the spectral weights

Use proper time $s$, write $u=s-s_{\rm on}>0$, and suppose the coupling rises as $F=\alpha u^r+O(u^{r+1})$ for an integer $r\ge1$. The previous finite Hamiltonian uses parent inertia $a^3$, bath frequencies $\omega_j/a$ and the interaction $y F\chi\sum_jg_jx_j$. Assume a smooth background and a factorized conformal bath vacuum at turn-on. Let $X_{q0}=\langle\chi_q^2\rangle$ at turn-on and $a_0=a(s_{\rm on})$. Include the native momentum weights $w_q$ once. Define two positive coefficients

$$
U_P=\frac{\alpha^2}{2(r+1)^2}
       \sum_{qj} w_q g_j^2\frac{1}{2\omega_j a_0^3},
\qquad
U_D=\frac{\alpha^2}{2(r+1)^2}
       \sum_{qj} w_q g_j^2\frac{X_{q0}}{a_0}.
$$

In the bare finite Hamiltonian, the excess **free** parent and daughter energies per comoving volume start as $U_Pu^{2r+2}$ and $U_Du^{2r+2}$. The interaction contribution is negative and larger in magnitude than their sum. In particular,

$$
\lim_{u\downarrow0}\frac{\delta E_{\rm total}}{\delta E_P+\delta E_D}
=-\frac{r}{r+2},
\qquad
\lim_{u\downarrow0}\frac{E_{\rm interaction}}{\delta E_P+\delta E_D}
=-\frac{2(r+1)}{r+2}.
$$

These ratios hold mode by mode and after positive momentum weighting. For the sine-to-the-fourth gate the ratios are $-2/3$ and $-5/3$. The statement concerns changes of energy relative to a prepared quantum state in the bare model. Both free sectors have a positive leading energy change.

To derive the result, use one bath mode and divide its induced mixed moments by $y g_j$. Denote the resulting moments of $\Pi x_j$, $\chi p_j$ and $\chi x_j$ by $R_x,Q_p,Q_x$. The canonical equations give

$$
R_x=-\frac{\alpha}{(r+1)2\omega_j}u^{r+1}+O(u^{r+2}),
\qquad Q_p=-\frac{\alpha X_{q0}}{r+1}u^{r+1}+O(u^{r+2}).
$$

Since $Q_x'=a^{-3}R_x+Q_p/a$, it follows that

$$
Q_x=-\frac{\alpha}{(r+1)(r+2)}
\left(\frac{1}{2\omega_j a_0^3}+\frac{X_{q0}}{a_0}\right)
u^{r+2}+O(u^{r+3}).
$$

The leading free-energy exchange rates are $-F a^{-3}\sum_jg_j^2R_x$ and $-F\sum_jg_j^2Q_p/a$ for the parent and daughter respectively. Integrating gives $U_P$ and $U_D$ above. At this order the interaction energy is $F\sum_jg_j^2 Q_x$. A time-varying background modifies the displayed leading coefficients only through its values at turn-on. Terms involving its derivatives contribute at later powers of $u$. These expansions take the short-time limit with the cutoff held fixed. At a large cutoff the interval where all phases $\omega_j u/a_0$ are small shrinks. [Finite Hamiltonian and spectral normalization](https://doi.org/10.5281/zenodo.24358059)

## Source entering the shell equations

For the bare Hamiltonian the leading coefficient of $a^3\delta p$ is $U_P+U_D/3$ at order $u^{2r+2}$. The free parent initially acquires kinetic energy, the bath has the radiation pressure relation, and the bare interaction has no explicit $a$ derivative. A differentiable scalar clock $F(\phi)$ makes the pressure and current satisfy $(a^3\delta\rho)'=-3ha^3\delta p+a^3\delta j\,\phi'$. Thus the leading coefficient of $a^3\delta j\,\phi'$ is $-(2r+2)r(U_P+U_D)/(r+2)$ at order $u^{2r+1}$.

For the static calibration the added Hamiltonian term is $y^2 a F^2\sum_j g_j^2\chi^2/(2\omega_j^2)$. Define $L=\tfrac12 a_0\alpha^2\sum_{qj}w_q g_j^2 X_{q0}/\omega_j^2>0$. Its explicit energy dominates the response of the free parent at turn-on. The leading comoving coefficients become

$$
a^3\delta\rho=L u^{2r}+O(u^{2r+1}),\qquad
a^3\delta p=-L u^{2r}/3+O(u^{2r+1}),\qquad
a^3\delta j\,\phi'=2rL u^{2r-1}+O(u^{2r}).
$$

The leading energy and scalar power thereby reverse signs between the bare and static prescriptions. Their leading daughter energy agrees. The exact native boundary combination $n(B-A)=-(\rho+p)/2$ eliminates the tension. For $r=4$ the bare source entering its matter term has comoving enthalpy $(U_P-U_D)u^{10}/3$ to leading order. For the fine quiet history at cutoff 60 this is $-0.1537332$ times the leading sum of free-energy changes. The statically calibrated source instead has positive leading enthalpy $2L u^8/3$. Consequently the explicit driving term in this boundary combination has opposite signs in these two finite models. The absolute boundary equations include the pre-existing matter and its response to the perturbed geometry. This sign test by itself does not establish a reversal in the total enthalpy or late expansion rate.

## A cutoff where the bare boundary driving term changes sign

The leading bare enthalpy can also change sign within a single calibration as the spectral band changes. For $\Lambda\ge8$ all 24 parent momentum nodes are below cutoff. Integration of the polynomial Dirac measure gives the signs of $U_P-U_D$ from two functions, with their common positive factor $\alpha^2/[2(r+1)^2]$ removed,

$$
P(\Lambda)=\sum_q \frac{w_q}{8\pi^2a_0^3}
  (\Lambda^3/3-q^2\Lambda+2q^3/3),\quad
D(\Lambda)=\sum_q\frac{w_q X_{q0}}{16\pi^2a_0}(\Lambda^2-q^2)^2.
$$

For the fine quiet history, the coefficients of $P(8+x)-D(8+x)$ in *ascending* powers of $x$ are approximately $(2.90071,1.93142,0.49904,-0.02122,-0.002192)$. There is precisely one sign change. Descartes’ rule together with the positive value at $x=0$ and negative leading coefficient gives a unique zero for $\Lambda>8$. It is at $\Lambda=21.20006$ on the fine grid and $21.20120$ on the coarse grid. Both use activation at the same scalar value. Direct quadrature of the unexpanded spectral integrals at each root gives relative parent–daughter differences below $2\times10^{-15}$. The leading matter enthalpy driving $n(B-A)$ is positive below this transition and negative above it. Exactly at the transition its $u^{10}$ coefficient vanishes, leaving higher orders to determine its sign. The threshold refers to the chosen band, preparation and activation event. It is a diagnostic of finite-cutoff model dependence.

## Independent numerical check

For each parent momentum and one bath frequency, the code advances a $4\times4$ canonical fundamental matrix and its first two response coefficients. Pairing each matrix with the initial covariance reconstructs the two free energies, the mixed interaction, pressure and scalar power. The spectral and momentum sums are taken afterwards. This tests the primary covariance calculation using a different numerical representation. Rescaling the response matrices by the analytically known time powers avoids a subtraction of two large vacuum energies when $u$ is small.

| Proper time after activation | Bare total / free energy change | Static pressure / density change |
|---:|---:|---:|
| 0.010 | −0.6686590 | −0.3339434 |
| 0.001 | −0.6667745 | −0.3337363 |
| 0.00003, refined | −0.6666696 | −0.3333466 |
| Short-time limit | −2/3 | −1/3 |

This calculation uses all 24 native parent momenta on the fine quiet wall at cutoff 60. At the smallest delay each of eight computed coefficients is within 0.015% of its predicted limit. In an additional comparison the fundamental matrix follows the exact delays recorded in the source tables used by the bulk simulation. Across the eight stress, energy and scalar-power entries, the largest fractional discrepancy is 2.74% two stored cells after activation, 0.0185% at eight cells, and 0.000079% at 32 cells. The earliest bare energies grow only as the tenth power of the delay. A ratio of these tiny numbers is sensitive to time resolution. [Onset ratios](ONSET.json) · [Comparison to injected matter](INJECTED.json)

~~~sh
OPENBLAS_NUM_THREADS=1 python theory/onset.py /path/to/prior/native
OPENBLAS_NUM_THREADS=1 python theory/onset.py /path/to/prior/native --tables /path/to/source/tables
~~~

## Why a table evaluated on the unperturbed wall gives a first variation

Let $Z$ denote the entire native discrete geometry and its pre-existing matter variables, with lapse and clock included. For smooth histories write the finite perturbation of its evolution schematically as $Z'=G(Z)+\lambda S[Z]+O(\lambda^2)$, with $\lambda=y^2$ and causal history dependence in the new source. The first variation $z=\partial_\lambda Z|_0$ obeys $z'=DG[Z_0]z+S[Z_0]$ with zero initial variation. Thus the explicit new source can be calculated on the archived unperturbed path. The geometry and free parent must still respond through $DG[Z_0]z$, and each pulsed history requires its own source table. A signed numerical derivative at zero coupling can use negative values of $\lambda$ as a mathematical continuation of the driving term. Only $\lambda\ge0$ corresponds to a real Yukawa coupling.

For a final readout specified by proper time, differentiating the fixed-coordinate output also requires its clock variation. If $s(T_0,0)=s_f$, then $\partial_\lambda Z|_{s_f}=z(T_0)-Z'_0(T_0)(\partial_\lambda s)_0/s'_0(T_0)$. Interpolation of each simulation onto the same proper-time readout implements this correction. Evaluating the same smooth forcing at a perturbed clock changes $\lambda S$ at order $\lambda^2$.

A leading pair energy of the form $\lambda D[Z_0]$ first feels this order-$\lambda$ change in geometry through a term of order $\lambda^2=y^4$. The present first variation therefore tests the response of the geometry, junctions and native parent. An absolute order-$y^4$ prediction for fermion emission would also require the higher-order fermion dynamics. Neither the local short-time law nor the tangent itself identifies a unique astronomical origin.
