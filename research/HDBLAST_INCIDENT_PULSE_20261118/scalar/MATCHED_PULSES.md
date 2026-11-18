# Compact scalar pulses with matching gravitational data on a B1 initial slice

A compact scalar pulse can satisfy both radial Einstein constraints and join the same quiet geometry at both ends of this slice. Three compensating scalar profiles suffice for any prescribed, sufficiently small smooth seed. All three controls can be confined to a chosen open interval farther from the wall. The result assumes a regular analytic quiet reference obeying the equations below and wall data in the certified box. It applies to the initial-value constraints. Its small-amplitude radius and the energy later transferred to particles require separate calculations.

Two explicit four-lobe examples give different numerical completions for the same outermost scalar lobe. A separate fixed-grid integration of the original absolute fields recovers all three matching quantities to within $3.5\times10^{-14}$ in both examples. The exact rank argument and the finite numerical examples use different checks. [Numerically matched coefficients](PROFILES.json) · [Interval certificate](RANK.json) · [Fixed-grid results](FIXED_GRID.json)

## Initial slice and endpoint map

Use the identity conformal chart of the B1 model on $Z\in[L,0]$, with $L=-0.20$ and the shell at zero. Set $A=B$ on the initial slice, $B_T=0$, $Q=A_Z$, and $u=A_T$. The quiet reference has $u_0=1$, scalar $\phi_0$, and $r=\phi_0'$. The inputs use `dc=0`, the zero-matter initial shell and the original B1 potential and tension. The model's numerical wall values are $\rho_b=7.838285538073243$ and $\phi_b=0.0006170897672752272$. They give $B(0)=\log\rho_b$, $Q(0)=\rho_b\sigma/6$, and $r(0)=-\rho_b\sigma_\phi/2$. [Source summary](https://github.com/maldonado-research/HDblast/blob/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/runs/main/main_ps0.5_G100_y1_bcut_dzf5e-4_summary.json) · [B1 constraint implementation](https://github.com/maldonado-research/HDblast/blob/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py)

Prescribe a real smooth compact scalar profile $F(Z)$ and set $\phi=\phi_0+F$, $p=\phi_T=-F'$. This selects a shell-directed perturbation of the scalar principal part at the initial time. Write $U=W_\phi^2/2-2W^2/3$, $W=1-\phi+\phi^3/3$. Integrating the two bulk constraints inward from the matched wall gives

$$
 B'=Q,\qquad Q'=u^2-Q^2-\frac{p^2+(r+F')^2}{6}-\frac{e^{2B}U(\phi_0+F)}3,
 \qquad u'=-\frac{p(r+F')}{3}.
$$

The quiet data solve these equations with $F=0$, along with $\phi_0'=r$ and $r'=e^{2B_0}U_\phi(\phi_0)-3Q_0r$. Define $\mathcal G(F)=(B-B_0,Q-Q_0,u-1)|_L$. If $\mathcal G(F)=0$, uniqueness for the radial ODE reproduces the quiet solution on an open collar next to $L$ where $F$ vanishes. The same uniqueness applies at the quiet wall. Both the scalar and the metric initial data then join the quiet background across these collars. An existing regular continuation of the quiet exterior can be used wherever it is defined. Exact wall and outer matching refer to the analytic construction. Rounded numerical coefficients have the finite mismatches reported below.

## Why three localized scalar controls suffice

Let $x=(\delta B,\delta Q,\delta u)$ denote variation with respect to $F$ at the quiet solution, and put $E=e^{2B_0}$. Values of $U$ and its derivatives here are taken at $\phi_0$. The linearized equations are $x'=Mx+H F+D F'$, where

$$
 M=\begin{pmatrix}0&1&0\\-2EU/3&-2Q_0&2\\0&0&0\end{pmatrix},\qquad
 H=(0,-EU_\phi/3,0)^\top,\qquad D=(0,-r,r)^\top/3.
$$

Subtracting $D F$ from $x$ gives a single forcing vector $K=MD+H-D'$. The quiet scalar equation reduces this to

$$
 K=\tfrac13(-r,(2-Q_0)r,-r')^\top.
$$

For the homogeneous transition matrix $\Phi$, the derivative of the outer mismatch is therefore $\mathcal G'(0)F=\int_0^L \Phi(L,z)K(z)F(z)\,dz$, since $F$ vanishes at the two endpoints. Write $N=\partial_z-M$. Differentiating the vector $h(z)=\Phi(L,z)K(z)$ gives $h'=\Phi NK$ and $h''=\Phi N^2K$.

The determinant of the three columns $[K,NK,N^2K]$ is bounded between **−1,589,246 and −1,532,659** throughout the following rational box for the wall data.

| Quantity | Lower | Upper |
|---|---:|---:|
| $E$ | 61.43 | 61.45 |
| $Q_0$ | 2.740 | 2.744 |
| $\phi_0$ | 0.00060 | 0.00064 |
| $r$ | 7.600 | 7.610 |

The certificate uses exact rational interval arithmetic on the quiet constraint and scalar equations. The archived decimal radius is the input for the numerical examples. For the exact local static wall at the same specified $\phi_b$ and tension parameters, the static first integral instead sets $\rho_b=(\sigma^2/36-\sigma_\phi^2/48+U/6)^{-1/2}\simeq7.83828553807424$. Direct rational checks on $E$, $Q_0^2$ and $r^2$ place this normalized wall in the same certified box. At the archived decimal wall inputs the determinant rounds to −1,560,778.14. Since $\Phi$ is invertible, $h,h',h''$ span three dimensions at the wall. Suppose that the span of $h$ over a nonempty open interval $I\subset (L,0)$ were a proper subspace. A fixed nonzero covector would annihilate $h$ on $I$. The quiet solution and coefficients are real analytic on the connected regular interval, so that scalar analytic function would vanish identically, including its first two derivatives at the wall. This contradicts the wall rank calculation. [Exact certificate and implementation](rank_certificate.py)

Choose three points of $I$ where the vectors $h$ are independent. Smooth normalized bumps $g_1,g_2,g_3$ supported in disjoint neighborhoods of these points give three independent integrated responses when their supports are sufficiently narrow. For any fixed smooth compact seed $f$ outside $I$, apply the finite-dimensional implicit function theorem to

$$
 (a,b_1,b_2,b_3)\ \longmapsto\ \mathcal G\!\left(a f+\sum_{j=1}^3 b_jg_j\right).
$$

This map is smooth near zero by ordinary differential equation dependence on finite parameters. At $a=0,b=0$ its value vanishes, and the derivative in $b$ is invertible. Hence the coefficients $b_j(a)$ exist smoothly for all sufficiently small $|a|$ and yield exact quiet matching. The argument allows different choices of controls for different intervals. It supplies no uniform lower bound on the usable amplitude as the interval shrinks, and no numerical error bound for a specified four-bump solution.

## Finite scalar-only examples

For $z_i=(-.17,-.125,-.08,-.035)$ and $w=.020$, set $F=\sum_i c_i f_i$ with $f_i=\exp(1-1/(1-x_i^2))$ for $|x_i|<1$, $x_i=(Z-z_i)/w$, and zero outside. These bumps have disjoint supports. The quiet collars are $[-.20,-.19]$ and $[-.015,0]$. Fix the outermost lobe coefficient at $c_1=-.0001$. Four continuation steps starting at the quiet state give the first solution; a second root is recovered from a separate initial guess. The recorded coefficients retain more digits than the table.

| Coefficient in increasing $Z$ | Continued from zero | Second completion |
|---|---:|---:|
| $c_1$ | −0.000100000000 | −0.000100000000 |
| $c_2$ | +0.000101859057 | −0.000236958015 |
| $c_3$ | −0.0000274874301 | +0.000298663941 |
| $c_4$ | +0.00000345303264 | −0.0000606359539 |

At the quiet state the three singular values of the $3\times4$ endpoint derivative are 0.708982, 0.187110 and 0.00472930 with the stated coordinates and unscaled amplitudes. The determinant for the last three controls is $3.17815\times10^{-4}$. A centered difference with amplitude $2\times10^{-7}$ agrees entry by entry with the variational calculation to $6.5\times10^{-9}$. Solving for three controls and repeating the calculation with a smaller radial step leaves endpoint discrepancies below $5\times10^{-17}$ in the difference formulation. Removing the wall-nearest lobe from the first solution leaves $|\Delta Q|\simeq1.49\times10^{-6}$. The isolated outermost lobe gives $|\Delta u|\simeq1.54\times10^{-5}$. [Shooting and controls](SHOOTING.json) · [Solver](shoot.py)

Independent classical RK4 evolves quiet and perturbed **absolute fields**. On 2,048 radial intervals the largest endpoint discrepancy is $2.1\times10^{-14}$ and $3.5\times10^{-14}$ for the two solutions. Centered differences on the resulting profiles directly evaluate the original Hamiltonian and momentum expressions. From 256 to 2,048 intervals, the maximum scaled momentum residual drops from $4.21\times10^{-5}$ to $2.88\times10^{-8}$ and from $1.04\times10^{-4}$ to $1.24\times10^{-7}$ respectively. These numerical tests do not enclose the error of time evolution. [Independent fixed-grid solver](fixed_grid.py)

The integrated momentum condition imposes $\int_L^0 F'r\,dz=-\int_L^0(F')^2\,dz$. The two positive squared-gradient terms are $3.19939\times10^{-6}$ and $2.40643\times10^{-5}$, with opposite background cross terms. These are unweighted coordinate integrals on one initial slice. Measuring physical energy reaching the shell requires the local stress current, the stated observer and an evolved solution.

## Reproduction

Use the unchanged B1 summary linked above, with SHA-256 `bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f`. The two independent checks use Python's standard library. The shooting and quadrature use NumPy and SciPy. The fixed coefficients in `PROFILES.json` are also usable in a separate evolution calculation. Re-solving the root on different machines can change its final floating-point digits.

~~~sh
python -I -B rank_certificate.py B1_summary.json
python -I -B fixed_grid.py B1_summary.json
python -I -B shoot.py B1_summary.json
~~~

The optional SciPy shooting program may report a subdivision warning for its independent momentum quadrature at the requested tolerance. The fixed-grid reconstruction and the analytic momentum identity provide separate checks of the reported matching.
