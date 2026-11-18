# A compact B1 scalar pulse with both quiet ends

The matter-free B1 wall admits a finite-amplitude incoming scalar perturbation whose initial geometry exactly rejoins the same quiet reference beyond its support. An interval calculation and the intermediate value theorem establish a nonzero amplitude between 0.00201 and 0.00203. Strict concavity makes it the only nonzero match in the amplitude range from zero through 0.003. A coordinate-invariant scalar measurement distinguishes every solution in that bracket from the corresponding quiet branch. The construction concerns initial data on one spatial slice. It leaves the wall's later energy gain to a time evolution.

## Initial slice and constraints

Take one smooth B1 interior and the identity conformal chart, with the shell at $Z=0$. The bulk metric is $e^{2B}(-dT^2+dZ^2)+e^{2A}d\mathbf{x}^2$. Keep the B1 potential, tuned tension and particle laws. Set the background seed displacement to zero, $dc=0$, with no particles or radiation initially. The quiet solution is $A_0=B_0(Z)+T$, $B=B_0(Z)$, $\phi=\phi_0(Z)$. Fix the planar coordinate normalization to this reference. Write $Q=B'_0$, $r_0=\phi'_0$, and take the initial slice $T=0$. Primes denote $\partial_Z$.

For a continuum reference using the recorded wall scalar and decimal tension parameters, impose the wall junctions $Q_b=\rho_b\sigma_b/6$, $r_b=-\rho_b\sigma_{\phi,b}/2$. The static first integral determines

$$
\rho_b=(\sigma_b^2/36-\sigma_{\phi,b}^2/48+U_b/6)^{-1/2}.
$$

Its numerical value is 7.83828553807424. The archived radius is 7.83828553807324. This rounding reconciliation lets the initial static constraint and wall junctions hold simultaneously. The certificate uses the specified decimal model parameters and this algebraic radius. The two smooth copies can be initialized by reflection as in B1. [B1 solver and model](https://github.com/maldonado-research/HDblast/tree/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source)

Choose $A=B_0$ and $A_T=1$ on the slice. Allow $B$ and $B_T$ to respond. Set $x=(Z+0.12)/0.02$, with $f=\exp[1-1/(1-x^2)]$ for $|x|<1$ and $f=0$ elsewhere. Prescribe $\phi=\phi_0+\eta f$, $p=\phi_T=-\eta f'$, $r=\phi_Z=r_0+\eta f'$. This selects the scalar characteristic moving toward increasing $Z$ at the initial instant. The gravitational solution and the scalar potential enter its subsequent propagation. The earlier initial-data trial fixed $A=B$ and $B_T=0$ on this slice. Here two different metric fields are held fixed.

Let $R=e^{B_0}$ and use the geometric planar mass $\mu=R^4e^{-2B}(1-Q^2)$. The quiet value is $\mu_0=R^2(1-Q^2)$. In the B1 normalization the radial Einstein equation gives

$$
\mu'=a\mu+s,\qquad
a=\frac{Q(p^2+r^2)-2pr}{3(1-Q^2)},\qquad
s=\frac23 R^4Q U(\phi).
$$

In particular the mass equation is linear once the scalar data, $A$ and $A_T$ have been prescribed. With $v=\eta f'$, $a-a_0=-2v(r_0+v)/[3(Q-1)]$. Integrate the difference $D=\mu-\mu_0$ from zero on the shell-facing collar,

$$
D'=a D + (a-a_0)\mu_0+\frac23 R^4Q\,[U(\phi_0+\eta f)-U(\phi_0)] .
$$

For $\mu/\mu_0>0$, reconstruction gives $B=B_0-\tfrac12\log(1+D/\mu_0)$. The momentum constraint then gives $B_T=(Q-B'+pr/3)/Q$. Substituting this expression into the original Hamiltonian imposes precisely the mass equation. One may also solve that Hamiltonian directly for $B'-Q$,

$$
B'-Q=
\frac{2[e^{2B}U(\phi)-e^{2B_0}U(\phi_0)]+p^2+r^2-r_0^2-2pr/Q}
{6(Q-1/Q)}.
$$

The full potential and the momentum term enter these formulas. The signed *change* in quasi-local mass relative to a background with a scalar gradient can cancel across the pulse. Vanishing mass difference at the two ends therefore permits a nontrivial scalar disturbance inside. [Planar mass normalization and Einstein identity](https://github.com/maldonado-research/HDblast/blob/d6c957f9e26816dcdb969a3246b714dfa8bf9700/research/HDBLAST_SCALAR_RESERVOIR_20261118/mass/MASS_BALANCE.md)

## An exact existence result for the specified data

Let $F(\eta)=D(-0.14)$ after integration from the wall toward decreasing $Z$. On the whole amplitude bracket $[0.001,0.003]$ the interval solver bounds $Q>2.7418$ and $\mu<-191.74$ throughout the construction from the wall to the inner edge of the support. Thus the mass equation and metric reconstruction stay regular, and $F$ is continuous there. The endpoint enclosures, rounded outward, are

| Amplitude | Certified interval for $F$ |
|---:|---:|
| 0.001 | [0.01252, 0.01263] |
| 0.00201 | [0.0001660, 0.0003764] |
| 0.00203 | [−0.0003324, −0.0001199] |
| 0.003 | [−0.03634, −0.03601] |

The narrower sign change places a root strictly between 0.00201 and 0.00203. For uniqueness, differentiate the linear ODE in amplitude twice. If $S$ denotes its forcing, the two variational equations in the same orientation are $D_1'=aD_1+a_1D+S_1$ and $D_2'=aD_2+2a_1D_1+a_2D+S_2$. Both variations vanish on the wall-facing collar. Interval integration bounds $-24848<F''(\eta)<-24416$ throughout $0\le\eta\le0.003$. Since $F(0)=0$, $F(\eta)/\eta$ is strictly decreasing for positive amplitudes in this range. The sign change therefore isolates exactly one nonzero root there. Other positive and negative amplitudes and changes of pulse shape are outside this uniqueness interval.

Within the pulse-free region $D'=a_0D$. If $D=0$ at the edge, it stays zero wherever the quiet solution continues. The smooth bump is flat at the edges. The resulting initial fields agree on both open collars, including all derivatives, so they can be extended by the original quiet data wherever that reference is defined. The exactly quiet solution preserves the shell junctions on its collar. The existence and uniqueness interval applies to this pulse shape and metric choice. [Interval method and replay](interval/README.md) · [Certified bounds](interval/BRACKET.txt)

For every amplitude in the bracket, the perturbed scalar at $Z=-0.13$ equals a quiet scalar value attained in $[-0.130,-0.129]$. The quiet field is strictly monotone there. The scalar invariant $X=g^{ab}\partial_a\phi\partial_b\phi=e^{-2B}(r^2-p^2)$ exceeds 0.290630 on the perturbed event and is below 0.282892 at the quiet event of equal field value. The difference exceeds 0.007738 throughout the entire amplitude bracket, and therefore also at the exact root. In the quiet spacetime $X_0(\phi_0)$ is independent of time. A mere change of slicing or radial label on that background cannot produce the measured gap on this branch.

## Numerical root and independent checks

Shooting gives $\eta\simeq0.002020949495$. Numerical integration in $D$ with DOP853 and an independent RK45 solution of the two original constraints for $B$ agree to $4\times10^{-14}$ in amplitude. At the numerical mass root, the maximum sampled $D$ within the pulse is 0.49675 in the geometric normalization. At $Z=-0.13$, $X$ exceeds its equal-field quiet value by 0.02328. The interval proof uses the wider certified bracket above.

The cancellation can be measured directly from the scalar stress. At this numerical root integrate $\mu'_{\rm pulse}-\mu'_{\rm quiet}$ in increasing $Z$ from $-0.14$ to $-0.10$. Write $K=e^{-2B}(p^2+r^2)/2$ and $P=e^{-2B}pr$. The integrand splits into $(2/3)R^4[Q\,\Delta U+Q\,\Delta K-\Delta P]$. The three integrals are respectively **+0.16815905**, **−0.13312650** and **−0.03503255** in geometric $\mu$ units for one smooth copy and the fixed planar orbit normalization. In particular the change in $K$ includes a cross term with the original scalar gradient and the change in the metric factor. The positive potential difference and negative changes in the other two terms cancel to below $10^{-12}$ here. These are spatial integrals comparing two initial slices. A time-integrated transfer of energy into shell matter requires evolution and a temporal flux measurement.

Separate finite differences of the original constraints give the following maxima after division by the sums of their absolute terms.

| Grid nodes | Hamiltonian | Momentum |
|---:|---:|---:|
| 401 | 1.020 × 10⁻⁵ | 8.389 × 10⁻⁵ |
| 801 | 1.470 × 10⁻⁶ | 1.123 × 10⁻⁵ |
| 1,601 | 1.014 × 10⁻⁷ | 7.724 × 10⁻⁷ |

These residuals test the recorded numerical root and the formulas used for reconstruction. The amplitude digits and derivative diagnostics are numerical estimates. [Mass integration and invariants](MASS_MATCH.json) · [Independent metric integration](METRIC_SHOOT.json)

For the floating-point calculations, save the linked public source summary as `B1_summary.json`, install NumPy and SciPy, and run

~~~sh
python -I -B mass_match.py B1_summary.json
python -I -B verify_lapse.py B1_summary.json
~~~

Both scripts authenticate the summary by SHA-256. The interval calculation needs a C++17 compiler and no third-party numerical library. The reference constants are explicit in its source.

~~~sh
c++ -std=c++17 -O2 -fno-fast-math -ffp-contract=off interval/bracket.cpp -o bracket
./bracket
~~~

[B1 input summary](https://raw.githubusercontent.com/maldonado-research/HDblast/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/runs/main/main_ps0.5_G100_y1_bcut_dzf5e-4_summary.json)
