# The scalar mean, signed quadrupole and angular variance

Fix the registered affine amplitude $\epsilon$ and the proper-clock label $q=-\tau/\rho_b$. Let $\Phi_\epsilon(q,n)$ be the scalar on the moved observer cut. Write $E_\epsilon[f]=\int f j_\epsilon\,d\mu_\kappa/\int j_\epsilon\,d\mu_\kappa$, where the positive $j_\epsilon$ is the cut's actual density. On the background $\Phi_0=\phi_b$ and $j_0=\rho_b^3\cosh^3 q$ is constant over the sphere. The following three quantities can be defined on every admitted cut.

$$
M_\epsilon=E_\epsilon[\Phi_\epsilon],\qquad
D_\epsilon=\frac{E_\epsilon[Y(\Phi_\epsilon-\phi_b)]}{E_\epsilon[Y^2]},\qquad
V_\epsilon=E_\epsilon[(\Phi_\epsilon-M_\epsilon)^2].
$$

The template $Y=n_1^2-n_2^2$ is transported from the original labels. Positivity of the metric gives $E_\epsilon[Y^2]>0$. Variance uses the mean of the *same* cut and is independent of how its labels are named. The template orientation matters for the signed quantity $D$. [Actual observer and positive measures](../HDBLAST_EFFECTIVE_OBSERVER_20261116/EFFECTIVE_OBSERVER_DOMAIN_AND_RESPONSE_BOUNDS.md) · [Earlier definition of the matched quadrupole](../HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/nonlinear/NONLINEAR_BRIDGE.md)

## Derivatives at zero amplitude

The field and clock identification gives $\Phi'_0=S(-q)Y$. The uniform distribution on the unit three-sphere is invariant under exchanging $n_1,n_2$, so $E_0[Y]=0$. To obtain the normalization, put $a=E_0[n_1^4]$ and $b=E_0[n_1^2n_2^2]$. A rotation through $\pi/4$ in this plane and reflection symmetry give $a=(2a+6b)/4$. Squaring $\sum_i n_i^2=1$ gives $4a+12b=1$. Thus $a=1/8$, $b=1/24$ and $E_0[Y^2]=1/6$. [Field identification](BRIDGE.md)

At zero amplitude all terms from the measure derivative in $M'_0$ cancel by normalization. In $D'_0$ they multiply the vanishing background contrast. This yields $M'_0=0$ and $D'_0=S(-q)$. Write $G_\epsilon=\Phi_\epsilon-M_\epsilon$. Twice differentiating $E_\epsilon[G_\epsilon^2]$ at fixed clock gives

$$
 V_0=V'_0=0,\qquad V''_0=2E_0[(G'_0)^2]=S(-q)^2/3.
$$

Terms containing a first or second density derivative, or the second derivative of the scalar, carry a factor $G_0=0$. This cancellation accounts for the moving measure. By the signed parameter symmetry $M$ and $V$ are even and $D$ is odd. The established $C^2$ dependence into continuous functions on the compact observer cylinder gives, uniformly in $q$,

$$
V_\epsilon(q)=\epsilon^2 S(-q)^2/6+o(\epsilon^2),\qquad
D_\epsilon(q)=\epsilon S(-q)+o(\epsilon^2).
$$

On a closed interval with $|S|\ge c>0$, for all sufficiently small nonzero $|\epsilon|$ the variance is positive and $D_\epsilon/\epsilon$ has the sign of $S$. For example, $V_\epsilon(q)>\epsilon^2(34748)^2/12$ throughout the window below once the uniform variance remainder divided by $\epsilon^2$ is smaller than $(34748)^2/12$. This is an existence statement for a subinterval of the admitted $|\epsilon|\le10^{-10^{500209}}$. Evaluating a usable amplitude threshold requires bounding a remainder. The two point witnesses in the replay have respectively negative and positive signs at $q/\ell=74999/50000$ and $149999/100000$. Uniform convergence and continuity thus give at least one sign-changing transition, with a zero of the finite-amplitude quadrupole strictly between them. Neither uniqueness of a zero nor the vanishing of the nonlinear variance there follows. [Exact signed-family covariance](SYMMETRY.md) · [Previous characteristic transfer and point witnesses](../HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/primary/PROOF.md)

## Bounding the whole observer interval

Write $x=q/\ell$ and $H=25/16$. A historical row with integer endpoints $L_i,U_i$ bounds $\ell S(+\ell x)$ between $L_i/2^{112}$ and $U_i/2^{112}$ throughout its **closed** time cell. The registered positive and negative time evolutions have opposite elapsed initial velocity. The characteristic transfer theorem bounds $|S(-q)-S(q)|<1/40$ in the units of $S$. Applying this *after* dividing by $\ell$ gives the interval

$$
 [a_i,b_i]=[L_i/(2^{112}\ell)-1/40,\; U_i/(2^{112}\ell)+1/40].
$$

The CSV contains 4,353 contiguous cells spanning $0\le x\le3$. Exactly 512 cells cover the window $[5/4,11/8]$. The smallest lower and largest upper bounds over this union give $-46653<S(-q)<-34748$ throughout the closed window. The upper end of the last cell and the lower end of the first are included. [Continuous future cells, §4](../HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/INERTIAL_SHELL_WAVEFORM.md) · [Past transfer, §§1–6](../HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/primary/PROOF.md)

For the whole interval define $\overline V_\epsilon=Q^{-1}\int_0^Q V_\epsilon(q)\,dq$ with $Q=\ell H$. The same average in physical proper duration $r=\rho_b q$ has the same value. The uniform Taylor remainder gives

$$
 \lim_{\epsilon\to0}\frac{\overline V_\epsilon}{\epsilon^2}
 =\frac{1}{6H}\int_0^H S(-\ell x)^2\,dx=C.
$$

For each cell intersected with $[0,H]$, multiply its overlap length by the lower and upper bounds for the square. The lower end is zero for an interval containing zero, and otherwise the smaller endpoint square. The upper end is the larger endpoint square. Sum these intervals and divide by $6H$. The 2,305 overlaps have total length $H$. Shared endpoints have zero integration measure. The calculation conservatively applies the global transfer bound even to the exact linear prefront. The resulting strict integer bounds are

$$
 35{,}373{,}971<C<43{,}294{,}110.
$$

These are bounds on the coefficient of $\epsilon^2$ in the averaged variance. Its second derivative is $2C$. The averaging interval is fixed across amplitudes and fits inside the certified observer domain. The data determine neither an injected energy nor a variance of the observed cosmic microwave background. [Exact rational output](RESULTS.json) · [Reproduction code](replay.py)
