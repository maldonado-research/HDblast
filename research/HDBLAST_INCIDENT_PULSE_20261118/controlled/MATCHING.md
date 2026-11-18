# Gluing an incident scalar pulse to the quiet B1 geometry at both ends

In the identity conformal chart of the tuned, planar B1 model, a compact scalar pulse can satisfy both initial Einstein constraints and join the quiet reference on an open collar at each end. The construction uses three smooth changes to the initial metric time derivative. They may all lie farther from the shell than the scalar pulse. A local existence result covers sufficiently small scalar amplitudes. Specific finite-amplitude examples have independently measured numerical matching errors.

The chart has metric $e^{2B}(-dT^2+dZ^2)+e^{2A}d\mathbf{x}^2$. The quiet reference is $A_0=B_0(Z)+T$, $\phi_0=\phi_0(Z)$, with $u_0=A_{0,T}=1$, $B_{0,T}=\phi_{0,T}=0$. Set $Q=A_Z$, $u=A_T$, $p=\phi_T$, $r=\phi_Z$, $v=B_T$. Empty matter on the shell initially gives $A_Z=B_Z=e^B\sigma/6$ and $\phi_Z=-e^B\sigma_\phi/2$. The source laws, units and the earlier one-pulse calculation are documented in the [original constraint experiment](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-6992402829) and [B1 solver](https://github.com/maldonado-research/HDblast/blob/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py).

## Constraint and matching equations

On the initial slice choose $A=B$, $\phi=\phi_0+F$ and $p=-F'$, where $F=\eta f(Z)$ is a specified smooth compact bump. The scalar principal characteristic of this perturbation initially travels towards the shell, at $dZ/dT=+1$. Primes denote $Z$ derivatives. With $r=r_0+F'$, the two B1 constraints give

$$
 B'=Q,\qquad Q'=u^2+uv-Q^2-(p^2+r^2)/6-e^{2B}U(\phi)/3,\qquad
 u'=Qv-pr/3.
$$

Integrate from the quiet wall to a cut $Z_L$ left of all supports. The endpoint residual is $(B-B_0,Q-Q_0,u-1)(Z_L)$. In the quiet neighborhood of the cut the ODEs, prescribed scalar and $v$ coincide with the reference; equality of these three values there gives equality of all initial fields on an open collar by ODE uniqueness. In the open collar at the wall the same argument preserves all three wall junctions and initial data compatibility inherited from the quiet solution. Reflecting the smooth-side data into the original second, symmetric copy preserves its matching convention.

## Existence for all sufficiently small amplitudes

Assume a regular analytic quiet background on a connected patch containing the wall, scalar support, control interval and cut. Linearize the endpoint map at $\eta=0,v=0$, and write $x=(\delta B,\delta Q,\delta u)$, $E=e^{2B_0}$, $Q=Q_0$, $r=r_0$, with $U,U_\phi$ evaluated at $\phi_0$. It obeys

$$
 x'=Lx+g\,\delta v,\qquad
 L=\begin{pmatrix}0&1&0\\-2EU/3&-2Q&2\\0&0&0\end{pmatrix},
 \quad g=(0,1,Q)^T .
$$

Let $\Phi(Z_L,Z)$ transport homogeneous variations from $Z$ to the cut. The response to a compact control is the signed integral of $\Phi g\,\delta v$ from the wall to the cut. For $D=d/dZ-L$, derivatives of $\Phi g$ are $\Phi Dg$ and $\Phi D^2g$. Directly,

$$
 Dg=(-1,0,Q')^T,\qquad D^2g=(0,-2EU/3-2Q',Q'')^T .
$$

Differentiate the quiet constraint $Q'=1-Q^2-r^2/6-EU/3$ and use the static scalar equation $r'=EU_\phi-3Qr$. This gives

$$
 \det(g,Dg,D^2g)=Q''+Q(2EU/3+2Q')=Qr^2-2ErU_\phi/3 .
$$

At the wall this determinant divided by $\rho_b^3$ is $\sigma\sigma_\phi^2/24+\sigma_\phi U_\phi/3$. The sign is strictly negative at the archived B1 wall field and tension parameters. An analytic function that is nonzero at the wall cannot vanish identically on any open subinterval of the regular connected patch. Thus any selected open control interval contains a point with nonzero determinant. A covector annihilating every compactly supported smooth-control response would annihilate $\Phi g$ throughout the interval and also its first two derivatives there. These three vectors are independent at that point, so the covector must vanish. This proves surjectivity to the three-dimensional endpoint space. Three fixed smooth compact controls can therefore be chosen to give an invertible derivative.

Use their three coefficients and the fixed scalar amplitude as parameters. Smooth ODE dependence and the finite-dimensional implicit function theorem give coefficients, smooth in $\eta$ near zero and zero at zero, satisfying **all three endpoint conditions exactly**. The control interval can be chosen strictly left of the pulse. The argument gives no explicit radius in $\eta$ and does not certify an exact root for any particular finite-amplitude numerical example or prescribed three basis shapes.

For the mathematical static background the matching and quiet first integral require $\rho_b^{-2}=\sigma^2/36-\sigma_\phi^2/48+U/6$ at the wall. The archived decimal field and radius only satisfy this identity to rounding. An exact compatible radius may be assigned from the indicated expression with the archived decimal B1 field and tension parameters. Nonvanishing of the determinant after division by $\rho_b^3$ is independent of that radius adjustment. The finite computations separately declare which convention they use. Neither construction establishes regularity through an arbitrary global bulk extension.

## The scalar disturbance can be distinguished from a coordinate change

Choose the controls behind the pulse and consider the regular quiet branch on which $r_0\ne0$ at the scalar support. The invariant $X=g^{ab}\partial_a\phi\partial_b\phi$ on the reference is a single-valued function $X_0(\phi_0)=e^{-2B_0}r_0^2$ on that monotone branch. Any coordinate description of the same quiet solution there obeys the same relation. If it held through the nonzero pulse as well, $p=-F'$ would give

$$
 F'=[e^{2B}X_0(\phi_0+F)-r_0^2]/(2r_0).
$$

This extra ODE, together with the constraint equations and $v=0$ between the wall and the compensators, is locally Lipschitz in $(F,B,Q,u)$. Its solution from the quiet data on the right of the pulse has $F=0$ by uniqueness. Therefore a nonzero compact pulse forces a change in the scalar-invariant relation somewhere on that branch. This excludes a mere relabeling of the specified quiet branch. The condition $p=-F'$ describes initial data and does not enforce exclusively one-way propagation during evolution.

## Concrete profiles

Each basis profile is $h=\exp[1-1/(1-s^2)]$ for $|s|<1$, with $s=(Z-c)/w$, and vanishes outside. The specified scalar pulse uses $c=-0.12,w=0.02$. For the delayed controls the three centers are $(-0.17,-0.22,-0.28)$ and the common width is $0.018$, with matching at $Z_L=-0.32$. At $\eta=+0.001$ the coefficients of $v=\sum_j k_jh_j$ are approximately $(-0.370528,0.678395,-0.297272)$; at $\eta=-0.001$ they are $(1.357949,-2.441326,1.048063)$. The scalar support nearest the wall is $Z=-0.10$ and the nearest control support is $Z=-0.152$. A measurement at the wall over $0.10<T<0.152$ can therefore receive the original scalar disturbance and its required metric adjustment before any part of those compact controls can causally arrive. No received radiation or energy is established by the initial-data calculation alone.

The interleaved comparison uses centers $(-0.05,-0.095,-0.158)$ with common width $0.022$ and cut $-0.20$. It solves the same matching problem with smaller metric-rate coefficients. In that case the nearest metric adjustment is at $Z=-0.028$.

The portable NumPy/SciPy module matched_pulse.py returns quiet, trial and difference arrays for $A,B,\phi,A_T,B_T,\phi_T$. It verifies the B1 summary SHA-256 and extends with the original quiet data outside the matched patch. Examples use the identity conformal chart, zero matter and $dc=0$; the archived production run used a different initial displacement and chart.

## Reproducing and assessing the finite roots

The adaptive calculation uses the archived radius $7.838285538073243$ and solves the reference Hamiltonian constraint. The independent standard-library calculation assigns the compatible radius $7.838285538074238$, evolves the static reference using $Q_0'=3-3Q_0^2-2e^{2B_0}U/3$, and applies fixed-step classical RK4 to the trial constraints and their variational equations. The endpoint measurements use the initial-value solution directly, before the generator fills the quiet extension. The replay also forms the original constraints using centered finite differences of the fields, and compares $X$ against the quiet $X_0$ at equal scalar field by interpolation.

| Basis | Scalar amplitude | Largest absolute endpoint residual, 1,600 intervals | Largest residual with controls omitted, 800 intervals | $X-X_0(\phi)$ at $Z=-0.11$ |
|---|---:|---:|---:|---:|
| Delayed | +0.001 | $6.3\times10^{-11}$ | $1.95\times10^{-4}$ | −0.014837 |
| Delayed | −0.001 | $6.4\times10^{-11}$ | $5.20\times10^{-4}$ | +0.014837 |
| Delayed | +0.0005 | $6.2\times10^{-11}$ | $1.10\times10^{-4}$ | −0.007418 |
| Delayed | −0.0005 | $6.7\times10^{-11}$ | $2.06\times10^{-4}$ | +0.007418 |
| Interleaved | +0.001 | $2.2\times10^{-12}$ | $1.95\times10^{-4}$ | −0.014827 |
| Interleaved | −0.001 | $1.9\times10^{-12}$ | $2.96\times10^{-4}$ | +0.014800 |

The quiet endpoint Jacobians in the native unscaled $(B,Q,u)$ basis have infinity-norm condition numbers about 1,372 (delayed) and 234 (interleaved). This measures sensitivity to the choice of metric-rate coefficients and coordinate normalization. The finite-difference normalized Hamiltonian and momentum residuals are at most $1.31\times10^{-5}$ at 1,600 intervals across these cases. The numerical results do not give rigorous enclosures for the finite roots. The source decimal gives $-256.106$ for the geometric control determinant at the wall.

Run the check from this directory:

```sh
python -I -B replay_matching.py
```

The full numerical table is in [MATCHING.json](MATCHING.json). [Input values and coefficient triples](MATCHING_INPUT.json) identify the exact original archived summary. To repeat the nonlinear shooting calculation with NumPy and SciPy, supply that summary and run

```sh
python -I -B shoot_controls.py B1_summary.json delayed -0.001 -0.0005 0.0005 0.001
```

[Initial-data generator](matched_pulse.py) · [Root solver](shoot_controls.py) · [Independent RK4 check](replay_matching.py)
