# Two-time localization of a physical n=2 perturbation

The result concerns the free, linear scalar perturbation on the **exactly stationary** closed-FRW hilltop, with its quantum state otherwise unrestricted. It is state independent within that model. It does not establish a probability for nonlinear reheating, for a gravitationally constrained quantum universe, or for a stochastic classical path.

## Normalization and transfer

Let x = Htop t, theta = phi/Mpl, and take a real unit-S³ harmonic with

\[
\frac{1}{2\pi^2}\int_{S^3}Y^2\,d\Omega=1,
\qquad -\nabla^2Y=n(n+2)Y.
\]

The perturbation is delta theta(x,Omega)=q(x)Y(Omega). At a stationary scalar background, both the scalar background derivative and the first potential derivative vanish. Consequently its first-order stress-energy perturbation vanishes and the inhomogeneous scalar mode is a physical decoupled scalar perturbation. Here n=2, whose Laplacian eigenvalue is 8 and whose real degeneracy is (n+1)^2=9. The unperturbed expanding closed de Sitter half-bounce has

\[
a(x)=a_b\cosh(h_*x),\qquad a_b=h_*^{-1},\qquad
\frac{H}{H_{\rm top}}=h_*\tanh(h_*x).
\]

The quadratic one-mode action, with hbar=1, is

\[
S_2=\frac{2\pi^2M_{\rm Pl}^2}{H_{\rm top}^2}\int dx\,
\frac{a^3}{2}\left[(q')^2-\left(\frac8{a^2}-\mu^2\right)q^2\right],
\quad \mu^2=-m^2/H_{\rm top}^2>0.
\]

Thus v0=q'(0), [q0,v0]=i kappa, and

\[
\kappa=\frac{H_{\rm top}^2}{2\pi^2M_{\rm Pl}^2a_b^3},\qquad
q''+3h_*\tanh(h_*x)q'+[8h_*^2\operatorname{sech}^2(h_*x)-\mu^2]q=0.
\]

Write q(T)=A(T)q0+B(T)v0, with A(0)=1,A'(0)=0,B(0)=0,B'(0)=1. Then

\[
[q0,q(T)]=i\kappa B(T)\equiv i\eta(T).
\]

The unequal-time commutator uses the bounce normalization; replacing it by a late-time instantaneous [q,q'] would give a different and incorrect coefficient. The Wronskian A B' - A' B = sech^3(h_*T) checks the canonical normalization.

The input branch parameters have mu² > 8h*². Therefore B(T)>0 and B'(T)>0 for T>0: starting at B=0,B'=1, the equation (a³B')'=a³(mu²-8/a²)B forbids a first loss of either positivity. All critical-time roots below are unique.

An independent closed form follows by setting y=h*T, beta=sqrt(mu²/h*²+9/4), nu=n+1/2:

\[
\begin{aligned}
A(T)&=\cosh^{-(n+2)}y\;
{}_2F_1\!\left(\frac{\beta-\nu}2,\frac{-\beta-\nu}2;\frac12;-\sinh^2y\right),\\
B(T)&=\frac{\sinh y}{h_*}\cosh^{-(n+2)}y\;
{}_2F_1\!\left(\frac{1+\beta-\nu}2,\frac{1-\beta-\nu}2;\frac32;-\sinh^2y\right).
\end{aligned}
\]

For n=2, B(T)=T+(mu²-11h*²)T³/6+O(T⁵). The code compares this formula to an independently integrated ODE.

## RMS constraint and its exact attainability

For any density operator with finite second moments and the usual real-line canonical representation,

\[
\Delta q0\,\Delta q(T)\ge\frac{|\eta(T)|}2.
\]

The same lower bound applies to the product of root mean square deviations about any two prescribed field values, since those root mean squares are at least the corresponding standard deviations. No energy bound or Gaussian assumption is needed for the inequality. It is sharp: a Gaussian in q0 with width sigma0 and a quadratic phase giving Cov(q0,q(T))=0 has Delta q(T)=|eta|/(2 sigma0), with both means freely adjustable. Explicitly, take

\[
\psi(q)\propto\exp\left[-\frac{(q-c_0)^2}{4\sigma_0^2}
+\frac{i}{\eta}\left\{(c_T-Ac_0)(q-c_0)-\frac A2(q-c_0)^2\right\}\right].
\]

Thus the existence of a state with RMS tolerances d0,dT is exactly equivalent to d0 dT >= |eta|/2. Arbitrary one-time squeezing does not evade this two-time product. For equal tolerances epsilon, the last possible time solves B(T)=2 epsilon²/kappa.

## A probability result without moment or energy assumptions

Let P=1{|q0-c0|<=d0} and Q=1{|q(T)-cT|<=dT}. For any one initial density operator rho, define p0=Tr(rho P), pT=Tr(rho Q). The projection inequality is

\[
p_0+p_T\le1+\|PQ\|.
\]

The eigenfunctions of Aq0+Bv0 have constant-modulus overlap with those of q0,

\[
|\langle q0=x\mid q(T)=y\rangle|^2=\frac1{2\pi|\eta|}.
\]

Consequently P Q is Hilbert-Schmidt, and its squared Hilbert-Schmidt norm is

\[
\operatorname{Tr}(PQP)=\|PQ\|_{\rm HS}^2
=\frac{(2d_0)(2d_T)}{2\pi|\eta|}.
\]

The trace is written as Tr(PQP), a positive trace-class operator; this avoids using an unqualified trace of a product of infinite-dimensional projections. Since operator norm is at most Hilbert-Schmidt norm,

\[
\boxed{\quad p_0+p_T\le1+\min\left\{1,
\sqrt{\frac{2d_0d_T}{\pi|\kappa B(T)|}}\right\}.\quad}
\]

In particular, both marginal probabilities can be at least p>1/2 only if

\[
\frac{d_0d_T}{|\kappa B|}\ge\frac\pi2(2p-1)^2.
\]

For p=.95 the necessary RHS is 1.2723450247; for p=.99 it is 1.5085927922. This is a rigorous necessary condition and does not assert sufficiency.

If the initial marginal is already specified to be at least pinit, a stronger useful form is the Landau-Pollak angle relation

\[
\arccos\sqrt{p_0}+\arccos\sqrt{p_T}\ge\arccos s,
\quad s=\|PQ\|.
\]

For pinit>=s² this yields

\[
p_T\le\left[s\sqrt{p_{\rm init}}+
\sqrt{(1-s^2)(1-p_{\rm init})}\right]^2.
\]

For a pure state this follows from the angular triangle inequality between its normalized P and Q projections; their overlap is at most s. Purification gives the same conclusion for mixed states because tensoring the projections with the identity preserves their operator overlap. Replacing s by the Hilbert-Schmidt upper bound gives a rigorous numerical ceiling when that upper bound is smaller than sqrt(pinit); otherwise the trivial ceiling 1 is used. At the departure reference times, s is tiny: with initial marginal probability >=.99, the final marginal is at most about .01, with the largest correction across the tabulated cases about 3.64e-7. This describes the same initial density operator and no intermediate measurement. It is not the conditional probability after measuring the initial interval.

The inequality holds for mixed states, non-Gaussian states, infinite-variance states, and states entangled with spectators, by applying it to the mode's reduced density operator. It needs the regular real-line canonical representation. Compact field domains or other boundary conditions require a new spectral analysis. It does not rely on a vacuum choice.

## Exact interval comparison with the prolate operator

Conjugating away the chirp from A and the phases from the two interval centers and scaling x=d0 u gives the exact operator

\[
(K_c f)(u)=\int_{-1}^1\frac{\sin[c(u-v)]}{\pi(u-v)}f(v)\,dv,
\quad c=\frac{d_0d_T}{|\eta|},
\quad K_c(u,u)=c/\pi.
\]

If lambda0(c) is its top eigenvalue, then **exactly**

\[
\|PQ\|=\sqrt{\lambda_0(c)},\qquad
\max_\rho\min(p_0,p_T)=\frac{1+\sqrt{\lambda_0(c)}}2.
\]

The second identity is attainable in the unrestricted density-operator class. If u lies in ran P and PQP u=s²u with s=sqrt(lambda0)>0, set v=Qu/s and psi=(u+v)/sqrt(2(1+s)). Then P v=s u and Q u=s v, so both probabilities equal (1+s)/2. This proves that the sharp two-time probability feasibility condition is lambda0(c)>=(2p-1)², independent of either interval's center and of A(T). The attaining vector can have infinite kinetic variance because of sharp interval edges. Smooth finite-energy vectors approximate it arbitrarily well, so the same value is a supremum without a fixed energy budget; attainment at equality is not asserted with an energy constraint.

Numerical roots using 96 and 160 Gauss-Legendre nodes agree to about 1e-14:

| desired probability at each time | exact-operator critical c (numerical) |
|---:|---:|
| .95 | 1.669665916406 |
| .99 | 2.706071984932 |

The exact operator identity is analytic; the decimal eigenvalues are numerical approximations with quadrature convergence checks, not interval-arithmetic certificates. The Hilbert-Schmidt inequality remains rigorous without these numerical eigenvalues. For c much smaller than one, lambda0(c)=2c/pi+O(c³), so the simple bound is asymptotically sharp.

## Full n=2 multiplet and spatial RMS balls

Choose all nine real n=2 harmonics with the same RMS normalization and orthogonal indices j. Then [q0_j,v0_k]=i kappa delta_jk, q_j(T)=Aq0_j+Bv0_j, and

\[
\frac1{2\pi^2}\int_{S^3}\delta\theta^2\,d\Omega
=\sum_{j=1}^9q_j^2.
\]

The spatial field RMS tolerance epsilon is therefore the Euclidean configuration ball |q|<=epsilon; it is invariant under a change of the real harmonic basis. No pointwise supremum bound is identified with it.

With vector RMS deviations r0=sqrt(<|q0-c0|²>), rT=sqrt(<|qT-cT|²>), coordinate Robertson inequalities and Cauchy-Schwarz give

\[
r_0r_T\ge\frac{9|\eta|}{2},
\]

with equality attainable by an isotropic product of the sheared Gaussians above. For ball projectors, the Fourier overlap kernel has squared modulus (2 pi |eta|)^(-9), and

\[
p_0+p_T\le1+\min\left\{1,
V_9\left(\frac{\epsilon_0\epsilon_T}{2\pi|\eta|}\right)^{9/2}\right\},
\quad V_9=\frac{32\pi^4}{945}=3.29850890274.
\]

Entanglement among the modes does not evade this bound. Necessary c=epsilon0 epsilonT/|eta| values for both probabilities at least .95 and .99 are respectively 4.70792824246 and 4.79786921144. These ball conditions are necessary, not sharp. The exact ball answer is a nine-dimensional truncated Fourier concentration problem; it is not the one-dimensional prolate lambda0 reused unchanged. One further upper comparison follows by enclosing the balls in the cubes [-epsilon,epsilon]^9: ||Pball Qball|| <= lambda0(c)^(9/2).

## What the numerical horizons mean

The code reconstructs Htop/Mpl=2 pi times the archived displacement scale, exactly as in the previous quantum note. This reconstruction is an input convention of that calibrated study, not a universal identity for arbitrary perturbation states. Here h* differs from one only by floating-point roundoff, giving kappa=3.93515010206e-6, 3.93515010206e-13 and 3.93515010206e-23 at detunings .001, 1e-10 and 1e-20.

`reference_homogeneous_departure_time_x_match` is the archived **homogeneous** at-rest departure time reaching a 1e-5 field displacement. It supplies a concrete comparison duration. The n=2 transfer is computed independently on the stationary hilltop for that duration; the calculation does not claim that an inhomogeneous quantum realization follows the homogeneous at-rest solution or its reheating outcome. Extending to the archived first minimum would go beyond the small-field stationary approximation; no such extension is used here.

At each epsilon=.01,.001,.0001, `sharp_one_mode_last_feasible_time` is the exact-operator two-endpoint feasibility horizon, evaluated numerically. Beyond it, no state can give both q0 and q(T) the specified Born probability inside their prescribed intervals. Up to it such a state does exist for these two marginal conditions, with no energy constraint; intermediate-time localization and gravitational backreaction are additional requirements. The reported initial-width requirements invert the same formula at the departure reference time. When that required initial width is enormous, it lies outside the intended small-field preparation range; it should not be interpreted as a physically acceptable stationary-background construction.

## Marginals, survival, and measurement disturbance

These are same-state Heisenberg marginal probabilities. In general P and Q do not commute, so there is no specified positive joint probability distribution for them. The result does **not** assert that at most about one half of trajectories survive near the hilltop; it asserts that the smaller of the two separate Born weights is at most the stated ceiling. A state may put nearly all its weight in either single interval while the other marginal is broad.

A sequential projective experiment has joint probability Tr(Q P rho P), which uses the state disturbed by the first measurement. It is a different experiment, not p0*pT and not min(p0,pT). For ideal sharp projections it also obeys Tr(Q P rho P)<=||PQ||² p0, by P Q P<=||PQ||² P, but that statement concerns this particular measurement protocol. Continuous monitoring, postselection, absorbing boundaries, and open-system decoherence define further distinct experiments. An endpoint localization obstruction is a necessary obstruction for any proposal that explicitly demands both the initial and final unmeasured marginals; it is not by itself a survival or reheating probability.

## Reproduction

From the preparation study directory:

```bash
python bounds/two_time.py --input prior/hilltop/HILLTOP.json --output bounds
```

The script writes `BOUNDS.json` and `BOUNDS.md`, compares direct ODE transfer against the independent closed form, checks the Wronskian at high precision, and compares prolate critical-c quadratures. No remote resource is edited.
