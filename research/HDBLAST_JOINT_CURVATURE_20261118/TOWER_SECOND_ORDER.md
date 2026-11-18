# A second-order residual for the fixed four-dimensional completion

The linear comparison that reproduces both the entire induced metric and the scalar makes a further prediction once its nonlinear couplings are specified. Keep the free massive canonical fields with positive kinetic terms, quadratic potentials, cosmological constant, $h_J=e^{2B}g_E$, $B=\sum b_j\psi_j$ and $\Phi=\phi_b+s$, $s=\sum d_j\psi_j$. Both $B$ and $s$ here are linear functions of the **full** fields. The shell identity can be tested with the intrinsic curvature and both the angular and temporal scalar gradients:
$$
 \mathcal C=R(h_J)-h_J^{ab}\partial_a s\,\partial_b s-F(\phi_b+s).
$$
The function $F$ comes from the two-sided shell Gauss and junction equations, with $F(\phi_b)=R_0$ and $F'(\phi_b)=K$. In the actual registered shell $\mathcal C=0$. The following calculation gives a necessary second-order prediction for the chosen four-dimensional completion. [Linear completion](COMPARISON.md) · [Background equations and graph junctions](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md)

## Eliminate the unknown second-order fields

Use $G_{ab}(g_E)+\Lambda g_{E,ab}=\kappa_4 T_{ab}$ with $\kappa_4>0$. Gradients and contractions on the right below use $g_E$ in signature $(-+++)$. The exact Einstein trace and uncoupled scalar equations are $R_E=R_0+\kappa_4\sum_j[(D\psi_j)^2+2M_j^2\psi_j^2]$ and $\Box_E\psi_j=M_j^2\psi_j$. Consequently
$$
\mathcal C=e^{-2B}\left[R_0+\kappa_4\sum_j\big((D\psi_j)^2+2M_j^2\psi_j^2\big)
-6\sum_j b_j M_j^2\psi_j-6(DB)^2-(Ds)^2\right]-F(\phi_b+s).
$$
The complete linear term vanishes mode by mode by $-(2R_0+6M_j^2)b_j=K d_j$. Hence terms linear in unknown second-order Cauchy corrections cancel too, including any allowed second-order metric solution. In the remaining quadratic terms one may use $h_0$ to contract first-order gradients. Curvature has already been eliminated by the field equations. This argument applies to a sufficiently differentiable on-shell family, or as a necessary identity of its on-shell perturbation jets. With countably many fields assume convergence strong enough to take the indicated derivatives and energy sums. The smooth first-order projections have the required graph norms. This test supplies no nonlinear existence theorem for an infinite-species Einstein system.

## Initial angular pattern

At $q=0$ write the prescribed first-order fields and their elapsed $q$ derivatives as $\alpha_j Y$ and $\beta_j Y$, with $Y=n_1^2-n_2^2$ on the unit $S^3$. The first-order readings $B,s$ and their gradients vanish on this initial cut because the radial data are supported a positive distance from the endpoint. The coefficient of $\epsilon^2$ in the residual is therefore
$$
 \mathcal C_{[2]}(0,n)=
 \frac{\kappa_4}{\rho_b^2}\left[
 \Big(\sum_j\alpha_j^2\Big)|\nabla_{S^3}Y|^2
 +\Big(\sum_j(2\rho_b^2M_j^2\alpha_j^2-\beta_j^2)\Big)Y^2\right].
$$
This convention uses the Taylor coefficient, half of the second derivative. Terms from a small change of clock, angle labels or area measure multiply vanishing background and first variations of the scalar residual.

At $n=(1,1,0,0)/\sqrt{2}$ the pattern has $Y=0$ and $|\nabla Y|^2=4$. Thus $\mathcal C_{[2]}=4\kappa_4\|\alpha\|_2^2/\rho_b^2>0$ there. A full-cut average isolates the same effect. Put $r=n_1^2+n_2^2$ and use $W=3-4r$ with the transported angle labels. With normalized round measure, $r$ is uniform on $[0,1]$ and the polar angle $\theta$ in the first plane is uniform independently. Then $Y=r\cos 2\theta$ and $|\nabla Y|^2=4(r-Y^2)$. Consequently
$$
 \langle W Y^2\rangle=0,\qquad
 \langle W |\nabla Y|^2\rangle=\tfrac23,\qquad
 \langle W \mathcal C_{[2]}\rangle=\frac{2\kappa_4}{3\rho_b^2}\sum_j\alpha_j^2>0.
$$
The weight is signed and has mean one. Integrating with the actual normalized cut area gives the same Taylor coefficient, since the residual and its first variation vanish pointwise at the background.

For the bounded-detector field normalization in the published tower,
$$
 \alpha_j=(5/2-\lambda_j)\langle E_j,U\rangle,\quad
 U=(Q(z/\ell),0),\quad
 \|\alpha\|^2=\|(A-5/2)U\|^2 .
$$
The profile is nonzero in $\ell<z<2\ell$. It belongs to $H^1_0$ of that interval. Cauchy–Schwarz, the Dirichlet Poincaré inequality and $V>5547$ on the collar imply $\|(A-5/2)U\|\ge(\pi^2/\ell^2+5547-5/2)\|U\|$. On the middle third $Q^2\ge e^{-1}>1/3$ by its explicit definition, so $\|U\|^2>\ell/9$ and $\|\alpha\|^2>9/\ell^3$. The weighted residual coefficient exceeds $6\kappa_4/(\rho_b^2\ell^3)$ in these conventions. The value is a model coefficient, before multiplication by $\epsilon^2$. Changing the canonical preparation normalization and compensating with detector gain can make the coefficient arbitrarily small over a class of completions. A physical use requires an energy, gain, gravitational coupling and precision budget. [Explicit packet](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/INERTIAL_SHELL_WAVEFORM.md) · [Operator and spectral preparation](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/TOWER.md) · [Angular moment replay](conformal_test.py)

The nonzero coefficient excludes this fixed completion at second order under the stated regularity and coupling hypotheses. Nonlinear changes to the conformal coupling or the operational scalar readout alter the second-order prediction even when their values and first derivatives at the background agree. The comparison therefore constrains specified dynamics and readouts. The dimension of the universe remains open.
