# Testing null convergence with the entire angular cut

One can make the null convergence test with two smooth, nonnegative angular exposure weights and a single direction field on the whole observer sphere. This removes the need to resolve the observer at a particular nodal label. It uses the same two lower time-window margins and the same physical geometry. [Pointwise curvature derivation](NULL_DIRECTIONS.md) · [Certified time-dependent coefficient](METRIC.md)

For $n=(n_1,n_2,n_3,n_4)$ on the initial unit sphere set
$$
 e(n)=(n_3,n_4,-n_1,-n_2),\quad
 Y=n_1^2-n_2^2,\quad Z=n_3^2-n_4^2,\quad w_\pm=1\pm Z.
$$
The vector $e$ is a smooth nonvanishing round-unit tangent at every $n$. Each $w_\pm$ is nonnegative with round mean one. Push $e$ forward by the actual observer cut map and normalize with its actual spatial metric to obtain $e_\epsilon$. Define $k_\epsilon=u_\epsilon+e_\epsilon$, where $u_\epsilon$ is the normal observer's unit timeline tangent. This gives a physical null direction at every observer. The two readings are proper-area normalized,
$$
 N_{\pm,\epsilon}(x)=
 \frac{\int_{S^3} w_\pm(n)\operatorname{Ric}[h_\epsilon](k_\epsilon,k_\epsilon)\,dA_\epsilon}
 {\int_{S^3} w_\pm(n)\,dA_\epsilon}.
$$
Both denominators are positive on every admitted smooth embedded cut. The background null curvature is zero pointwise, so density, denominator, event and frame derivatives do not enter the first readings. The indices $\pm$ refer to two positive angular weightings of a single field of null directions.

To evaluate them, write $F=U(q)Y(n)$, $a=\rho_b\cosh q$. The background null contraction of its Hessian consists of terms proportional to $Y$, terms proportional to $D_eY=2(n_1n_3-n_2n_4)$, and the remaining term $(U/a^2)\operatorname{Hess}_\kappa Y(e,e)$. Indeed the $qq$ derivative is $U''Y$, the $qi$ components are $(U'-\tanh q\,U)D_iY$, and the purely spatial connection correction is proportional to $\kappa_{ij}Y$. A reversal of the time direction changes only the mixed coefficient. With the round probability measure,
$$
 \langle w_\pm Y\rangle=\langle w_\pm D_eY\rangle=0,
 \qquad \operatorname{Hess}_\kappa Y(e,e)=2(Z-Y).
$$
The first two identities follow by coordinate interchanges or oddness. The sphere Hessian follows by differentiating $Y(n\cos s+e\sin s)$ twice at zero. The standard fourth moments on the unit three-sphere give $\langle Z^2\rangle=1/6$ and $\langle ZY\rangle=0$. Thus $\langle w_\pm\operatorname{Hess}_\kappa Y(e,e)\rangle=\pm1/3$. Applying $J'_{ab}k^ak^b=-2\nabla^2F(k,k)$ yields
$$
 N'_{+,0}=-\frac{2U}{3\rho_b^2\cosh^2q},\qquad
 N'_{-,0}=+\frac{2U}{3\rho_b^2\cosh^2q}.
$$
Upon multiplying by $3\mu\rho_b^4/(2\ell)$ these give precisely $\mp M(x)\operatorname{sech}^2(\ell x)$. Hence each angular reading may also be averaged in time on either of the two registered windows. Its absolute first coefficient in this common positive normalization exceeds $0.03912$ on taking the first time mean and $0.23770$ on taking the second. Uniform parameter differentiability on the entire compact observer cylinder gives opposite signs for the two actual readings at sufficiently small nonzero amplitude on both time windows. A common rescaling of angular weights cancels in the proper-area ratio. Multiplying a curvature reading by a negative or unknown-sign instrumental gain would require calibrating its sign before physical interpretation.

For any exact four-dimensional Einstein metric with matter satisfying the pointwise null energy condition, every integrand in each reading is nonnegative. Therefore any *negative* reading, including its ordinary positive time average, contradicts that comparison under the specified angular and clock calibration. One does not need to infer a negative pointwise value by subtracting two large positive signal channels. A finite-width time or angular exposure here averages curvature along the observer congruence. It carries no conclusion about the complete-geodesic averaged null energy condition.
