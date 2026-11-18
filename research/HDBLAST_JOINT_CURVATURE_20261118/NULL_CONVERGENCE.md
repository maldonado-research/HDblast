# A finite, positive angular average detects negative intrinsic null convergence

Does the induced four-dimensional geometry satisfy the null convergence condition $R_{ab}k^a k^b\ge0$ in every null direction? The background de Sitter shell saturates this inequality. The registered packet departs from it with opposite signs in two directional averages. Both can be formed with nonnegative weights over the entire angular cut. This also gives a test of any four-dimensional Einstein model with positive gravitational coupling for which the induced metric is the physical metric and the total matter obeys the null energy condition.

## A measurement on a cut

Use the transported angular label $n\in S^3\subset\mathbb R^4$. For each label draw a unit tangent direction $v$ from the uniform two-sphere in $T_nS^3$ of the reference unit metric $\kappa$. Push $v$ to the actual proper-clock cut along its label map and divide by its actual spatial length to obtain the unit spatial vector $e_\epsilon(v)$. The shell observers are orthogonal to these cuts. Their unit timelike velocity $u_\epsilon$ thus gives a true null direction $k_\epsilon=u_\epsilon+e_\epsilon(v)$. This convention fixes its scaling. Use the actual normalized proper-area measure over the cut. Set
$$
Y=n_1^2-n_2^2,\qquad D=v_1^2-v_2^2,\qquad Z=D+Y/3,\qquad w_\pm=1\pm Z .
$$
These are two sampling protocols. For an orthonormal pair $(n,v)$, the bound $v_i^2+n_i^2\le1$ gives $|Z|\le1$. Both weights are nonnegative. Given $n$, the tangent-direction average satisfies $\mathbb E[v_i^2\mid n]=(1-n_i^2)/3$, hence $\mathbb E[Z\mid n]=0$. Each protocol has mass one even with the actual perturbed cut measure. Define $N_\pm(\epsilon,q)$ to be the corresponding average of $\mathrm{Ric}(h_\epsilon)(k_\epsilon,k_\epsilon)$.

## Linear response and sign

The background curvature is proportional to the metric. Replace its Ricci tensor under this contraction by $\mathrm{Ric}(h_\epsilon)-3h_\epsilon/\rho_b^2$, which vanishes as a tensor at zero amplitude. Thus smooth changes of frame, sampling measure, observer positions and gauge make no extra first-order contributions. The physical metric derivative is conformal with factor $2U(q)Y$, up to a coordinate change. The conformal Ricci formula in four dimensions leaves only $-2\operatorname{Hess}_{h_0}(UY)(k,k)$ in the null contraction. [Intrinsic metric and gauge identification](METRIC.md)

Write $a(q)=\rho_b\cosh q$. The purely temporal Hessian and the expansion part of the spatial Hessian multiply $Y$ and average to zero under both weights. The time–space terms are odd in $v$ and also average to zero. The remaining spatial part uses $\operatorname{Hess}_{\kappa}Y(v,v)=2(D-Y)$. At the background both $n$ and $v$ marginally have the uniform law on $S^3$. Thus $\mathbb E D^2=\mathbb E Y^2=1/6$, $\mathbb E YD=-1/18$ and $\mathbb E ZD=4/27$. The full angular response is
$$
 N'_\pm(0,q)=\mp\frac{16}{27}\,\frac{U(q)}{a(q)^2}.
$$
The fixed axis labels can be chosen using the specified packet polarization. Interchanging labels 1 and 2 interchanges the two protocols. The averaging uses an ensemble of local null directions at each timelike-observer cut. It does not integrate along a null geodesic.

Using $M=\mu\rho_b^2 U/\ell$, the normalized exposure is $\mp(16/27)M\operatorname{sech}^2q$. The first metric window has mean $M>0.03912$ and the second has mean $M>0.23770$ and the factor $\operatorname{sech}^2q$ changes those enclosures by less than $2\times10^{-12}$ in $M$ units. The independent calculator integrates the enclosed accelerations cell by cell and recovers the same metric window endpoints as the retarded-kernel calculation. It gives exact rational intervals for both signs. [Full-sky replay](null_sky.py) · [Signed null exposures](NULL_SKY.json)

Both exposures vanish at zero amplitude. Their clock averages have nonzero first derivatives of opposite sign. Smooth dependence of the actual induced metric, curvature and observer cut on the family parameter implies that for **every sufficiently small nonzero amplitude**, at least one of the two nonnegative-weight averages of the actual null Ricci curvature is negative. The protocol giving the negative sign swaps when the amplitude sign changes. A negative average guarantees a local negative contraction somewhere in its support. Two orthogonal spatial axes at a critical zero of $Y$ give a pointwise version. A smooth tangent vector field gives another protocol using one direction at each observer. [Two point directions](NULL_DIRECTIONS.md) · [One direction per observer](ANGULAR_EXPOSURE.md) The proofs supply no numerical amplitude radius or instrumental sensitivity. [Observer regularity, §2](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_EFFECTIVE_OBSERVER_20261116/EFFECTIVE_OBSERVER_DOMAIN_AND_RESPONSE_BOUNDS.md)

## How large is the directional signal in the same curvature units?

The first variations of both readings have units of inverse length squared. For the plus-weighted null exposure and the quadrupole coefficient of the intrinsic Ricci trace, their ratio of means over one clock window is
$$
\frac{\overline{N'_+}}{\overline{(\delta R)_Y}}
=-\frac{8\ell^2}{81}\,
\frac{\overline{M\operatorname{sech}^2 q}}{\overline A}.
$$
Indeed $(\delta R)_Y=2p_bB_b S=6B_b A/(\rho_b^3\ell)$, $N'_+=-16\ell M\operatorname{sech}^2q/(27\mu\rho_b^4)$, and the endpoint normalization satisfies $\mu\rho_b B_b=1$. The exact enclosures give a ratio between **−9.677 and −9.602 times $10^{-16}$** in the first window, and between **1.59167 and 1.60653 times $10^{-14}$** in the second. The minus-weighted result has the opposite sign. The small factor arises from the registered packet width $\ell=10^{-6}$. The independently normalized plot makes the shapes visible and uses different vertical scales.

The reference in the denominator is the perturbation of the Ricci trace projected on $Y$, with the background removed. A common amplitude and common curvature unit cancel in this first-order comparison. The relative response of the two curvature measurements, null normalization, time windows and angular registration still require calibration. These figures neither specify actual instrumental sensitivity nor compare the signal with the entire background curvature. [Rational endpoints](NULL_SKY.json)

## What the sign excludes

In four-dimensional Einstein gravity with positive coupling, $R_{kk}=\kappa_4 T_{kk}$ for true null vectors even in the presence of a cosmological constant. Matter satisfying the pointwise null energy condition gives a nonnegative value to **each** of these observations, irrespective of the initial states or number of fields. Such a model cannot match the registered intrinsic metric or both of the null-curvature exposures at the stated small amplitudes. This is a conditional exclusion for the metric to which the observers couple.

In the five-dimensional model, the Gauss equation instead gives $R^{(4)}_{kk}=\tfrac23(k\Phi)^2-E_{kk}$, where $E$ is the normal projection of the bulk Weyl tensor. The bulk Einstein–scalar equations give $R^{(5)}_{KK}=(K\phi)^2\ge0$ for bulk-null $K$. At a shell event with negative intrinsic null Ricci contraction the projected Weyl contribution exceeds the nonnegative local scalar term. The induced-metric sign does not require negative bulk scalar energy. It is also distinct from the expansion of a light beam, which depends on initial data and shear. A four-dimensional theory with a different physical metric or effective stress can realize negative intrinsic null convergence. Earlier brane models already exhibited this general effect. [Gauss projection, Eq. B](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/expansion/GAUSS_IDENTITY_STANDARD_APPLICATION.md) · [Physical four-dimensional alternatives](COMPARISON.md) · [Relevant precedents](LITERATURE.md)
