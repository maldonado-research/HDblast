# Two second-order measurements from intrinsic curvature and the shell scalar

The four-dimensional Ricci scalar of the induced shell metric and the value of the shell scalar field form a local consistency test. In the specified five-dimensional Einstein–scalar model, subtracting a known function of the scalar leaves its **signed spacetime gradient squared**. Two angular weightings of this subtraction isolate, at second order, how fast the scalar changes around the sphere and along the clock. They need neither a numerical derivative of the observed signal nor the unknown second-order metric and scalar perturbations.

## The exact junction constraint

Write $h$ for the Lorentzian four-dimensional shell metric, $n$ for one outward spacelike unit normal and $\Phi$ for the restriction of the bulk scalar. The spacelike-normal Gauss equation gives $R(h)=R_5-2\operatorname{Ric}_5(n,n)+(\operatorname{tr}K)^2-K_{ab}K^{ab}$. Substitute $\operatorname{Ric}_5=d\phi\otimes d\phi+(2/3)Vg$, $K_{ab}=\sigma h_{ab}/6$ and $n(\phi)=-\sigma'/2$ to obtain

$$
J=R(h)-F(\Phi)=h^{ab}D_a\Phi D_b\Phi,\qquad
F=2V-\frac{(\sigma')^2}{4}+\frac{\sigma^2}{3}.
$$

These quantities are evaluated on the physical timelike shell. The scalar curvature of a single three-dimensional observer cut is a different observable. In the background $F(\phi_b)=12/\rho_b^2$ by the radial Hamiltonian constraint. Since $p_b=-\sigma'_b/2$, $a_b=\sigma_b/6$, and $B_b=V'_b/p_b-4a_b+\sigma''_b/2$, one finds $F'(\phi_b)=2p_b B_b<0$. In particular the joint linear variation of $R(h)$ and $\Phi$ on every admitted observer is $R_1=F'(\phi_b)\Phi_1$. [Model and physical junctions, §1](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/inputs/REGISTERED_EXACT_FAMILY_BINDING.md) · [Background and strict sign, §1](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_SHELL_RESPONSE_ERROR_CERTIFICATE_20261015/archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md)

## The proper-clock and angular derivatives

Use coordinates supplied by the intrinsic observers on the admitted slab. At amplitude zero $h_0=\rho_b^2(-dq^2+\cosh^2q\,\kappa)$ with the unit-sphere metric $\kappa$, and $\Phi_0=\phi_b$ is constant. Put $s(q)=S(-q)$ and $Y=n_1^2-n_2^2$. The registered first derivative is $\Phi_1=s(q)Y$. Because $d\Phi_0=0$, both $J_0$ and $J'_0$ vanish. Taking two amplitude derivatives of the exact constraint gives

$$
J''_0=\frac{2}{\rho_b^2}\left[-s'(q)^2Y^2+s(q)^2\operatorname{sech}^2q\,|\nabla_\kappa Y|^2\right].
$$

Here $s'$ denotes a derivative in $q$. The factor of two comes from the convention that $J''_0$ is the second derivative, so the coefficient of $\epsilon^2$ is $J''_0/2$. The previous parameter and observer regularity is strong enough to differentiate the scalar with respect to $q$ and the sphere labels. Equivalently, the exact constraint transfers the required regularity of the field and observer metric to this combination of curvature and scalar. [Observer and field identification](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_SHELL_RESPONSE_20261118/BRIDGE.md) · [Common proper-clock domain and regularity](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_EFFECTIVE_OBSERVER_20261116/EFFECTIVE_OBSERVER_DOMAIN_AND_RESPONSE_BOUNDS.md)

A pointwise check follows at the transported label $(1,1,0,0)/\sqrt{2}$ where $Y=0$ and $|\nabla Y|^2=4$. There $J''_0=8s^2\operatorname{sech}^2q/\rho_b^2$. This label lies on a regular part of the nodal set. At the angular extremum $(1,0,0,0)$ the same quantity is $-2s'^2/\rho_b^2$. A whole-sphere implementation avoids isolating either point.

Let $E_\epsilon$ be the normalized actual area average on a cut, and transport the following weights with the observer labels:
$$
w=1-\frac{20}{9}Y^2,\qquad v=\frac32w-1.
$$
As $J_0=J'_0=0$, $(E_\epsilon[wJ])''_0=E_0[wJ''_0]$, and likewise for $v$. Write $U=n_1^2+n_2^2$. For uniform area on the unit three-sphere, $U$ is uniform on $[0,1]$ and $Y=U\cos(2\theta)$ with uniform independent angle $\theta$. This follows from the angular coordinates $(\sqrt{U}\cos\theta,\sqrt{U}\sin\theta,\sqrt{1-U}\cos\varphi,\sqrt{1-U}\sin\varphi)$ whose area element is $dU\,d\theta\,d\varphi/2$. Thus $E_0[Y^2]=1/6$, $E_0[Y^4]=3/40$, $E_0[UY^2]=1/8$, and $|\nabla_\kappa Y|^2=4(U-Y^2)$. In particular $E_0[wY^2]=0$ and $E_0[w|\nabla Y|^2]=8/9$. It follows that

$$
P_s(q)=(E_\epsilon[wJ])''_0=\frac{16}{9\rho_b^2}\operatorname{sech}^2q\,s^2,\qquad
P_t(q)=(E_\epsilon[vJ])''_0=\frac{s'^2}{3\rho_b^2}.
$$

Both are nonnegative. The unweighted second derivative is $(E_\epsilon[J])''_0=3P_s/2-P_t$. The angular weights are signed and require the same frame and full-sphere area convention as the registered quadrupole.

## One common finite exposure

Average these quantities over $W=[5/4,11/8]$ in $x=q/\ell$. The linear relation between $x$, $q$ and background-normalized observer proper duration makes these normalized means identical. Use $A=(\ell/G)s$ for the scalar coefficient in the archived normalization. The 512 original continuous cells, their past-transfer correction and the positive dilution factor give
$$
16.1465<\overline{A^2}<16.2298 .
$$
For $k=0,\ldots,16$ take the equally spaced points $x_k=5/4+k/128$. Intersect the two closed cell bounds at each endpoint. If $g_k$ is the distance between adjacent endpoint intervals (zero if they overlap), Cauchy–Schwarz gives $\int_{x_k}^{x_{k+1}}A_x'^2\,dx\ge g_k^2/(x_{k+1}-x_k)$. Summing proves $\overline{A_x'^2}>111.40$. The endpoint intervals constrain the integrals of squared derivatives even though no derivative samples are supplied.

In the model normalization the exact arithmetic gives
$$
166{,}449<\overline{P_s}<209{,}840,\qquad
\overline{P_t}>2.1532\,10^{17},\qquad
\frac{\overline{P_t}}{\overline{P_s}}>1.2869\,10^{12}.
$$
The ratio uses $(3/16)\ell^{-2}\overline{A_x'^2}/\overline{\operatorname{sech}^2(\ell x)A^2}$ so the constant gain and $\rho_b$ cancel. The unweighted mean satisfies $\overline{(E_\epsilon[J])''_0}<-2.1532\,10^{17}$. This is a statement about a time-and-area average of a Lorentzian norm. Its sign may differ from the sign at individual angular positions. The short dimensionless packet duration $\ell=10^{-6}$ enters through the squared time derivative. These numbers have no independent conversion to observed curvature units. [Arithmetic check](powers.py) · [Exact bounds](POWERS.json)

Since the linear contributions vanish and the second-derivative margins are strict, the corresponding time-and-area averages of $wJ$ and $vJ$ are positive and the average of $J$ is negative for all sufficiently small nonzero amplitudes of the registered family. The admitted regularity supplies a uniform Taylor remainder on this compact slab. Its amplitude radius has not been evaluated. The upper and lower numerical magnitudes above refer to the derivatives at zero. Applying this signal comparison would require measurement of both four-dimensional intrinsic curvature and the shell scalar, a fixed normalization and properly matched spacetime exposures.
