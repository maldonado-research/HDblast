# Curvature of the registered shell on its observer clock

On the timelike shell the spacelike unit normal, pure tension, and scalar boundary row are $K_{ab}=\sigma h_{ab}/6$ and $n\phi=-\sigma'/2$. Trace reversal of the five-dimensional action gives $R^{(5)}_{AB}=\partial_A\phi\partial_B\phi+(2/3)U g_{AB}$. Contracting the Gauss equation yields the exact intrinsic four-dimensional identity
$$
R(h)= |D\Phi|_h^2+F(\Phi),\qquad F=2U-\sigma'^2/4+\sigma^2/3 .
$$
Here $\Phi$ is the actual pulled-back scalar, $D$ the spacetime derivative for $h$, and the norm is Lorentzian. The sign of the normal drops out. At the regular background $F(\phi_b)=12/\rho_b^2$, $p_b=-\sigma'_b/2$ and $a_b=\sigma_b/6$. The static scalar equation gives $B_b=U'_b/p_b-4a_b+\sigma''_b/2$. Consequently $F'_b=2p_b B_b$. Both $R(h_0)$ and $\Phi_0$ are constant on the entire observer cylinder, which removes first-order label and clock displacement terms. The first curvature perturbation is therefore $F'_b S(-q)Y$ on the same proper-clock observer labels as the scalar. [Physical shell and first scalar coefficient](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_SHELL_RESPONSE_20261118/BRIDGE.md) · [Background and junction conventions](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/inputs/previous/archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md)

The curve uses the registered $U=W'^2/2-2W^2/3$, $\sigma=2W+\delta s$, $s=1+c\phi$, $\delta=1/1000$ and $W=1/3+(\phi-1)^2+(\phi-1)^3/3$. Algebraically,
$$
F=\delta\big[(4/3)Ws-cW'\big]+\delta^2\big[s^2/3-c^2/4\big] .
$$
This expression cancels the two undetuned contributions before interval evaluation. Using the outward rational shell box,
$$
-0.000929224<F'_b<-0.000929203,\qquad
0.00306451<F''_b<0.00306458 .
$$
Together with the 512 continuous time cells in $W_1=[5/4,11/8]$ and the 128 in $W_2=[49/32,25/16]$, the first curvature quadrupole's *mean coefficient* lies in $(37.1998,41.1312)$ on $W_1$ and $(-15.0782,-13.5978)$ on $W_2$. These values multiply the dimensionless perturbation amplitude. They use calibrated model curvatures; they are not astronomical samples.

## A test that removes the unknown second perturbation

At zero background gradient, $\epsilon^{-2}(R-F(\Phi))$ converges to $|D(SY)|_{h_0}^2$. This requires first parameter differentiability of the scalar into a space controlling its first spacetime and angular derivatives, and continuous metrics and cut densities. In proper observer coordinates $q$, $h_0=\rho_b^2[-dq^2+\cosh^2 q\,\kappa]$. For $u=n_1^2+n_2^2$ one has $|\nabla_\kappa Y|^2=4(u-Y^2)$. Take the transported signed weight $w=3-4u$. Uniform $S^3$ measure gives $u$ uniform on $[0,1]$ and $Y=u\cos(2\theta)$ with a uniform angle independent of $u$. Thus $E_0[wY^2]=3/6-4/8=0$ and $E_0[wu]=3/2-4/3=1/6$. The time derivative, however large, drops out. The normalized proper-area expectation $E_\epsilon$ gives
$$
Q(q)=\lim_{\epsilon\to0}\epsilon^{-2} E_\epsilon[w(R-F(\Phi))]
=\frac{2S(-q)^2}{3\rho_b^2\cosh^2 q}.
$$
The density derivative carries a zero background contrast. Because the integrand starts at order two, changes of density, clock registration and observer labels at first order do not change $Q$. Spatial derivatives are evaluated on the observer's smooth shell spacetime. This identity does not turn the Lorentzian norm into a positive pointwise norm. The positive coefficient comes specifically from its signed angular projection. [Earlier variance calculation](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_SHELL_RESPONSE_20261118/VARIANCE.md)

The earlier leading proper-area variance is $V_2=S^2/6$. Hence $Q=4 V_2/(\rho_b^2\cosh^2 q)$, with the factor in $(0.00023695958,0.00023696162)$ throughout the admitted clock interval. This is an identity even at a zero of $S$. A normalized version uses the measured cut volume radius $r_\epsilon=(\operatorname{Vol}(\text{cut})/(2\pi^2))^{1/3}$, whose background value is $\rho_b\cosh q$. Uniformly on the compact interval,
$$
r_\epsilon^2 E_\epsilon[w(R-F(\Phi))]-4\operatorname{Var}_\epsilon(\Phi)=o(\epsilon^2).
$$
The ratio of the window averages of the first and second observables (before multiplying variance by four) tends to four whenever the window has positive leading variance. This avoids dividing instantaneous measurements close to a zero.

| Finite window in $x=q/\ell$ | Mean $Q$, outward bounds in model units |
| --- | --- |
| $[5/4,11/8]$ | $(63455.81,77574.04)$ |
| $[49/32,25/16]$ | $(8750.36,10755.12)$ |
| Entire admitted $[0,25/16]$ | $(8382.20,10259.05)$ |

The replay integrates continuous closed cell bounds for the physical scalar waveform, transfers the two time orientations with absolute error $1/40$, and uses $1-q^2\le\operatorname{sech}^2q\le1$. The output gives exact rational endpoints. The coefficient of $\epsilon^2$ differs by a factor of two from a second derivative. The finite-amplitude convergence does not supply a numerical allowed radius in $\epsilon$. [Source of the shell box](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_SHELL_RESPONSE_ERROR_CERTIFICATE_20261015/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json) · [Reproduction](curvature.py) · [Exact results](CURVATURE.json)
