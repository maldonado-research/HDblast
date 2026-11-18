# Distance at fixed observed redshift for freely falling emitters

Choose receivers and emitters from the same shell congruence. Each worldline is a timelike geodesic of the intrinsic shell metric, with initial tangent normal to the originally specified quiet cut. The receiving observer has initial label $n_*=(0,0,1,0)$. Reception is at *fixed physical past proper duration* $\rho_b q_0$ from that cut, with $q_0=49\ell/32$ and $\ell=10^{-6}$. Define $q_1=99\ell/64$, $s_i=\sinh q_i$, $c_i=\cosh q_i$, and specify the observed redshift by
$$
 1+z_*=c_0/c_1 .
$$
This is a small blueshift on the chosen contracting portion of the background. The comparison fixes the received direction and this physical redshift. It selects emitting events along a continuous distribution of freely falling clocks. It does not track a galaxy of fixed label or compare emission at a prescribed proper time.

## The emitter motion needs to be included

On the material shell the physical first metric derivative is $2Fh_0$ up to a tangential Lie derivative, where $F=U(q)Y(n)$, $Y=n_1^2-n_2^2$, and
$$
 h_0=\rho_b^2[-dq^2+\cosh^2q\,\kappa],\qquad M(x)=\mu\rho_b^2 U(\ell x)/\ell .
$$
Here $q=-t$ is the **conformal-representative clock coordinate**. It equals the dimensionless observer proper clock at zero perturbation, and along the critical receiver through first order. For other perturbed emitters these two clocks generally differ. Established diffeomorphism covariance permits working with the representative $e^{2\epsilon F}h_0$ for first derivatives. The chart change has vanishing initial data and first time derivative, and the original shell packet starts strictly away from the material cut. Consequently the initially normal physical congruence and its label axes have the same initial data in this representative. The full metric at finite perturbation need not be conformal to the background. [Metric tangent and exact relation to the scalar](https://github.com/maldonado-research/HDblast/blob/72257d4ace4a37a241adf721f9f78494f2db89bc/research/HDBLAST_JOINT_CURVATURE_20261118/METRIC.md) · [Gauge data, §2](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_FIRST_GAUGE_WAVE_20261016/operator/GAUGE_OPERATOR_AND_CONSTRUCTION.md)

Write $c(q)=\cosh q$ and $I(q)=\int_0^q U(y)\,dy$. For convenience orient the timelike emitter tangent and photon tangent both toward the past and take their positive frequency to be minus their metric contraction. This produces the same emission-to-reception frequency ratio as reversing both. If $dn^i/dq=\epsilon V^i+o(\epsilon)$, then the geodesic momentum equation to first order gives
$$
 \partial_q(c^2 V_i)=-D_iF,\qquad V_i(0)=0,qquad
 V_i=-\frac{I(q)}{c^2}D_iY .
$$
In fact the spatial covariant momentum per unit mass is $\epsilon\rho_b c^2 V_i$ and its derivative in past proper duration is $-\epsilon D_iF$. The spatial connection quadratic in the peculiar velocity has no linear term. Here indices on $V$ and $DY$ use $\kappa$. In conformal time $d\eta=dq/c$ the emitter's peculiar velocity in the past frame, relative to comoving cylinder coordinates, is $\epsilon\boldsymbol\beta$, with $\boldsymbol\beta=-I(q)DY/c(q)$. Evaluating this first-order field at the unperturbed emission event also gives the Eulerian velocity there since the background peculiar velocity is zero.

At $n_*$ one has $Y=DY=0$. Both the induced first angular velocity and the metric lapse variation vanish along the unperturbed receiver's entire past history. Thus its fixed-proper-duration reception event has no first coordinate displacement. A rotation of the initial sky basis would leave the isotropic background redshift and distance unchanged; the oriented basis at zero amplitude fixes the coefficient below.

## Propagation, frequency, and the selected endpoint

Let $v=(v_1,v_2,v_4)$ be a unit vertex sky direction along the label axes $n_1,n_2,n_4$, and put $d_v=v_1^2-v_2^2$. The photon is traced backwards from the receiving observer. In the conformal cylinder its angular path is
$$
 n(q)=n_*\cos\theta+v\sin\theta,\qquad
 \theta=\arctan(\sinh q)-\arctan(s_0).
$$
The conformal perturbation preserves this path as an unparametrized null curve to linear order. At the emission point its unit angular propagation direction is $m=dn/d\theta$. Along this ray $Y=d_v\sin^2\theta$ and $D_mY=2d_v\sin\theta\cos\theta$.

For the past photon normalized to unit frequency at reception, the physical frequency measured in the moving emitter frame at fixed endpoint coordinate $q$ is proportional, to first order, to $e^{-\epsilon F}(1-\epsilon\boldsymbol\beta\cdot m)/(\rho_b c)$. This follows by evaluating the timelike emitter tangent against the conformally rescaled null tangent in $e^{2\epsilon F}\rho_b^2c^2(-d\eta^2+\kappa)$. The cylinder null tangent rescales inversely with the square of the conformal factor and the emitter's time component inversely with its first power. The angular direction can be taken at zero amplitude in the first-order Doppler term. The reception normalization has $F=\boldsymbol\beta=0$. Consequently
$$
 \frac{(\partial_\epsilon z)_{q,0}}{1+z_0(q)}
   =d_v B(q),\qquad
 B(q)=-U(q)\sin^2\theta+
       \frac{2I(q)}{c(q)}\sin\theta\cos\theta,
 \qquad 1+z_0(q)=\frac{c_0}{c(q)} .
$$
The term involving $I$ is the freely falling emitter's accumulated velocity. Declaring an emitter to remain at rest in the conformal coordinates would remove this term and change the source specification.

Define $\Delta=s_1-s_0>0$ and $L(q)=\rho_b(\sinh q-s_0)/c_0$. At zero amplitude the area distance to the cut of coordinate time $q$ is $\rho_b c(q)\sin\theta=L(q)$. At fixed $q$ its first fractional change is $F(q,n(q))=d_vU(q)\sin^2\theta$. The unparametrized rays and observed angles are unchanged at this order, and the physical transverse length at emission scales with the metric conformal factor. Changing from a constant-$q$ to a direction-dependent radial cut contributes to transverse area only through the nominal distance evaluated at the new endpoint: the extra angular Jacobi vector along the ray is null and orthogonal to the vertex cone. Projecting source area onto a screen orthogonal to the emitter also preserves this quotient-screen area.

Since $z_0'(q_1)=-(1+z_*)\tanh q_1\ne0$, the required endpoint shift in coordinate $q$ at fixed observed redshift is
$$
 (\partial_\epsilon q_e)_{z_*,0}=d_v B(q_1)/\tanh q_1 .
$$
The desired response is therefore
$$
 \frac{(\partial_\epsilon D_A)_{z_*,0}}{L(q_1)}
 =d_v\left[U(q_1)\sin^2\theta_1+
       \frac{c_1^2}{\Delta s_1} B(q_1)\right]
 =d_v\frac{1+s_0s_1}{c_0^2c_1^2s_1}
       \big[2c_1 I(q_1)-\Delta U(q_1)\big].
$$
For the simplification use $\sin\theta_1=\Delta/(c_0c_1)$, $\cos\theta_1=(1+s_0s_1)/(c_0c_1)$ and $c_1^2-s_1^2=1$. For a check using the affine calculation, setting $\lambda=L(q_1)$ instead shifts $q_e$ by $-2\int_0^{L(q_1)}F\,d\lambda/L'(q_1)$ at this receiver. Combining the affine result with the physical redshift endpoint term cancels that integral and restores the expression above.

The formula for the $v_1\ge4/5$ cap has the prefactor $d_v$ replaced by $18/25$. For $v_2\ge4/5$ the prefactor is $-18/25$. This follows by differentiating the square root of cap area divided by vertex solid angle at its isotropic background value. After removing that angular factor the axial coefficient is
$$
 \frac{\ell^2}{\mu\rho_b^2}\frac{1+s_0s_1}{c_0^2c_1^2s_1}
 \left[2c_1\int_0^{99/64} M(x)\,dx-\frac{\Delta}{\ell}M(99/64)\right].
$$
Both the endpoint and the history of the angular curvature coefficient are required. The history of the receiver and the emitter starts on the original quiet cut, rather than at first photon reception. All derivatives use the chosen signed perturbation parameter. The shift and its sign need independent source and velocity information for any other tracer population.

## Physical scope of the endpoint selection

The emitting congruence fills a neighborhood of this background ray segment. The observer map is a small smooth deformation of the embedded background proper-clock cylinder, with nonvanishing cut Jacobian and a timelike transverse tangent. On a compact interior neighborhood of the target endpoints it has a smooth local inverse for sufficiently small parameters. The full sky can be covered by finitely many such charts with a common amplitude neighborhood. The photon geodesic, the redshift as the ratio of two invariant frequency contractions, and the selected point then vary differentiably. The implicit function theorem applies uniformly to the local branch near $q_1$. The small nonzero redshift derivative makes the problem sensitive, but does not invalidate local existence. The endpoint has a margin of $\ell/64$ to the upper observer clock boundary and stays in an admissible neighborhood for sufficiently small amplitude. A numerical amplitude radius has not been obtained. A globally unique astrophysical redshift-distance map is not required for this local branch. [Registered observer flow and regularity, §2](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_EFFECTIVE_OBSERVER_20261116/EFFECTIVE_OBSERVER_DOMAIN_AND_RESPONSE_BOUNDS.md)

Positive or negative distance responses at fixed redshift do not by themselves test the null energy condition. The affine inequality $D_A\le\lambda$ in Einstein gravity with null-energy-respecting matter needs the actual affine length along each ray and a caustic-free pencil. At fixed redshift that length varies with the emitting congruence and metric. The displayed endpoint calculation concerns clocks specified by the local shell model. The continuum calculation extends it to a bounded interval of blueshifts. [Continuum and local limit](LOCAL_LIMIT.md) It does not identify emitters, a received astronomical source, or the global expansion history.
