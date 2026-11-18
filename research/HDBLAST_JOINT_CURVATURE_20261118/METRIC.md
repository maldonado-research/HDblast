# An intrinsic curvature observable on the same shell as the scalar

The registered local Einstein–scalar solution predicts an anisotropic shell Ricci tensor in the same two clock windows as its scalar quadrupole. The two exposure means of the measured curvature coefficient have the same sign even though the scalar means have opposite signs. Both follow from the same initial packet. At one initial sphere label the curvature additionally gives opposite signs of null convergence in two physical directions. [Null contractions and four-dimensional interpretation](NULL_DIRECTIONS.md)

## Induced metric on the material shell

Let $h$ denote the shell's intrinsic four-dimensional Lorentz metric, $q=-t$ the dimensionless elapsed observer time at zero amplitude, and $x=q/\ell$, $\ell=10^{-6}$. The background metric is $h_0=\rho_b^2[-dq^2+\cosh^2q\,\kappa]$, with $\kappa$ the unit three-sphere metric. The registered copy-even perturbation is a straight material shell in longitudinal gauge. Its metric derivative is $2U(q)Yh_0$ with $Y=n_1^2-n_2^2$. Indeed the longitudinal reconstruction has $h^{(1)}_{ab}=2\Psi h_{0,ab}$, $\Psi=P/\rho^2$ and $P=\alpha\cosh^{-3/2}t\,wY$ with $\alpha_b=1$. The graph correction $a_b\zeta$ in the induced conformal factor vanishes in this representative. The relative-harmonic bridge has zero normal displacement at the shell and introduces a tangential Lie derivative. The actual evolved family has this first jet by the registered linear uniqueness argument. Its intrinsic normal observers and cuts have a common parameter domain. [Induced metric and embedding, equations (3)–(4)](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md) · [Master conversion, §1](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/INERTIAL_SHELL_WAVEFORM.md) · [Actual metric and scalar tangent](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_SHELL_RESPONSE_20261118/BRIDGE.md)

Set $F=UY$, $C(x)=\cosh^{3/2}(\ell x)$. Then $U(\ell x)=w(-\ell x,0)/(\rho_b^2 C(x))$. For $J_{ab}=\operatorname{Ric}_{ab}[h]-R[h]h_{ab}/4$ the background tensor vanishes identically. A tangential coordinate shift, displacement of the observer event, or first variation of its frame multiplies this zero background or its zero covariant derivative. The first conformal variation in four dimensions gives
$$
 J'_{ab}=-2\left(\nabla_a\nabla_b F-\tfrac14 h_{0,ab}\Box_{h_0}F\right).
$$
Here a prime on curvature denotes an amplitude derivative at zero. This follows from $\delta\operatorname{Ric}_{ab}=-2\nabla_a\nabla_bF-h_{0,ab}\Box F$, $\delta R=-2R_0F-6\Box F$ and $R_0=12/\rho_b^2$. Computing with the vanishing-background tensor automatically accounts for the trace term and the variation of null or orthonormal frames.

Project $J$ to a proper observer cut and take its spatial trace-free part $\Pi_{ij}$ with that cut's metric. At zero amplitude $\Pi=0$. With $D$ the sphere derivative write $Z_{ij}=D_iD_jY+(8/3)\kappa_{ij}Y$. Since $\Delta_\kappa Y=-8Y$ and the spatial-to-temporal Christoffel term in the Hessian is isotropic,
$$
 \Pi'_{ij}=-2U Z_{ij}.
$$
Indices here are the original transported sphere labels in a background coordinate frame, not unit physical components. For a full cut, contract with $Z$ using the background $\kappa$ and divide the resulting mean by $-2\langle |Z|^2\rangle$ to recover $U$. The sphere Bochner identity gives $\langle |Z|^2\rangle=[(2/3)8^2-2\cdot8]\langle Y^2\rangle=40/9$. Physical unit components instead carry the expected inverse factor $\rho_b^2\cosh^2 q$.

There are other curvature readouts with different information. The intrinsic four-dimensional Weyl tensor is zero at first order because $h'$ is conformal plus a diffeomorphism of de Sitter. The trace of the intrinsic Ricci tensor obeys $\delta R=2p_bB_b S(-q)Y$. To see this, differentiate the prior exact Gauss identity $R=(D\Phi)^2-\sigma'^2/4+2V+\sigma^2/3$ at constant background $\Phi$, using $\sigma=6a$, $\sigma'=-2p$, $V'=p_y+4ap$, $B=p_y/p+\sigma''/2$. The projected traceless tensor is different from that algebraic trace. With the same conventions, the spatial trace-free part of the projected Einstein equation gives $\Pi'_{ij}=-(E^{(5)}_{ij})'_{\rm spatial\,TF}$, where $E^{(5)}_{ab}=C^{(5)}(n,e_a,n,e_b)$ uses the shell's spacelike normal. All tangential scalar-gradient products vanish to first order. [Exact Gauss relations (B)–(C), definitions and background](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_INTRINSIC_SHELL_BUDGET_20261101/expansion/GAUSS_IDENTITY_STANDARD_APPLICATION.md) · [Actual field equations and geometric boundary rows](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/PHYSICAL_RECONSTRUCTION_AND_BOUNDARY_JOIN.md)

## Causal response of the geometry to the registered scalar

Use the physical normalized scalar $A(x)=(\ell/G)S(-\ell x)$, where $G=3/(\rho_b^3p_b)$, and $u=C A$. Define
$$
 z(x)=\frac{\mu}{\ell} w(-\ell x,0),\qquad
 M(x)=\frac{\mu\rho_b^2}{\ell}U(\ell x)=z(x)/C(x).
$$
The dynamic endpoint law gives $w_z-\gamma w=\mu[w_{qq}+(5/2-(35/4)\tanh^2q)w]$. The registered scalar equals $G(w_z-\gamma w)/C$, so the exact linear boundary relation is
$$
 z''+\ell^2[5/2-(35/4)\tanh^2(\ell x)]z=u(x).
$$
Primes here denote $x$ derivatives. Both endpoint Cauchy data are zero. Finite propagation at the initially separated boundary gives $z=u=0$ through $x=1$. The pair $A,M$ are different measured quantities constrained by this causal second-order relation. Thus independently measuring the geometry tests a dynamical prediction of the model. Its two signals are not statistically independent random fields. [Endpoint equation and exact quiet interval, §§1–4](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Scalar window normalization](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md)

For $X=25/16$ let $z_0(x)=\int_1^x(x-s)u(s)\,ds$ after the front and $z_0=0$ earlier. Put $d=X-1$, $\vartheta=(5/2)\ell^2 d^2/2$ and $B=\|u\|_\infty d^2/[2(1-\vartheta)]$, using the cell enclosure for the supremum. The coefficient in brackets is positive and at most $5/2$ here. Volterra's identity gives $\|z\|_\infty\le B$ and $\|z-z_0\|_\infty\le\vartheta B$. With $1-10^{-11}<C^{-1}\le1$, it follows that $\|M-z_0\|_\infty<7.0\times10^{-12}$. A separate deterministic sensitivity check keeps the endpoint law, known zero data and all normalizations fixed. An additional absolute input error $1/2000$ in $A$ throughout $[0,X]$, even before the front, propagates to an absolute error below $0.000611$ in $M$. Apply the same bound on the entire interval with length $X$ and conjugated input error $(1/2000)/(1-10^{-11})$. This is a sensitivity to the assumed scalar input, not an estimated experimental error in a curvature measurement.

For a normalized box on $[a,b]$ the kernel against $u(s)$ is $[(b-s)_+^2-(a-s)_+^2]/[2(b-a)]$. The replay integrates it exactly on each closed cell. At a fixed $x$ its primitive kernel $(x-s)_+$ is nonnegative and increases with $x$. Testing separately with the two endpoints bounds $z_0(x)$ throughout each window with arbitrary input signs; the lower and upper endpoints of the input cell intervals are multiplied with the correct signed bounds. The cells authenticate their grid, gain and past-transfer source. The displayed bounds round outward. [Rational replay](metric.py) · [Exact rational outputs](METRIC.json)

| Elapsed window $x$ | Mean of $M$ | Enclosure valid for every $M(x)$ in the window |
| --- | --- | --- |
| $[5/4,11/8]$ | $(0.03912,0.03932)$ | $(0.01214,0.07710)$ |
| $[49/32,25/16]$ | $(0.23770,0.23863)$ | $(0.22142,0.25472)$ |

The late-to-early ratio of the two metric means is strictly between $6.0456$ and $6.0996$. The previously enclosed scalar ratio is strictly between $-0.36754$ and $-0.36459$. Each ratio is invariant under a fixed nonzero gain on its own channel. The different signs preclude a constant proportional readout of the scalar as the measured geometry. These are first coefficients and extrapolations in the model’s clock and units. [Previous scalar means and assumptions](https://doi.org/10.5281/zenodo.24351909)
