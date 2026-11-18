# Four-dimensional gravity and the joint linear signal

Use signature $(-+++)$, physical shell radius $\rho_b$ and background metric $h_0=\rho_b^2\gamma$. Its Ricci scalar is $R_0=12/\rho_b^2$. Write $s=\delta\Phi$ for the measured scalar perturbation, $r=\delta R$ for the intrinsic Ricci scalar perturbation, and $K=2p_b B_b\ne0$. The registered shell predicts $r=Ks$. Both background scalars are constant, so an infinitesimal change of shell coordinates leaves both first variations invariant. [Graph, master and junctions, §§1–2](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_REGULAR_SHELL_CONTACT_20261013/junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md)

## Minimally coupled canonical matter

Take the physical metric as the Einstein metric. Use the four-dimensional field equations $G_{ab}+\Lambda g_{ab}=\kappa_4 T_{ab}$ with $\kappa_4>0$ and canonical scalar fields. Suppose all background fields are constant at a stationary point of their total potential. Absorb its constant value into $\Lambda=3/\rho_b^2$. The first derivative of every kinetic term vanishes and the potential has zero first variation. The trace equation then gives $r=0$, even if the scalar readout has an arbitrary nonzero first variation. For positive canonical kinetic terms with no other matter on an exact Einstein de Sitter background, stationarity follows directly: every null vector $k$ obeys $0=R_{ab}k^ak^b=\kappa_4\sum_j(k^a\partial_a\psi_j)^2$. Every background scalar gradient therefore annihilates all null vectors and vanishes. The source-free scalar equations then put the potential at a stationary point. The trace argument covers finite sets and convergent towers for which the summed stress is differentiable at this background. Free vacuum gravitational waves may have nonzero linear Weyl curvature, but obey the same Ricci trace. This exclusion uses Einstein gravity with stationary scalar matter and no further stress with a varying trace.

## A metric conformally coupled only to the measured scalar

Suppose an observer instead sees $h_J=e^{2a(\Phi-\phi_b)}g_E$ to first order, with an arbitrary constant $a$ and Einstein-frame $\delta R_E=0$. The conformal transformation of the four-dimensional Ricci scalar gives $r=-6a(\Box_{h_0}+4/\rho_b^2)s$. Demanding $r=Ks$ with $a\ne0$ requires the *observed scalar combination* to satisfy a constant-mass Klein–Gordon equation on both test intervals. A smooth function of the one readout has the same first-order restriction. Multiple species and arbitrary preparation do not change this argument.

A quantitative test uses the prior two nonnegative B-spline weights $\psi_A,\psi_B$. Write $s=s_1Y$, $r=r_1Y$, $x=q/\ell$, $C=\cosh^{3/2}(\ell x)$, $u=C\ell s_1/G$, $v=C\ell r_1/(GK)$. For each weight define the dimensionless measured ratio
$$
J_i=\frac{\int\psi_i v}{\int\psi_i u},\qquad
k_i=\frac{-\int\psi_i'' u-\ell^2\int k(x)\psi_i u}{\int\psi_i u},
\quad k(x)=\tfrac{35}{4}\operatorname{sech}^2(\ell x)-\tfrac94 .
$$
Here $q$ is the dimensionless elapsed de Sitter clock, $\ell=10^{-6}$, and $G$ is the archived conversion between scalar trace and normalized model amplitude. The denominators retain one sign on the two admitted test supports. Conjugation of $\Box_{h_0}$ in the quadrupole sector followed by integration by parts yields $J_i=t(k_i+4\ell^2)$ for this conformal comparator, with $t=-6a/(K\rho_b^2\ell^2)$. The shell instead has $J_A=J_B=1$. For $X=k_A+4\ell^2>Z=k_B+4\ell^2>0$, minimizing $\max(|tX-1|,|tZ-1|)$ over all real $t$ gives $(X-Z)/(X+Z)$, attained at $t=2/(X+Z)$. This expression increases with $X$ and decreases with $Z$.

The pinned interval certificate gives a discrepancy above **10.246%** in at least one ratio for an exact match of the scalar profile. If the competing scalar profile departs by at most 0.0005 uniformly in normalized $u$ units on both supports, the widened mass intervals still give a discrepancy above **9.192%**. Each ratio uses that model's own scalar denominator. These bounds concern deterministic ratios of two weighted measurements with a specified clock, relative gain and common angular template. Translating errors from individual noisy channels requires propagating their uncertainties through the ratios. [Weak tests and normalizations](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md) · [Rational check](conformal_test.py) · [Exact endpoints](CONFORMAL_TEST.json)

## A distinct conformal coupling reproduces the full induced metric

The previously constructed four-dimensional tower has Einstein-frame canonical fields $\psi_j$, physical masses $M_j^2=(\lambda_j-13/2)/\rho_b^2>0$, bounded coefficients $d_j=G\mu e_j(0)$ and first readout $s=\sum_j d_j\psi_j$. The eigenfunctions $e_j$ come from the finite causal collar. The registered fields are $\psi_j=(5/2-\lambda_j)\chi_j$ at first order, where $\chi_j=\cosh^{-3/2}(q)y_jY$ and $y_j$ are coefficients of the radial master. Prescribe a second, fixed coupling vector and let test clocks and rulers see the following physical metric:
$$
 h_J=e^{2B}g_E,\qquad B=\sum_j b_j\psi_j,\qquad
 b_j=\frac{e_j(0)}{\rho_b^2(5/2-\lambda_j)}
     =-\frac{K d_j}{6(M_j^2+4/\rho_b^2)} .
$$
The equality follows from $G=3/(\rho_b^3p_b)$ and $\mu=1/(\rho_b B_b)$. Both coefficient vectors are square-summable. The local Einstein-frame action consists of the Einstein–Hilbert term with $\Lambda$ and independent canonical massive scalars. Test matter couples to $h_J$ and its stress is omitted from the background and field calculation. The observer's scalar reading is $\phi_b+\sum d_j\psi_j$. This specifies the metric used by the observer. Merely relabeling the Einstein metric would leave the physical measurement unspecified.

Take the homogeneous linear Einstein-frame metric perturbation to be zero, which solves its field equation. The field equations then give $r=-6\sum_j(M_j^2+4/\rho_b^2)b_j\psi_j=Ks$. More strongly,
$$
 B=\rho_b^{-2}\sum_j e_j(0)\chi_j=\rho_b^{-2}P_b .
$$
The radial master has $P=\alpha\cosh^{-3/2}(q)wY$ and $\alpha(0)=1$ at the shell. Endpoint evaluation of the spectral expansion gives the last equality. For the selected even graph the full shell metric variation is $2(P_b/\rho_b^2)h_0$. Thus the observer's metric matches its complete first variation in the reconstruction chart, along with the scalar, using the same Cauchy data and a single fixed pair of readout vectors. Common pullback identifies the corresponding geometric first variations. The smooth compact packet belongs to graph domains of every order, which permits these sums and derivatives on the admitted compact clock interval. This is a classical four-dimensional comparison with deliberately chosen couplings. A physical theory of those couplings or of preparation probabilities is additional input. [Masses, data, convergence and scalar readout](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/TOWER.md) · [Shell normalization, §1](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/INERTIAL_SHELL_WAVEFORM.md)
