# How changing a four-dimensional readout affects the quadratic test

The shell obeys the exact Gauss trace constraint $R(h)=F(\Phi)+|D\Phi|_h^2$. The specified four-dimensional model in the [comparison](COMPARISON.md) reproduces the entire linear induced metric and scalar. Its linear conformal exponent and linear scalar readout leave a strictly positive quadratic discrepancy at a quiet cut. How dependent is that conclusion on the prescribed detector and scalar couplings? The following construction makes that dependence explicit.

Take finitely many canonical fields $\psi_i$ on an Einstein-frame metric $g$ with positive gravitational coupling $\kappa$, cosmological constant $\Lambda$, and quadratic potential $\tfrac12\sum_i M_i^2\psi_i^2$. The unperturbed fields vanish and $R_0=4\Lambda$. Set the physical metric to $h=e^{2B(\psi)}g$ and the reported scalar to $\varphi(\psi)$, with $B(0)=0$ and $\varphi(0)=\phi_b$. Write $b_i=B_{,i}(0)$ and $d_i=\varphi_{,i}(0)$. The Einstein trace, the convention $\Box_g\psi_i=M_i^2\psi_i$, and the four-dimensional conformal curvature formula give the exact identity
$$
R(h)-|D\varphi|_h^2-F(\varphi)
=P(\psi)-F(\varphi)
+e^{-2B}\sum_{ij}
\big[\kappa\delta_{ij}-6B_{,ij}-6B_{,i}B_{,j}-\varphi_{,i}\varphi_{,j}\big]\,
g^{ab}\partial_a\psi_i\partial_b\psi_j ,
$$
where
$$
P(\psi)=e^{-2B}\left[R_0+2\kappa\sum_i M_i^2\psi_i^2-6\sum_i B_{,i} M_i^2\psi_i\right].
$$
These formulas also display the effect of intrinsic metric perturbations through the contractions. No fixed metric or free choice of the Einstein trace after specifying the fields has been assumed.

For couplings satisfying $-(2R_0+6M_i^2)b_i=F'_b d_i$, specify the Hessian at the origin by
$$
 B_{,ij}(0)=[\kappa\delta_{ij}-6b_i b_j-d_i d_j]/6 .
$$
One concrete choice is the quadratic polynomial with this Hessian and the given first derivative. The shell coefficient $F'_b$ is nonzero. The inverse function theorem thus supplies a local operational scalar $\varphi=(F|_{\rm local})^{-1}(P(\psi))$. It has $\varphi(0)=\phi_b$ and $\varphi_{,i}(0)=d_i$. The field-value term in the identity vanishes by construction, and the bracket multiplying each gradient product is zero at the background. For smooth fields of amplitude $O(\epsilon)$ and first derivatives $O(\epsilon)$, the residual is $O(\epsilon^3)$ on a compact patch where the functions are smooth. Thus the *same* linear couplings permit a four-dimensional theory and readout satisfying the trace constraint through quadratic order for all smooth finite-field preparations. A two-field rational example checks both mixed and diagonal terms. [Numerical algebra check](second_order_freedom.py) · [Exact output](FREEDOM.json)

The finite-field construction respects the prescribed linear readouts on each truncation. The earlier countable tower attains exact agreement at linear order on the full interval. Its quadratic forms here can be interpreted termwise on smooth spectrally decaying data, when all mass-weighted sums converge. A nonlinear infinite-field solution, a fit to the individual second-order metric and scalar waveforms, and the ultraviolet physics of a detector implementing these couplings require additional work. Reproducing one joint consistency relation is a weaker requirement than predicting the two observables separately. The source of these cancellations is an expressly chosen nonlinear readout. Whether a physical observer has such a coupling has to be fixed independently.
