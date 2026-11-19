# Initial scalar stress and force at fixed occupation

The initial scalar energy, pressure and force can agree for two states that subsequently produce different fermion counts. The physical mode basis fixes which part of the initial correlation enters each source. The criterion below applies to one minimally coupled real scalar with action density $-[(\partial\chi)^2+m^2(\phi)\chi^2]/2$ on a spatially flat, homogeneous shell. The scalar force convention is $j=(\partial_\phi m^2)\langle\chi^2\rangle/2$. The comparison uses spatially homogeneous and isotropic states with zero one-point functions. Initial sources are evaluated with the Yukawa switch off.

## Preparation and one-point sources

Write $k=q/a_i$, $E=\sqrt{k^2+m_i^2}$ and $d\mu=q^2\,dq/(2\pi^2a_i^3)$. At preparation the positive-frequency mode and its proper-time derivative are $u_i=(2a_i^3E)^{-1/2}$ and $\dot u_i=-i E u_i$. Apart from the momentum delta functions, the two oscillator moments are $n=\langle a^\dagger_{\bf q}a_{\bf q}\rangle$ and $c=\langle a_{\bf q}a_{-{\bf q}}\rangle=R+iI$ for each magnitude $q$. At fixed occupation, compare this preparation to $c=0$ on the same history. Differences are restricted to the populated band, so common state-independent subtraction terms cancel.

Contraction of $\chi_{\bf q}=u_i a_{\bf q}+u_i^*a^\dagger_{-{\bf q}}$ gives these covariance differences. The last expression uses the symmetric field–velocity product.

$$
\begin{aligned}
\Delta\langle\chi^2\rangle &= \int d\mu\,R/E, &
\Delta\langle\dot\chi^2\rangle &=-\int d\mu\,E R,\\
\Delta\langle (\nabla\chi)^2/a_i^2\rangle &=\int d\mu\,k^2R/E, &
\Delta\langle\{\chi,\dot\chi\}/2\rangle &=\int d\mu\,I.
\end{aligned}
$$

Substituting into the minimally coupled stress tensor gives the following sources above the common initial vacuum. The isotropic pressure is half the kinetic variance minus one sixth of the physical gradient variance minus half the mass-weighted field variance.

$$
\begin{aligned}
\rho &=\int d\mu\,n E,\\
p &=\int d\mu\,\left[\frac{n k^2}{3 E}
          -\left(m_i^2+\frac{2k^2}{3}\right)\frac{R}{E}\right],\\
j &=\frac{\partial_\phi m_i^2}{2}\int d\mu\,\frac{n+R}{E}.
\end{aligned}
$$

Mean momentum flux and anisotropic stress vanish by isotropy. Define $M_0=\int d\mu\,R/E$ and $M_2=\int d\mu\,k^2 R/E$. Then $\Delta p=-m_i^2 M_0-2 M_2/3$, $\Delta j=(\partial_\phi m_i^2)M_0/2$, and $\Delta\rho=0$ mode by mode. Here $m_i^2=100^2(\phi_i-0.5)^2$ in archived units, with nonzero first derivative on the initial slice. Thus matching the two real moments to zero is necessary and sufficient to match all three sources at fixed $n$. Either physical $k^2$ or comoving $q^2$ gives the same zero condition on a specified slice.

## Applying the criterion on both meshes

The four shell arms on each mesh have the same prefix. The two meshes differ slightly at preparation. The archived fine prefix gives $a_i=0.5316569642$, $m_i=46.09314802$ and $a_i m_i=24.50574314$. The coarse prefix gives $0.5317148519$, $46.09141720$ and $24.50749107$, respectively. A witness with $R=0$ at every momentum satisfies the two source constraints exactly on both meshes. For a profile with nonzero $R$, evaluate both moments separately on each mesh.

An equivalent four-moment formulation avoids subtracting nearly identical mesh weights. Let $u=\sqrt{q^2+\lambda_f}$, $v=\sqrt{q^2+\lambda_c}$, $d=\lambda_c-\lambda_f$ and $\lambda=(a_i m_i)^2$. For either mesh the two zero moments are equivalently orthogonality of $R$, with measure $q^2dq$, to the respective frequency and inverse frequency. Thus for both meshes together they amount to orthogonality to each of

$$
u,\qquad 1/u,\qquad 1/[uv(u+v)],\qquad 1/[uv(u+v)^3].
$$

The change of basis follows from $1/u-1/v=d/[uv(u+v)]$ and $v-u-d(1/u+1/v)/4=-d^3/[4uv(u+v)^3]$. An admissible covariance must also satisfy $|c|^2\le n(n+1)$ pointwise. A classical mixture of product coherent states with zero means can realize the smaller disk $|c|\le n$ for each opposite-momentum pair.

## What the shell sources specify

At fixed geometry and initial bulk data, equality of $\rho,p,j$ gives equality of the parent contributions to the mean metric and scalar junction conditions for the two preparations. Both still contribute the full stress of the initially populated band relative to the archived empty state. The quantum covariance contains a mixed field–velocity component that the initial stress and scalar force do not measure.

For free evolution near preparation, $\ddot\chi_q+3H\dot\chi_q+E_q(t)^2\chi_q=0$. On a common shell with the two moments matched, this gives

$$
\Delta\dot\rho_i=0,\qquad
\Delta\dot j_i=(\partial_\phi m_i^2)\int d\mu\, I,\qquad
\Delta\dot p_i=-2\int d\mu\,(m_i^2+2k^2/3) I .
$$

The first identity also follows from $\dot\rho+3H(\rho+p)=j\dot\phi$. Equality of the initial sources therefore permits subsequent source differences, including for purely imaginary $c$. Determining the geometry at later times requires coupled evolution with a specified initial quantum state.

[Initial oscillator and mass conventions](https://github.com/maldonado-research/HDblast/tree/0c79f123c00918f16647df9c868b4ca889336dc5/research/HDBLAST_TAGGED_PARENT_20261119/primary) · [Signed shell junctions](https://doi.org/10.5281/zenodo.22922928) · [Preregistered follow-up](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-6996927941)
