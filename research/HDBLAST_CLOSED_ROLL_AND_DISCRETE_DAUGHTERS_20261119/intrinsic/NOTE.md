# Can the rolling slices be replaced by a round cavity of equal volume?

An instantaneous spatial frequency answers a geometric question. Freeze the physical three-metric and the daughter kinetic factor on one proper-time slice. Find the harmonic frequencies a free spectator would have if those coefficients were held fixed. On the actual evolving spacetime both vary with time and the spatial eigenfunctions can mix. The frequency calculation alone assigns no particles to the field.

Write the spatial metric as $e^{2L}d\chi^2+e^{2M}\sin^2\chi\,d\Omega_2^2$. These are the full metric exponents, including the round background scale. The minimally coupled real matter scalar has Einstein-frame kinetic weight $Z=f(\Theta_*)/f(\Theta)$. Its normalization at the hilltop cancels from the eigenproblem. At fixed spherical angular degree $\ell$ a test function is $v(\chi)Y_{\ell m}$. The weak radial eigenproblem has the two bilinear forms

$$
K=\int d\chi\,Z\left[e^{2M-L}\sin^2\chi\,v'g'+\ell(\ell+1)e^L vg\right],\quad
M_0=\int d\chi\,Ze^{L+2M}\sin^2\chi\, vg .
$$

We solve $K v=\omega^2 M_0 v$ in the regular Gegenbauer basis $\sin^\ell\chi\,C_{n-\ell}^{\ell+1}(\cos\chi)$. Each radial eigenvalue occurs $2\ell+1$ times. In the uniform control they combine into the $(n+1)^2$ modes of the round sphere at $\omega_n=\sqrt{n(n+2)}/a$. We repeat the calculation with $Z=1$ on exactly the same geometry and compare both with the round sphere defined by the actual proper volume. At fixed $n$ the reported fractional width is the range over $\ell$ divided by the equal-volume round frequency. Increasing the radial basis from degree 20 to 28, angular quadrature from 96 to 160 points, and spacetime resolution from degree 32 to 48 tests the low eigenvalues.

On the initial constrained bounce the metric perturbation starts at order $A^2$, while $Z$ can have a first-order spatial gradient. For $q=A U_2$ and $z_1=(d\ln Z/d\Theta)_*$, the fractional frequency change is

$$
\frac{\delta\omega_{n\ell}}{\omega_n}
= -\frac{2 z_1 A}{n(n+2)}
\left(1-\frac{2\ell(\ell+1)}{n(n+2)}\right)+O(A^2).
$$

First-order perturbation of $K-\omega_n^2 M_0$ gives the integral of $z_1 U_2(|\nabla Y|^2-n(n+2)|Y|^2)$. Integration by parts moves half a Laplacian onto $U_2$, giving $-4z_1\langle U_2\rangle_{n\ell}$. The Gegenbauer recurrence for multiplication by $\cos\chi$ gives the parenthetical expectation above. The first-order weighted centroid vanishes. The code independently evaluates this formula using both signs of a small constraint-satisfying seed.

At time 0.47 on the negative-minimum branch, before first passage for either tested sign, the matter-weighted $n=7$ multiplet has a width of **0.203%** or **0.604%** of its round frequency for $A=+0.001$ or $-0.001$. Holding $Z$ constant on those same geometries gives **0.286%** and **0.367%**. Inhomogeneity and its sign therefore matter even at fixed proper-volume radius. These widths alone are insufficient for predicting resonant transfer, which depends on time evolution, spatial overlap and the state of the daughter. [Numerical spectra and convergence](FREQUENCIES.json) · [Eigenvalue calculation](frequencies.py)
