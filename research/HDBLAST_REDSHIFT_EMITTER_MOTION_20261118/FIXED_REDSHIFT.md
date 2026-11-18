# Independent calculation at fixed emission-to-reception frequency ratio

Take the receiving observer and the two launch caps of the affine comparison. Work to first order in its canonical amplitude. The shell perturbation is $2Fh_0$ in a longitudinal representative, $F=U(q)Y(n)$, $Y=n_1^2-n_2^2$. Write $\eta=\arctan(\sinh q)$, $a_{\!s}(\eta)=\rho_b\cosh q$, so $h_0=a_{\!s}^2(-d\eta^2+\kappa)$ and increasing $q$ goes into the past. A conformal representative with this first derivative suffices to calculate these invariant derivatives, provided the source and receiver worldlines are transformed with the metric. The actual shell need be conformal to the background only to this order.

The receiver event is anchored to its normal observer at initial label $n_*=(0,0,1,0)$, at elapsed proper clock $\rho_bq_0$, with $q_0=\ell(49/32)$. Along this worldline both $Y$ and $D Y$ vanish. Thus the first-order source/receiver prescriptions below agree there and the receiver’s clock, direction triad and peculiar velocity need no first-order correction in the chosen representative. Null directions are labelled on its physical rest sky. The background ray at $\eta$ has sphere position $n=n_*\cos\theta+v\sin\theta$, with $\theta=\eta-\eta_0$ and round-unit $v$. The angular tangent $e$ has $D_e Y=(v_1^2-v_2^2)\sin 2\theta$.

## Two ways to supply sources

For a free source initialized normal to the quiet cut, write $V=dn/d\eta$ for its first-order peculiar velocity relative to conformal coordinate rest, with the tangent taken *toward the past*. To first order the timelike geodesic equation is
$$
 V'+(a_{\!s}'/a_{\!s})V=-D F,\qquad
 V(q,n)=-\operatorname{sech}q\int_0^qU(q')\,dq'\,D Y(n) .
$$
The prime in this equation differentiates with respect to conformal time. These equations are evaluated on each background timelike source curve. Its first-order displacement changes the evaluated $V$ or $F$ only at higher order. An actual observer/source intersection is defined using a congruence smooth over the whole sky and through the needed clock interval.

As a separate model set $V=0$ for supported emitters in this longitudinal conformal representative. One covariant way to implement this to first order is to attach to each initial label $n$ its physical rest vector of norm $|D Y(n)|_\kappa$ and initial components proportional to the raised round gradient. Fermi–Walker transport this vector $w$ along the worldline. At elapsed source proper time $r$ impose acceleration $\epsilon U(r/\rho_b)w/[\rho_b\cosh(r/\rho_b)]$ and initial rest at the quiet cut. The background transported vector is the corresponding round gradient in the orthonormal rest frame. Conformal coordinate rest has first acceleration $D F/a_{\!s}$ in that frame; hence this controller gives $V=0$ at first order. Proper acceleration is unchanged when a unit timelike tangent is reversed. This example requires the prescribed force template and labels and is not a prediction for undriven matter. Different preparations introduce different Doppler terms.

## Frequency ratio and endpoint shift

For the past-directed photon normalized at reception, the positive endpoint frequency can be calculated using either sign convention for the unit emitter tangent and taking the absolute contraction. Conformal transport of the null wave vector gives at fixed source coordinate $q$
$$
 \left.\delta\ln(1+z)\right|_q= F_o-F_s+V_o\cdot e_o-V_s\cdot e_s.
$$
For instance a past unit tangent is proportional to $\partial_\eta+V$, and its positive frequency is proportional at first order to $a_{\!s}^{-1}e^{-\epsilon F}(1-\epsilon V\cdot e)$. The receiver terms vanish at $n_*$. The unperturbed redshift is $1+z=c_0/c$, with $c_i=\cosh q_i$ and $c=\cosh q$. For positive $q_1>q_0$, this is a **blueshift**. Set $q_1=\ell(99/64)$ and hold $z_*=c_0/c_1-1$ fixed, rather than imposing a time or affine endpoint in the perturbed geometry. The nonzero background slope $d_q\ln(1+z)=-\tanh q$ gives the local endpoint derivative
$$
 \delta q_s=\left.\delta\ln(1+z)\right|_{q_1}/\tanh q_1 .
$$
Positivity of the interior clock margin and uniform nonvanishing of this background slope near $q_1$ imply a local branch for sufficiently small amplitudes, by the implicit function theorem uniformly on the compact sky. A localized source catalogue without a continuous congruence need not contain emitters at these endpoints.

Conformal propagation from this receiver gives $D_A|_{q_s}=e^{\epsilon F_s}a_{\!s}(q_s)\sin\theta$ through first order. The screen area is invariant under changes of the orthonormal source frame: the screen is the metric quotient of the null orthogonal space by the ray tangent. Direction-dependent endpoint shifts along the ray likewise do not tilt its measured transverse area. For $\epsilon=0$, $D_0(q)=\rho_b(\sinh q-\sinh q_0)/c_0$. Put $s_i=\sinh q_i$, $\Delta s=s_1-s_0$ and $L=D_0(q_1)$. At fixed redshift,
$$
 \frac{\delta D_A}{L}= F_s+K_z\left.\delta\ln(1+z)\right|_{q_1},\qquad
 K_z=\frac{c_1^2}{s_1\Delta s}.
$$
In particular $\sin\theta=\Delta s/(c_0c_1)$ and $\cos\theta=(1+s_0s_1)/(c_0c_1)$. Write $b=99/64$, $M(x)=\mu\rho_b^2U(\ell x)/\ell$ and $T=(1+s_0s_1)/(s_1 c_1^2 c_0^2)$. The signed axial coefficients for the two source protocols are then
$$
 \frac{\mu\rho_b^2}{\ell}\frac{\delta D_A(v)}{L}
    =(v_1^2-v_2^2)Z,\qquad
 Z_{\rm supported}=-T\Delta s\,M(b),\qquad
 Z_{\rm free}=T\left[2\ell c_1\int_0^b M(x)\,dx-\Delta s\,M(b)\right].
$$
For a cap $v_1\ge4/5$, the angular factor averages to $18/25$; for $v_2\ge4/5$ its mean is $-18/25$. Since the unperturbed distance is isotropic, the variation for $\sqrt{\langle D_A^2\rangle_{\rm cap}}$ is the cap average of the pencil variations. The physical distance $L$ here is the common **background** value. The actual terminal affine lengths vary with the metric, the emitter congruence and the direction in the sky.

## Cellwise calculation and robustness

For $a=49/32$, the replay uses the archived $u(t)$ intervals in all 2240 cells from the front $x=1$ to $b$. Each future interval is widened by $3\times10^{-6}$ in the normalized units. Exact integration against $b-t$ gives the endpoint of the zero-data flat double integral $z_0(b)$, while integration against $(b-t)^2/2$ gives $\int_1^bz_0(x)\,dx$. They enclose their true metric analogues after adding the archived uniform Volterra/geometry error and $(b-1)$ times that error respectively. The source is exactly quiet before 1. This gives $0.23776519<M(b)<0.23868027$ and $0.03184748<\int_0^b M<0.03197999$. For $q\le q_1$ use $1\le\cosh q<C=1+10^{-11}$, $q\le\sinh q\le Cq$ and $(b-a)\le\Delta s/\ell\le C(b-a)$. All enclosure operations in the replay use rational arithmetic.

The positive cap has $\delta D_{\rm cap}/L$ between $7.27866\times10^{-10}$ and $8.06483\times10^{-10}$ for the free preparation and between $-4.99243\times10^{-11}$ and $-4.50939\times10^{-11}$ for the specified supported emitters. These are dimensionless first variations per unit canonical amplitude. The baseline $z_*$ is about $-2.40478515625\times10^{-14}$. The cap area coefficient is twice its distance coefficient. At fixed redshift the cap response can change sign when only the source acceleration changes.

There is substantial numerical conditioning in matching redshifts. The exact background relation is $d\ln D_0/d\ln(1+z)=-K_z$ with $K_z\simeq4.13737\times10^{13}$. For the cap half-difference $(D_+-D_-)/(2L)$, a common first-order shift of the log-redshift threshold cancels. A differential offset $e=(e_+-e_-)/2$ contributes $-K_z e$ at first order about the isotropic background. Unmodeled source Doppler errors can have the same magnitude with a sign depending on whether the error refers to the selected redshift or the predicted source frequency. At fixed other inputs, differential threshold uncertainty below approximately $1.75\times10^{-23}|\epsilon|$ (free) or $1.08\times10^{-24}|\epsilon|$ (supported) is sufficient to preserve these lower sign margins at first order. The numbers describe model normalization and sensitivities. They do not constitute an instrument forecast or a finite-amplitude error budget.

The null energy condition bounds $D_A$ by the actual frequency-normalized affine length before caustics. At a specified redshift that length itself varies, so an increase over the reference $L(z_*)$ alone does not establish a violation. An observational interpretation needs an identified emitter population and its velocities or acceleration law, redshift calibration, source-size or other distance information, and a range free of caustics. [Independent numerical replay](redshift_check.py) · [Exact endpoints](REDSHIFT_CHECK.json)
