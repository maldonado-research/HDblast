# The distance response near the receiving clock

At the vertex, both the redshift and angular-diameter distance equal zero for an emitter comoving with the receiver. The quantity with a local physical meaning is the slope of distance against redshift. In this model the positive $v_1$ cap of freely falling emitters has a strictly positive first fractional change in that slope. It comes from shear of the freely falling clocks' velocity field. The undisturbed shell contracts toward the future at this reception event, so nearby emitters have a blueshift.

## The directional expansion gives the local slope

Fix the actual receiving event at prescribed elapsed proper duration on the worldline labelled $n_*=(0,0,1,0)$. Write $u$ for the **future** unit timelike field of the freely falling emitters and receiver. Write $k=-u+e$ for the initial **past** null tangent with a unit spatial viewing direction $e$. Along the affine null geodesic, the positive measured frequency in this past-ray convention is $E=u\cdot k$, with $E(0)=1$. The redshift to an event on this ray is $z=E-1$. The geodesic equation gives
$$
 \left.\frac{dz}{d\lambda}\right|_0
  = k^ak^b\nabla_a u_b
  = e^ae^b\nabla_a u_b =: H_e .
$$
Here $u^a\nabla_a u_b=0$ and $u^b\nabla_a u_b=0$ remove the temporal contractions. The background $H_e=H_0=-\tanh q_0/\rho_b$ is negative. A regular vertex pencil has $D_A(0)=0$ and $D_A'(0)=1$. In particular $z=H_e\lambda+O(\lambda^2)$ and $D_A=\lambda+O(\lambda^3)$. For nonzero $H_e$ the slope of $D_A$ against signed redshift at the vertex is $1/H_e$. Its fractional first variation is $-\delta H_e/H_0$. One can equivalently work with the positive slope of distance against blueshift magnitude.

To evaluate $\delta H_e$, use the conformal representative of the shell metric tangent, $h'=2UYh_0$. The physical initially normal congruence and the specified observer and axes must transform together under a chart change. In the dimensionless elapsed coordinate $q$, the freely falling covariant momentum in the spatial sphere directions with **future** orientation is
$$
 u_i=\epsilon\rho_b I(q)D_iY+o(\epsilon),\qquad
 I(q)=\int_0^q U(y)\,dy .
$$
The sign is opposite to the momentum of a past-oriented timelike emitter. This follows from the Euler–Lagrange momentum equation and zero initial angular velocity in the quiet cut. At the critical receiver $Y=DY=0$. There is no first perturbation of its lapse, position, proper-time registration or peculiar velocity. The conformal metric connection variation also vanishes there since $d(UY)=0$. The angular gradient of the velocity survives:
$$
 \delta H_e=\frac{I(q_0)}{\rho_b\cosh^2 q_0}
    \operatorname{Hess}_{\kappa}Y(v,v)
   =\frac{2I(q_0)}{\rho_b\cosh^2 q_0}(v_1^2-v_2^2).
$$
The unit viewing angle $v$ refers to the original orthonormal axes $n_1,n_2,n_4$ of this observer. A first variation of a normalized screen orientation leaves the isotropic background $H_0$ invariant. The trace of the Hessian at this label is $-8Y=0$, so the expansion scalar has zero first variation there. The surviving quadrupole is the symmetric trace-free velocity gradient. This distinguishes the freely falling congruence from a supported congruence at rest in the specified conformal coordinates, whose peculiar velocity gradient is zero. Its four-acceleration and the conformal lapse and expansion corrections vanish at the critical label to first order as well.

Consequently the local fractional slope variation along a direction is
$$
 -\frac{\delta H_e}{H_0}
 =\frac{2I(q_0)}{\cosh q_0\sinh q_0}(v_1^2-v_2^2).
$$
The same expression follows by sending the emission coordinate to $q_0$ from above in the finite-distance result. The cap $v_1\ge4/5$ multiplies this coefficient by $18/25$; $v_2\ge4/5$ gives the opposite result. At zero perturbation the distance slopes are the same in every direction, so differentiating a cap root-mean-square area distance gives the cap mean of the fractional slope at first order. No individual distance divided by its background value is defined exactly at the vertex. [Physical congruence and fixed-redshift derivation](FREE_FALL.md)

## Bounds over the blueshift interval

Put $A=49/32$, $B=99/64$, $q_0=\ell A$, $q=\ell x$, and let $z_0(q)=\cosh q_0/\cosh q-1$. The range of emission times is $A<x\le B$. Let $J(x)=\int_0^x M(y)\,dy$ in the metric normalization of the preceding geometric result. The quiet interval allows $J(A)$ to be bounded by integrating the source cells against the nonnegative kernel $(A-y)^2/2$ from 1 to $A$, then adding the previous uniform geometric reconstruction error integrated over this range. For $x\in[A,B]$, $M(x)>0$ and its previously computed upper bound controls $J(x)-J(A)$. With $c=\cosh q$, $s=\sinh q$, and $K=1+10^{-11}$, the factor $\ell(1+\sinh q_0\sinh q)/(\cosh^2q_0\cosh^2q\sinh q)$ lies between $1/(K^5 B)$ and $(1+K^2\ell^2 AB)/A$. Indeed $1\le\cosh q<K$, $q\le\sinh q<Kq$ throughout the interval. The difference $(\sinh q-\sinh q_0)/\ell$ lies between zero and $K(B-A)$. These inequalities do not divide by the emission–reception separation.

| Quantity for the cap $v_1\ge4/5$ | Outward interval |
| --- | --- |
| Receiver history $J(A)$ | $(0.02825082,0.02836952)$ |
| Fractional first variation of the free-fall cap distance for every $A<x\le B$ | $(6.37512,8.74944)\times10^{-10}$ |
| Fractional first variation of the free-fall cap's local reciprocal-expansion slope | $(6.92819,7.67301)\times10^{-10}$ |
| Fractional first variation of the supported cap distance for $A<x\le B$ | $(-5.38217\times10^{-11},0)$ |

All variations use the signed canonical amplitude. The supported fractional slope change tends to zero at the vertex. The strictly positive finite-separation lower bound for the free-fall cap is allowed by its nonzero limiting shear. The two tracer choices constrain different velocity fields even in the same spacetime. [Exact source and interval calculation](continuum.py) · [Exact rational outputs](CONTINUUM.json)

For the freely falling congruence the common reception condition $z_\epsilon(v,0)=0$ holds for every direction and parameter. On a small compact collar containing the background endpoints, the background affine redshift derivative is strictly negative. Smooth dependence of the redshift and its affine derivative on the chosen metric and velocity field preserves this inequality uniformly over the compact sky for sufficiently small amplitudes. A margin beyond $q=\ell B$ permits inversion on the same local blueshift branch through $z_0(\ell B)$ for every direction. Near the vertex $z_\epsilon(v,\lambda)/\lambda=\int_0^1\partial_\lambda z_\epsilon(v,t\lambda)\,dt$ extends smoothly and stays away from zero; the Jacobi equation gives the same regular extension of $D_A/\lambda$ with value one. These extensions justify the limit of the distance ratio on approaching zero redshift. For any other tracer prescription, the analogous uniform statement through the vertex needs the same reception frequency and a sufficiently regular, monotone local branch. A comparison with a noncomoving emitter velocity at the receiving event need not satisfy this common anchoring condition. Coordinate corrections to redshift or observer angular velocity used en route depend on gauge. The redshift, slope, and normalized area distances compare fixed physical frequency ratios, clocks and directions. They have no dependence on a choice of conformal coordinates when the same physical tracer prescription is mapped with the metric.
