# Spherical mass, wall pressure and an exact Maxwell bulk

Use a two-dimensional Lorentzian metric $h_{ab}$ and area radius $r$, so $ds^2=h_{ab}dx^a dx^b+r^2 d\Omega_3^2$. Set $\Lambda=-6/\ell^2$. Three independent canonical Abelian fields have identical spherical amplitudes $\mathcal A$ in left-invariant directions, $A^{(I)}=\mathcal A\sigma_I/\sqrt3$, with $d\sigma_I=-\epsilon_{IJK}\sigma_J\wedge\sigma_K$ (sum over all $J,K$). All integrals below factor out the unit sphere volume $2\pi^2$. The two-dimensional bulk Lagrangian with measure $\sqrt{-h}$ is $-r(D\mathcal A)^2/2-2\mathcal A^2/r$. The Maxwell equation and base stress are

$$
D_a(r D^a\mathcal A)=4\mathcal A/r,\qquad
T_{ab}=\frac{D_a\mathcal A D_b\mathcal A}{r^2}
-\frac12 h_{ab}\left[\frac{(D\mathcal A)^2}{r^2}+B^2\right],
\qquad B=2\mathcal A/r^2.
$$

The pressure in each sphere direction is $p_\Omega=-(D\mathcal A)^2/(6r^2)+B^2/6$. The spherical reduction keeps one co-closed angular mode per gauge field. The three stresses sum isotropically.

## The enclosed mass seen by a moving wall

Define $F=(Dr)^2$ and the Schwarzschild–AdS mass parameter $\mu=r^2(1+r^2/\ell^2-F)$. The base Einstein equations give

$$
D_aD_b r=\left[\frac{1-F}{r}-\frac{\Lambda r}{3}\right]h_{ab}
-\frac{\kappa_5^2 r}{3}(T_{ab}-T^c{}_c h_{ab}),\qquad
D_a\mu=\frac{2\kappa_5^2 r^3}{3}(T_a{}^b-T^c{}_c\delta_a^b)D_b r .
$$

Here $T^c{}_c=-B^2$ is the **two-dimensional** trace. At a wall of radius $a$, let $U$ be its unit tangent and $n$ the outward normal of the retained interior directed toward the wall. Write $X=U(\mathcal A)/a$, $Y=n(\mathcal A)/a$, $N=n(r)$, and $q=-T_{Un}=-XY$ for incident net power per unit proper area from that side. The energy density and normal pressure are

$$
\rho_b=\tfrac12(X^2+Y^2+B^2),\qquad P=\tfrac12(X^2+Y^2-B^2),
\qquad |P|\leq\rho_b,\qquad P+B^2/2\geq |q|.
$$

This last bound follows from $X^2+Y^2\geq2|XY|$. A magnetic field can produce negative normal pressure while respecting the bulk null energy condition. The two derivatives of enclosed mass are

$$
\dot\mu=-\frac{2\kappa_5^2 a^3}{3}\big(P\dot a+qN\big),\qquad
n(\mu)=\frac{2\kappa_5^2 a^3}{3}\big(\rho_b N+q\dot a\big).
$$

The flux determines only part of $\dot\mu$ when the wall moves. For a radially outgoing null stream $P=q>0$ the familiar Vaidya loss proportional to $(N+\dot a)q$ is recovered. For two identical retained interiors the positive angular junction gives $N/a=\kappa_5^2(\lambda+\rho_s)/6$, where $\lambda$ is tension and $\rho_s,p_s$ refer to all material on the wall. The temporal junction and the Hessian identity give

$$
\frac{\ddot a}{a}=\frac{\kappa_5^4}{36}(\lambda+\rho_s)(\lambda-2\rho_s-3p_s)
-\ell^{-2}-\frac{\mu}{a^4}-\frac{\kappa_5^2P}{3}.
$$

For distinct interiors set $b_i=N_i/a$. The corresponding two-sided expression is $(\ddot a/a+\ell^{-2})\sum_i 1/b_i=\kappa_5^2(\lambda-2\rho_s-3p_s)/3-\sum_i(\mu_i/a^4+\kappa_5^2P_i/3)/b_i$, with $\sum_i b_i=\kappa_5^2(\lambda+\rho_s)/3$ and the same $\ell$ on both sides.

A continuous tangential gauge field with surface action $K(\psi)[a\dot{\mathcal A}^2-4\mathcal A^2/a]/2$ has surface density $K(X^2+B^2)/2$, pressure one-third that density, and equation $(K a\dot{\mathcal A})^{\displaystyle\cdot}+4K\mathcal A/a=-a\sum_i n_i(\mathcal A)$. Its energy ledger is

$$
\dot\rho_{\rm em}+4H\rho_{\rm em}=\sum_i q_i-\tfrac12\dot K(X^2-B^2).
$$

The wall scalar gains the last term with opposite sign. Any modeled conversion of scalar energy to wall radiation is internal to $\rho_s$. Then $\dot\rho_s+3H(\rho_s+p_s)=\sum_i q_i$. No single local inequality in this derivation determines a capture fraction or a bounce.

## Characteristic data and checks

Choose $h=-2S\,du\,dv$, $S=e^{2\sigma}>0$, with $r_u<0<r_v$ in a future-oriented untrapped region. Subscripts denote partial derivatives. Along each null line the constraints are $r_{uu}-2\sigma_u r_u=-\kappa_5^2\mathcal A_u^2/(3r)$ and its $v$ counterpart. The evolution equations are

$$
r_{uv}=-S\left(\frac r{\ell^2}+\frac{\mu}{r^3}-\frac{2\kappa_5^2\mathcal A^2}{3r^3}\right),\qquad
2r\mathcal A_{uv}+r_u\mathcal A_v+r_v\mathcal A_u=-4S\mathcal A/r ,
$$

$$
\sigma_{uv}=S\left[-\frac1{2\ell^2}+\frac{3\mu}{2r^4}-\frac{5\kappa_5^2\mathcal A^2}{3r^4}\right]
-\frac{\kappa_5^2\mathcal A_u\mathcal A_v}{6r^2}.
$$

In particular $\mu_u=(2\kappa_5^2/3)[2\mathcal A^2r_u/r-r\mathcal A_u^2r_v/S]\leq0$, and the interchange $u\leftrightarrow v$ gives $\mu_v\geq0$. Positivity of enclosed mass needs a nonnegative reference mass and an outward integration in the untrapped region. A compatible characteristic initial value prescription gives $\mathcal A$ on intersecting null curves, fixes the two null coordinate freedoms and the corner geometry/mass, then integrates the two constraints along those curves. At a timelike excision tube one prescribes one incoming wave characteristic and receives the other. The wall gauge equation and gravitational jumps relate the one-sided normal gradients to its fields and motion. Independently prescribing both wave characteristics at a physical emitting boundary and all values at the wall would generally overdetermine a forward evolution.

## Exact nonvacuum benchmark and its wall limitation

In five-dimensional flat-spatial FLRW coordinates use
$ds^2=-dt^2+c(t)^2(d\chi^2+\chi^2d\Omega_3^2)$, $r=c\chi$, and $\mathcal A=k\chi^2$. In Cartesian comoving spatial coordinates these are three constant spatial two-forms of definite duality with equal magnitude. Maxwell's equation holds because $F^{tM}=0$ and $\partial_j(\sqrt{-g}F^{ji})=0$. Summing their stress gives zero in all four spatial directions and $\rho=4k^2/c^4$. Define $C=2\kappa_5^2 k^2/3$. A full self-gravitating bulk solution on $0<2(t-t_0)/\ell<\pi$ is

$$
c^2=\ell\sqrt C\sin[2(t-t_0)/\ell],\qquad
(\dot c/c)^2=C/c^4-\ell^{-2},\qquad \mu=C\chi^4.
$$

The homogeneous spatial magnetic stress behaves as pressureless matter in four spatial dimensions. This local bulk benchmark may be restricted to an untrapped ball. Its big-bang and crunch endpoints already enter through the chosen FLRW solution. A comoving wall between two copies of the interior has $\dot{\mathcal A}=0$, $n_i(\mathcal A)=2k\chi/c$. Its Maxwell jump would require $K=-a$ for $k\chi\ne0$. A positive surface kinetic term therefore excludes *this comoving gluing*.

The reproducible check constructs all five-dimensional Cartesian stress components directly and differentiates the analytic scale factor using five-point finite differences. It tests the double-null equations, both mass derivatives in 105 moving frames including contraction and expansion, the pressure inequality and a vacuum Schwarzschild–AdS control. The largest absolute error in the five-dimensional Einstein equation is below $3.2\times10^{-9}$ in the sampled units. [Calculation](identities.py) · [Check results](IDENTITIES.json)

A separate inverse annulus reconstruction uses the selected wall history to evaluate the characteristic actually needed at a distinct source surface. Its data are archived under the `collar` key in the same numerical product. [Source reconstruction and convergence](COLLAR.md)
