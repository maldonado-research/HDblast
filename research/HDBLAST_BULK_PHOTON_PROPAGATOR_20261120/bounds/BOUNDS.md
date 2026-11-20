# What selects a photon coupling, and how accurately does it behave as a four-dimensional field?

On the tuned wall, the five-dimensional scalar is $\phi(y)=-\tanh y$, the shell retains $y\ge b$, and $A'(y)=-W(\phi)/3$ with $W=1-\phi+\phi^3/3$. All lengths below are in units of $L_0$. The original Einstein–scalar action leaves the photon sector unspecified. Once one adds a Maxwell action $-\int \sqrt{-g_5}\,G(\phi) F_{MN}F^{MN}/(4g_5^2)$, one can ask what further physical requirements would determine its kinetic function $G$. We use a fixed shell and match the even, source-free vector across two identical copies of its retained space. Norms and kernels below are per copy, with the common factor $1/g_5^2$ omitted.

## An exact uplift of a null electromagnetic wave selects one weight

Consider the separated ansatz $ds^2=dy^2+e^{2A(y)}\hat g_{\mu\nu}(x)dx^\mu dx^\nu$, $\phi=\phi(y)$, $F_{\mu y}=0$, $F_{\mu\nu}=\hat F_{\mu\nu}(x)$. We require a four-dimensional source-free Einstein–Maxwell solution with nonzero electromagnetic stress and $\hat F_{\mu\nu}\hat F^{\mu\nu}=0$ to lift while retaining the original warp factor and scalar profile. Write $\hat G_{\mu\nu}=\lambda\hat T_{\mu\nu}$ with $\lambda>0$ and $\hat T_{\mu\nu}$ the Maxwell stress with no gauge prefactor. The five-dimensional Einstein equation after subtracting the original flat-slice background in the tangential components reads

$$
\hat G_{\mu\nu}=\frac{\kappa_5^2}{g_5^2}G(\phi(y))e^{-2A(y)}\hat T_{\mu\nu}.
$$

Separation therefore requires $G(\phi(y))=C e^{2A(y)}$ wherever the uplift is imposed, with $\lambda=\kappa_5^2 C/g_5^2$. This also suffices locally. The $yy$ geometric change is $-e^{-2A}\hat R/2=0$ by four-dimensional tracelessness. Its Maxwell source is proportional to $\hat F^2=0$, as is the additional scalar source $G_{,\phi}e^{-4A}\hat F^2/(4g_5^2)$. The mixed stress vanishes. The tangential Maxwell equation reduces to $\hat\nabla_\mu\hat F^{\mu\nu}=0$, and the normal Maxwell equation and source-free normal boundary condition hold. Since $K_{\mu\nu}$ remains proportional to the induced metric and the scalar normal derivative stays fixed, the original tuned tension meets the metric and scalar junction conditions. The construction holds on a static cut of the tuned branch. A moving wall changes the normal in the electromagnetic boundary condition.

Freitas, Alencar and Landim (2020) already obtained the separated tangential factor and the choice $G\propto e^{2A}$ for a constant Maxwell zero mode in section 5.2, equation 62 and the following paragraph. Alencar and colleagues (2026) also discuss local trace conditions. For the null field here, both the intrinsic Maxwell trace and its radial pressure vanish. We use that earlier separation condition on the specific HDblast wall, check the full null ansatz, and calculate the contrasts and continuum response it predicts. [Freitas et al. (2020)](https://arxiv.org/abs/2001.01267) · [Alencar et al. (2026), conditions 119–121](https://arxiv.org/html/2601.05190#S3.E120)

For the registered cubic wall this gives an explicit function on the open scalar branch,

$$
G(\phi)=C(1+\phi)^{5/9}(1-\phi)^{-1/9}\exp(-\phi^2/9),\qquad -1<\phi<1 ,
$$

in the global convention $A(0)=0$. Equivalently $d\log G/d\phi=2W/[3(1-\phi^2)]$ follows from dividing $2A_y$ by the scalar gradient.

A nonempty example has four-metric $-2\,du\,dv+dx_1^2+dx_2^2+H\,du^2$, wave $\hat F=E(u)\,du\wedge dx_1$ and $H=-\lambda E(u)^2(x_1^2+x_2^2)/2$. Then $\hat R_{uu}=-\nabla_\perp^2 H/2=\lambda E(u)^2$. This is a local exact null wave. It supplies neither a production process nor a cosmological photon state. The uplift uses the extra requirement of separation with unperturbed profiles and an independent Maxwell sector. Other five-dimensional solutions with responding scalar or warp factor need not select this coupling.

Use **one global convention** $A(0)=0$ for all cuts, set $C=1$, and let $N_b=\int_b^\infty e^{2A(y)}dy$. The archived gravitational coefficient and the projected electromagnetic coefficient obey $f_b=2e^{-2A_b}N_b$ and $I_b^2=2N_b/g_5^2$. Independently resetting $A=0$ on each new cut inside the definition of $G$ would change the microscopic coupling. Taking the endpoint values of the two archived first descents gives

| First minimum | $N$ at launch | $N$ at first minimum | $\log(I_f/I_i)$ | Marginal coefficient at first minimum |
|---|---:|---:|---:|---:|
| $+0.5$ | 1.037904 | 1.676865 | +0.239861 | 1.144121 |
| $-0.5$ | 1.033639 | 0.588743 | −0.281425 | 1.566618 |

All four reconstructed $f_b$ agree with the archived rolling-leaf values to $10^{-9}$ or better. These are contrasts between static cuts associated with the indicated field values. They do not give a separated solution during the roll. The selected weight behaves as $(1+\phi)^{5/9}$ toward the limiting throat value $\phi=-1$. It is smooth on the open branch and has unbounded first derivative in $\phi$ at that endpoint. A completion smooth through that scalar value would require additional work. In an orthonormal frame built from a fixed four-dimensional observer, the wave stress at a point of nonzero four-dimensional flux scales as $G e^{-4A}\propto e^{-2A}$ and can grow without bound down the throat. Vanishing scalar invariants of null electromagnetic stress do not establish global regularity or an admissible incoming scattering state there.

## An endpoint restriction independent of the strength of localization

For comparison take $G=(1-\phi^2)^r=\operatorname{sech}^{2r} y$, $r>0$. The constant Maxwell mode has finite single-copy norm

$$
N_b=B(1/2,r)\,I_{(1+\phi_b)/2}(r,r),
$$

where $I_x$ is the regularized incomplete beta function. This follows by changing the variable in the integral of $G$ to $(1-\tanh y)/2$. More generally, if a positive integrable $G(y)$ is even about $y=0$, the positive-$\phi_b$ launch already contains more than half its entire weight. Moving to any larger $\phi_b$ can increase its integral by a factor smaller than two. Thus on that roll **every such weight** obeys $0<\Delta\log I<\log(2)/2=0.346574$. This bound uses symmetry of the weight. The warp-selected weight above is asymmetric.

On the other archived roll, both cut positions satisfy $0<b_i<b_f$. Its contrast is negative, and the magnitude increases strictly with $r$. To see this, write $N_f/N_i$ as the tail probability for $y\ge b_f$ in the density proportional to $\operatorname{sech}^{2r}y$ on $y\ge b_i$. Increasing $r$ multiplies that density by a strictly decreasing function of $y$ and hence strictly reduces the probability. The ratio tends to one as $r\to0^+$ and to zero as $r\to\infty$ on these distinct cuts. There is one exponent with $\Delta\log I=-3$. It is $r=16.186541$ for the archived launches. The positive roll at that same exponent gives only $+0.341182$. At $r=2$ the two contrasts reproduce $+0.260028$ and $-0.579973$ in the earlier study. Contrast alone leaves production and leakage undetermined. [Archived numerical inputs](../prior/origin/ORIGIN.json)

## Upper and lower bounds on the static current response

At Euclidean four-momentum $Q>0$, denote by $K(Q)$ the action kernel for a vector transverse to that momentum with unit boundary value. Use the local momentum frame with $A(b)$ subtracted in the metric. Define $w=G$, $p=w e^{2(A-A_b)}$, $W_t(y)=\int_y^\infty w$, and $N=W_t(b)$. The measure and stiffness are positive. The spectral problem in this coordinate is

$$K(Q)=\inf_{u(b)=1}\int_b^\infty [p(u')^2+Q^2 w u^2]\,dy.$$

The infimum is over finite-action profiles with a natural open-throat endpoint. Define the relative loss in this **kernel** by $\Delta=1-K/(NQ^2)$. A constant test function gives $\Delta\ge0$. The response to an externally imposed current is proportional to $1/K$, with a corresponding different denominator for its relative error.

For finite $R>b$, put $N_R=\int_b^R w$, $J_R=\int_b^R W_t^2/p$, $h_R(y)=\int_b^{\min(y,R)}W_t/p$ and $H_R=h_R(R)$. Integration by parts gives $\int w h_R=J_R$ and $\int p(h'_R)^2=J_R$. A trial profile $1-d h_R$, optimized over the real constant $d$, together with $\int w h_R^2\le H_R J_R$, gives the lower bound below. For the upper bound, first discard the positive tail of the functional. On the remaining interval use $N_R=\int_b^R wu-\int_b^R (\int_y^R w)u'$ and Cauchy–Schwarz. Substituting $W_t$ for the truncated tail only weakens its denominator. The result is

$$
\frac{Q^2J_R}{N(1+Q^2H_R)}
\ \leq\ \Delta\ \leq\
\frac{N-N_R}{N}+\frac{N_R}{N}\frac{Q^2J_R}{N_R+Q^2J_R}.
$$

Each side can use a different $R$. If $J=\int_b^\infty W_t^2/p$ is finite, $\Delta\le Q^2J/(N+Q^2J)$ at *all* momenta. For example, at $r=2$ the coefficients $J/N$ on $\phi_b=0,+0.5,-0.5$ are 0.277516, 0.811375 and 0.144069. The four-dimensional kernel is guaranteed within 1% of $NQ^2$ for $QL_0\le0.19078,0.11158,0.26479$ respectively. The lack of a bulk mass gap permits this level of accuracy because the current's overlap with arbitrarily slow bulk modes can be small. Static accuracy does not establish a bound on radiation from a moving wall.

## Three different low-momentum laws on the same geometry

Write $k=5/9$, so $A(y)=-ky+A_\infty+o(1)$, $A_\infty=(2/9)\log2-1/18$. For the symmetric family $w\sim C_g e^{-2ry}$ with $C_g=2^{2r}$. The selected weight has the same form with $r=k$ and $C_g=e^{2A_\infty}$. For both, $W_t\sim w/(2r)$, hence $W_t^2/p$ is asymptotic to a positive constant times $e^{2(k-r)y}$ and $H_R$ grows as $e^{2kR}$. In conformal radial distance $dz=e^{A_b-A}dy$ the canonical zero mode $(G e^{A-A_b})^{1/2}$ falls as $z^{-1/2-r/k}$. Its radial potential tends to zero as an inverse square, leaving a continuum down to zero mass on the retained infinite leaf.

- When $r>k$, $J$ is finite and $\Delta/Q^2\to J/N$. For the lower limit, first take $Q\to0$ at fixed $R$ in the trial bound, then $R\to\infty$.
- When $0<r<k$, $\Delta$ is bounded above and below by positive constants times $Q^{2r/k}$ for sufficiently small $Q$. Choose $R=b+\log(1/Q)/k$ in both bounds. Then $Q^2H_R$ stays bounded, $N-N_R$ is of order $Q^{2r/k}$ and $Q^2 J_R$ is of that same order. The kernel has no finite coefficient for an ordinary $Q^4$ correction.
- When $r=k$, with $T=\log(1/Q)$ and momenta in the stated units,

$$
\lim_{Q\to0}\frac{\Delta}{Q^2 T}
=\frac{C_g\,e^{2(A_b-A_\infty)}}{4k^3N}.
$$

For the upper bound take $R=b+T/k$. Then $J_R/T\to C_g e^{2(A_b-A_\infty)}/(4k^3)$ and the omitted norm is only of order $Q^2$. For the lower bound take $R=b+(T-\tfrac12\log T)/k$. The same $J_R/T$ limit holds, while $Q^2 H_R\to0$. The coefficient for the selected weight reduces to **$1/(2k^3 f_b)$**. Thus the conditional null uplift fixes both the endpoint normalization and the leading logarithmic departure from the four-dimensional kernel using the archived geometry. The logarithm's reference scale affects only lower-order-in-$T$ terms. The coefficient concerns the normalized kernel on a fixed leaf and carries length squared in physical units.

## Independent checks

The program compares beta-function norms with direct proper-distance integrations. For each of the 15 symmetric-family leaf/exponent pairs, and at each of the four selected-weight endpoints, it brackets the kernel numerically by independently integrating from a distant Dirichlet or Neumann boundary. With $B=-e^{A-A_b}u'/(Qu)$ and $s=(B-1)/(B+1)$ the regular bounded variable obeys $s'=2Qe^{A_b-A}s-(\log w+A)'(1-s^2)/2$. The distant values $+1$ and $-1$ represent Dirichlet and Neumann data respectively. Their computed losses at $Q=0.01,0.1$ satisfy the analytic enclosures within the specified ODE tolerance. The code also applies the analytic enclosures at $Q=10^{-40}$ to check the marginal coefficients in a regime where directly subtracting two almost equal kernels would be inaccurate. The bounds are analytical; reported decimal values use floating-point quadrature.

[Executable calculation](bounds.py) · [All results](BOUNDS.json)
