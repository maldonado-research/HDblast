# When does a current on the wall see a four-dimensional photon?

Place a weak transverse test current on a motionless wall. The response contains the familiar massless photon and a continuum of modes that vary along the additional spatial dimension. Counting how much each part contributes to the static field of the current answers a practical question. At which spatial wavelengths does a four-dimensional description work for this specified gauge interaction?

## Boundary action and observable

This calculation uses the tuned flat leaves with dimensionless proper coordinate $y$, wall at $b=-\operatorname{arctanh}\phi_b$, and bulk scalar $\phi(y)=-\tanh y$. The fixed-background metric is $ds^2=dy^2+e^{2A(y)}ds^2_{4,\mathrm{flat}}$, with the wall normalization $A(b)=0$. The warp derivative is $A'=-(1-\phi+\phi^3/3)/3$. Following the earlier localized Maxwell extension, take the positive weight $w=G=(1-\phi^2)^2=\operatorname{sech}^4y$ and $p=w e^{2A}$. Fixing this gauge weight is an extra microscopic assumption. We omit the common gauge coupling and take one copy of the bulk. On two identical copies the kinetic action doubles for a fixed even boundary field, so the pole and continuum shares are unchanged.

For one spacelike Fourier component of the transverse current, $Q>0$ denotes proper momentum in units of $L_0^{-1}$ at the wall. The finite-action radial solution at fixed unit wall amplitude obeys

$$
(pu_Q')'=Q^2 w u_Q,\qquad u_Q(b)=1,\qquad
D(Q)=-p(b)u_Q'(b).
$$

Integrating the equation by parts gives $D(Q)=\int_b^\infty[p(u_Q')^2+Q^2wu_Q^2]dy$. This is the coefficient in the quadratic boundary action, often called the Dirichlet-to-Neumann map. Varying the wall amplitude against a weak boundary current gives the boundary propagator $P(Q)=1/D(Q)$. Both statements assume a free transverse Maxwell sector and no additional wall kinetic or current-contact term.

The constant normalizable mode has squared norm $N=\int_b^\infty w\,dy=(1+\phi_b)^2(2-\phi_b)/3$. Its contribution to $P$ is $P_0=1/(NQ^2)$. The fraction of the full propagator supplied by all the other modes is

$$
\eta(Q)=\frac{P(Q)-P_0(Q)}{P(Q)}=1-\frac{D(Q)}{NQ^2}.
$$

This fraction pertains to the static current response at fixed momentum. A pulse of real radiation needs a timelike source spectrum and an outgoing boundary condition.

## Analytic control at every momentum

Write $T(y)=\int_y^\infty w$. For $Q^2=s$ near zero, $u_Q=1+s v+\cdots$. At fixed $y$, $p v'=-T$ and $v(b)=0$. Equivalently the exact identity $D/s=\int w u_Q$ gives the first derivative at zero by taking a one-sided limit. Since $0<u_Q\le1$, the integral of $w(1-u_Q)/s$ is bounded by $C$ below. Fubini's theorem and dominated convergence then give

$$
C=\int_b^\infty \frac{T(y)^2}{p(y)}\,dy,\qquad
D(Q)=NQ^2-CQ^4+o(Q^4).
$$

The integral converges. In the throat, $A'\to-5/9$ and $T^2/p$ decays as $e^{-(4-10/9)y}$. This argument does not need analyticity of the full response at a gapless threshold.

The self-adjoint radial Maxwell operator for the source-free even channel has nonnegative quadratic form $\int p(u')^2$ and Neumann boundary condition at the wall. Its boundary spectral measure is positive. Separating the constant mode gives

$$
P(Q)=\frac1{NQ^2}+R(Q^2),\qquad
R(s)=\int_{0^+}^\infty\frac{d\mu(m)}{s+m^2},\qquad
R(0)=\frac C{N^2}.
$$

The localized zero mode is the only zero-energy state in this channel. In conformal distance the radial potential tends to zero as an inverse square, so the continuum extends down to zero. Since $R(s)\le R(0)$ and $sR(s)$ is increasing, the extra-mode fraction increases with momentum and satisfies the all-momentum upper bound

$$
0\le\eta(Q)\le\frac{Q^2 C/N}{1+Q^2 C/N}.
$$

There is a related guarantee for the static potential measured by charges at wall separation $r>0$. Each massive mode contributes its Yukawa factor $e^{-mr}$ relative to the massless $1/r$ interaction. Since $(mr)^2e^{-mr}\le4e^{-2}$,

$$
0\le \frac{V(r)-V_0(r)}{V_0(r)}
=N\int_{0^+}^\infty e^{-mr}d\mu(m)\le
\frac{4}{e^2}\frac{C}{Nr^2}.
$$

Here $V$ and $V_0$ refer to the same charges. This is a sufficient separation for a given accuracy. A momentum reciprocal is only a length scale and does not give an exact position-space crossover.

## Results on the three leaves

| Wall field $\phi_b$ | $N$ | $C/N$ | Extra share at $Q=0.1$ | Extra share at $Q=1$ | $Q$ at 1% | $Q$ at 10% | $Q$ at 50% | Sufficient $r/L_0$ for at most 1% static-potential correction |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.666667 | 0.277516 | 0.276% | 19.074% | 0.191589 | 0.659755 | 2.466208 | 3.8760 |
| +0.5 | 1.125000 | 0.811375 | 0.803% | 40.924% | 0.111730 | 0.375467 | 1.237343 | 6.6275 |
| −0.5 | 0.208333 | 0.144069 | 0.144% | 11.004% | 0.266789 | 0.941579 | 3.858281 | 2.7927 |

The results are dimensionless in the radial unit $L_0$. For a distance guarantee use the unrounded threshold in the numerical results. At long wavelength the massless mode controls the response on every leaf. As the wall moves from one candidate position to another, the range where the continuum is important changes substantially, even with the same five-dimensional action. The three calculations use separately stationary walls. Their differences alone do not track the excitation of continuum photons by a moving wall.

## Numerical checks

The script propagates the dimensionless logarithmic derivative $r_Q=-u'_Q/(Qe^{-A}u_Q)$ backward from optical depth $Q\int_b^{y_{\rm end}} e^{-A}dy=28$ with $r_Q(y_{\rm end})=1$. It independently propagates the full second-order equation twice from depth 19, once with a Dirichlet endpoint and once with a Neumann endpoint. At exact arithmetic those finite-domain results bracket the infinite-domain kinetic coefficient. A Neumann end minimizes over unrestricted endpoint data, whereas the Dirichlet solution can be extended by zero into the full domain. Across all 36 reported momenta the largest disagreement between the main result and either finite-domain solver is below $5.3\times10^{-13}$ relative. The differences include integration error and do not constitute interval-arithmetic certificates.

An independent quadrature recovers the analytic zero-mode norm. The $Q=0.003$ results recover the low-momentum quadratic coefficient to within 0.001%. At the other end, a local expansion of the logarithmic derivative gives $D(Q)/w(b)=Q+[A'(b)+(\log w)'(b)]/2+O(1/Q)$, consistent with direct propagation at $Q=512$. The program checks spectral positivity, the all-momentum bound and monotonicity on every leaf. [Solver](euclidean.py) · [Numerical results](EUCLIDEAN.json)

## A gauge weight selected by separation of the null gravitational source

A second test uses the weight proportional to the warp squared that is selected by the restricted null-field uplift. The static transverse probe is still calculated on each fixed flat leaf. To keep a single microscopic law across leaves, define $A_{\rm ref}(0)=0$ and set

$$
G_{\rm glob}(\phi)=e^{2A_{\rm ref}(y)}
=(1+\phi)^{5/9}(1-\phi)^{-1/9}e^{-\phi^2/9}.
$$

On a single leaf it is convenient to divide $G$ by its wall value. Then $w=e^{2A}$ and $p=e^{4A}$ with $A=A_{\rm ref}(y)-A_{\rm ref}(b)$. This overall normalization cancels in $\eta$. When comparing the coupling between leaves, the factor $G_{\rm glob}(\phi_b)$ must be restored. With this fixed scalar law, the log change in $I$ starting at $\phi_b=0$ is $+0.240890$ at $+0.5$ and $-0.282456$ at $-0.5$. These are differences among stationary geometries. The closed cosmological launches have their own initial field displacements.

The pole survives and the continuum remains gapless. Its low-momentum scaling now changes. In wall-normalized units $T/w\to 1/(2k)$ at large $y$, where $k=5/9$. The earlier correction integral $C$ diverges because $T^2/p\to1/(4k^2)$. We obtain

$$
\frac{D(Q)}{Q^2}=N-\frac{Q^2}{4k^3}\log(1/Q)+O(Q^2),\qquad
\eta(Q)=\frac{Q^2}{4k^3N}\log(1/Q)+O(Q^2).
$$

Here $Q$ is in inverse radial units. Changing the reference scale inside the logarithm shifts the finite $O(Q^2)$ term.

To see the logarithm and evaluate it without cancellation, define $t=Qe^{-A}$, $M=-pu'_Q/(Q^2 u_Q)$, $U=T/w$ and $B=(T-M)/Q^2$. Direct differentiation gives
$$
U'=-1-2A'U,\qquad B'=-(U-t^2 B)^2 .
$$
The decaying solution has $0\le M\le T$, by its integrated equation and positivity. At the point $y_Q$ where $t=1$, $0\le B(y_Q)\le U(y_Q)$ is bounded independently of $Q$. On $b\le y\le y_Q$ we have $-B'\le U^2$ and $B(y)\le \mathrm{const.}\,(1+y_Q-y)$. The difference $U^2-(U-t^2B)^2$ is bounded by $2Ut^2B$. Its integral remains bounded, since $t^2$ decreases exponentially with $y_Q-y$ throughout this interval. Also $U^2$ approaches $1/(4k^2)$ exponentially and $y_Q=k^{-1}\log(1/Q)+O(1)$. Thus $B(b)=(4k^3)^{-1}\log(1/Q)+O(1)$. The normalization $w(b)=1$ is used in this derivation.

| $\phi_b$ | $N$ with $w(b)=1$ | $N$ for the fixed global weight | Coefficient of $Q^2\log(1/Q)$ in $\eta$ | Extra share at $Q=0.1$ | Extra share at $Q=1$ | $Q$ at 1% | $Q$ at 10% |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 1.035771 | 1.035771 | 1.407647 | 2.842% | 42.721% | 0.051659 | 0.241892 |
| +0.5 | 1.274341 | 1.676865 | 1.144121 | 2.694% | 45.610% | 0.054324 | 0.241816 |
| −0.5 | 0.930667 | 0.588743 | 1.566618 | 2.976% | 41.949% | 0.049812 | 0.237872 |

The direct solver and the subtracted system agree on their overlapping momenta. Between $Q=10^{-7}$ and $10^{-8}$, the numerical slopes of $\eta/Q^2$ against $\log(1/Q)$ agree with the predicted coefficients to within $2\times10^{-11}$ relative. Thus the same background can give different static force responses under different gauge laws. Neither low-momentum result by itself fixes how many photons emerge from a time-dependent source.
