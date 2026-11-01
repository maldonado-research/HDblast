# Independent audit: first-only shell curvature defect and Buchert coefficient

This audit responds to the parallel shell-expansion derivation. It uses
standard integrated scalar-curvature variation and averaged hypersurface
identities, applied to the retained shell first field. It is not a new
generic theorem or a claim of new physical dynamics.

Let the shell proper-time slice metric be g_e, g_0=A²kappa, with A=rho_b
cosh(t), proper time tau=rho_b t and H=A'/A. All primes in this note are
proper-time derivatives. Let V=integral dV, S=integral R3 dV and

    C(g)=S/V-6(2pi²/V)^(2/3).

At the round metric Ric=2g_0/A². On the closed three-sphere the integrated
scalar-curvature first variation is

    delta S=integral [-Ric+(R3/2)g_0]:delta g dV
           =(1/A²)integral tr(g_0^{-1}delta g)dV,
    delta V=(1/2)integral tr(g_0^{-1}delta g)dV.

Consequently delta C=delta S/V−2delta V/(A²V)=0 for every symmetric metric
variation. Thus along a C2 family C'' depends only on g_1, not on g_2.
Diffeomorphism invariance also kills every Hessian contribution to g_1:
C(g_0)=0, DC(g_0)=0, and differentiating DC_g[Lie_Xg]=0 gives
D²C(g_0)[v,Lie_Xg_0]=0 for every v. This handles cross terms as well as
the Hessian square.

In the anchored Gaussian shell representation, the first metric is

    g_1=2A²[EY kappa-J Hess_kappa(Y)],
    E=u-HP, P'=u, P(0)=0,
    J'=P/A², J(0)=0.

The second slot is Lie_{-J gradY}g_0 and drops out of C''. We can use the
conformal path g_e=exp(2eEY)g_0. The conformal scalar-curvature law and an
integration by parts give

    V_e/V_0=<exp(3eEY)>,
    S_e/V_0=A^{-2}<exp(eEY)[6+2e²E²|DY|²]>.

With <Y>=0, <Y²>=1/6 and <|DY|²>=4/3,

    V''_0/V_0=3E²/2,
    S''_0/V_0=19E²/(3A²),
    (<R3>)''_0=−8E²/(3A²),
    (6/Rvolume²)''_0=−6E²/A²,
    C''_0=10E²/(3A²).

The individual conformal-path derivatives in the middle are only a
convenient computational path. Their difference is independent of the
second completion. Other paths with the same first variation can change
the two individual derivatives together.

For the Buchert quantity use extrinsic expansion K=(1/2)g^{-1}g',
theta=tr K and shear=K−theta I/3, with the squared norm of shear containing
**no factor one half**. Define

    Q_B=(2/3)(<theta²>−<theta>²)−<|shear|²>.

The background theta_0=3H is spatially constant and shear_0=0. Therefore
all second metric, second observer and measure variations cancel in
Q_B''; only the first theta and shear occur. Indeed for B=g_0^{-1}g_1,
K_1=B'/2 and

    theta_1=(3E'+8J')Y,
    shear_1=−J'[HessY+(8/3)Y I],
    <|HessY+(8/3)YI|²>=8−(64/3)(1/6)=40/9.

Thus

    Q_B''=2E'^2+(32/3)E'J'+(16/3)J'^2.

Since H'=1/A², put W=u'−Hu and j=P/A²; E'=W−j and J'=j give

    Q_B''=2W²+(20/3)Wj−(10/3)j².

No universal sign for this polynomial is asserted. Numerical signs must
use the actual jointly enclosed u,u',P. A free comparator alone does not
supply them.

The exact Lorentzian Gauss identity in this convention is
R4=R3+2theta'+theta²+|K|². Moving-volume averaging satisfies
<theta'>=(<theta>)'−<theta²>+<theta>². With H_D=<theta>/3, substitution
therefore gives

    6(H_D'+2H_D²)=<R4>−<R3>+Q_B.

At the background this gives R4=12/rho_b², checking the curvature sign.
This relation does not make each of <R4> or <R3> separately first-only;
only the critical defect C and vanishing-background variance/shear
quantity have the cancellations established above.

The Gaussian first-field relation also cancels shell-preserving gauge
amplitudes in the old proper-clock source: its X=u+b_t,
Z=u+tanh(t)b give Cclock=integral X d tau=P+rho_b b (b(0)=0), hence
Z−H Cclock=u−HP. The angular Hessian mode is a spatial pullback. This
uses the same physical initial section and proper-clock prescription.

`check_curvature_buchert.py` checks these reductions independently with
exact symbolic arithmetic and no source/producer imports. Normal and
optimized reports agree byte for byte. Conditional exact-family existence,
regularity, and the source-bound x<3/2 clock domain remain as recorded in
`DERIVATION.md`; the audit supplies no effective finite parameter radius.
