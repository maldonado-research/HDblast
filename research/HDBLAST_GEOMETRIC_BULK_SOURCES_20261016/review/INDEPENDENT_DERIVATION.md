# Independent quadratic-geometry derivation

16 October 2026. Internal mathematical referee derivation. The producer
implementations are reviewed separately; the identities in this section do
not follow from finite tests.

Write g(e)=g0+e q and phi(e)=phi0+e eta. Every Q is the second parameter
derivative at zero, twice the e^2 coefficient. The reference connection is
the Levi-Civita connection of g0. Set

    P_abc = (nabla_b q_ca+nabla_c q_ba-nabla_a q_bc)/2,
    P^a_bc = g0^ad P_dbc,
    J^a_bc = -2 q^a_d P^d_bc,
    C1^a = g0^bc P^a_bc.

The exact connection difference is K(e)^a_bc=e g(e)^ad P_dbc. The inverse
identity g(e)^-1=g0^-1-e q^sharp+e^2 q^a_c q^cb+O(e^3), obtained by
multiplying by g(e), proves K=eP+e^2 J/2+O(e^3). The exact coordinate Ricci
formula, regrouped using the reference connection without commuting any
derivatives, is

    Ric(g)_ab = Ric(g0)_ab + nabla_c K^c_ab-nabla_b K^c_ac
                +K^c_cd K^d_ab-K^c_bd K^d_ac.

Taking two parameter derivatives proves, for arbitrary smooth q,

    QRic_ab = nabla_c J^c_ab-nabla_b J^c_ac
               +2(P^c_cd P^d_ab-P^c_bd P^d_ac).

This establishes the indexed formula universally. No background Einstein
equation, gauge condition, angular restriction, or normal-coordinate choice
is needed. The derivative of J uses the tensor product rule. In computing
nabla_e nabla_d q_ab, the inner derivative index d is a covariant tensor
index: the term -Gamma^f_ed nabla_f q_ab is required. A noncoordinate frame
does not allow its omission. All time-index raising uses the Lorentz metric.

For E_ab=Ric_ab-phi_a phi_b-(2/3)V(phi)g_ab, ordinary product differentiation
gives

    QE_ab = QRic_ab-2 eta_a eta_b
             -(2/3)[V'' eta^2 (g0)_ab+2 V' eta q_ab].

The coefficient of the kinetic term is two. The coefficient of the potential
cross term inside brackets is two. There is no V times a second metric
variation on the affine path. The corresponding scalar residual S=Box_g phi
-V'(phi) has general affine Hessian

    QS = 2 q^a_d q^db nabla_a nabla_b phi0
         -2 q^ab nabla_a nabla_b eta
         +2 q^ab P^c_ab (phi0)_c -V''' eta^2
         +2 q^c_d C1^d (phi0)_c-2 C1^c eta_c.

This follows directly by expanding the inverse metric and the connection in
g^ab(partial_ab phi-Gamma^c_ab partial_c phi). The final two terms vanish
for the actual first harmonic field; rectangular independent jet inputs do
not automatically satisfy that premise.

## The two quadratic gauge covectors

The varying-metric lowering is especially useful. Exactly,

    C_b(e)=g(e)_ba C^a(e)=e g(e)^cd P_bcd.

Thus QClower_actual,b=(C_b)''=-2q^cd P_bcd for arbitrary q. By contrast,

    g0_ba QCupper^a = -2q_bd C1^d-2q^cd P_bcd.

Their difference is 2q_bd C1^d. Merely using the same lower index does not
make them the same source on general jets. For the registered reduction
Ered=E-nabla^g_(a C_b), with symmetrization by one half, the product rule gives

    QEred_general,ab = QE_ab-nabla^0_(a QClower_actual,b)
                       +2 P^c_ab (C1)_c.

This general identity needs only q and its covariant two-jet. If C1=0 as a
field, it reduces to QEred=QE-nabla^0_(a QC_b), and actual and background
lowering agree together with their covariant derivatives. C1=0 at a single
point alone does not imply equality of those two derivatives. The root note
GAUGE_COVECTOR_BOOKKEEPING.md uses precisely these conventions and is correct.

## Independent coordinate implementation and its scope

independent_coordinate_oracle.py constructs actual rational coordinate
two-jets of g0+e q on (y,t,x1,x2,x3). The sphere embedding is

    n=(2x1,2x2,2x3,1-|x|^2)/(1+|x|^2),
    kappa_ij=4 delta_ij/(1+|x|^2)^2,
    Y=n0^2-n1^2.

The background is dy^2+rho(y)^2[-dt^2+c(t)^2 kappa]. At the evaluation
point rho=2, rho_y=3/7, rho_yy=5/11, c=5/4, c_t=3/4, c_tt=5/4; these are
genuine jets of a smooth warped metric with c=cosh(t+t0). All seven metric
amplitudes and their ordinary y,t derivatives through order two are populated
independently. Their coordinate representatives are

    qyy=A Y, qyt=rho B Y, qyi=rho c C partial_iY,
    qtt=rho^2 D Y, qti=rho^2 c E partial_iY,
    qij=rho^2 c^2(L Y kappa_ij+M Hess_kappa(Y)_ij), eta=sY.

A small value/gradient/Hessian algebra expands the coordinate inverse and
coordinate Christoffel formula through e^2, then forms Ricci, the matter
residual, the upper and actual-lowered gauge vector, the full varying-
connection reduction and the scalar residual. A separate indexed path
constructs covariant q jets and applies the identities above. This is an
independent implementation witness for curved, nonnormal, nonharmonic
families. It is not claimed as a universal proof. Negative controls remove
the derivative-slot connection, reverse time raising, omit the general
2PC1 correction, or confuse the two lowerings.

## Angular completeness and moving shell

An all-point projection must reconstruct every reduced polynomial component
modulo n0^2+n1^2+n2^2+n3^2-1. A zero remainder proves equality on the whole
sphere. An invertible rational matrix of monomial coefficients proves
independence of the selected scalar, vector, and tensor features; it is not
collocation. Formal independence of all radial/time two-jet symbols then
proves the source identities for all such jets, rather than just substituted
manufactured values. Both chiral TT tensors must remain in the tensor basis.

The fixed-material-shell gauge pullback is a separate requirement from this
bulk calculation. A spacetime gauge transformation must pull back the shell
embedding as well as the bulk fields. In a fixed-shell representative the
normal first gauge vector vanishes on the shell, so eta|shell=vY when the
common scalar variation is eta=(v+pf)Y and f|shell=0. Tangential derivatives
of f|shell vanish too. Bulk QRic/QE must not be substituted for physical shell
Hessians; the shell normal, induced metric, and embedding variations remain
part of the existing shell construction.

The optional intrinsic scalar I=gamma^mu nu d_mu(phi|shell)d_nu(phi|shell)
has I0=I1=0 because phi0 is constant along the background shell. Therefore
I''=2 gamma0^mu nu d_mu eta_material d_nu eta_material and is invariant under
joint field/embedding gauge changes and shell relabellings. With eta_material
=vY its normalized sphere mean is rho_b^-2[-vt^2/3+8v^2/(3c^2)]. This is a
kinematic quadratic diagnostic, not a forcing or second-solution result.
