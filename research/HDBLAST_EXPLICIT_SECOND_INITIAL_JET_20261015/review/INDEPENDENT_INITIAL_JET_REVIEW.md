# Independent review of the compact second Cauchy jet

Status: the explicit second **initial constraint jet** is valid in the physical
Cauchy convention below. This separate review changes no frozen capsule and
does not establish a second evolution, gauge completion, nonlinear remainder,
or an exact nonlinear initial-data family.

## Conventions and the actual first constraints

Use spatial dimension four and

\[
 \mathcal H=R[\gamma]+K^2-|K|^2-\Pi^2-|d\phi|^2-2U(\phi),
 \qquad \mathcal M_i=D_j(K^j{}_i-\delta^j_iK)-\Pi\,\partial_i\phi.
\]

Here `K=(2N)^-1 gamma_t` when the shift vanishes, and `Pi=N^-1 phi_t`.
These are a covariant tensor and scalar velocity, not canonical momentum densities.
Changing the sign of K requires the corresponding change to the momentum
convention and scalar momentum. At `t=0`,

\[
 \gamma_0=dy^2+\rho^2\kappa_{S^3},\quad K_0=0,\quad\Pi_0=0,
\]
\[
 \gamma_{1,yy}=-4uY,\quad\gamma_{1,yi}=0,\quad
 \gamma_{1,ij}=2uY\gamma_{0,ij},\quad\phi_1=vY,
\]
\[
 K_{1,yy}=-2u_tY/\rho,\quad K_{1,ij}=u_tY\gamma_{0,ij}/\rho,
 \qquad\Pi_1=v_tY/\rho.
\]

The retained aligned reconstruction gives

\[
 v=-\frac3p(u_y+2au),\qquad
 v_t=-\frac3p(u_{yt}+2au_t),
 \quad a=\rho_y/\rho,\quad p=\phi_{0,y}\ne0.
\]

This follows directly for the specified actual packet: with `C0=alpha/rho^2`,
`B=3C0/(rho p)`, and `k=alpha_z/alpha`, one has
`C0_z/C0=k+2rho a` and `d_y=-rho^-1 d_z`. Substitution into
`u=C0 Q`, `u_t=C0 e^-z Q` reproduces the stated `v` and `v_t` exactly.

The background identities used to check the first constraints are

\[
 p_y+4ap=U',\qquad a_y=-\rho^{-2}-p^2/3.
\]

The background Hamiltonian itself additionally uses the retained radial first
integral `a^2=rho^-2+p^2/12-U/6`. Background momentum is zero because
`K0=Pi0=0`. The first integral is not needed again in the first-variation
cancellations below.

Direct linearization gives

\[
 \mathcal H_1=\big[-6u_{yy}-36au_y
 -(24a_y+48a^2+12\rho^{-2}+4p^2)u-2pv_y-2U'v\big]Y,
\]
\[
 \mathcal M_{1,y}=\rho^{-1}[-3(u_{yt}+2au_t)-pv_t]Y,
 \qquad \mathcal M_{1,i}=0.
\]

Substituting the reconstruction and first background identity reduces the
Hamiltonian row to `-4u(3a_y+3rho^-2+p^2)Y=0`. The momentum row vanishes
by the time derivative of the reconstruction. Thus the actual first data
used by the construction satisfy both constraints; arbitrary independent
interval selections need not do so. The first spatial reconstruction alone
is sufficient for this initial Hamiltonian calculation; the first evolution
equation is not needed here.

## The Hamiltonian source and its nonzero witness

Let `S_H` denote the second parameter derivative along the physical Cauchy affine
path `(gamma0+eps gamma1, eps K1, phi0+eps phi1, eps Pi1)`. It is

\[
 S_H=A_HY^2+\rho^{-2}(36u^2-2v^2)|\nabla_\kappa Y|^2,
\]
\[
\begin{aligned}
 A_H={}&-24u_y^2-288auu_y-24uu_{yy}
 -(192a_y+384a^2+32p^2)u^2\\
 &-\rho^{-2}(336u^2+12u_t^2+2v_t^2)
 -2v_y^2-16puv_y-2U''v^2.
\end{aligned}
\]

This agrees with my earlier independent affine-spacetime calculation. The
two Hamiltonian Hessians agree because background K and Pi are zero.

One derivation uses the radial foliation of the spatial affine metric. Set
`A=1-4eps psi`, `C=1+2eps psi`, `L=sqrt(A)`, `sigma=log(C)/2`, and
`b=(a+sigma_y)/L`, with `psi=uY`. The exact radial Gauss identity is

\[
 R^{(4)}=R[\rho^2C\kappa]-12b^2-6L^{-1}b_y
                  -2L^{-1}\Delta_{\rho^2C\kappa}L.
\]

Using the three-dimensional conformal curvature and Laplacian identities,
its Hessian is

\[
\begin{aligned}
 R''={}&-24\psi_y^2-288a\psi\psi_y-24\psi\psi_{yy}
 -(192a_y+384a^2)\psi^2\\
 &+\rho^{-2}(48\psi^2+48\psi\Delta_\kappa\psi
                       +36|\nabla_\kappa\psi|^2).
\end{aligned}
\]

Since `Delta Y=-8Y`, the coefficient `48-384=-336` follows. The remaining
terms follow directly by expanding the inverse spatial metric in
`|dphi|^2`; the kinetic Hessians are `-12u_t^2Y^2/rho^2` and
`-2v_t^2Y^2/rho^2`. The independent checker expands the exact Gauss formula
in rational epsilon series, rather than importing the proposed source.

At `t=0,z*=3ell/2,n=(1,0,0,0)`, `Y=1`, `|grad Y|=0`, and
`Q=1,Q_z=0,Q_zz=-32/ell^2`. Write `C0=alpha/rho^2`,
`B=3alpha/(rho^3p)`, `k=alpha_z/alpha`. Then

\[
 u=C_0,\quad u_y=-C_{0,z}/\rho,\quad
 u_{yy}=(C_{0,zz}-32C_0/\ell^2)/\rho^2+aC_{0,z}/\rho,
 \quad u_t=C_0e^{-z_*},
\]
\[
 v=Bk,\quad v_t=B(k-1)e^{-z_*},\quad
 v_y=32B/(\rho\ell^2)-(B_zk+Bk_z)/\rho.
\]

The authenticated actual coefficient report gives
`2.3e9<|v_y|<3e9`. It also gives the conservative bounds
`|u|<1e-4, |u_y|<1, |u_yy|<1e6, |u_t|<1, |v|,|v_t|<4`,
`|a|,|a_y|,|p|<1`, and `|U''|<4`. All terms of `A_H` other than
`-2v_y^2` have combined absolute value below `5e6`. Exact rational
arithmetic therefore proves

\[
 \boxed{S_H(t=0,z_*,n)<-10^{19}.}
\]

The witness excludes naive zero second initial data. It is not an estimate
of a solved second field.

## Local Hamiltonian and momentum inverses

For a scalar harmonic `Y_l` with `lambda=l(l+2)`, choose

\[
 \gamma_{2,yy}=\gamma_{2,yi}=0,\quad
 \gamma_{2,ij}=f_l\gamma_{0,ij}Y_l,\qquad
 \phi_2=-3f_l'Y_l/(2p).
\]

The spatial curvature variation is
`-3f_l''-12af_l'+(2lambda-6)f_l/rho^2`. The scalar variation contributes
`-2p phi2_y-2U'phi2=3f_l''+12af_l'`, using `p_y+4ap=U'`. Hence

\[
 L\mathcal H[\gamma_2,\phi_2]
       =2(\lambda-3)f_lY_l/\rho^2.
\]

For `lambda=0,8,24`, the respective denominators are `-6,10,42`, all
nonzero. If `S_H=sum S_H,l Y_l`, choose

\[
 f_l=-\rho^2S_{H,l}/[2(\lambda_l-3)].
\]

No division by the excluded `lambda=3` sector occurs:
`Y^2` has only degrees `0,2,4`, and
`|gradY|^2=(Delta(Y^2)+16Y^2)/2` has those same degrees.

For the **physical Cauchy affine** momentum source, independent coordinate
divergence gives

\[
 S_{M,y}=G Y^2,\quad S_{M,i}=F\,\partial_i(Y^2),
\]
\[
 G=(12uu_{yt}-48auu_t-6u_yu_t-2v_tv_y)/\rho,
 \qquad F=(18uu_t-vv_t)/\rho.
\]

The proposed correction uses the full four-dimensional spatial trace:

\[
 K_2=q\gamma_0,\quad q=FY^2/3,\qquad
 \Pi_2=(G-F_y)Y^2/p.
\]

Indeed `L M[K2,Pi2]=-3dq-p Pi2 dy=-S_M`. The factor three is spatial
dimension minus one. This K2 includes its `yy` component; it is not merely
a trace on the three-sphere. Because `K0=Pi0=0`, the Hamiltonian correction
does not enter linear momentum, and the K2/Pi2 correction does not enter
linear Hamiltonian. Thus both inverses can be combined without interference.

## Affine spacetime data are a different momentum convention

Along `g0+eps h,phi0+eps chi` with zero shift, the lapse at `t=0` is
`N=rho sqrt(1+2eps psi)`. This induces the physical Cauchy second derivatives

\[
 K_{2,\mathrm{affine\ spacetime}}=-2\psi K_1,\qquad
 \Pi_{2,\mathrm{affine\ spacetime}}=-2\psi\Pi_1.
\]

On the first-constraint data, their linear momentum contribution is
`+6u_yu_tY^2/rho` in the radial row and zero in the angular rows. Thus

\[
 G_{\mathrm{affine\ spacetime}}
 =(12uu_{yt}-48auu_t-2v_tv_y)/\rho,
\]

whereas the physical Cauchy G above includes `-6u_yu_t/rho`. Both formulas are
correct in their respective conventions. A report must not interchange them.

To realize the completed **physical Cauchy** second jet in spacetime coordinates
with zero shift, the initial coordinate time derivatives must obey

\[
 \gamma_{2,t}=2\rho(K_2+2\psi K_1),\qquad
 \phi_{2,t}=\rho(\Pi_2+2\psi\Pi_1).
\]

The second lapse does not contribute here because the background time
derivatives vanish. Additional lapse/shift/gauge choices are not solved by
the constraint inversion.

Throughout, the family convention is `X=X0+eps X1+eps^2 X2/2`. The sources
are Hessians, and the second equations are `L X2=-S`, with no extra half.
If one instead calls the coefficient of `eps^2` the second field, every
second coefficient and source equation must be changed consistently.

## Support and verification scope

The retained bump and all its derivatives vanish outside `[ell,2ell]` and
are flat at its endpoints. Every source and every correction above is a
local expression in those first data and smooth background coefficients;
`p` stays bounded away from zero. Thus the second initial data are compact
and flat at the same support endpoints, and vanish in a neighborhood of the
physical shell. The construction supplies the initial physical constraints;
it permits nonzero second scalar displacement and velocity and therefore
changes any prior convention that selected those data to be zero.

`check_initial_jet_review.py` independently verifies 75 exact seeded check
groups: Gauss-formula Hamiltonian Hessians, actual first constraints,
physical Cauchy momentum Hessians, the affine-spacetime lapse conversion,
Hamiltonian inverses for all three eigenvalues, momentum inversion, and the
finite angular sectors. It authenticates the actual coefficient report and
reproduces the strict witness. Normal and optimized Python agree. These
finite checks supplement the all-point algebra above; they are not a proof
of evolved constraints, a gauge solve, or a nonlinear initial family.

I also read the complete producer note `EXPLICIT_COMPACT_SECOND_INITIAL_JET.md`
and replayed its checker independently: 48 exact checks and seven intentional
errors pass in each Python mode. The seven errors are rejected for their
specified reasons and leave no output report. No mathematical objection
remains to publishing the restricted initial-constraint result with the
physical Cauchy conventions and stated limitations retained.
