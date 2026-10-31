# A conditional proper-time Gaussian-screen curvature coefficient

For the exact family selected in the [exact-family TT binding](../HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md), launch future unit-normal geodesics from its fixed initial physical hypersurface and retain the initial radial and angular labels. At the fixed proper time and initial radial label defined below, the selected magnetic screen projection has second parameter derivative

\[
\boxed{\quad
h_{\rm geo}\in
\left[
-\frac{2099308763896503245797815094104589}
       {5192296858534827628530496329220096},
-\frac{37909098171726239036888826682085}
       {2596148429267413814265248164610048}
\right]
\subset[-0.404312161091,-0.014602053466].\quad}                 \tag{1}
\]

The decimal endpoints are rounded outward. The exact interval comes from the [complete numerical report](numeric/NORMAL.json), including the boosted metric response, the full scalar quadratic curvature, and the change of observer and screen direction. This is a new evaluation of the specified proper-time observation. Its sign is established by this complete sum.

All coordinates and elapsed times below use the retained model normalization, with `ell=10^-6`. No conversion to an SI length or time scale is supplied.

## Conditional statement and retained family

Retain the same five-dimensional Einstein–scalar model, regular branch, canonical initial curve, relative-harmonic completion, and physical shell conditions as the exact-family binding. In its conventions,

```
Ric_AB = phi_A phi_B+(2/3)V(phi)g_AB,
Box_g phi = V'(phi),
K_out,ab = sigma(phi)h_ab/6,
n_out(phi) = -sigma'(phi)/2,

W(phi)=1-phi+phi^3/3,
V=(W')^2/2-(2/3)W^2,
sigma=2W+delta(1+c_sigma phi),
delta=1/1000,
c_sigma=96714321425692/161839258930731.
```

The background is the selected near-`phi=+1` member,

\[
g_0=dy^2+\rho(y)^2[-dt^2+\cosh(t)^2\kappa_{S^3}],
\qquad a=\rho_y/\rho,
\qquad dy=-\rho\,dz.
\]

Use `x=t/ell`, `r=z/ell`; r increases inward. Keep the parameter convention
`g_e=g0+e q+e² k/2+o(e²)`. Thus k and the coefficient in (1) are second parameter derivatives themselves.

The analytic premises remain those of the retained localized constraint theorem, prescribed-slab shell evolution theorem, actual master/gauge solutions and their differentiated error certificates, and exact-family TT identification. In particular, retain the existence of the selected exact family with nonlinear regularity `X14` and C2 parameter dependence in `X10`, which controls the lower piecewise spacetime C4 topology used here. The selected initial correction is supported in
`ell/2<z<5ell/2`, and its source-bound common coordinate slab is

\[
0\le t\le3\ell/2,\qquad 0\le z\le6\ell.
\]

The exact-family proof identifies q with the complete actual common-gauge first solution and the two selected TT components of k with the equal scalar response `U_T`. It proves their zero TT initial displacement and velocity, although other second initial sectors need not equal those of an earlier finite-angular construction. The same family and identification are used here. None of these analytic existence or regularity premises is replaced by source authentication or finite algebra checks.

Under these premises there is a common, possibly ineffective, open parameter interval on which the Gaussian observation below exists, has a nondegenerate screen, and defines a C2 scalar functional `J_geo(e)`. The result is

\[
J_{\rm geo}(0)=J'_{\rm geo}(0)=0,
\qquad J''_{\rm geo}(0)=h_{\rm geo},                         \tag{2}
\]

with (1). Consequently, after decreasing that parameter interval if necessary, every sufficiently small nonzero e satisfies

\[
J_{\rm geo}(e)\le {u\over4}e^2<0,
\qquad u=-0.014602053466.
                                                               \tag{3}
\]

No effective lower bound on the allowed parameter radius is claimed.

## Fixed anchors, proper clock, and screen axes

Let `Sigma={t=0}` be the initial physical hypersurface in the selected common chart. Its initial embedding, radial labels `z_i`, angular labels, outward orientation, round probability measure, and fixed tensors `Tplus,Tminus` are supplied anchor data. Let `u_e` be the future physical unit normal to Sigma. Define

\[
F_e(\tau,z_i,\Omega)
 =\exp^{g_e}_{i(z_i,\Omega)}(\tau u_e),
\qquad U_e=\partial_\tau F_e.
\]

The geodesic equation and initial unit normalization imply `g_e(U_e,U_e)=-1`; therefore tau is elapsed proper time on every trajectory. Before crossings or focal points, Gauss' lemma gives
`F_e^*g_e=-d tau²+gamma_e`. On a tau slice use the outward spatial unit normal to `z_i=constant`, namely `-grad_gamma(z_i)/|grad_gamma(z_i)|`. The minus sign is required because the initial radial label increases inward.

The screen is tangent to `tau,z_i=constant`. Its axes are obtained from the angular connecting vectors `partial_Omega F_e` by the positive symmetric inverse Gram square root relative to the fixed initial round frame. This defines a Gaussian-label screen convention; it imposes no parallel-transport condition on those axes.

Use the retained curvature sign
`R^a_bcd=partial_c Gamma^a_db-partial_d Gamma^a_cb+...` and outward sphere orientation `EPS[0,1,2,3]=+1`. The magnetic screen tensor is

\[
H_{AB}(e)=\tfrac12\epsilon(s_e)_{(A}{}^{CD}
                    C(g_e)_{B)U_eCD},
\]

with normalized symmetrization and normalized screen components `Hhat_e`. At the fixed anchor pair `(tau_dagger,z_i_dagger)` set

\[
J_{\rm geo}(e)=3\,\operatorname{mean}_{S^3}
       \langle\widehat H_e,T_+-T_-\rangle_\kappa.             \tag{4}
\]

The integration measure and templates in (4) are fixed initial data. Their probability norms are `||Tplus||²=||Tminus||²=1/6`, with zero cross inner product. Thus (4) extracts the coefficient of `Tplus-Tminus`. Reversing the specified orientation reverses the coefficient and the associated scalar functional.

## Background shooting and the common family domain

Write `e0=rho^-1 partial_t`, `ey=partial_y`, and

\[
U_0=\cosh\chi\,e_0+\sinh\chi\,e_y,
\qquad n_0=\sinh\chi\,e_0+\cosh\chi\,e_y.
\]

The background geodesic launched at rest from `z_i` obeys, with coordinate t as independent variable,

```
dz/dt=-tanh(chi),
dchi/dt=-rho a,
dtau/dt=rho sech(chi),
z(0)=z_i, chi(0)=tau(0)=0.
```

Define `z_i_dagger` to be the unique actual-background shooting root for
`z(5ell/4;z_i)=ell/8`, and define `tau_dagger` by the last integral along that root. Freeze these exact implicit quantities for the whole parameter family. Neither an interval endpoint nor a floating-point representative is used as launch or clock data. The background event is exactly `(x,r)=(5/4,1/8)`.

The [background proof](background/DERIVATION.md) and [rational receipt](background/BACKGROUND_NORMAL.json) enclose the root, rapidity, clock, and Jacobi fields for the same retained actual branch. They use `14<rho a<15`, `-216<(rho a)_z<0`, elementary integral inequalities, and strict shooting monotonicity. If
`j=partial_zi z|t`, then the background identities include

\[
\det {\partial(t,z)\over\partial(\tau,z_i)}
 ={j\cosh\chi\over\rho}>0,
\qquad
\gamma_{z_i z_i}=\rho^2j^2\cosh^2\chi>0,
\qquad
\gamma_{AB}=(\rho\cosh t)^2\kappa_{AB}.
\]

Positive `j` gives uniqueness of the radial launch at fixed coordinate time, while positive `tau_t` gives uniqueness along each trajectory. Angular labels are unchanged in the background. This proves injectivity on the bounded background tube, in addition to the local determinant and screen-Gram margins.

The auxiliary background certificate covers `tau<=1/5000` and `x<=8/5`. For the selected exact family, apply it only on the smaller tube

\[
0\le\tau\le1/5500,\qquad
\ell/16\le z_i\le3\ell/16,\qquad\Omega\in S^3.              \tag{5}
\]

The [family-domain join](background/FAMILY_DOMAIN_JOIN.md) proves

\[
\tau_\dagger<131/800000<1/5500,
\qquad
x\le {62500000000000\over44343749987229}<10/7<3/2,
\]

throughout (5), with enclosing radial range `1/32<=r<=1/4`. Thus the entire chosen observer tube lies strictly inside the already supplied family slab and inside the actual first-field domain. The larger background overcover supplies no extension of the exact family through `x=8/5`.

ODE dependence on the parameter and initial labels preserves these strict clock, inverse-map, Gram, spatial, and upper-time margins for a common sufficiently small parameter interval. Compactness and C1 stability of the background embedding preserve injectivity. At the initial face this is a forward, one-sided construction; no negative-time extension is needed. The retained C2 family in spacetime C4 supplies the derivatives required for the geodesic flow, connecting vectors, curvature evaluation, and angular integration in (4). In particular `J_geo` is C2.

Every exact initial datum equals the background near the whole label band in (5): that band lies below both the first-seed support `[ell,2ell]` and the full correction support beginning at `ell/2`. The prescribed initial lapse and shift also equal their background values there. Hence the fixed initial embedding and initial unit normal have exactly zero first variations, including label derivatives. This establishes the zero initial conditions used in the first-position estimate below.

## Complete curvature reduction

The background Weyl tensor vanishes. For the complete scalar first field, the full first symmetric screen magnetic tensor vanishes for both radial basis normals `e0` and `ey`, and hence for every radial boost. Codazzi reduces their scalar terms to the symmetric curls of `Y kappa` and `Hess Y`, plus isotropic gradient terms; each is zero. This proves both vanishings in (2) before any sphere integration. It does not assert that the full first Weyl tensor is zero.

For the second variation compare the Gaussian frame with the metric-dependent registered frame boosted through the fixed background rapidity. Its complete first magnetic tensor is identically zero. Consequently the moving-event correction, radial boost correction, and relative in-screen rotation act on that zero tensor. Second frame and event variations multiply the background Weyl and vanish. The surviving relative changes are the two angular tilts of the observer and screen normals acting on the first Weyl tensor.

The [full geometric derivation](geometry/DERIVATION.md) retains the baseline connection, normal, epsilon raising/volume, and positive-Gram changes. Put `R=rho cosh(t)`, `theta=tanh(t)`, `c=cosh(chi_dagger)`, `s=sinh(chi_dagger)`. With ordinary physical derivative suffixes, define

```
F_t = ((D+4M-L)M_t-2M L_t)/(4rho)-E(D+L)/(2R)
      +(C/R-M_y/2)[R(C_t-theta C)/(2rho)-R E_y/2+B/2],

eta = R(C_t-theta C)/(2rho)-R E_y/2-B/2,
F_y = ((-A+4M-L)M_y-2M L_y)/4+C(A-L)/(2R)
      +(E/R-M_t/(2rho))eta,

mathcal A_t = -(M_tt+theta M_t)/(2rho²)
              +E_t/(rho R)-(D+L)/(2R²),
mathcal A_y = -(M_yy+a M_y)/2+C_y/R-theta M_t/(2rho²)
              +theta E/(rho R)+(L-A)/(2R²),
mathcal K = -(M_yt+theta M_y)/(2rho)
            +(C_t+theta C)/(2rho R)+E_y/(2R)-B/(2R²).
```

The capital letters A through M are the retained seven scalar metric amplitudes. The symbols `mathcal A_t,mathcal A_y` above are curvature combinations, distinct from derivatives of the metric amplitude A.

Let `xi^t Y,xi^y Y` be the first coordinate position variations at fixed proper time and initial labels, and put `X=rho xi^t`, `Z=xi^y`. The complete coefficient is

\[
h_{\rm geo}
 ={3(cU_{Tx}-sU_{Tr})\over2\ell\rho^2\cosh t}
 +{4\over R}\left[
 cF_t+sF_y-X(c\mathcal A_t+s\mathcal K)
             -Z(s\mathcal A_y+c\mathcal K)\right].             \tag{6}
\]

For completeness, the off-screen changes of the normalized angular axes are
`delta e_a=p D_aY U0+q D_aY n0`, where `p=(cX-sZ)/R` and `q=(-sX+cZ)/R`. Orthonormality forces `delta U=p gradY` and `delta n=-q gradY`. Varying all four Weyl slots gives
`delta H=-[p(g+e_chi)+q k_chi]Aax`, where
`Aax=epsilon_(a^{cd}D_cY Hess(Y)_b)d=2(Tplus-Tminus)`.
The two factors of two come from that angular identity and from the first-frame/first-curvature cross term in a second derivative. The Schouten terms cancel in the three combinations
`mathcal A_t=g+e_t`, `mathcal A_y=e_y-g`, and `mathcal K=k`.
Substituting p,q yields (6). Thus angular trajectory coefficients and second geodesic variations are unnecessary for this coefficient, although the congruence proof still requires the label-map margins.

The TT term in (6) follows from `H_TT=(U0 u)curl(T)/(2R)`, the fixed eigenvalues `+3,-3`, and `partial_y=-(ell rho)^-1 partial_r`. Scalar and vector linear second sectors have zero selected TT projection by the retained spherical representation separation. This retains the exact-family TT binding while changing the observation prescription.

## Actual first-position and TT enclosures

The [first-position lemma](geometry/FIRST_POSITION_LEMMA.md) uses
`P=rho delta V^t`, `Q=delta V^y`, and proper-time differentiation:

```
X'=a sinh(chi) X+P,             Z'=Q,
P'=-a sinh(chi)P-2a cosh(chi)Q-2a_y cosh(chi)sinh(chi)Z-S_P,
Q'=-2a cosh(chi)P-(2a²+a_y)cosh(chi)²Z-S_Q,

S_P=aB-[cosh(chi)² D_t-sinh(chi)² A_t]/(2rho)
    -cosh(chi)sinh(chi)D_y-sinh(chi)²B_y,
S_Q=cosh(chi)²[B_t/rho-a(A+D)-D_y/2]
    +cosh(chi)sinh(chi)[A_t/rho-2aB]+sinh(chi)²A_y/2.
```

Here `A_t,A_y` are metric derivatives. All four initial values are zero by the exact open collar established above. The t,y subsystem uses only A,B,D and their first derivatives; background radial motion does not couple its first variation to angular trajectory variables.

The [frozen campaign](numeric/FREEZE.json) uses 64 time cells covering `0<=x<=5/4` and the fixed radial box `1249/10000<=r<=1251/10000`. The background monotonicity and `rho a<20` imply a center-path displacement less than `20(5ell/4)²/2`, which places the whole center trajectory inside that box. Each cell is inside the inherited actual-field domain. Sixty-four calls to the complete actual first-jet oracle include its differentiated global errors; one further apex call supplies the metric two-jet for (6).

For a cell of width `Delta x`, `d tau/dx=ell rho sech(chi)<=ell rho_max`. If `S_j` bounds both absolute forcing components on that entire cell, define

\[
\mathcal B=\sum_j\ell\,\Delta x_j\rho_{j,\max}S_j.
\]

With T an upper bound on `tau_dagger`, define uniform drift bounds

```
M=max(sup[|a sinh(chi)|+2|a cosh(chi)|],sup 2|a cosh(chi)|),
K=max(sup 2|a_y cosh(chi)sinh(chi)|,
      sup |(2a²+a_y)cosh(chi)²|),
b=sup |a sinh(chi)|.
```

Integration and `sup|Z|<=T W`, for `W=max(sup|P|,sup|Q|)`, give

\[
W\le{\mathcal B\over1-MT-KT^2/2},
\qquad |Z|\le TW,
\qquad |X|\le{TW\over1-bT}.                                  \tag{7}
\]

The report verifies both denominators are positive. Useful loose rational consequences of its exact entries are
`T<0.000162407`, `W<0.001808806`, and `|X|,|Z|<0.000000293762`.
The source integral uses whole-cell ranges, rather than point quadrature or comparator values.

The [new TT calculation](tt/DERIVATION.md) encloses the combined actual jet in (6) using the retained authenticated physical source boxes. It separately reconstructs the radial derivative and charges the full scalar feedback. For `T_x=5/4`, `R_x=1/8`, `L=T_x-R_x=9/8`, write

```
Fplus(v)=F_T(v,T_x+R_x-v),
Fminus(v)=F_T(v,|v-L|),
V_x=-G_T(L)+(1/2)integral(Fplus+Fminus),
V_r=+G_T(L)+(1/2)integral(Fplus-Fminus).
```

Thus

\[
cV_x-sV_r=-e^\chi G_T(L)
 +\tfrac12\int_0^{T_x}
       [e^{-\chi}F_{\rm plus}+e^\chi F_{\rm minus}]\,dv.
\]

The opposite endpoint signs retain the inhomogeneous Neumann history and physical reflection. The inherited C1 contraction estimate bounds each actual derivative error by the complete resolvent allowance d, so the combined error is at most `exp(|chi|)d`. The rapidity interval is charged in both exponential weights and the gain. This bounds the actual TT solution, including every scalar feedback iterate. The [TT receipt](tt/TT_NORMAL.json) uses authenticated retained source ranges; it does not rerun their source-oracle campaign.

## Complete sum, verification, and small-parameter consequence

The primary calculation inserts the exact apex intervals, (7), and the combined TT enclosure into (6) using rational outward arithmetic. For orientation, the boosted scalar quadratic term lies inside `[3.4815e-18,1.4205e-17]` and the entire first-geodesic tilt lies inside `[-2.6178e-17,2.6178e-17]`; these deliberately loose decimal ranges are consequences of the exact report. Neither contribution is dropped. The complete sum gives (1).

The verification has distinct roles. The [geometric checks](geometry/FINAL_FORMULA_NORMAL.json) verify the scalar reductions, Lorentz slots, sphere normalization, and clock factors, with an exact full-curvature fixture and actual incorrect-formula controls. The [independent geometric review](independent/INDEPENDENT_REVIEW.md) uses separately constructed coordinate curvature and analytic reductions. The background and family-domain checks supply the shooting and observer margins. The TT check independently reconstructs weighted endpoint integrals and feedback from the retained box inputs. Finally the independent numerical reconstruction recomputes the first-variation source caps, apex map, positive denominators, and complete coefficient arithmetic from the producer-recorded actual-jet boxes; it does not rerun the 65 oracle calls. These checks share the explicitly retained analytic and actual-field premises; they are not an external specialist review or a proof-assistant verification.

Because `J_geo` is C2 and satisfies (2), Taylor's theorem gives

\[
J_{\rm geo}(e)=\tfrac12h_{\rm geo}e^2+o(e^2).
\]

Let `u=-7301026733/500000000000`, the negative outward upper endpoint displayed in (1). The little-o remainder is bounded in absolute value by `(-u)e²/4` on some sufficiently small punctured parameter interval. Since `h_geo<=u`, this proves (3), equivalently

\[
J_{\rm geo}(e)\le
-{7301026733\over2000000000000}\,e^2<0.
\]

This is an ineffective small-amplitude conclusion for the selected exact family and anchored proper-time observation. The scalar is a fixed sphere projection; its negativity does not assert pointwise negativity at every angle. The construction depends on the chosen initial labels, clock value, screen axes, templates, and orientation. It yields no observer-independence, detector, radiation-energy, or cosmological-origin claim.
