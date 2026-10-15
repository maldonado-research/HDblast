# Review of the completed initial spacetime one-jet in the background frame

15 October 2026. Independent mathematical review of the proposed frame conversion and coarse bound integration. No frozen gauge formula or producer file was changed. This report concerns the second parameter derivative and its full first spacetime derivative **at the initial slice**. It does not bound derivatives of those time velocities, a positive-time solution, or a Sobolev evolution norm.

**Judgment: accepted.** The proposed frame connection, all five listed time-derivative component families, and the scalar conversion are correct. The metric derivative cap \(6.4\times10^{29}\) and scalar derivative cap \(3.06\times10^{39}\) follow from the cited physical bounds and elementary first-data bounds. Two exact cancellations make the dominant metric component particularly transparent.

## 1. Fixed background frame and component convention

Use the background metric

\[
g_0=dy^2+\rho^2[-dt^2+\cosh^2t\,\kappa_{S^3}],\qquad
a=\rho_y/\rho.
\]

At \(t=0\), let \(e_0=\rho^{-1}\partial_t\), \(e_y=\partial_y\), and \(e_a=\rho^{-1}E_a\), where the \(E_a\) form a sphere orthonormal frame. The indices \(a,b,c\) are sphere-frame indices. The derivative and tensor slots are all evaluated in this fixed background frame; it is not a parameter-dependent orthonormalization of the perturbed metric. The bounds mean the supremum of the absolute values of individual frame components, including the derivative index.

The time connection at the initial slice is

\[
\nabla_0 e_0=a e_y,\qquad \nabla_0 e_y=a e_0,
\qquad \nabla_0 e_a=0.
\]

All radial frame derivatives vanish. In angular directions,

\[
\nabla_a e_0=0,\qquad \nabla_a e_y=a e_a,\qquad
\nabla_a e_b=(\hbox{sphere connection})-a\delta_{ab}e_y.
\]

These signs follow from \(\Gamma^y_{tt}=a\rho^2\), \(\Gamma^t_{ty}=a\), and \(\Gamma^y_{ij}=-a\rho^2\kappa_{ij}\). The de Sitter time expansion contributes nothing at \(t=0\), since \(\partial_t\cosh t=0\) there.

For \(k=(g)_2\), with \(U=uY\), \(T=u_tY\), and \(f=\sum_l f_lW_l\), its initial frame components are

\[
k_{00}=-2U^2,\qquad k_{ab}=f\delta_{ab},\qquad
k_{yy}=k_{ya}=k_{0y}=k_{0a}=0.
\]

## 2. All initial time-direction covariant derivatives

The complete nonzero component families are

\[
\boxed{(\nabla_0 k)_{00}
=-\frac{2\dot N_2}{\rho^2}-\frac{4UT}{\rho}
=-8q+\frac{44UT}{\rho},}                                      \tag{1}
\]
\[
\boxed{(\nabla_0 k)_{0y}
=\frac{\dot\beta_2^y-8U\dot\beta_1^y}{\rho^2}+2aU^2
=-26UU_y-\frac32f_y-3af,}                                    \tag{2}
\]
\[
\boxed{(\nabla_0 k)_{0a}
=\frac{\dot\beta_2^a}{\rho}
=\frac{22U E_aU-\tfrac12E_af}{\rho},}                         \tag{3}
\]
\[
\boxed{(\nabla_0 k)_{yy}=2q-\frac{8UT}{\rho},\qquad
(\nabla_0 k)_{ab}=\left(2q+\frac{4UT}{\rho}\right)\delta_{ab}.} \tag{4}
\]

The \(ya\) components vanish, and symmetry gives the reversed tensor slots. In (3), \(\dot\beta_2^a\) means components of the angular shift vector in the unit-sphere frame \(E_a\); a coordinate angular component must first be converted into that frame. There is no additional factor \(\rho\) hidden in that notation.

For (2), the connection correction is
\(-k(\nabla_0e_0,e_y)-k(e_0,\nabla_0e_y)=-ak_{yy}-ak_{00}=2aU^2\).
It cancels the coordinate-velocity term \(-2aU^2\). Omitting this correction would give an incorrect tensor derivative, even though its numerical effect here is small. Formula (1) also retains the lapse cross term \(-4N_1\dot N_1\) in the second coordinate derivative of \(g_{tt}\).

The scalar conversion is

\[
\boxed{\nabla_0\phi_2=\Pi_2+2U\Pi_1.}                         \tag{5}
\]

It is the derivative along the fixed background \(e_0\). The physical coefficient \(\Pi_2=(n_\epsilon\phi_\epsilon)_2\) alone is not that fixed-frame derivative.

## 3. Spatial covariant derivatives are also complete

For the time-time component,

\[
(\nabla_yk)_{00}=-4UU_y,\qquad
(\nabla_ak)_{00}=-4U\rho^{-1}E_aU.
\]

The spatial-spatial component derivatives coincide with the producer's full covariant derivatives of \(\gamma_2\):

\[
(\nabla_yk)_{ab}=f_y\delta_{ab},\quad
(\nabla_ak)_{bc}=\rho^{-1}(E_af)\delta_{bc},\quad
(\nabla_ak)_{yb}=(\nabla_ak)_{by}=-af\delta_{ab}.
\]

All spatial derivatives with precisely one time tensor index vanish. In particular there is no missing angular connection term involving \(k_{00}\), because \(\nabla_a e_0=0\) at this slice. This enumerates the full first derivative tensor, including components whose field value is zero.

## 4. Coarse exact bound check

The inspected producer report supplies strict caps

\[
|f|<1.06\times10^{24},\quad |f_y|<4.25\times10^{29},\quad
\rho^{-1}|D_\kappa f|<1.34\times10^{22},\quad |q|<7.96\times10^7,
\]
\[
|\phi_2|<5.13\times10^{33},\qquad
\max_I|\nabla_I\phi_2|<3.06\times10^{39},\qquad |\Pi_2|<2\times10^{11}.
\]

The authenticated background intervals imply, writing \(\varkappa=\alpha_z/\alpha\) for the reconstruction coefficient named `k` in the input report,

\[
129\le\rho\le131,\quad |a|<14/125,\quad
|C|<6\times10^{-5},\quad |C_z|<13/1000,\quad
|B|<12/1000,\quad |\varkappa|<270.
\]

Using the existing bump bounds \(|q_\ell|\le1\), \(|q_{\ell,z}|\le16\times10^6\), \(u=Cq_\ell\), \(u_t=e^{-z}u\), and \(|D_\kappa Y|\le2\), the convenient overestimates

\[
|U|,|T|<10^{-4},\qquad |U_y|<8,\qquad |D_\kappa U|<2\times10^{-4}
\]

follow. For example \((13/1000+(6\times10^{-5})(16\times10^6))/129<8\). Also \(|d|=|3/(\rho p)|<194\) and

\[
|\Pi_1|\le\frac{|B|(16\times10^6+270)+194\times10^{-4}}{129}<1500.
\]

The dominant new metric bound is therefore

\[
|(\nabla_0k)_{0y}|
<26(10^{-4})(8)+\frac32(4.25\times10^{29})
+3\frac{14}{125}(1.06\times10^{24})
=\frac{398437722600000000000000000000013}{625}
<6.4\times10^{29}.                                           \tag{6}
\]

The other time metric families are bounded respectively by
\(8(7.96\times10^7)+44\times10^{-8}/129\),
\(44\times10^{-8}/129+(1.34\times10^{22})/2\), and
\(2(7.96\times10^7)+8\times10^{-8}/129\). They are all far below (6). The added spatial \(00\) derivatives are bounded by \(0.0032\) radially and \(8\times10^{-8}/129\) angularly. The producer's spatial metric derivative cap is already smaller than (6).

Finally \(|\nabla_0\phi_2|<2\times10^{11}+0.3\), so the existing scalar spatial derivative cap also bounds every scalar spacetime derivative component. Exact rational comparisons, rather than displayed floating values, verify all the inequalities in this section.

The resulting valid summary is

\[
\boxed{\|k\|_{\mathrm{frame},\infty}<1.06\times10^{24},\quad
\|\nabla^{g_0}k\|_{\mathrm{frame},\infty}<6.4\times10^{29},}
\]
\[
\boxed{\|\phi_2\|_\infty<5.13\times10^{33},\quad
\|d\phi_2\|_{\mathrm{frame},\infty}<3.06\times10^{39},}
\]

all on \(t=0\), with the derivatives including the time direction. Here the norm is a component maximum, not a sum of field and derivative bounds or a Lorentzian contraction.

## 5. Scope and inspected evidence

It is accurate to call these bounds a **complete initial spacetime one-jet bound after the specified harmonic initialization**, meaning \((k,\nabla^{g_0}k,\phi_2,d\phi_2)|_{t=0}\). They do not bound \(\nabla\dot k\), provide a spacetime \(C^1\) bound on a positive-time slab, prove reduced or physical evolution, or establish a Sobolev well-posedness theorem. The physical data remain unchanged. Their transformed first evolution and second forcing still require the gauge-wave construction and subsequent evolution work.

Inspected numerical inputs:

- `producer/SECOND_INITIAL_NORMS_NORMAL.json`, SHA-256 `e18b348a86a2702a6b520c0cc56e9205918cb725d0677eb970e9db68a129598d`.
- `producer/inputs/BACKGROUND_RECONSTRUCTION_JETS.json`, SHA-256 `3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22`.

The two new all-point Laurent checks in `review/check_gauge_adm_identities.py`, named `NORMAL_MIXED_FRAME_CANCELLATION` and `NORMAL_LAPSE_FRAME_CANCELLATION`, implement exactly (2) and (1). I read both and found their signs, frame powers and factorial conventions correct. This report is an internal mathematical review, not external peer review.
