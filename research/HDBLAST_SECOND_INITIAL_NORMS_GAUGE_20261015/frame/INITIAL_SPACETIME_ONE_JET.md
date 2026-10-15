# Completed initial spacetime second 1-jet in the background frame

15 October 2026. Conditional on the authenticated actual regular background,
the specified explicit physical second initial jet and its local relative
harmonic completion. On \(t=0\), \(\ell\le z\le2\ell\),
\(\ell=10^{-6}\), and the whole sphere, the following strict component
bounds hold:

| Initial spacetime coefficient | Value | All first spacetime covariant directions |
|---|---:|---:|
| \(k=g_2\) | \(1.06\times10^{24}\) | \(6.4\times10^{29}\) |
| \(\phi_2\) | \(5.13\times10^{33}\) | \(3.06\times10^{39}\) |

These are second parameter derivatives, multiplying \(\epsilon^2/2\).
The norm is the largest absolute **individual component**, including the
derivative index, in the fixed background orthonormal frame
\(e_0=\rho^{-1}\partial_t\), \(e_y=\partial_y\),
\(e_a=\rho^{-1}E_a\) at initial time. Here \(e_0\) is timelike; the
component maximum uses Euclidean absolute values of its coefficients and is
not a positive Lorentzian energy. The normal metric derivative majorant is
approximately \(6.36314\times10^{29}\); the other displayed convenient
caps are inherited from the tighter physical spatial bounds.

## Fixed initial representative and frame conversion

Use the physical conventions \(K=\tfrac12\mathcal L_n\gamma\),
\(\Pi=n\phi\), with future Cauchy normal \(n\). They are not ADM
momentum densities. As in the local harmonic completion, set

\[
U=uY,\quad T=u_tY,\quad q=(F/3)Y^2,\quad
f=\sum_{l=0,2,4}f_lW_l,\qquad N_1=\rho U,\quad N_2=0,
\quad\beta_1=\beta_2=0
\]

for initial lapse and shift **values**. The notation \(q\) in this note
is the pure-trace \(K_2\) coefficient, not the first master bump.
Their time derivatives are the local algebraic values enforcing the initial
relative harmonic condition. In particular they are not set to zero.

The second metric's initial background-frame values are

\[
k_{00}=-2U^2,\quad k_{ab}=f\delta_{ab},\qquad
k_{yy}=k_{0i}=k_{ya}=0.                                       \tag{1}
\]

The authenticated harmonic-completion derivation supplies the full coordinate
velocities, including the first-shift cross term in \((\dot g_{ty})_2\)
and the first-lapse cross term in \((\dot g_{tt})_2\). The background
normal-frame connection at initial time is

\[
\nabla_0e_0=ae_y,\qquad\nabla_0e_y=ae_0,\qquad\nabla_0e_a=0.
\]

Conversion therefore gives every nonzero normal derivative class:

\[
\begin{aligned}
(\nabla_0k)_{00}&=-8q+44UT/\rho,\\
(\nabla_0k)_{0y}&=-26UU_y-\tfrac32f_y-3af,\\
(\nabla_0k)_{0a}&=22U D_aU/\rho-\tfrac12D_af/\rho,\\
(\nabla_0k)_{yy}&=2q-8UT/\rho,\\
(\nabla_0k)_{ab}&=(2q+4UT/\rho)\delta_{ab},\\
\nabla_0\phi_2&=\Pi_2+2U\Pi_1.
\end{aligned}                                                 \tag{2}
\]

Here \(D_a\) is a unit-sphere derivative, and
\(\Pi_1=v_tY/\rho\). The mixed \(0y\) cancellation matters:
\((\dot g_{ty})_2/\rho^2\) includes \(-2aU^2\), whereas the
frame connection contributes \(-a k_{00}=+2aU^2\). Dropping either
term changes (2). One cannot obtain the normal derivative of \(k_{00}\)
by differentiating only its initial value \(-2U^2\) while pretending
that the second lapse remains zero off the initial section.

Spatial covariant derivatives of the \(\gamma_2\) block are already
bounded by the physical norm certificate, including its radial-angular
projector connection terms. At \(t=0\), \(\nabla_Ie_0=0\) for each
spatial index. The additional spatial metric components are therefore only

\[
(\nabla_Ik)_{00}=-4U\partial_IU,
\]

and all spatial derivatives of the initially zero \(k_{0i}\) vanish.
Together with (2), this covers all indices of the initial spacetime 1-jet.

## Exact quantitative majorants

The producer's physical report supplies exact rational upper bounds
\(M_\gamma,M_{\gamma,y},M_{\gamma,\mathrm{ang}},M_{\nabla\gamma},
M_\phi,M_{\nabla\phi},M_q,M_\Pi\) for the corresponding amplitudes
and component classes. No displayed float is reused. From the actual
background and bump inputs let

\[
C_*\ge|C|,\ C_{z,*}\ge|C_z|,\ B_*\ge|B|,\ k_*\ge|k|,
\quad r_*\le\rho,\ a_*\ge|a|,\quad Q_1=16/\ell.
\]

The first initial data satisfy

\[
|U|,|T|\le C_*,\quad
|U_y|\le L_y=(C_{z,*}+C_*Q_1)/r_*,\quad
|D U|\le L_s=2C_*,\quad
|\Pi_1|\le P_1=B_*(Q_1+k_*+1)/r_*.
\]

The normal derivative majorants used in the exact checker are

\[
\begin{array}{ll}
00:&8M_q+44C_*^2/r_*,\\
0y:&26C_*L_y+3M_{\gamma,y}/2+3a_*M_\gamma,\\
0a:&22C_*L_s/r_*+M_{\gamma,\mathrm{ang}}/2,\\
yy:&2M_q+8C_*^2/r_*,\\
ab:&2M_q+4C_*^2/r_*,\\
\phi:&M_\Pi+2C_*P_1.
\end{array}
\]

The metric value bound is \(\max(M_\gamma,2C_*^2)\). The full
spatial derivative bound additionally takes the maximum of
\(M_{\nabla\gamma},4C_*L_y,4C_*L_s/r_*\). The final full metric
derivative bound is the maximum of this spatial bound and the five metric
rows above. The full scalar derivative bound is the maximum of its inherited
spatial bound and the last row. All operations here use exact Fractions.

## Reproduction, authentication and scope

`INPUT_PINS.json` authenticates copies of the physical norm report, actual
background and bump reports, and the gauge completion report and derivation.
It also checks the physical report's background/bump dependency hashes. Paths
are relative to `--input-root`, defaulting to this checker directory. No prior
scientific module is imported or executed.

```sh
python check_initial_spacetime_frame.py --run-controls --output INITIAL_SPACETIME_FRAME_NORMAL.json
python -O check_initial_spacetime_frame.py --run-controls --output INITIAL_SPACETIME_FRAME_OPTIMIZED.json
```

Forty-two exact checks at three nonspecial rational point jets independently
convert the coordinate velocities with the full frame connection and compare
all component classes in (2). The all-point reason for these identities is
the coordinate/frame calculation above and the authenticated exact local gauge
prescription. Four controls omit the normal clock factor, omit the frame
connection, omit the mixed \(f_y\) term, or halve the scalar lapse cross
term. Each must fail for its intended reason and leave no output report.
Both ordinary and optimized Python give the same mathematical bounds.

These estimates hold on the initial section only. The compact profiles and
local gauge completion preserve the original annulus and flat endpoint jets.
The bounds therefore also hold for their smooth zero extensions. They do not
bound a future time slab, spatial derivatives of the gauge time velocities,
second spacetime derivatives, a high Sobolev norm or an exact nonlinear
remainder. In particular they supply no validated second evolution or shell
boundary response. Evolving harmonic sources must still be transformed in one
common representative; the old longitudinal forcing cannot be attached
unchanged merely because the initial gauge was completed.
