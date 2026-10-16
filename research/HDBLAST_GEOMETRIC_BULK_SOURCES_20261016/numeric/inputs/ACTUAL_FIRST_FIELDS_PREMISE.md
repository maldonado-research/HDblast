# Explicit first relative-harmonic metric/scalar pair

Unpublished follow-on, 16 October 2026. This note is separate from the frozen
first-gauge-wave capsule. It derives the actual transformed first fields and a
direct quadratic-source recipe; it does not enclose or evaluate that forcing,
construct a second gauge vector, or solve a second evolution.

## 1. Conventions and the actual first pair

Use

\[
g_0=dy^2+\rho(y)^2[-dt^2+c(t)^2\kappa],\quad c=\cosh t,
\quad H=\tanh t,\quad a=\rho_y/\rho,\quad p=\phi_{0,y},
\quad dy=-\rho\,dz.
\]

The quadrupole is \(Y=n_1^2-n_2^2\), with
\(\Delta_\kappa Y=-8Y\). The longitudinal first pair is

\[
h_{yy}=-4uY,\quad h_{y\mu}=0,\quad
h_{\mu\nu}=2uY(g_0)_{\mu\nu},\qquad \chi=vY.
\tag{1}
\]

Here \(\mu,\nu\) run over \(t\) and the three sphere coordinates. With the
**actual registered** master \(w\), not the displacement comparator, put

\[
R=\alpha/\rho^2,\quad B=3\alpha/(\rho^3p),\quad
k=\alpha_z/\alpha,\quad \vartheta=c^{-3/2},\qquad
u=R\vartheta w,\quad v=B\vartheta(w_z+kw).
\tag{2}
\]

The raw contravariant gauge vector certified in the preceding work is

\[
\xi^y=fY,\qquad \xi^t=bY,\qquad
\xi^i=dD_\kappa^iY.                                      \tag{3}
\]

No factor of \(\rho\) has been removed from these amplitudes. The plus sign
in \(h'=h+\mathcal L_\xi g_0\), \(\chi'=\chi+\mathcal L_\xi\phi_0\)
is the pullback convention used by the gauge-wave certificate. It is paired
with \((\Box_0+\operatorname{Ric}_0)\xi=-DC[h]\).

Direct use of
\((\mathcal L_\xi g_0)_{AB}=\xi^C\partial_C(g_0)_{AB}
+(g_0)_{CB}\partial_A\xi^C+(g_0)_{AC}\partial_B\xi^C\)
gives all fifteen independent metric components:

\[
\begin{aligned}
h'_{yy}&=(-4u+2f_y)Y,\\
h'_{yt}&=(f_t-\rho^2b_y)Y,\\
h'_{yi}&=(f+\rho^2c^2d_y)D_iY,\\
h'_{tt}&=-2\rho^2(u+af+b_t)Y,\\
h'_{ti}&=\rho^2(c^2d_t-b)D_iY,\\
h'_{ij}&=2\rho^2c^2[(u+af+Hb)Y\kappa_{ij}
                         +dD_iD_jY],\\
\chi'&=(v+pf)Y.
\end{aligned}                                                    \tag{4}
\]

In the angular calculation, \(\mathcal L_{D Y}\kappa=2\operatorname{Hess}Y\).
In particular the clock term in the angular trace is \(Hb\), while the
time-time component contains \(b_t\). The scalar shift has sign \(+pf\).

## 2. Background-frame and scaled formulas

Take the background orthonormal frame
\(e_y=\partial_y\), \(e_0=\rho^{-1}\partial_t\),
\(e_a=(\rho c)^{-1}E_a\), where \(E_a\) is \(\kappa\)-orthonormal.
The timelike frame vector still has squared length \(-1\). Then

\[
\begin{aligned}
h'_{yy}&=(-4u+2f_y)Y,\\
h'_{y0}&=(f_t/\rho-\rho b_y)Y,\\
h'_{ya}&=[f/(\rho c)+\rho c\,d_y]D_aY,\\
h'_{00}&=-2(u+af+b_t)Y,\\
h'_{0a}&=(cd_t-b/c)D_aY,\\
h'_{ab}&=2[(u+af+Hb)Y\delta_{ab}+dD_aD_bY].
\end{aligned}                                                    \tag{5}
\]

Let \(x=t/\ell\), \(r=z/\ell\), \(\ell=10^{-6}\), and use the same symbols
for the raw amplitudes evaluated at \((\ell x,\ell r)\). Put
\(W(x,r)=w(\ell x,\ell r)\). Thus

\[
u=R\vartheta W,\qquad
v=B\vartheta(\ell^{-1}W_r+kW),\quad
\partial_t=\ell^{-1}\partial_x,\quad
\partial_y=-(\ell\rho)^{-1}\partial_r.                         \tag{6}
\]

The six frame sectors in (5) become, in their displayed order,

\[
\begin{aligned}
&[-4u-2f_r/(\ell\rho)]Y,\\
&[f_x/(\ell\rho)+b_r/\ell]Y,\\
&[f/(\rho c)-c d_r/\ell]D_aY,\\
&-2[u+af+b_x/\ell]Y,\\
&[cd_x/\ell-b/c]D_aY,\\
&2[(u+af+Hb)Y\delta_{ab}+dD_aD_bY].
\end{aligned}                                                    \tag{7}
\]

The frame factors in (5)--(7) are part of the conversion: coordinate-amplitude
bounds for \(f,b,d\) alone are not metric-component bounds. Also derivatives
of frame components are not covariant derivatives until the frame connection
is included.

For a future covariant evaluator, the complete frame connection is

\[
\begin{gathered}
\nabla_y e_y=\nabla_y e_0=\nabla_y e_a=0,\qquad
\nabla_0 e_y=ae_0,\quad \nabla_a e_y=ae_a,\quad
\nabla_0e_0=ae_y,\quad \nabla_0e_a=0,\\
\nabla_ae_0=(H/\rho)e_a,\qquad
\nabla_ae_b=(\rho c)^{-1}\omega^d{}_{ab}e_d
                      +(H/\rho)\delta_{ab}e_0-a\delta_{ab}e_y,
\end{gathered}
\]

where \(\omega\) is the unit-sphere frame connection. At each angular point
one may use a sphere-normal frame. For the actual quadratic harmonic, the
angular derivative identities, with all indices in that sphere frame, are

\[
\begin{aligned}
D_cD_aD_bY&=-\delta_{ca}D_bY-\delta_{cb}D_aY-2\delta_{ab}D_cY,\\
D_dD_cD_aD_bY&=-\delta_{ca}D_dD_bY-\delta_{cb}D_dD_aY
                                      -2\delta_{ab}D_dD_cY.
\end{aligned}
\]

For example these follow by restricting the ambient quadratic form
\(A=\operatorname{diag}(1,-1,0,0)\):
\(D_aD_bY=2\langle E_a,AE_b\rangle-2Y\delta_{ab}\).
Uniform individual sphere-frame component caps for derivative orders zero
through four are therefore \(1,2,4,8,16\). They also bound absolute angular
directional components; they are not summed tensor norms.

For any scalar amplitude \(F(x,r)\), write \(L=(\ell\rho)^{-1}\). The radial
conversion needed through order three is

\[
\begin{aligned}
F_y&=-LF_r,\\
F_{yy}&=L^2F_{rr}+aLF_r,\\
F_{yyy}&=-L^3F_{rrr}-3aL^2F_{rr}+(a_y-a^2)LF_r.
\end{aligned}                                                    \tag{8}
\]

Time derivatives commute with these radial conversions and contribute the
corresponding power of \(\ell^{-1}\). These formulas use
\(\rho_y=a\rho\); treating \(\rho\) as constant would lose the last terms.

## 3. Shell and initial rows

The certified physical rows are \(f(t,0)=0\), \(b_y(t,0)=d_y(t,0)=0\).
Tangential/time differentiation of \(f|_b=0\) gives \(f_t|_b=0\) and
\(D_i(fY)|_b=0\). Equation (4) therefore proves

\[
h'_{yt}|_b=h'_{yi}|_b=0,\qquad \chi'|_b=vY.                    \tag{9}
\]

The vector is tangent to the shell, so the shell is fixed as a set; \(b,d\)
can change its clock/angular parameterization. There is no assertion that
the tangential transformed metric equals the longitudinal metric at the
shell. No extra Dirichlet row for \(b,d\) is imposed.

Initially \(f=b=d=f_t=b_t=d_t=0\) identically in the radial coordinate.
Thus \((h',\chi')|_{t=0}=(h,\chi)|_{t=0}\), and
\(\chi'_t|_{t=0}=\chi_t|_{t=0}\). The gauge-wave equations give

\[
f_{tt}|_0=-6\rho^2(u_y+4au)|_0,\qquad b_{tt}|_0=d_{tt}|_0=0.
\tag{10}
\]

Consequently the sole change among initial coordinate metric velocities is
\(\partial_t h'_{yt}-\partial_t h_{yt}=f_{tt}Y\). Equality of the initial
metric does not mean equality of its entire first spacetime jet.

The longitudinal gauge defect is
\(DC[h]_y=-6(u_y+4au)Y\), \(DC[h]_\mu=0\). The solved vector equation and
\(DC[\mathcal L_\xi g_0]=(\Box_0+\operatorname{Ric}_0)\xi\) establish
\(DC[h']=0\). This is the actual relative-harmonic representative.

## 4. A direct full quadratic-source recipe

The retained natural residuals, with fixed actual potential \(\mathcal V\),
are

\[
\mathscr E_{AB}=R_{AB}-\phi_A\phi_B-\tfrac23\mathcal V(\phi)g_{AB},
\qquad \mathscr S=\Box_g\phi-\mathcal V'(\phi).                 \tag{11}
\]

In this section put \(q=h'\) and \(\eta=\chi'\); \(q\) is a metric tensor,
not the compact packet. All indices, contractions and covariant derivatives
below use \(g_0\). Define

\[
P^A{}_{BC}=\tfrac12g_0^{AD}
 (\nabla_Bq_{CD}+\nabla_Cq_{BD}-\nabla_Dq_{BC}),\quad
J^A{}_{BC}=-2q^A{}_D P^D{}_{BC},\quad
C^A=g_0^{BC}P^A{}_{BC}.                                      \tag{12}
\]

Along the **affine pair path** \((g_0+\epsilon q,\phi_0+\epsilon\eta)\),
\(P=\Gamma'\), \(J=\Gamma''\). Hence

\[
\begin{aligned}
Q^E_{AB}={}&\nabla_CJ^C{}_{AB}-\nabla_BJ^C{}_{AC}
 +2(P^C{}_{CD}P^D{}_{AB}-P^C{}_{BD}P^D{}_{AC})\\
&-2\eta_A\eta_B
 -\tfrac23[\mathcal V''\eta^2(g_0)_{AB}
                         +2\mathcal V'\eta q_{AB}].
\end{aligned}                                                     \tag{13}
\]

\[
\begin{aligned}
Q^\phi={}&2q^A{}_Dq^{DB}\nabla_A\nabla_B\phi_0
          -2q^{AB}\nabla_A\nabla_B\eta
          +2q^{AB}P^C{}_{AB}\partial_C\phi_0
          -\mathcal V'''\eta^2\\
&+2q^C{}_D C^D\partial_C\phi_0-2C^C\partial_C\eta.
\end{aligned}                                                     \tag{14}
\]

For the actual pair \(C=DC[q]=0\), so the last line of (14) vanishes.
This simplification is justified by the first relative-harmonic solve, rather
than by deleting the general trace terms before that solve.

Equations (13),(14) are second parameter derivatives, so the coefficient of
\(\epsilon^2\) is one half of the displayed source. For example the scalar
stress term is \(-2\eta_A\eta_B\); there is no second metric term in the
potential product along this affine path. They retain every tensor component
and impose no second-order longitudinal ansatz.

For completeness the quadratic relative-harmonic vector of this affine path is

\[
Q_C^A=g_0^{BC}J^A{}_{BC}-2q^{BC}P^A{}_{BC}
     =-2q^A{}_D C^D-2q^{BC}P^A{}_{BC}.
\tag{15}
\]

For the actual harmonic first pair it is simply
\(Q_C^A=-2q^{BC}P^A{}_{BC}\). This is a useful quadratic gauge row; it is
not a prescribed second vector or a complete boundary condition.

The retained reduction is explicitly
\(\mathscr E^{\mathrm{red}}_{AB}=\mathscr E_{AB}-\nabla_{(A}C_{B)}\),
with symmetrization containing one half, \(C_B=g_{BA}C^A\), and fixed
reference connection \(\Gamma[g_0]\). This is the convention in
`FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`, section 3, equation (6a), and
`GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md`, section 1. Since
\(C_0=C_1=0\), differentiation of lowering and covariant differentiation
adds no surviving first-variation cross terms. Thus its direct Hessian is

\[
Q^{E,\mathrm{red}}_{AB}=Q^E_{AB}
                -\nabla^0_{(A}(Q_C)_{B)},\qquad
(Q_C)_B=(g_0)_{BA}Q_C^A.                                      \tag{16}
\]

The second reduced equation has forcing
\(-Q^E_{AB}+\nabla^0_{(A}(Q_C)_{B)}\), while the scalar forcing is
\(-Q^\phi\). This matches the registered sign in equation (6a).
Formula (16) needs the metric two-jet, hence at most the gauge three-jet;
it does not invoke the fourth-derivative covariance shortcut.

## 5. Required derivatives and the next finite computation

Direct use of (13) needs \(q\) through two derivatives and \(\eta\) through
one; (14) needs \(q\) through one and \(\eta\) through two. Formula (4)
therefore requires \(f,b,d\) through order three in the Einstein source,
and through order two in the scalar source. In (2), \(u\) through order two
uses master order two, while \(v\) through order two uses master order three.
Thus the actual master C3 and actual gauge-vector C3 certificates suffice
to define and continuously bound the direct full quadratic natural and
registered reduced source. A first spatial/time derivative of those source
terms for a later higher-regularity estimate would in general need one more
master/gauge derivative; no such extra certificate is asserted here.

This statement of sufficiency is not an evaluated bound. The next finite
calculation is to combine the authenticated full collar coefficient jets,
actual master mixed jets and actual gauge mixed jets in (4), then in the
covariant recipes (12)--(15), including the frame connection and sphere
derivatives. A coefficient-wise Taylor-jet or exact interval evaluator would
produce actual relative-harmonic source bounds. The present checker verifies
the finite transformation/source algebra only; it does not perform that
enclosure or replay the prior analytic certificates.

There is a separate covariance identity for the unreduced natural residual:
if \(\Phi_1' = \Phi_1+\mathcal L_\xi\Phi_0\) and
\(G=2\mathcal L_\xi\Phi_1+\mathcal L_\xi^2\Phi_0\), then formally on shell
\(Q(\Phi_1',\Phi_1')=Q(\Phi_1,\Phi_1)-E'[G]\).
Separately forming or bounding the second-order operator \(E'[G]\) can use
four derivatives of \(\xi\), since \(G\) already contains two. No C4 gauge
certificate is claimed. This does not prevent the direct C3 calculation
above; use of the covariance shortcut needs extra regularity or a proved
cancellation.

## 6. Portable draft checks and scope

`check_first_harmonic_fields.py` is a standalone standard-library draft. It
checks all metric components against the coordinate Lie derivative at exact
generic local metric jets, the scalar pullback, the normalized radial
conversion, and the full quadratic connection, Ricci, matter and scalar
recipes against direct affine expansions. These exact finite witnesses supplement the derivation;
they are not an all-point evaluation of the registered master/gauge fields.
Controls must fail for the named reason in both Python modes. The checker
does not read scientific reports and stores no absolute paths or timestamps.

From this folder, choose new output paths and run

```sh
python check_first_harmonic_fields.py --run-controls --output /tmp/first-fields-normal.json
python -O check_first_harmonic_fields.py --run-controls --output /tmp/first-fields-optimized.json
```

Existing output paths are rejected. The retained draft runs each pass 385
exact finite comparisons and five mutation-specific controls, with
byte-identical normal/optimized reports. `DRAFT_HANDOFF_PINS.json` records
artifact hashes and the prior documents used to fix conventions; these
reference hashes are provenance, not unperformed report replays.

No files of the frozen public capsule are modified. There is no claim here
of a second evolution, a nonlinear family, an invariant force, a boundary
compatibility solve, or a numerical full relative-harmonic forcing enclosure.
