# Quantitative values and first spatial derivatives of the actual second initial jet

15 October 2026. Conditional on the retained actual regular background and its
authenticated initial-constraint construction. This note gives continuous
upper bounds for the specified **physical Cauchy second jet**, including its
full first spatial covariant derivatives. It does not compute its evolution,
harmonic-gauge completion, a Sobolev energy, or a nonlinear amplitude interval.

On the original closed packet annulus \(\ell\le z\le2\ell\),
\(\ell=10^{-6}\), and the whole unit sphere, the following strict bounds hold:

| Second initial field | Absolute adapted-frame component | Full first spatial covariant derivative component |
|---|---:|---:|
| \(\gamma_2\) | \(1.06\times10^{24}\) | \(4.25\times10^{29}\) |
| \(K_2\) | \(7.96\times10^7\) | \(1.96\times10^{13}\) |
| \(\phi_2\) | \(5.13\times10^{33}\) | \(3.06\times10^{39}\) |
| \(\Pi_2\) | \(2.00\times10^{11}\) | \(5.02\times10^{16}\) |

The numbers are coarse upper bounds in registered model units. They bound the
second **parameter derivatives** multiplying \(\epsilon^2/2\). They are not
finite-amplitude physical fields, response amplitudes, lower bounds, or ADM
momentum-density estimates. Every maximum is over the absolute values of
individual tensor components in the adapted positive background spatial frame.
Both identical interior copies use the same initial data and bounds.

## 1. Specified second initial data

Retain

\[
\gamma_0=dy^2+\rho^2\kappa_{S^3},\quad K_0=\Pi_0=0,\quad
a=\rho_y/\rho,\quad p=\phi_{0,y}\ne0,\quad dy=-\rho\,dz.
\]

Here \(K=\tfrac12\mathcal L_n\gamma\), \(\Pi=n\phi\), and
\(n\) is the future Cauchy normal. These are the same physical Cauchy
variables as the preceding explicit constraint checkpoint, rather than metric
or scalar momentum densities. The potential, common tension, shell parameter
and actual regular branch are unchanged. The copied source derivation
`inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md` specifies that model and the
full initial constraint proof.

Write \(E=e^{-z}\), \(h=U'/p\),
\(\alpha=(\rho/\rho_b)^{3/2}p/p_b\), and

\[
C=\alpha/\rho^2,\quad B=3\alpha/(\rho^3p),\quad
k=\alpha_z/\alpha=\rho(5a/2-h),\quad d=3/(\rho p).
\]

For the exact registered bump

\[
q_\ell(z)=\begin{cases}
\exp\!\left(4-\dfrac{\ell^2}{(z-\ell)(2\ell-z)}\right),&\ell<z<2\ell,\\
0,&\text{otherwise},
\end{cases}
\]

the actual initial amplitudes obey

\[
u=Cq_\ell,\quad u_t=Eu,\qquad
v=B(q_{\ell,z}+kq_\ell),\quad v_t=E(v-du).                     \tag{1}
\]

These identities use the actual velocity data \(w_t=e^{-z}q_\ell\)
at \(t=0\); no evolved-wave approximation is made. They also obey

\[
v=-\frac3p(u_y+2au),\quad p_y=p(h-4a),\quad
a_y=-\rho^{-2}-p^2/3.                                         \tag{2}
\]

For \(Y=n_1^2-n_2^2\), \(N=n_1^2+n_2^2\), use

\[
Z_2=N-1/2,\qquad H_4=Y^2-N/2+1/12.
\]

The Hamiltonian quadratic source of the explicit physical-Cauchy-affine path
is \(Q_H=A Y^2+D|\nabla Y|^2\), where the previous expression was

\[
\begin{aligned}
A={}&-24u_y^2-288auu_y-24uu_{yy}
 -(192a_y+384a^2+32p^2)u^2\\
&-\rho^{-2}(336u^2+12u_t^2+2v_t^2)
 -2v_y^2-16puv_y-2U''v^2,\qquad
D=\rho^{-2}(36u^2-2v^2).
\end{aligned}
\]

Set

\[
S_0=A/6+4D/3,\quad S_2=A/2+2D,\quad S_4=A-4D,
\qquad (\lambda_0,\lambda_2,\lambda_4)=(0,8,24).
\]

The actual second initial fields to be bounded are

\[
\begin{aligned}
f_l&=-\frac{\rho^2S_l}{2(\lambda_l-3)},\qquad
\zeta_l=-\frac{3f_{l,y}}{2p},\\
H&=f_0+f_2Z_2+f_4H_4,\qquad
\Phi=\zeta_0+\zeta_2Z_2+\zeta_4H_4,\\
(\gamma_2)_{yy}&=(\gamma_2)_{yi}=0,\qquad
(\gamma_2)_{ij}=H(\gamma_0)_{ij},\qquad \phi_2=\Phi,\\
K_2&=\mathcal K Y^2\gamma_0,\qquad
\Pi_2=\mathcal P Y^2.                                         \tag{3}
\end{aligned}
\]

The momentum amplitudes in the earlier source derivation are

\[
F=\frac{18uu_t-vv_t}{\rho},\quad
G=\frac{12uu_{yt}-48auu_t-6u_yu_t-2v_tv_y}{\rho},\qquad
\mathcal K=F/3,\quad\mathcal P=(G-F_y)/p.                       \tag{4}
\]

These choices solve all five second differentiated initial constraints.
The present calculation bounds this particular choice, not the alternative
weighted-inverse correction prescribing zero second scalar/density data.

## 2. Two useful exact cancellations before interval evaluation

Bounding the unsimplified formulas separately would lose substantial
correlation. Equations (2) give the exact Hamiltonian reduction

\[
\boxed{\begin{aligned}
A={}&-2v_y^2-8puv_y+8p(h+2a)uv
 +(16p^2-192\rho^{-2})u^2\\
&-(2U''+8p^2/3)v^2-\rho^{-2}(12u_t^2+2v_t^2).
\end{aligned}}                                                \tag{5}
\]

All explicit \(u_y,u_{yy},a_y\) terms disappear. For example,
\(u_y=-pv/3-2au\) and

\[
u_{yy}=-\frac{p(h-6a)}3v-\frac p3v_y+(4a^2-2a_y)u;
\]

substituting these into the old \(A\), then using \(a_y\) in (2),
gives (5). Its derivatives are identities for the actual smooth background
and packet as well, rather than merely pointwise substitutions in a bound.

For the momentum correction, \(E_y=E/\rho\) and
\(d_y=d(3a-h)\). Differentiating \(F\) in (4) with (1) shows that
the large \(vv_y\) terms cancel in \(G-F_y\). Using (2) yields

\[
\boxed{\begin{aligned}
\mathcal K&=\frac{E}{3\rho}(18u^2-v^2+duv),\\
\mathcal P&=\frac{E}{\rho p}\left[
duv_y+\big(d(h-\rho^{-1})+10p\big)uv
 +(30a-6\rho^{-1})u^2+(2\rho^{-1}-a)v^2\right].
\end{aligned}}                                                \tag{6}
\]

The producer evaluates (5) and (6). At three nonspecial rational jets it
checks (5) against the full previous Hamiltonian source, and (6) against
\((G-F_y)/p\), imposing exactly the background/first-data identities used
above. These finite tests supplement the explicit analytic cancellations;
they are not the all-point justification by themselves.

In particular, the \(\Pi_2\) bound is a bound for the same particular
physical second initial jet as before. No second data choice or affine-path
convention has been changed to obtain the cancellation.

## 3. Authenticated actual coefficient and bump jets

The background report `inputs/BACKGROUND_RECONSTRUCTION_JETS.json` supplies
exact rational intervals for normalized jets
\(D_z^j f/j!\), \(0\le j\le4\), of
\(\rho,a,p,C,B,k,h,U''\) on the whole closed packet annulus. It
comes from interval differentiation of the actual static ODE with the
retained strict collar and shell boxes, not sampled background values.
Its report includes a separate formal total-derivative audit and the
\(-\rho^{-1}\partial_z\) conversion. The actual \(p\) interval stays
strictly negative and \(\rho\) strictly positive.

The independent bump certificate supplies the dimensionless derivative bounds
for \(Q(r)=q_\ell(\ell r)\):

\[
\|Q\|\le1,\quad\|Q'\|\le16,\quad\|Q''\|\le256,\quad
\|Q'''\|\le46656/7,\quad
\boxed{\|Q''''\|\le16777216/49.}
\]

These bounds include the flat zero extension at both endpoints. Physical
derivatives obey \(\|q_\ell^{(j)}\|\le\ell^{-j}\|Q^{(j)}\|\).
The fourth derivative is needed for \(\phi_{2,y}\): \(A\) contains
\(v_y^2\), and its second radial derivative includes \(v_yv_{yyy}\).
The scalar reconstruction makes \(v_{yyy}\) contain \(q_\ell''''\).

There is no inference from a value bound to a derivative bound. Every needed
derivative is propagated as a Taylor coefficient. The retained polynomial
orders at the end of the calculation are \(A:2\), \(f_l:2\),
\(\zeta_l:1\), and \(\mathcal P:2\), in the normalized z-jet
representation. Thus \(\zeta_{l,y}\) and all other first radial
derivatives explicitly exist in the certificate.

## 4. Exact outward interval Taylor computation

At an arbitrary point of the annulus, let \(J_f[j]=D_z^jf/j!\).
Products, reciprocals, and derivatives use

\[
J_{fg}[j]=\sum_{k=0}^jJ_f[k]J_g[j-k],\quad
J_{1/f}[0]=1/J_f[0],\quad
J_{1/f}[j]=-J_{1/f}[0]\sum_{k=1}^jJ_f[k]J_{1/f}[j-k],
\]

\[
J_{f_z}[j]=(j+1)J_f[j+1],\qquad
J_{f_y}=J_{-1/\rho}\,J_{f_z}.
\]

Each interval operation rounds outward to denominator \(2^{144}\).
Addition, multiplication and division use exact rational endpoints before
rounding; reciprocal intervals may not cross zero. For a squared value
whose interval crosses zero, its lower endpoint is tightened to zero.
Higher Taylor coefficients still use the full product recurrence. This
tightening is valid for the constant coefficient only and never discards
the correlations needed by differentiation.

Initialize the bump value interval as \([0,1]\), and its jth normalized
coefficient as \([-Q_j,Q_j]/(j!\ell^j)\). For the exponential,
\(1-2\ell\le E\le1\) and \(D_z^jE=(-1)^jE\), with the same
factorial normalization. The background inputs plus (1), (5), and (6)
then provide interval jets for every profile in (3). Multiplication discards
cross-field correlations conservatively. The full exact intervals and
majorants are in `SECOND_INITIAL_NORMS_NORMAL.json`.

The dyadic denominators, exact rational numerators and strict decimal caps
are stored separately from display-only floating values. Every displayed
cap is obtained by an exact upward rounding of a rational upper bound;
no floating value decides acceptance.

## 5. Whole-sphere and full covariant derivative bounds

Take \(e_r=\partial_y\), \(e_i=\rho^{-1}E_i\), where
\(E_i\) is any sphere orthonormal frame. Along y, this adapted frame is
parallel: radial covariant component derivatives equal ordinary radial
derivatives of the amplitudes. In angular directions the radial-angular
connection must be retained.

Since \(0\le N\le1\) and \(0\le Y^2\le N^2\),

\[
|Z_2|\le\tfrac12,\quad -\tfrac5{12}\le H_4\le\tfrac7{12}.
\]

Also
\(|\nabla Z_2|^2=4N(1-N)\le1\) and
\(|\nabla Y^2|^2=16Y^2(N-Y^2)\le4N^2\le4\), because
\(Y^2(N-Y^2)\le N^2/4\). Hence

\[
|\nabla Y^2|\le2,\qquad
|\nabla H_4|\le|\nabla Y^2|+\tfrac12|\nabla N|\le\tfrac52.
\]

For any three amplitude bounds \(b_l\), define
\(T(b)=b_0+b_2/2+7b_4/12\) and \(L(b)=b_2+5b_4/2\).
Then \(|H|\le T(|f|)\), \(|H_y|\le T(|f_y|)\), and
\(|\nabla_{S^3}H|\le L(|f|)\). The same relations hold for
\(\Phi\) and \(\zeta\).

The nonzero adapted-frame derivative components of \(\gamma_2\) have
the forms

\[
\begin{aligned}
(\nabla_{e_r}\gamma_2)_{ij}&=H_y\delta_{ij},\\
(\nabla_{e_a}\gamma_2)_{bc}&=\rho^{-1}(E_aH)\delta_{bc},\\
(\nabla_{e_a}\gamma_2)_{rb}&=(\nabla_{e_a}\gamma_2)_{br}
                         =-aH\delta_{ab}.
\end{aligned}                                                   \tag{7}
\]

Consequently the full component maximum is bounded by

\[
\max\{T(|f_y|),\ \rho_{\min}^{-1}L(|f|),\ |a|_{\max}T(|f|)\}.
\]

The scalar derivative bound is
\(\max\{T(|\zeta_y|),\rho_{\min}^{-1}L(|\zeta|)\}\).
For \(K_2=\mathcal K Y^2\gamma_0\), metric compatibility gives
\(\nabla K_2=d(\mathcal K Y^2)\otimes\gamma_0\); its component
bound is
\(\max\{|\mathcal K_y|,2|\mathcal K|/\rho_{\min}\}\).
The identical formula with \(\mathcal P\) bounds \(\nabla\Pi_2\).
These estimates include all first spatial directions and both mixed tensor
index orders in (7). They are not merely coordinate radial derivatives.

The source profiles remain smooth and flat at the packet endpoints. Their
zero extensions and first covariant derivatives are zero outside the closed
annulus, so the same bounds hold for the entire initial section wherever
those extended correction fields are defined. Division by \(p\) is used
only on a neighborhood of the packet where its gap is certified, not at a
possible remote vertex or outside the controlled branch.

## 6. Reproduction and precise scope

The producer imports only the Python standard library and reads three
authenticated input files. It does not import or execute previous scientific
producer modules. From this directory run

```sh
python check_second_initial_norms.py --run-controls --output SECOND_INITIAL_NORMS_NORMAL.json
python -O check_second_initial_norms.py --run-controls --output SECOND_INITIAL_NORMS_OPTIMIZED.json
```

The normal and optimized runs return the same scientific bounds. Six
mutation-specific controls flip a Hamiltonian reduction cross term, drop the
momentum cancellation term, reverse proper radial conversion, drop the mixed
projector connection, confuse the second Taylor derivative with its half
coefficient, or omit the physical bump scale. Each must fail for its expected
reason and leave no output report. No correctness check uses Python `assert`.

The prior constraint theorem identifies these explicit profiles as a valid
second initial constraint jet. The new result is a quantitative value and
first spatial derivative certificate for those profiles, conditional on the
retained actual background input. The large upper bounds reflect a narrow
packet, derivatives and conservative interval products. They give no allowed
finite \(\epsilon\), no positivity neighborhood for an actual finite
perturbation, no harmonic-gauge initial completion, no Hs or Einstein energy,
no evolved U2, no invariant force and no nonlinear remainder estimate.
