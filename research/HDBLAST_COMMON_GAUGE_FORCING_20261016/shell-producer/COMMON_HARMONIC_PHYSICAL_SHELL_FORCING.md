# Full physical shell forcing in the common first harmonic gauge

16 October 2026. Mathematical producer; quantitative report and independent
review accompany this note. This is an affine quadratic source for the
actual registered first pair, not a second solution or an invariant force.

## 1. Registered signs, one tension and source normalization

The retained `SCALAR_SHELL_BOUNDARY_SYSTEM.md`, sections 1--2, states that
increasing proper \(y\) is outward on both identical interiors and that
\(a_b=\sigma_0/6\), \(p_b=-\sigma'_0/2\). The same signs follow from the
untraced Israel and scalar residuals in `MATERIAL_SECOND_ORDER_JUNCTIONS.md`,
section 2:

\[
I_{ab}=\sum_i(K^{\rm out}_{i,ab}-K_i^{\rm out}\gamma_{ab})
                  +\sigma\gamma_{ab},\qquad
P=\sum_i n_i^{\rm out}\phi+\sigma'.
\]

Trace reversal in four intrinsic dimensions gives
\(\sum_iK^{\rm out}_{i,ab}-\sigma\gamma_{ab}/3=0\). There is one common
tension. For the registered copy-even state the convenient one-copy residuals
are

\[
\mathcal B_{ab}=K^{\rm out}_{ab}-\sigma(\varphi)\gamma_{ab}/6,
\qquad \mathcal B_\phi=n^{\rm out}\phi+\sigma'(\varphi)/2.
\tag{1}
\]

The two-copy outward residuals are exactly \(\mathcal J=2\mathcal B\).
The inward normals are \(n_i^{\rm in}=-n_i^{\rm out}\), so the corresponding
two-copy equations are
\(\sum_iK^{\rm in}_{i,ab}+\sigma\gamma_{ab}/3=0\),
\(\sum_i n_i^{\rm in}\phi-\sigma'=0\); their residuals are
\(-2\mathcal B\). These are equivalent orientations, not a tension change.

The actual fixed registered tension has \(\sigma''_b=4\varphi_b\) and
\(\sigma'''_b=4\), including the detuning term whose first derivative
enters \(\sigma'_b=-2p_b\). Their numerical intervals are retained in the
authenticated coordinate-shell source report. No tension parameter varies.

Taylor coefficients here are parameter derivatives:
\(\Phi_\epsilon=\Phi_0+\epsilon\Phi_1+\epsilon^2\Phi_2/2\).
Thus the affine Hessians \(Q_K=D^2\mathcal B\),
\(Q_p=D^2\mathcal B_\phi\) enter the second equations as
\(D\mathcal B[\Phi_2]=-Q_K\), \(D\mathcal B_\phi[\Phi_2]=-Q_p\).
The two-copy outward Hessians are \(2Q_K,2Q_p\), and inward Hessians their
negatives. All adapters are stated explicitly below.

These ten metric rows are the invertible trace-reversed Israel equations.
For the original untraced two-copy Israel residual \(I\), its actual affine
Hessian is \(Q_I=2[Q_K-(\operatorname{tr}_{\gamma_0}Q_K)\gamma_0]\).
There are no varying-trace cross terms left here because the background and
actual first residuals vanish. No physical metric row is lost by using (1).

## 2. Actual first harmonic trace and fixed-shell affine path

Let \(q=h+\mathcal L_\xi g_0\), \(\eta=\chi+\mathcal L_\xi\phi_0\), where
the raw actual gauge amplitudes are
\(\xi^y=fY\), \(\xi^t=bY\), \(\xi^i=dD_\kappa^iY\),
\(Y=n_1^2-n_2^2\), and \(c=\cosh t\), \(T=\tanh t\).
The physical shell rows \(f|_b=0\), \(b_y|_b=d_y|_b=0\) imply
\(q_{ya}|_b=0\). With \(H=q_{yy}|_b\) and \(m_{ab}=q_{ab}|_b\), the
background-shell orthonormal components are

\[
\begin{aligned}
H&=JY,\quad J=-4u+2f_y,\qquad \eta=vY,\\
m_{00}&=-2X Y,\quad X=u+b_t,\\
m_{0a}&=V D_aY,\quad V=cd_t-b/c,\\
m_{ab}&=2[Z Y\delta_{ab}+dD_aD_bY],\quad Z=u+Tb.
\end{aligned}                                                     \tag{2}
\]

Here \(e_0=\rho_b^{-1}\partial_t\), \(e_a=(\rho_bc)^{-1}E_a\);
\(\gamma\) has frame matrix \(\operatorname{diag}(-1,1,1,1)\).
The tangential tensor is generally anisotropic. The scalar trace is unchanged
by the first normal gauge displacement because \(f|_b=0\).

For the affine path use the fixed coordinate shell, \(g_{ya}=0\) there for
every parameter value, \(g_{yy}=1+\epsilon H\),
\(\gamma_\epsilon=\gamma_0+\epsilon m\), and
\(\phi_\epsilon=\phi_0+\epsilon\eta\). Its lapse is
\(N=\sqrt{1+\epsilon H}\) and outward unit normal is \(N^{-1}\partial_y\).
Although \(g_{ya}\) can be nonzero away from the shell, its tangential
derivatives vanish on the shell, so the exact boundary expressions are

\[
K^{\rm out}_{ab}=\frac{\partial_y\gamma_{\epsilon,ab}}{2N},
\qquad n^{\rm out}\phi_\epsilon=
                    \frac{p+\epsilon\eta_y}{N}.               \tag{3}
\]

Differentiating gives the full tensor affine Hessians

\[
\begin{aligned}
Q_{K,ab}&=\tfrac34aH^2\gamma_{ab}-\tfrac12H(m_y)_{ab}
 -\tfrac16[\sigma''\eta^2\gamma_{ab}+2\sigma'\eta m_{ab}],\\
Q_p&=\tfrac34pH^2-H\eta_y+\tfrac12\sigma'''\eta^2.
\end{aligned}                                                     \tag{4}
\]

All indices in (4) are lower tangential indices. In particular the product
\(H m_y\) is not a trace replacement. Equation (3) accounts for the
normal normalization; using the coordinate vector \(\partial_y\) as the
perturbed normal would delete the first terms.

## 3. On-shell reduction eliminates the normal derivatives

The actual first pair satisfies the physical first rows in this same gauge:

\[
m_y=2am+(aH+\sigma'\eta/3)\gamma,
\qquad \eta_y=pH/2-\sigma''\eta/2.                            \tag{5}
\]

These follow either by differentiating (1), or by covariance of the physical
rows under the solved vector tangent to the fixed shell. Substitution into
(4), only at the shell, gives

\[
\boxed{Q_{K,ab}=A\gamma_{ab}-L m_{ab},\quad
A=\frac a4H^2-\frac{\sigma'}6H\eta-\frac{\sigma''}6\eta^2,
\quad L=aH+\frac{\sigma'}3\eta,}                              \tag{6}
\]

\[
\boxed{Q_p=\frac p4H^2+\frac{\sigma''}2H\eta
                               +\frac{\sigma'''}2\eta^2.}    \tag{7}
\]

Thus every actual physical source value needs only \(H,m,\eta\).
No normal derivative of the first fields remains. This is an on-shell
boundary simplification, not an identity on neighboring parallel leaves;
one must not differentiate (5) normally to estimate a source derivative.

Writing \(A=A_0Y^2\), \(L=L_0Y\), equations (2),(6) expose all ten rows:

\[
\begin{aligned}
Q_{K,00}&=(-A_0+2L_0X)Y^2,\\
Q_{K,0a}&=-L_0V\,Y D_aY,\\
Q_{K,ab}&=(A_0-2L_0Z)Y^2\delta_{ab}
                         -2L_0d\,Y D_aD_bY,\\
Q_p&=P_0Y^2.
\end{aligned}                                                     \tag{8}
\]

No angular or time-angular tensor row has been dropped. The longitudinal
specialization \(f=b=d=0\) gives
\(-Q_K=(\sigma''v^2/6-12au^2)Y^2\gamma\) and
\(-Q_p=(-4pu^2+2\sigma''uv-2v^2)Y^2\), exactly the old one-copy
coordinate coefficients. This comparison fixes the factorial and sign.

## 4. Second coordinate physical rows and proper derivatives

With the second mixed metric rows \(k_{ya}|_b=0\) imposed separately, the
coordinate physical equations, multiplied in the metric row by two, are

\[
k_{ab,y}-2ak_{ab}-(a k_{yy}+\sigma' w/3)\gamma_{ab}
                                      =-2Q_{K,ab},\qquad
w_y+\sigma''w/2+\sigma'k_{yy}/4=-Q_p.                         \tag{9}
\]

The metric source in (9) is also minus the two-copy outward metric residual
Hessian. The scalar row in (9) is the one-copy row, so its forcing is half
minus the two-copy scalar Hessian. This adapter prevents normal-sum
normalization from being silently transferred between equations.

Let \(x=t/\ell\), \(r=z/\ell\), \(\ell=10^{-6}\). The required field
conversion is
\(f_y=-f_r/(\ell\rho_b)\), \(f_{yt}=-f_{xr}/(\ell^2\rho_b)\),
\(b_t=b_x/\ell\), \(b_{tt}=b_{xx}/\ell^2\), and similarly for \(d\).
The proper shell clock derivative is \(D_0=\rho_b^{-1}\partial_t\).
Using (5) gives valid bounds for the *first field* proper-normal derivatives
without higher normal reconstruction. The harmonic row additionally gives
\(H_y=-4aH+2a\tau+4\sigma'\eta/3\), where
\(\tau=\operatorname{tr}_{\gamma}m\). In a radially parallel background frame,
\(\nabla_y m=(aH+\sigma'\eta/3)\gamma\); the \(2am\) term in (5) belongs
to the derivative of the coordinate warp factors. No normal derivative of
the shell-only source itself is assigned an intrinsic meaning.

The report encloses source values, proper shell time derivatives, and
full first tangential covariant derivatives, with explicit angular and
de Sitter connection terms. These are individual frame-component bounds,
not summed tensor norms or a Sobolev theorem.

## 5. Quantitative inputs and verification

The portable checker consumes authenticated copies of the actual global
coefficient report, actual master C3 report, actual gauge C3 report, and
actual coordinate-shell report containing tension and first shell trace
enclosures. The shell parameters must agree across their independent pins.
It never substitutes the free-wave comparator for the registered master.
The finite algebra checks compare all ten anisotropic rows and the scalar
row to an independent affine Christoffel/unit-normal expansion, as well as
the longitudinal specialization and one/two-copy adapters. These finite
generic-jet checks supplement (3)--(8); the continuous enclosure instead uses
the analytic inequalities below. The coordinate-shell input manifest is also
retained and authenticated: its common master inputs must equal the C3
master report's input pins. The independently retained shell boxes for
\(\rho_b,a_b,p_b\) must match exactly, as must
\(\sigma''_b=4(1+\eta_b^{\rm background})\).

For reproducibility write \(y(x)=w(\ell x,0)\) and
\(\mathcal I(x)=\ell(w_z-\gamma_b w_b)(\ell x,0)\), distinguishing
\(\gamma_b\) from the induced metric. The master report certifies suprema
of \(y,y',\mathcal I,\mathcal I'\). Since \(\alpha_b=1\),

\[
u=\rho_b^{-2}\vartheta y,\qquad
v=B_b\vartheta[\ell^{-1}\mathcal I+(\gamma_b+k_b)y],\qquad
B_b=\frac3{\rho_b^3p_b},\quad \vartheta=\cosh(t)^{-3/2}.
\tag{10}
\]

The exact old shell report also encloses \(u,v\); the checker takes the
smaller of the two valid absolute bounds. With \(t_* =3\ell\),
\(0<\vartheta\le1\), \(|\vartheta_t|\le3t_*/2\), it uses

\[
\begin{aligned}
|u_t|&\le\rho_{\min}^{-2}
       (\|y'\|/\ell+3t_*\|y\|/2),\\
|v_t|&\le |B_b|[\|\mathcal I'\|/\ell^2
 +( |\gamma_b|+|k_b|)\|y'\|/\ell\\
&\hspace{39mm}+\tfrac32t_*
 (\|\mathcal I\|/\ell+(|\gamma_b|+|k_b|)\|y\|)].
\end{aligned}                                                     \tag{11}
\]

Every gauge bound is taken from the actual continuous C3 report; an order
\(i+j\) scaled mixed derivative divides by \(\ell^{i+j}\). Equations
(2), (6)--(8), differentiated in time, then give exact rational amplitude
majorants. Normal conversion uses \(dy=-\rho\,dz\), including the factor
\(\rho_b^{-1}\) in \(f_{yt}\). These value bounds need gauge order one;
their first tangential derivatives need gauge order two and master order two.
The available C3 certificates are therefore more than the shell calculation
requires. No master or gauge order four is used.

For \(S=n_1^2+n_2^2\), the sphere identities are
\(|D Y|^2=4(S-Y^2)\),
\(|\operatorname{Hess}Y|^2=8-8S+24Y^2\). Hence

\[
|Y|\le1,\quad |DY|\le2,\quad |YDY|\le1,\quad
|(\operatorname{Hess}Y)_{ab}|\le4,\quad
|(D\operatorname{Hess}Y)_{cab}|\le8.
\tag{12}
\]

The last inequality follows directly from
\(D_cD_aD_bY=-\delta_{ca}D_bY-\delta_{cb}D_aY-2\delta_{ab}D_cY\).
Consequently \(|D(Y^2)|\le2\),
\(|D_c(YD_aY)|\le8\), and
\(|D_c(YD_aD_bY)|\le16\). The checker retains the exact trace contractions

\[
\begin{aligned}
\tau&=(8u+2b_t+6Tb-16d)Y,\\
\operatorname{tr}[(\gamma^{-1}m)^2]
&=(4X^2+12Z^2-64Zd)Y^2
                  -2V^2|DY|^2+4d^2|\operatorname{Hess}Y|^2.
\end{aligned}                                                     \tag{13}
\]

The negative time-angular term is a Lorentzian contraction. Its absolute
bound uses the triangle inequality; it is never treated as a positive norm.
The report exports these caps and their time derivatives for the separate
gauge-boundary calculation.

All first covariant source derivatives use the *induced background shell*
connection. In its frame, \(D_0 e_0=D_0e_a=0\),
\(D_a e_0=(T/\rho_b)e_a\), and
\(D_a e_b=(\rho_bc)^{-1}\omega^d{}_{ab}e_d
+(T/\rho_b)\delta_{ab}e_0\). Thus ordinary sphere derivatives divide by
\(\rho_bc\), and the time-angular boost terms are also included. In the
mixed source derivative the isotropic terms cancel exactly:
\(Q_{K,ca}+\delta_{ca}Q_{K,00}\) contains no \(A\) contribution, and
\(X-Z=b_t-Tb\). This avoids a large spurious mixed-component bound.

In the retained outward one-copy residual convention the resulting strict
uniform caps on \([0,3\ell]\times S^3\) are

| Source | Each frame component | Each first intrinsic tangential covariant derivative |
|---|---:|---:|
| Full ten-component \(Q_K\) | \(1.4512\,10^9\) | \(1.3485\,10^{15}\) |
| Scalar \(Q_p\) | \(4.3539\,10^9\) | \(4.0456\,10^{15}\) |

The time-angular metric source is below \(1.8984\,10^{-11}\), and each
off-diagonal purely angular source is below \(1.1391\,10^{-16}\). All ten
row-by-row fractions, their time/angular derivative bounds, and the two-copy
and coordinate-forcing adapters are retained in the JSON report. The metric
forcing of (9) has maximum component below \(2.9024\,10^9\); its scalar
forcing has the \(Q_p\) cap above. These are upper bounds, not lower bounds
or sign statements for those anisotropic rows.

`check_common_harmonic_shell.py` uses only the standard library, hard-pins
the manifest and each input, rejects existing output paths, and accepts
portable `--input-root`, `--pins`, and required `--output`. For example,

```sh
python check_common_harmonic_shell.py --run-controls --output /tmp/common-shell-normal.json
python -O check_common_harmonic_shell.py --run-controls --output /tmp/common-shell-optimized.json
```

Both modes must produce byte-identical reports. Controls require the
mutation-specific rejection reason and absence of an output file. The
reported arithmetic is exact `Fraction` arithmetic with outward dyadic
rounding on a \(2^{96}\) grid; displayed decimal caps are strictly larger
than the certified positive rational upper bounds. An exactly zero normal
field row is labeled as an exact zero rather than a strict positive cap.

The GH/gauge boundary rows and compatibility analysis are separate work.
These physical source rows do not prescribe a second gauge vector, establish
a compatible reduced solve, propagate constraints, or produce a second
metric/scalar response or a geometric force.
