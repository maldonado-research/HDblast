# Registered longitudinal bulk sources: derivative reduction and finite envelopes

15 October 2026. Conditional mathematical continuation of the retained actual
aligned quadrupole. No second coefficient or gauge-invariant force is solved.

The new calculation gives coarse, continuous upper bounds for the complete
quadratic **longitudinal** Einstein and scalar residuals by combining certified
master jets and actual reconstruction coefficients. On the complete triangle
and sphere defined below, each individual adapted-frame Einstein component
and the scalar obey

\[
\boxed{|S^E_{AB}|<1.2\times10^{21},\qquad |S^\phi|<10^{13}.}
\]

The exact rational majorants are approximately \(1.07293\times10^{21}\)
and \(5.84252\times10^{12}\). They bound one-copy source coefficients in
registered model units; both identical copies satisfy the same bounds.
They do not give a source waveform, source signs, a normal-sum junction
coefficient, or a norm of a solved spacetime response. The generic
jet checks are a separate verification of the tensor formulas. Neither those
checks nor a master norm alone would establish the actual-source bound.

## Conventions and full-component formulas

Use the retained equations

\[
\mathscr E_{AB}=R_{AB}-\phi_A\phi_B-\tfrac23U g_{AB},\qquad
\mathscr S=\Box_g\phi-U',
\]

on \(g_0=dy^2+\bar h\), \(\bar h=\rho^2\gamma\),
\(\gamma=-dt^2+\cosh^2t\,\kappa_{S^3}\). The letters \(a=\rho_y/\rho\)
and \(p=\phi_{0,y}\) retain their background meaning. Let
\(\psi=uY\), \(\chi=vY\), \(Y=n_1^2-n_2^2\), with

\[
h_{yy}=-4\psi,\quad h_{ya}=0,\quad h_{ab}=2\psi\bar h_{ab}.
\]

Here \(a,b\) as indices range over the four slice coordinates. Write
\(D\), \(\Box_{\bar h}\), and \(|D\psi|^2\) for the slice connection,
wave operator, and Lorentzian contraction. The Hessians below are second
parameter derivatives along the **affine** path
\((g_0+\epsilon h,\phi_0+\epsilon\chi)\); the second Taylor coefficient
of a solution obeys \(D\mathscr E[k,w_2]=-S^E\),
\(D\mathscr S[k,w_2]=-S^\phi\). In particular each displayed source includes
the factor two relative to an \(\epsilon^2\) coefficient.

Define

\[
K=8\psi_y^2+104a\psi\psi_y+8\psi\psi_{yy}
       +48(a_y+4a^2)\psi^2.
\]

Then the full, unprojected sources are

\[
\begin{aligned}
S^E_{yy}={}&-8\psi_y^2-32a\psi\psi_y+16\psi\psi_{yy}
 +16|D\psi|^2-8\psi\Box_{\bar h}\psi
 -2\chi_y^2-\tfrac23U''\chi^2+\tfrac{16}3U'\psi\chi,\\
S^E_{ya}={}&12\psi\partial_a\psi_y-48a\psi\partial_a\psi
 -2\chi_y\partial_a\chi,\\
S^E_{ab}={}&24\psi D_aD_b\psi+12\partial_a\psi\partial_b\psi
 -2\partial_a\chi\partial_b\chi\\
&+\bar h_{ab}\left[4\psi\Box_{\bar h}\psi+4|D\psi|^2-K
 -\tfrac23U''\chi^2-\tfrac83U'\psi\chi\right].
\end{aligned}                                                     \tag{1}
\]

All fifteen independent Einstein components are present. This extends no
second-order scalar-master ansatz to the unknown second metric. The tensor
source can subsequently be decomposed into its finite angular sectors.
The formulas agree with the retained degree-four source rows; their
unprojected form is useful for a single whole-sphere envelope.

For completeness, before using the background or first scalar equation,

\[
\begin{aligned}
S^\phi={}&8\psi(\chi_{yy}+4a\chi_y)-4\psi\Box_{\bar h}\chi
 +12\psi_y\chi_y+48p\psi\psi_y\\
&+32(p_y+4ap)\psi^2-U'''\chi^2.
\end{aligned}                                                     \tag{2}
\]

Equations (1) and (2) hold for arbitrary smooth \(\psi,\chi\) in this
representative, with no master substitution. The background identity
\(p_y+4ap=U'\) and actual first scalar equation

\[
\chi_{yy}+4a\chi_y+\Box_{\bar h}\chi-U''\chi
                  +6p\psi_y+4U'\psi=0
\]

give the cancellation

\[
\boxed{S^\phi=-12\psi\Box_{\bar h}\chi
               +8U''\psi\chi+12\psi_y\chi_y-U'''\chi^2.}          \tag{3}
\]

The potentially large \(p\psi\psi_y\) and \(U'\psi^2\) contributions
cancel exactly, rather than being bounded individually. For the actual
quadrupole the **entire scalar source** is \(s^\phi(t,y)Y^2\), with

\[
\boxed{s^\phi=12u\rho^{-2}(v_{tt}+3\tanh t\,v_t
                  +8\operatorname{sech}^2t\,v)
               +8U''uv+12u_yv_y-U'''v^2.}                        \tag{4}
\]

Thus its degree-zero/two/four coefficients use the retained decomposition
\(Y^2=1/6+Z_2/2+H_4\). Equation (4) supplies a coefficient bound, so
the same absolute bound applies to the full scalar field since \(|Y|\le1\).

## Why the scalar still requires third master jets

Write \(k=\alpha_z/\alpha\), \(c(t)=\cosh^{-3/2}t\),
\(C=\alpha/\rho^2\), \(B=3\alpha/(\rho^3p)\). The physical
reconstruction is exactly

\[
u=Ccw,\qquad v=Bc(w_z+kw),\qquad \partial_y=-\rho^{-1}\partial_z.
\tag{5}
\]

Equation (1) uses \(u\) through order two and \(v\) through order one.
It therefore needs master derivatives only through order **two**. Equation
(4) instead has \(v_{tt}\), whose highest derivative is \(Bc\,w_{ttz}\).
The master equation implies

\[
w_{ttz}=w_{zzz}-\mathcal W_z w-\mathcal W w_z,
\]

which exchanges one third derivative for another. It does not bound it by
the two-jet. At the principal-jet level, fixing all derivatives through
order two leaves \(w_{zzz}\) free, while the differentiated equation only
fixes the difference \(w_{ttz}-w_{zzz}\). The coefficient of that free
third jet in (4) is \(12uBc/\rho^2\), nonzero wherever \(u\ne0\).
Thus this substitution does not algebraically lower the scalar-source
requirement below C3 for the actual smooth first field.

There is a distinct Noether observation. For the trace-reversed residual,
the exact identity is

\[
\nabla^A\mathscr E_{AB}-\tfrac12\partial_B
 (g^{AC}\mathscr E_{AC})+\mathscr S\,\partial_B\phi=0.
\]

On a background and first variation satisfying **all** their equations,
its quadratic coefficient is

\[
\nabla_0^A S^E_{AB}-\tfrac12\partial_B
 (g_0^{AC}S^E_{AC})+S^\phi\partial_B\phi_0=0.                      \tag{6}
\]

Because \(p\ne0\), the scalar source can be recovered from the radial
row of (6). This still differentiates the Einstein source and therefore
still uses third master derivatives. A sufficiently regular solution of
the **full** inhomogeneous Einstein equations would imply the scalar
equation by this identity; this does not by itself establish a reduced
hyperbolic solve or propagate its constraints and shell rows. Noether
compatibility is an analytic consequence of the actual on-shell first
field, not something proved by generic finite point tests.

## Exact envelope calculation

The domain is the backward-closed triangle
\(0\le t\le3\ell\), \(z\ge0\), \(t+z\le8\ell\),
\(\ell=10^{-6}\), times the entire unit sphere. Physical order-\(j\)
master derivatives equal \(\ell^{-j}\) times the certified scaled
derivatives of \(w(\ell x,\ell r)\). Losing these powers would invalidate
the bulk source enclosure.

The actual background report supplies rational intervals for
\(C,C_z,C_{zz},B,B_z,k,k_z,\rho,a,a_z,U',U'',U'''\). Its \(h\)
means \(U'/p\), so the correct identity is
\(k=\rho(5a/2-h)\). In particular no substitution of \(h=p_y/p\)
is made. The checker differentiates (5), uses
\(|c|\le1\), \(|c'/c|\le3T/2\), \(|c''/c|\le3/2\), and

\[
u_{yy}=\rho^{-2}u_{zz}+a\rho^{-1}u_z.
\]

This gives explicit global bounds for
\(u,u_t,u_{tt},u_y,u_{yy},u_{yt},v,v_t,v_{tt},v_y\), recorded in
`ACTUAL_LONGITUDINAL_BULK_SOURCE_BOUNDS.json`. Correlations discarded in
these inequalities make the final bounds deliberately coarse.

For reference the simple strict amplitude bounds checked in the receipt are

| Amplitude | Strict absolute upper bound |
|---|---:|
| \(u\) | \(1/16000\) |
| \(u_t\) | \(1000\) |
| \(u_{tt}\) | \(1.6\times10^{10}\) |
| \(u_y\) | \(8\) |
| \(u_{yy}\) | \(10^6\) |
| \(u_{yt}\) | \(1.25\times10^8\) |
| \(v\) | \(190000\) |
| \(v_t\) | \(3\times10^{12}\) |
| \(v_{tt}\) | \(8\times10^{19}\) |
| \(v_y\) | \(2.4\times10^{10}\) |

At every point take the adapted background frame
\(n=\partial_y\), \(e_0=\rho^{-1}\partial_t\),
\(e_i=(\rho\cosh t)^{-1}E_i\), where \(E_i\) is any orthonormal
spherical frame. The bounds apply individually to every component in
such a frame; arbitrary Lorentz boosts are not part of the norm. On S3,
\(|Y|\le1\), \(|\nabla Y|\le2\), and
\(\|\operatorname{Hess}Y\|_{\rm op}\le4\), from the restricted
quadratic form with ambient operator norm one. For example,

\[
\begin{aligned}
|D\psi|^2_{\rm absolute}&\le\rho^{-2}(U_t^2+4U_0^2),\\
|D_a\psi|&\le\rho^{-1}\max(U_t,2U_0),\\
|D_aD_b\psi|&\le\rho^{-2}\max\{U_{tt},
 2(U_t+TU_0),4U_0+TU_t\},\\
|\Box_{\bar h}\psi|&\le\rho^{-2}(U_{tt}+3TU_t+8U_0).
\end{aligned}
\]

Here capital \(U\) symbols denote already certified absolute amplitude
bounds in this paragraph only, not the bulk potential. The first line
bounds the absolute value of a Lorentzian contraction, not a positive
Lorentzian norm. Triangle inequalities in (1) and (4) then give the
actual registered source bounds in the JSON receipt.

## Verification and remaining work

`check_longitudinal_sources.py` uses only Python's standard library. It
constructs three rational warped background two-jets with arbitrary nonzero
tangential metric derivatives and unrestricted first-field two-jets. A
direct inverse-metric/Christoffel/Ricci epsilon expansion independently
checks all fifteen curvature Hessian components in (1); the direct scalar
wave expansion checks (2), and exact substitution checks (3). These 51
generic tests supplement the analytic derivation; they are not an all-point
proof for the actual waveform. The actual global enclosure instead follows
from the authenticated master/background bounds and the explicit inequalities.
Mutations exercise the Ricci factorial, mixed-row term, scalar cancellation,
and physical derivative scale under ordinary and optimized Python.

From a capsule root containing `potential/`, `jets/`, and `bulk/`, run

```sh
python bulk/check_longitudinal_sources.py --run-controls --output bulk/ACTUAL_LONGITUDINAL_BULK_SOURCE_BOUNDS.json
python -O bulk/check_longitudinal_sources.py --run-controls --output bulk/ACTUAL_LONGITUDINAL_BULK_SOURCE_BOUNDS_OPTIMIZED.json
```

`bulk/INPUT_PINS.json` authenticates the two exact report inputs, including
the master report's own hash link to the background report. The default paths
are relative to the checker location, and no prior scientific producer is
imported or executed. The master and background checkers should be replayed
first to audit the analytic and numerical dependencies of these final
envelopes. Use `--identities-only` only to request the generic finite-jet
checks separately; that mode explicitly reports no actual source enclosure.

The source is in the same **first longitudinal representative** as the
retained coordinate shell sources. It is not the source of the first
relative-harmonic representative: that requires solving the first gauge
vector, transforming both metric and scalar fields, and then evaluating
the corresponding Hessian and quadratic harmonic-vector terms. Source
components and these bounds are coordinate dependent. Material boundary
terms or first-gauge transformed bulk terms cannot be interchanged without
the full covariance and boundary conversion.

The quantitative bounds do not compute a particular second initial
constraint correction, transformed physical/gauge/constraint boundary
forcing in a common solve, a continuous residual certificate for the
second coupled system, U2, an invariant force, a nonlinear remainder or an
effective finite-amplitude family. They supply actual first-field bulk
forcing envelopes for that future calculation.
