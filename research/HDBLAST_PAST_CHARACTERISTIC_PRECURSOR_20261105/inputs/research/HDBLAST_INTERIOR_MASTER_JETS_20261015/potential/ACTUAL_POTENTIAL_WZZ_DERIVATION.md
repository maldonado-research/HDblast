# Retained actual potential: a rational second spatial derivative certificate

15 October 2026. Research continuation, separate from the immutable
coordinate-shell package. This is a conditional coefficient enclosure for the
retained actual regular background. It is not a numerical background trajectory,
an interior wave approximation or a nonlinear existence claim.

The certified bounds, valid on the rectangle `0<=t<=3ell, 0<=z<=8ell` and
therefore on the backward-closed triangle `z+t<=8ell`, are

\[
\ell=10^{-6},\qquad
|\mathcal W|<102000,\quad |\mathcal W_z|<30000000,\quad
|\mathcal W_t|\le\frac{21}{400000},\quad
\boxed{|\mathcal W_{zz}|<30000000000.}
\]

All decisive comparisons use exact fractions. The exact interval produced for
`Wzz` is in `ACTUAL_POTENTIAL_JET_BOUNDS.json`; its decimal display is roughly
`[-12264849902.54,20065346881.93]`. No rounded endpoint is used to prove the
simple bound.

## Actual input identity and first exit

`INPUT_PROVENANCE.json` records byte-identical copies with source locations,
sizes and SHA-256 hashes. The checker pins the actual coordinate-shell input
`BACKGROUND_DERIVATIVE_BOUNDS.json` to
`fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c`.
It was also compared byte for byte with the original inertial-waveform archive
copy. The registered analytic branch theorem, exact polynomial certificate,
prior shell enclosure and actual-background audit are copied and authenticated.

The retained model has `eta=phi-1` and

\[
W_{\rm super}=\tfrac13+\eta^2+\tfrac13\eta^3,\quad
U=\tfrac12(W_{{\rm super},\eta})^2-\tfrac23W_{\rm super}^2,
\]
\[
c=\frac{96714321425692}{161839258930731},\quad\delta=\frac1{1000},
\qquad\sigma=2W_{\rm super}+\delta(1+c+c\eta).
\]

These are the equations in the authenticated
`PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md`, Section 1.
The retained branch theorem is the mathematical input specifying this actual
solution. The checker verifies its model identity by reconstructing the shell
Hamiltonian polynomial cancellation and matching its coefficients with the
registered exact certificate. No similarly named historical model is used.

Inward depth satisfies `dy=-rho dz`, and the ODE is

\[
\rho_z=-\rho^2a,\quad\eta_z=-\rho p,\quad
p_z=-\rho(U_1-4ap),\quad a_z=\rho^{-1}+\rho p^2/3.
\]

Here `U_j` means the j-th eta derivative. The outer rectangle is

\[
129\le\rho\le131,\quad-9/10^5\le\eta\le-8/10^5,\quad
-14/10^5\le p\le-12/10^5,\quad111/1000\le a\le112/1000.
\]

The new checker recomputes all four ODE images and requires exact endpoint
equality with the actual retained audit. Starting from the authenticated prior
shell box, it independently verifies
`I_shell+[0,8ell]F(I_outer)` lies strictly inside `I_outer` in every component.
The usual first-exit argument and smooth ODE continuation prove the entire
collar exists in the rectangle. Also `131*8ell<9`, the retained proper endpoint
lower bound. All denominators stay nonzero. This proof is interval coverage
of the whole collar; no sampled trajectory or numerical integration is used.

## Derivation through Wzz

Write

\[
h=U_1/p,\quad L=U_2-h^2+4ah,\quad
Q=2h^2-10ah+35a^2/4-U_2-p^2/6,
\quad\mathcal W=\rho^2Q-35\tanh^2(t)/4.
\]

The ODE gives

\[
h_z=-\rho L,\qquad
h_{zz}=-\rho_zL-\rho[U_3\eta_z+(-2h+4a)h_z+4a_zh],
\]
\[
\rho_{zz}=2\rho^3a^2-\rho^2a_z,
\quad a_{zz}=-\rho_z/\rho^2+\rho_zp^2/3+2\rho pp_z/3,
\]
\[
p_{zz}=-\rho_z(U_1-4ap)-\rho(U_2\eta_z-4a_zp-4ap_z).
\]

Differentiate the retained first derivative exactly:

\[
Q_z=(4h-10a)h_z+(-10h+35a/2)a_z+\rho pU_3-pp_z/3,
\]
\[
\begin{aligned}
Q_{zz}={}&4h_z^2+(4h-10a)h_{zz}-20a_zh_z+\tfrac{35}{2}a_z^2
 +(-10h+\tfrac{35}{2}a)a_{zz}\\
&+(\rho_zp+\rho p_z)U_3+\rho pU_4\eta_z
 -(p_z^2+pp_{zz})/3.
\end{aligned}
\]

The crucial product-rule factor is then

\[
\mathcal W_z=-2\rho^3aQ+\rho^2Q_z,
\]
\[
\boxed{\mathcal W_{zz}
=(6\rho^4a^2-2\rho^3a_z)Q-4\rho^3aQ_z+\rho^2Q_{zz}.}
\]

All temporal dependence is explicit, so
`|W_t|=(35/2)|tanh(t)|sech(t)^2 <= (35/2)T`.
The elementary estimates `|tanh(t)|<=t` and `sech(t)^2<=1` avoid any numerical
transcendental evaluation. For `W` itself the checker encloses
`rho^2 Q - [0,(35/4)T^2]`.

The new program evaluates all displayed expressions with rational intervals.
It first requires exact equality with the retained `h,h_z,Q,Q_z,W_z` endpoint
pairs. Thus the extension is anchored to the actual audited potential, rather
than merely to a polynomial with the same name.

## Independent formal product-rule check and negative controls

In addition to interval evaluation, the checker contains a separate sparse
Laurent polynomial differentiator over seven formal variables
`rho,a,p,h,U2,U3,U4`. Its total derivative substitutes the background ODE,
`U2_z=-rho*p*U3` and `U3_z=-rho*p*U4`. This does not sample points.
It proves exact zero polynomial defects for `hzz,azz,pzz,Qz,Qzz,Wz,Wzz`.
The formal `pzz` check uses the ratio form `p_z=-rho*p*(h-4a)`, while the
interval `pzz` uses the original `U1` expression. They agree because `h=U1/p`.

`run_potential_controls.py` verifies normal and optimized Python produce
identical complete report bytes. Both modes must reject each of:

- a changed authenticated input byte;
- halving the Wzz cross coefficient from `-4` to `-2`;
- deleting the Qzz mixed term `-20*a_z*h_z`;
- claiming the false bound `|Wzz|<1`.

The rejection checks require the intended exception text, a nonzero status,
and absence of an output report. They remain active under `python -O` because
the checker does not use `assert` for validation.

## Background coefficient jets for first-field reconstruction

The same report gives exact intervals useful for converting certified master
jets into the first longitudinal fields. Put

\[
\alpha=(\rho/\rho_b)^{3/2}p/p_b,\qquad
k=\alpha_z/\alpha=\rho(5a/2-h),\qquad q(t)=\cosh^{-3/2}t,
\]
\[
C=\alpha/\rho^2,\quad B=3\alpha/(\rho^3p),\qquad
u=Cqw,\quad v=Bq(w_z+kw).
\]

The report contains `k,k_z,k_zz`, `C,C_z,C_zz` and `B,B_z,B_zz`.
For example

\[
C_z=lC,\quad C_{zz}=(l_z+l^2)C,\quad l=\rho(9a/2-h),
\]
\[
B_z=mB,\quad B_{zz}=(m_z+m^2)B,\quad m=3\rho a/2.
\]

Since `alpha(0)=1` and `|k|<=K`, `KZ<1` and the integral logarithm identity
give `1-KZ <= alpha <= 1/(1-KZ)`. This uses the elementary exponential
majorant and avoids fractional powers or rounded shell constants.
The temporal factor obeys
`|q|<=1`, `|q_t|<=3T/2`, `|q_tt|<=3/2+9T^2/4`.
Derivative conversion is
`partial_y=-rho^-1 partial_z` and
`partial_y^2=rho^-2 partial_z^2+(a/rho)partial_z`.

These are coefficient bounds only. Actual master jets, common gauge, angular
factors, complete tensor/scalar residuals and their physical signs must still
be supplied by the interior/source continuation before any bulk-source claim.

## Reproduction

Run from this directory, using only the Python standard library:

```sh
python -I -S -B check_actual_potential.py
python -I -S -B run_potential_controls.py
```

`ACTUAL_POTENTIAL_JET_BOUNDS.json` and `POTENTIAL_CONTROL_RESULTS.json` are
the exact numerical and verification receipts, respectively.
Both scripts accept `--output PATH`. For replay from a read-only source
capsule, supply output paths in a writable replay directory. The controls
runner creates its scratch directory in the system temporary directory and
writes its control receipt only to the selected `--output` path.
