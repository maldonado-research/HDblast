# Methods for a rigorous HDBLAST shell-waveform certificate

This note independently derives an error estimate for the HDBLAST linear scalar sector. No theorem is transferred from another research project. It supplies a conditional estimate for an existing sufficiently regular solution, together with a plan for certifying its inputs. It does not establish a nonlinear gravitational solution, a tame inverse estimate, or the remaining gauge, support and evolution hypotheses.

## 1. Specify the model, data and observation

On a collar `0 <= z <= Z`, write the target linear equation and inertial boundary condition as

\[
w_{tt}-w_{zz}+W(t,z)w=0,\qquad
-w_z(t,0)+\mu w_{tt}(t,0)+\beta(t)w(t,0)=0,
\quad \mu>0.
\]

The potential is allowed to depend on time. The inertial mass `mu` is constant. Fix the initial bulk displacement and velocity, their boundary traces and compatibility conditions, and the map from `w` to the requested shell observable before producing a numerical certificate.

Let `e = w - w_app`, and define its **actual** residuals

\[
r=e_{tt}-e_{zz}+We,\qquad
q=-e_z(t,0)+\mu e_{tt}(t,0)+\beta e(t,0).
\]

If `w` solves the displayed homogeneous problem, these are the negatives of the corresponding approximate-solution residuals. Norm bounds are insensitive to that sign, but the energy identity below uses this definition. The initial error is an independent input; a small equation residual does not constrain a homogeneous initial-state error.

Assume the error and its needed traces are sufficiently regular, `W` and `beta` are differentiable in time to the extent needed below, and the far endpoint contributes no flux. One sufficient condition is that both target and approximate waves vanish near `Z` for the entire certified interval. Initial support inside `[ell,2ell]` and a separately justified unit propagation speed keep the outer front below `2ell+T`. For `ell=10^-6`, `T=3ell`, and `Z=8ell`, this is strictly inside the collar. This is a conditional support statement for the linear sector; it requires no assertion about the background beyond that collar.

## 2. Include the boundary kinetic energy and both coefficient derivatives

Define the quadratic form

\[
Q_t[e]=\int_0^Z (|e_z|^2+W|e|^2)\,dz
       +\beta(t)|e(t,0)|^2
\]

and energy

\[
E(t)=\frac12\left[
\int_0^Z|e_t|^2\,dz+\mu|e_t(t,0)|^2+Q_t[e]
\right].
\]

For real fields, multiply the equation by `e_t` and integrate by parts. For complex fields take real parts. Using the boundary residual gives the exact identity

\[
\boxed{
E'(t)=\operatorname{Re}\int_0^Z r\overline{e_t}\,dz
+\operatorname{Re}\bigl(q\overline{e_t(t,0)}\bigr)
+\frac12\int_0^Z W_t|e|^2\,dz
+\frac12\beta'(t)|e(t,0)|^2.
}
\]

The boundary kinetic term cancels `mu e_tt e_t` from the boundary work. If the far endpoint does not vanish, its flux `Re(e_z(t,Z) overline(e_t(t,Z)))` must be retained and bounded. A finite-domain computation cannot discard that term merely because its numerical residual is small. This is the scalar-sector energy, not a gravitational energy.

## 3. Prove coercivity on the full certified domain

Suppose actual coefficient bounds give

\[
W(t,z)\ge w_0>0,\qquad
\beta(t)\ge b_*>-\sqrt{w_0}
\]

on the full collar and time interval, and `e(t,Z)=0`. The one-dimensional trace estimate is

\[
|e(t,0)|^2\le2\|e\|_2\|e_z\|_2
\le\frac{1}{\sqrt{w_0}}
\int_0^Z(|e_z|^2+w_0|e|^2)\,dz.
\]

Consequently, with

\[
\rho=\min\{1,1+b_*/\sqrt{w_0}\}>0,
\qquad C_{\rm tr}=\frac1{\sqrt{w_0}+b_*},
\]

one has

\[
Q_t[e]\ge\rho\int_0^Z(|e_z|^2+w_0|e|^2)\,dz,
\qquad |e(t,0)|^2\le C_{\rm tr}Q_t[e].
\]

These are sufficient coercivity hypotheses. If they fail, one must establish a different coercive energy and its trace constants before claiming an error bound.

Put `A(t)=sqrt(2E(t))`, and let

\[
F(t)=\|r(t,\cdot)\|_{L^2(0,Z)}+\frac{|q(t)|}{\sqrt\mu},
\]

\[
G(t)=\frac12\left[
\frac{\|W_t(t,\cdot)\|_\infty}{\rho w_0}
+C_{\rm tr}|\beta'(t)|\right].
\]

Cauchy–Schwarz and the preceding form bounds imply `A' <= F+GA` wherever `A>0`; the integrated inequality extends through zeros by regularization. Thus, with `H(t)=int_0^t G(s) ds`,

\[
\boxed{
A(t)\le e^{H(t)}\left[A(0)+\int_0^t e^{-H(s)}F(s)\,ds\right],
\qquad
|e(t,0)|\le\sqrt{C_{\rm tr}}A(t).
}
\]

This estimate converts a complete bulk-and-boundary residual bound plus an initial-state enclosure into a shell-trace error. Further transformations to the physical observable need their own bounded error map.

## 4. Conditional specialization to the wider collar bounds

If the HDBLAST coefficient audit supplies, throughout the stated collar and interval,

\[
W\ge10000,\quad \beta\ge200,\quad\mu>1/500,
\quad |W_t|<21/400000,
\quad |\beta'|<21/160000000,
\]

then `rho=1`, `C_tr=1/300`, and

\[
G\le\frac12\left[
\frac{21}{400000\cdot10000}
+\frac{21}{160000000\cdot300}\right]
=\frac{91}{32000000000}.
\]

For `T=3*10^-6`, `GT <= 273/32000000000000000`, so `exp(GT)<1+10^-14`. The elementary estimate `exp(x)<=1/(1-x)` for `0<=x<1` verifies this using rational comparisons. Since `1/sqrt(mu)<sqrt(500)<23`, a conservative bound is

\[
\boxed{
A(t)\le(1+10^{-14})\left[
A(0)+\int_0^t(\|r(s,\cdot)\|_2+23|q(s)|)\,ds
\right],
\quad |e(t,0)|\le A(t)/\sqrt{300}.
}
\]

This specialization depends on the coefficient, regularity and support premises; it does not supply them. Small coefficient growth does not ensure a small waveform error unless the initial error and both residual integrals are small.

## 5. Certificate design

The following practices can be implemented directly on the HDBLAST equations without importing any external result:

- Use rational parameter cells or outward interval enclosures covering the entire target domain. A sampled plot does not prove between-sample positivity. A cell-center margin propagates only after a HDBLAST-specific Lipschitz or derivative bound has been established.
- Certify polynomial signs over whole cells, for example with exact Bernstein coefficients; enclose every nonpolynomial remainder and spatial tail separately. Failure of a sufficient sign test is inconclusive.
- Record the model, units, domain, coefficient bounds, initial error, bulk residual, boundary residual, observation error map, and all assumed support/regularity hypotheses in the result. Preserve negative or unresolved margins instead of silently promoting them to success.
- Audit the boundary sign and time-derivative terms using manufactured waves by a calculation independent of the production discretization. Keep immutable input hashes, arithmetic precision/directions, executed checks and code revision in a reproduction receipt.

For the nonlinear remainder, the missing step is a proved map from the nonlinear defect to a forcing norm controlled by an appropriate solution norm, without unaccounted derivative loss. Only then can a bootstrap or comparison argument close a remainder estimate. The displayed linear energy alone does not establish that map, the localized inverse, full nonlinear support, a gauge construction or nonlinear evolution.

Positive spectral-measure methods would require a separately constructed fixed nonnegative self-adjoint generator and verified common Laplace-moment bounds. Time-dependent real-time waveform samples do not supply those premises. They therefore add no theorem to this shell-waveform certificate.
