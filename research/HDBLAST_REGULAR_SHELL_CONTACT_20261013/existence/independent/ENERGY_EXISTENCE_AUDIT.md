# Independent scalar energy and constructive-existence audit

Date: 2026-10-14 UTC. Scope: the canonical scalar equation only, with Minkowski signature `(-,+,...)`, outward spacelike normal, and `Box_eta psi - V(R^2-T^2) psi = 0` inside `R < sqrt(b^2+T^2)`, `b>0`. The even boundary condition is `n psi - c Box_Sigma psi + d psi = 0`, `c>0`; the odd condition is `n psi + d psi = 0`. No claim about a complete coupled gravitational or material system is made.

## 1. Direct fixed-ball energy, valid on every compact regular time slab

Let the spatial dimension be `m`, `a(t)=sqrt(b^2+t^2)`, `h=t/a^2`, `x=X/a(t)`, and `u(t,x)=psi(t,a(t)x)`. With `D=partial_t-h x.grad_x`, the bulk equation is

`D^2 u - a^(-2) Delta_x u + V(a^2|x|^2-t^2)u = 0`.

Its expanded principal part is `u_tt-2h x.grad u_t +(h^2 x_i x_j-a^(-2)delta_ij)u_ij`; its remaining first derivative is `(h^2-h')x.grad u`. Define

`A=a^m`, `B_i=A h x_i`, `C_ij=a^(m-2)delta_ij-A h^2 x_i x_j`.

Then

`partial_t(A u_t-B.grad u)-div(B u_t+C grad u)+A V u=0`.

The eigenvalues of `C` are `a^(m-2)` tangentially and `a^(m-2)-a^m h^2|x|^2` radially. The latter is at least `b^2 a^(m-4)>0`. All these bounds are uniform on a compact finite slab; they are not uniform as the slab tends to infinite Minkowski time.

Write `q=u|S^(m-1)`. At the boundary,

`n psi=(t/b)q_t+(b/a^2)u_r`,

`(B u_t+C grad u).nu = b a^(m-2) n psi`.

The induced boundary metric is `-(b/a)^2 dt^2+a^2 sigma`. Put

`W=b a^(m-2)`, `K=a^m/b`, `L=b a^(m-4)`.

The even boundary equation becomes

`W n psi + c partial_t(K q_t) - c L Delta_S q + d W q=0`.

Consequently the principal bulk flux cancels against the surface kinetic and tangential-gradient energy. For a constant `Lambda>0`, define

`E=1/2 int_B [A u_t^2+C_ij u_i u_j+Lambda A u^2]`

`  +1/2 int_S [c K q_t^2+c L |grad_S q|^2+d W q^2]`.

In the odd sector delete the two terms proportional to `c`. Negative `d` does not invalidate positivity: the trace interpolation inequality `||q||_L2(S)^2 <= epsilon ||grad u||_L2(B)^2+C_epsilon ||u||_L2(B)^2` and sufficiently large `Lambda` give a positive energy equivalent to the natural norm on a fixed compact slab. This argument covers arbitrary fixed real `d`. The bulk potential need only be bounded for this estimate; no derivative of `V` is needed.

For smooth solutions the exact differentiated identity is

`E' = int_B [-A'/2 u_t^2+B'_i u_i u_t+C'_ij/2 u_i u_j`

`                  +Lambda A'/2 u^2+A(Lambda-V)u u_t]`

`       +int_S [-c K'/2 q_t^2+c L'/2 |grad_S q|^2+d W'/2 q^2]`.

Thus `|E'| <= C E`, and `E(t)<=E(t0)exp(C|t-t0|)` for homogeneous evolution. This is a finite-time coercive estimate, not positivity of the unshifted physical Hamiltonian or a claim of spectral stability.

## 2. Hemisphere transformation and its actual domain

For spatial dimension four, use

`T=b sin(tau)/(cos(tau)+cos(chi))`,

`R=b sin(chi)/(cos(tau)+cos(chi))`, `Q=b/(cos(tau)+cos(chi))`,

and `psi=Q^(-3/2)U`, with `0<=chi<=pi/2`. The bulk equation is exactly

`U_tautau-Delta_S4 U+[9/4+Q^2 V]U=0`.

At the equator, `Q_b=b sec(tau)`,

`mu=c/Q_b=c cos(tau)/b`,

`beta=(d-3/(2b))Q_b-mu[(3/2)sec(tau)^2+(3/4)tan(tau)^2]`,

and the even condition is

`partial_chi U+partial_tau(mu U_tau)-mu Delta_S3 U+beta U=0`.

For the odd condition replace `mu` by zero and use `beta_odd=(d-3/(2b))Q_b`.

The exact invariant is

`s=R^2-T^2=b^2[cos(tau)-cos(chi)]/[cos(tau)+cos(chi)]`.

On `|tau|<=tau_*<pi/2`, its range is `[-b^2 tan(tau_*/2)^2,b^2]`. Hence regularity of the actual background on `-s_-<s<=b^2` permits this rectangle when `tau_*<2 arctan(sqrt(s_-)/b)`. Bounding its boundary Minkowski time by `sqrt(s_-)` would be unnecessarily restrictive: the shell always has `s=b^2`.

Nevertheless, geometric coverage must be stated accurately. A hemisphere rectangle covers every original Minkowski slice only up to `|T|<=T0=b tan(tau_*/2)`. Its boundary reaches `|T|=b tan(tau_*)`. A theorem on this rectangle does not by itself prove existence on the entire original regular Minkowski slab. The local contact construction below requires only the covered slab.

## 3. A complete finite-harmonic construction without a broad regularity claim

Set `x=pi/2-chi`, let `Y_l` satisfy `-Delta_S3 Y_l=lambda_l Y_l`, `lambda_l=l(l+2)`, and write `U=f(tau,x)Y_l`. The rescaled function `v=cos(x)^(3/2)f` satisfies

`v_tautau-v_xx+W_l(tau,x)v=0`,

`W_l=(lambda_l+3/4)sec(x)^2+Q^2 V`.

At `x=0`, with `z(tau)=v(tau,0)` and `k=mu lambda_l+beta`, the even condition is

`mu z''+mu' z'-v_x+kz=0`.

Choose smooth initial position `F(x)` supported in a collar `[a,d_0]`, with `a>0`, `inf supp F=a`, and zero initial velocity `G=0`. Choose `d_0+tau_*<pi/2`. Then all corner compatibility conditions hold because the data vanish on a neighborhood of the boundary. Smoothly extend the coefficient past an interior cutoff greater than `d_0+tau_*`; ordinary finite propagation makes this extension irrelevant, and no axis regularity issue arises.

Here is an explicit causal construction. The incoming Riemann invariant `r=v_tau+v_x` obeys

`r(t,x)=F'(x+t)+G(x+t)-int_0^t W_l(s,x+t-s)v(s,x+t-s) ds`.

Since `v_x(t,0)=r(t,0)-z'(t)`, the boundary unknown solves

`mu z''+(mu'+1)z'+kz`

` =F'(t)+G(t)-int_0^t W_l(s,t-s)v(s,t-s) ds`.

Recover the bulk field from

`v(t,x)=F(x-t)+int_0^t r(s,x-t+s)ds`, if `x>=t`,

`v(t,x)=z(t-x)+int_(t-x)^t r(s,x-t+s)ds`, if `x<t`.

These equations form a Volterra system. Put `p=z'` and `z=z(0)+int p`; on a short interval the bulk difference is bounded by `t||delta p||+C t^2||delta v||`, while the ODE difference is bounded by `C t||delta p||+C t^2||delta v||`. This is a contraction for sufficiently short time in the supremum norm for `(v,p)`. On a compact interval `mu` is bounded below by a positive constant, and the linear Volterra estimates give continuation and uniqueness by Gronwall. Differentiating the integral equations gives smoothness successively because every coefficient and initial datum is smooth. The initial-boundary corner is identically zero on a neighborhood, so differentiating across the characteristic issuing from the corner causes no compatibility jump.

The odd condition is `-v_x+beta_odd z=0`; the same construction uses the first-order ODE `z'+beta_odd z=r(t,0)`. It is at least as regular.

A finite sum over `l` is therefore a smooth genuine solution. Its support remains inside `x<=d_0+|tau|`, so extension by zero near the axis is smooth on the full hemisphere. No uniform-in-`l` estimates or convergence of an arbitrary infinite spherical-harmonic series have been supplied by this construction; those are not needed for existence of this finite-harmonic family.

## 4. Finite propagation and first boundary activation

The causal integral system proves `z(t)=0` for `0<=t<a`: the initial forcing vanishes on the entire backward cone of these boundary points, and the homogeneous Volterra system has only the zero solution there. The same system gives the ordinary bulk support bound `x<=d_0+t`. A positive dynamic boundary mass does not create a signal before the incoming light cone arrives.

All coefficients in the second-order problem are even in `tau` in the appropriate sense: `W_l`, `mu`, `beta` are even and `mu'` is odd. With zero initial velocity, uniqueness therefore makes the solution even in time.

Suppose the trace were zero throughout `(0,tau0)` for some `a<tau0<=tau_*`. Evenness gives `v(tau,0)=0` on `|tau|<tau0`, and the boundary equation gives `v_x(tau,0)=0` there. Regard

`v_xx-v_tautau=W_l v`

as a one-dimensional wave equation with `x` as its evolution variable. Standard cone energy uniqueness with bounded smooth lower-order coefficient gives `v=0` in `x>=0`, `|tau|+x<tau0`. Evaluating at `tau=0` yields `F(x)=0` for `0<x<tau0`, contradicting `inf supp F=a`.

This sideways argument is legitimate here because after harmonic separation there is exactly one tangential time variable: it is an ordinary `1+1` hyperbolic Cauchy problem. It is not an appeal to general well-posedness of a higher-dimensional timelike Cauchy problem and needs no analyticity.


Thus the trace is nonzero at some point in every interval `(a,a+epsilon)` contained in the construction interval. A smooth compactly supported datum is flat at its support edge, so the trace is still zero at the exact first arrival `tau=a`; no nonzero value at that exact time is claimed. In physical shell time the onset has infimum `T=b tan(a)`. If the claim also requires onset before the fully covered ambient slice endpoint `T0=b tan(tau_*/2)`, impose `a<tau_*/2`; the proposed stronger `d_0<tau_*/2` ensures this.

## 5. What a broader Galerkin theorem would still require

For the even hemisphere problem the correct fixed spaces are

`H=L2(hemisphere) direct_sum L2(equator)`,

`V={(u,q):u in H1(hemisphere), q in H1(equator), trace(u)=q}`.

The mass is `diag(1,mu(tau))`. A boundary velocity in `H` is an independent component at the finite-energy level: the bulk `L2` velocity has no ordinary trace. A proof must not silently identify it with a pointwise boundary trace.

Uniform positive mass, smooth bounded coefficients, and coercivity after a scalar shift give the Galerkin energy bounds and a weak limit. To assert uniqueness and strong continuity one still supplies the time-regularization or duality energy argument. To assert smoothness one supplies, in addition, an elliptic Ventcel estimate and a compatibility hierarchy. A representative common operator domain is the trace-compatible `H2(bulk) direct_sum H2(boundary)` space for a smooth boundary and uniformly positive surface elliptic coefficient, with its graph norm justified by elliptic regularity.

Successive initial derivatives are recursively determined by the bulk and surface equations. Their traces must agree to every order for a smooth solution. Arbitrary Galerkin projections of smooth initial data need not produce uniformly bounded projected initial accelerations or higher derivatives; differentiating the finite-dimensional ODE and citing the zeroth-order energy bound does not fix that defect. One may instead use graph-norm-compatible approximants or a proved common-domain nonautonomous evolution theorem. The finite-harmonic Volterra construction above avoids this separate theorem rather than claiming it for free.

The original-time fixed-ball weak form has mixed space-time terms. Its formal cancellation is sound, but converting it into a general abstract evolution theorem requires care, especially in the odd sector where there is no surface `L2` velocity. A formal multiplier identity alone does not prove existence.

## Verification performed

Independent symbolic simplification with SymPy passed for the fixed-ball radial lower bound and flux coefficient, the radial Liouville potential `(lambda_l+3/4)sec(x)^2`, and the induced boundary conformal coefficient `(3/2)sec(tau)^2+(3/4)tan(tau)^2`. The Volterra construction and activation argument were checked analytically. No external network or publication operation was used.
