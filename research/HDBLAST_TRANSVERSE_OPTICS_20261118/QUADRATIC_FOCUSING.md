# Quadratic bulk volume focusing along the registered tangent rays

18 November 2026. This derivation concerns the one-sided smooth bulk of either registered Einstein–scalar copy, at a shell-tangent ray. It uses the actual local family and its regularity, rather than an independently chosen second metric variation. The numerical conclusions are consequences of the retained analytic and waveform enclosures. No numerical finite-amplitude remainder is evaluated.

## Result and normalization

Use the same vertex and physical affine normalization as the archived shell optical calculation: initial sphere label `n_*=(0,0,1,0)`, vertex clock `q0=ell*(49/32)`, terminal reference clock `q1=ell*(99/64)`, and

    ell=10^-6,
    L=rho_b [sinh(q1)-sinh(q0)]/cosh(q0),
    |h_e(u_e,k_e)|=1 at the vertex.

Both `L` and the normalization are held fixed as the canonical amplitude `e` varies. Let `J_e` be the **three-dimensional bulk-screen** Jacobi matrix along this shell-tangent null ray, with `J_e(0)=0`, `J_e'(0)=I_3`. Write

    D3_e=(det J_e)^(1/3),
    V3_e=det J_e,
    Q(v)=v1^2-v2^2,
    chi(v)=lim_(e->0) [1-D3_e(v,L)/L]/e^2.

The central result is

    chi(v) = (1/(3L)) integral_0^L lambda(L-lambda)
                   { p1(v,lambda)^2 + 6 f_v'(lambda)^2 } d lambda,       (1)

where `p1=partial_e(k_e Phi_e)|0` and `f_v=delta D_shell(v,lambda)/lambda`. The prime on `f_v` means physical affine differentiation; `delta` means the amplitude derivative. In particular `chi>=0`. For the registered profile, the exact rational lower bounds proved below are

    chi(v) > 2.2755*10^-51

at every direction in either tangent cap `v1>=4/5` or `v2>=4/5`. On either axis `Q=+1` or `Q=-1`, `chi>2.9024*10^-50`. The average of `chi` over either tangent cap exceeds `1.5917*10^-50`.

These are coefficients per squared canonical amplitude. The actual expansions are

    D3_e/L = 1-e^2 chi(v)+o(e^2),
    V3_e/L^3 = 1-3e^2 chi(v)+o(e^2).                           (2)

Thus the normalized optical volume deficit has coefficient `3chi`, not `chi`; its amplitude second derivative is `-6chi`. Both signs of sufficiently small nonzero amplitude give strict bulk volume focusing along the indicated tangent pencils. An amplitude neighborhood for the numerical lower bound is not evaluated.

A tangent cap here is a cap of the shell observer's `S^2` tangent sky. It is not a positive-measure `S^3` bulk sky aperture. Each pencil has a three-dimensional screen, including the one-sided normal Jacobi direction. A claim about a finite transverse aperture additionally requires the corresponding inward-ray retention and variation construction. None is assumed in the numerical replay.

## Geometric input and one-sided interpretation

The archived local physical solution satisfies on each smooth bulk copy

    Ric5 = dphi tensor dphi + (2/3) U(phi) g,
    K_ab = c h_ab,  c=sigma(Phi)/6.

The shell scalar is `Phi=phi|shell`. These are the actual pure-tension Israel and bulk equations, with variable scalar-dependent `c`; its constancy is unnecessary. They are recorded in `HDBLAST_SUBSIDIARY_CLOSURE_20261116/PHYSICAL_RECONSTRUCTION_AND_BOUNDARY_JOIN.md`, and in the retained `PR145_UNREDUCED_SOLUTION.md` copied by the effective-observer capsule. For every tangent null vector,

    Ric5(k,k)=(k Phi)^2.                                      (3)

The shell is exactly umbilic, so a shell null geodesic with affine tangent `k` obeys `nabla5_k k=nabla4_k k+/-K(k,k)n=0`. It is consequently also a bulk geodesic along the boundary. This statement uses only the smooth metric and its boundary jets on one copy; it does not posit a bulk crossing law at the doubled distributional shell.

Choose two shell screen vectors `E_1,E_2`, parallel along the ray modulo `k`, and use the unit normal `n` as the third screen direction. Their bulk normal derivative contribution is `K(k,E_i)=0`; moreover `nabla_k n=c k` in the stated `K` convention, so the normal class is screen-parallel. A multiple of `k` can be added to obtain ordinary parallel representatives without changing any screen curvature contraction.

Define the optical tidal matrix by the sign convention in which the Jacobi equation is

    J''=-T J,    tr T=Ric5(k,k).

Gauss eliminates the extrinsic-curvature products in the tangent screen block because `h(k,k)=h(k,E_i)=0`. Codazzi eliminates its mixed normal/tangent entries: the only derivatives of `c` multiply those same vanishing contractions. Hence, exactly at every amplitude,

    T5 = diag( T4, (k Phi)^2-Ric4(k,k) ),                     (4)

where `T4` is the shell's two-dimensional optical tidal matrix. No electric-Weyl term is omitted: it is encoded by the difference in the last entry. This derivation also shows that the result is independent of reversing the normal on the second identical copy.

The one-sided `J_e` is the unique screen Jacobi solution using these actual boundary curvature jets. It is independent of a smooth auxiliary metric extension across the boundary. Where obtained as the derivative of an inward family of bulk rays, it is precisely that family's normal-angle derivative. We do not infer a two-sided family of physical bulk rays through a distributional shell from this ODE alone.

## First strain and the missing linear bulk trace

At zero amplitude the shell metric is de Sitter,

    h0=rho_b^2[-dq^2+cosh(q)^2 kappa].

Its null tidal matrix vanishes. The background `Phi0` is constant along the shell, so `k0 Phi0=0`. Equation (4) therefore gives `T5,0=0`, `J0=lambda I_3`, and `D3,0=lambda`.

The actual shell metric tangent is `2Fh0` plus a tangential Lie derivative, with `F=U(q)Y(n)`, `Y=n1^2-n2^2`. The first intrinsic Weyl tensor vanishes. Since the constant-curvature background tidal matrix is zero, evaluating-event, ray and frame variations add no first tidal forcing. The first tangent optical matrix is isotropic, and the physically normalized vertex Jacobi data are fixed. Thus

    delta J4 = d1(v,lambda) I_2,
    d1(v,lambda)=delta D_shell(v,lambda).

By (3), `delta Ric5(k,k)=0`. Taking the trace of the first bulk Jacobi equation yields `tr delta J5=0`; the zero initial data fix the integration constants. Consequently

    delta J5 = diag(d1,d1,-2d1).                              (5)

The third direction compensates the shell's linear area change. It cannot be dropped when using the five-dimensional null-energy condition.

Let `B=J' J^-1`, `theta=tr B`, and `Sigma=B-(theta/3)I_3`. Put `f_v=d1/lambda`, extending continuously to the vertex. Differentiating at zero gives

    delta B = diag(f_v',f_v',-2f_v'),
    delta theta=0,
    Sigma1=diag(f_v',f_v',-2f_v'),
    tr(Sigma1^2)=6(f_v')^2.                                  (6)

Here the shear norm is the full positive Frobenius square `tr(Sigma^2)`. This explicit convention avoids the factor-of-two convention sometimes used for two-dimensional optical shear.

The scalar tangent supplies an additional nonnegative term. The archived field tangent on the shell is `delta Phi=S(-q)Y`. Because the background restriction `Phi0` is constant, its pullback and ray/frame variations vanish. Along the background ray,

    p1 = d/dlambda [ S(-q(lambda)) Q(v)
                        lambda^2/(rho_b^2 cosh(q(lambda))^2) ].           (7)

The bound below discards this term. It requires neither an enclosure of `S'` nor a numerical differentiation of waveform cells.

## Focusing equation and the quadratic limit

The one-sided Jacobi equation itself proves zero twist: its symmetric tidal matrix conserves `W=J^T J'-J'^T J`, and the vertex data give `W=0`. While `J` is invertible, multiplication by its inverse and transpose shows that `B=J'J^-1` is symmetric. This argument does not assume an open cone of physical rays across the shell. The three-screen optical Riccati equation gives

    theta'=-theta^2/3-tr(Sigma^2)-Ric5(k,k),
    D3''=-(1/3){tr(Sigma^2)+(k Phi)^2}D3.                     (8)

The actual family is `C2` in the retained lower bulk/trace spaces; the physical metric and scalar have more joint differentiability than the curvature, geodesic, parallel-frame and Jacobi ODEs used here require. See the archived `HDBLAST_EVOLVED_PARAMETER_20261116/README.md` and `EVOLVED_PARAMETER_AND_PHYSICAL_CHAIN.md`. In fact the limit in (1) only needs continuous first parameter derivatives of those ODE data. It does not require estimating the unknown second metric tangent.

For completeness, put `J=lambda H`. Its regular Volterra equation gives `H=I+O(lambda^2)` and the same estimate for first parameter differences on this compact segment. The regular part `B-I/lambda` and its first parameter derivative are `O(lambda)`. Thus

    Sigma_e=e Sigma1+o(e),
    k_e Phi_e=e p1+o(e),
    D3_e=lambda+o(1)

hold with integrable uniform vertex control. The archived common interior clock margin and continuous dependence keep the tangent rays in the physical shell neighborhood. Because the background matrix is `lambda I`, they remain free of conjugate points for sufficiently small amplitudes. Twice integrating (8), dividing by `e^2`, and taking the limit yields (1). This also proves (2). The qualitative remainder is enough for strict focusing at sufficiently small nonzero amplitude, but supplies no numerical amplitude threshold.

Equation (8) separately gives the exact ceiling `D3_e(lambda)<=lambda` on any actual positive, caustic-free segment, from the canonical scalar null-energy condition. The new content here is a strictly positive, evaluated quadratic coefficient from the actual registered waveform.

## A lower bound from the affine profile, without its derivatives

Define the constants and background ray geometry

    K=ell/(mu rho_b^4)=ell B_b/rho_b^3,
    Omega=K L^2,
    theta=atan(sinh q)-atan(sinh q0),
    sin(theta)=lambda/(rho_b cosh q).

Along the great circle from `n_*`, `Y=Q(v)sin(theta)^2`. The archived fixed-affine shell Jacobi formula is

    d1(v,l)=l[F(0)+F(l)]-2 integral_0^l F(lambda) d lambda.

Here `F(0)=0`. With `H(x)=mu rho_b^2 U(ell x) sech(ell x)^2/ell`, it follows, for **every prefix** `0<l<=L`, that

    f_v(l)=Q(v) K l^2 C(l),
    C(l)=H(x(l))-2 integral_0^1 t^2 H(x(lt)) dt.              (9)

Every prefix stays within the archived second window. If its uniform bounds are `H_-<=H<=H_+`, positive integration gives

    H_--(2/3)H_+ <= C(l) <= H_+-(2/3)H_-.

The authenticated `METRIC.json` endpoints imply the simpler strict bounds

    c_-=0.051613 < C(l) < 0.107100=c_+.                       (10)

The use of all prefixes matters: the value of `f_v(L)` alone cannot give this weighted shear lower bound, because the reciprocal endpoint weight is not integrable on `[0,L]`.

Use the strictly interior affine subinterval `[L/8,7L/8]`. For `Q!=0`, (9)–(10) imply

    |f_v(7L/8)-f_v(L/8)| > |Q| Omega d,
    d=(49c_--c_+)/64=2421937/64000000 >0.03784.                (11)

Weighted Cauchy–Schwarz, with weight `lambda(L-lambda)`, gives

    integral_(L/8)^(7L/8) lambda(L-lambda)(f_v')^2 d lambda
       >= L [f_v(7L/8)-f_v(L/8)]^2/log(49).

Discarding the scalar term and the remaining nonnegative shear integral in (1),

    chi(v) > 2 Q(v)^2 Omega^2 d^2/log(49)
           > (1/2) Q(v)^2 Omega^2 d^2,              Q!=0.    (12)

The final inequality follows from `log(49)<4`; the replay verifies it by the positive Taylor lower sum `sum_(k=0)^7 4^k/k!>49` for `exp(4)`. No transcendental floating-point comparison is used in the certificate. At `Q=0`, (1) is still valid and gives `chi=0`, because both (6) and (7) are proportional to `Q`. More generally the exact directional dependence is `chi(v)=Q(v)^2 chi_axis`.

For the physical scale, monotonicity of `cosh q` on this interval gives

    L/rho_b >= ell*(99/64-49/32)=ell/64,
    Omega >= ell^3*(1/64)^2*(B_b)_lower/(rho_b)_upper
           >6.3666728160*10^-24.                            (13)

All shell coefficients in (13) are the archived exact rational enclosures. Substitution in (12) proves the axis bound stated above. On `v1>=4/5`,

    Q>=v1^2-(1-v1^2)>=7/25;

on `v2>=4/5` the corresponding negative bound holds. This proves the uniform pointwise cap bound. For completeness, normalized `S^2` cap averaging makes `x=v1` uniform on `[4/5,1]`, with

    E[v2^2|x]=(1-x^2)/2,
    E[v2^4|x]=3(1-x^2)^2/8,
    E[Q^2]=E[(19/8)x^4-(7/4)x^2+3/8]=10283/18750.

Thus either tangent-cap average of `chi` obeys the stated stronger lower bound. This is an average over tangent directions, not an assertion that the corresponding set has positive measure in the full bulk sky.

## Reproduction and scope

Run `python quadratic_focusing.py` in this directory. The script authenticates the retained sibling `METRIC.json` and `inputs/SHELL_BOX.json`, recomputes every rational inequality in (10)–(13), checks the exact cap second moment, and prints the `QUADRATIC.json` report. Its computations certify finite arithmetic consequences. The continuum Gauss/Codazzi, screen transport, regularity and ODE arguments are the written proof above.

The source normalization is the half-action Einstein–scalar convention in which (3) has coefficient one. A different scalar normalization changes the scalar square, while the shear-only bound remains tied to the registered shell metric and canonical amplitude. The construction uses no higher-dimensional energy-condition violation. It distinguishes intrinsic two-screen defocusing from the compensating third screen and the positive quadratic bulk volume focusing. It claims neither an experimental detection nor a nonlinear global or cosmological conclusion.
