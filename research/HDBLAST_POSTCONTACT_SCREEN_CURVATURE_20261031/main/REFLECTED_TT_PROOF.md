# Frozen first experiment after shell contact

The point is P=(5/4,1/8), with shell history ending at 9/8. Its physical
backward set is C={0<=s<=5/4,max(0,s-9/8)<=v<=11/8-s}. This is backward
closed under direct and Neumann reflected characteristics, and is contained
in the existing physical comparison domain D+. Its first-field footprint is
within the authenticated domain 0<=s<=3,v>=0,s+v<=8.

The homogeneous-Neumann Green operator acts on the complete physical-side
source f=A u_x+B u_r+D u+F, extended by f(s,|v|). In particular reflecting
u_r alone while retaining an even B would change the equation. Endpoint
differentiation gives

    (K_N f)_x=1/2 integral_0^x [f(s,r+x-s)+f(s,|r-x+s|)] ds,
    (K_N f)_r=1/2 integral_0^x [f(s,r+x-s)-f(s,|r-x+s|)] ds.

There is no extra sign on the lower reflected endpoint. These formulas give
(K_N f)_r(x,0)=0. For the boundary convention u_r=G,

    L G=-1_(x>r) integral_0^(x-r) G(s)ds,
    (L G)_x=-G(x-r), (L G)_r=+G(x-r) for x>r.

The initial data are zero, and the retained compatible-corner premises
include G(0)=0; thus the lift is C1 across x=r. No new corner regularity
theorem is claimed here. The physical positive-y direction is opposite +r;
the registered boundary_forcing adapter already includes that conversion.

For T=5/4 put X(u)=max(2 sup|u|/T,sup|u_x|,sup|u_r|). Direct absolute
integration gives X(K_N f)<=T sup|f|: reflection changes multiplicity but not
the total integration length 2(x-s). Also X(LG)<=2 sup|G|. With a,b,d the
absolute coefficient caps, sigma=a+b+Td/2 and q=T sigma, the actual equation
U=LG+K_N F+K_N La U yields

    X(U)<=(2g+Tf)/(1-q),
    |U_x(P)-(LG+K_N F)_x(P)|<=q(2g+Tf)/(1-q),  q<1.

This is the complete feedback tail, not one omitted iterate. It uses only
C1 U and continuous sources. It is sufficient for the selected screen map,
which uses only U_x. The seed V=LG+K_N F is introduced only for this actual-U
enclosure; it is not the registered finite W1=LG+K_N(La LG+F). No W2 or C2
error improvement is claimed.

The two characteristic integrals are enclosed using whole rectangular source
boxes containing every point of each segment, times the exact segment length.
The lower characteristic is split at s=9/8. These are interval Riemann bounds,
not point quadrature and not finite differences. Uniform source and shell caps
retain the entire reflected history. All source derivatives are ordinary
physical derivatives in the original source map, which already supplies the
scaled PDE forcing -2 ell^2 rho^2 QEred. The shell source conversion is the
original 2 ell rho QKab. No extra ell factors are inserted in integration.

The original fixed first-metric quadratic screen map is evaluated at this
same point with actual reflected first-field errors. The final interval adds
the whole feedback error before applying the positive factor
3/(2 ell rho^2 cosh(t)), and then adds 4Fquad/(rho cosh(t)). All observer,
radial-screen, connection, and epsilon corrections in that original map are
retained. The result is a conditional coordinate coefficient with a specified
slicing observer. It does not establish a nonlinear shell-contact solution or
a geometrically transported gauge-invariant observable.
