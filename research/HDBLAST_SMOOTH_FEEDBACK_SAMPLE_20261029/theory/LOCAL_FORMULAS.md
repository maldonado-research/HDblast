# Local evaluation formulas for the consistent smooth candidate feedback

This note concerns the **first feedback of a nominal candidate** in the frozen
constant-matrix 23-mode system. It does not solve the actual second-order field
equations or give a numerical C² error. In particular affine bulk/shell sources
remain omitted by this candidate and their errors remain part of the inherited
conditional error ledger.

Let the candidate use one smoothed displacement u(z), its **own derivative**
u'(z), and the scaled x-velocity v(z). Let A,B,D be the fixed full23 matrices
in the declared Gram normalization. Matrix multiplication occurs in the displayed
order; no commutation, independent parity-equivariance, or diagonal hypothesis
is used. Work at 0<x<r, a=r+x, q=r−x>0; then the characteristic diamond lies on
the physical half-line and the boundary lift is zero. The formulas do not cover
post-contact or shell points without the reflected terms.

Define normalized primitives V_j of v with V_0=v and dV_j/dz=V_{j−1}, and

    F(z) = [(A+B)(u'(z)+v(z))+D(u(z)+V_1(z))]/2,
    G(z) = [(B−A)(u'(z)−v(z))+D(u(z)−V_1(z))]/2.

Then A B0_x+B B0_r+D B0 = F(r+x)+G(r−x) for the free wave B0. Change
variables in the unreflected Duhamel diamond to get the exact local expression

    B1 = (1/4) integral_q^a [(z−q)F(z)+(a−z)G(z)] dz.

Differentiating its endpoint weights at fixed r or x gives

    B1_x = x[F(a)+G(q)]/2 + integral_q^a(F+G)/4,
    B1_r = x[F(a)−G(q)]/2 + integral_q^a(G−F)/4.

No global primitive subtraction is needed. The following integration-by-parts
form makes interval evaluation substantially tighter because only one local
integral of u is required. Write u_a=u(a), u_q=u(q), p_a=u'(a), p_q=u'(q),
v_a=v(a), v_q=v(q), and

    Iu = integral_q^a u(z) dz,
    Iv = integral_q^a v(z) dz,
    Jv = integral_q^a (z−r)v(z) dz,
    Kv = integral_q^a [x²−(z−r)²]v(z) dz.

The three feedback values are exactly

    B1 = { x A(u_a+u_q)+x B(u_a−u_q)+(xD−A)Iu
                 +x A Iv+B Jv+(D/2)Kv }/4,

    B1_x = { x[(A+B)p_a+(B−A)p_q+D(u_a+u_q)
                         +(A+B)v_a+(A−B)v_q+D Iv]
                 +B(u_a−u_q)+A Iv+D Iu }/4,

    B1_r = { x[(A+B)p_a−(B−A)p_q+D(u_a−u_q)
                         +(A+B)v_a+(B−A)v_q]
                 −A(u_a−u_q)−B Iv+D Jv }/4.

For example, F+G=B u'+A v+D u, G−F=−A u'−B v−D V_1,
and the remaining v primitive cancels using

    integral_q^a (z−r)V_1(z) dz = Kv/2,
    x[V_1(a)+V_1(q)]−integral_q^a V_1(z)dz = Jv.

The velocity moments can be evaluated exactly from its smoothed P1 hinges:

    Iv = V_1(a)−V_1(q),
    Jv = x[V_1(a)+V_1(q)]−[V_2(a)−V_2(q)],
    Kv = 2x[V_2(a)+V_2(q)]−2[V_3(a)−V_3(q)].

Those rational expressions may contain cancellation, but are exact arithmetic
on a finite polynomial representation. They do not introduce interval width.
Only Iu and the endpoint u,u' quantities need enclosures. One may bound
u_s(z) via authenticated nominal u0(z) plus the global smoothing remainder
alpha_u=M2 delta_u²/18. Integrating a nominal interval enclosure on each local
panel and adding (a−q)alpha_u encloses Iu without nested convolution.
Endpoint derivatives similarly use beta_u=35 M2 delta_u/128, or a separately
validated convolution quadrature. This inherits the prior global nominal
second-derivative bound and smooth-construction premises; it does not prove them.

All quantities are formal coefficient amplitudes in scaled coordinates. A
computed nonzero nominal candidate value, if later transferred through a valid
conditional error bound, would still concern only that formal coefficient and
would not establish a physical force, cosmological event, or nonlinear family.
