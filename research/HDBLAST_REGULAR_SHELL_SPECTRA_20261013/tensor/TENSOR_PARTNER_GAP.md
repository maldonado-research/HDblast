# Tensor partner operator and an analytic mass threshold

Working derivation, 13 October 2026 UTC. This is a conditional result for the
registered static backgrounds, not a statement about the origin of the Big
Bang or full nonlinear stability. Factorization and de Sitter mass thresholds
are established methods; the contribution assessed here is their exact
application to the already certified nonconstant branch.

## Assumptions and physical reduction

Use the Einstein–scalar action and regular doubled-interior static solution
of the October 11 capsule. Write p=eta', a=rho'/rho, and

    a' = -rho^(-2) - p^2/3,
    a^2 = rho^(-2) + p^2/12 - U/6.

At fixed positive detuning, 0<y<=T is finite, rho~y at the regular cone,
p=O(y), a~1/y, and a(T)>0. The standard transverse-traceless reduction is

    h'' + 4 a h' + mu_T^2 rho^(-2) h = 0,
    m_T^2 = mu_T^2 H_b^2.

The September 27 full linearized tensor calculation already verifies this
equation and the one-copy Neumann junction. For two identical glued interiors,
continuity of the induced transverse-traceless metric and the sum junction
give h_+=h_- and h_+'+h_-'=0. Thus the even sector has h_e'(T)=0 and the odd
sector has h_o(T)=0. These conditions concern the pure tensor reduction;
scalar-derived exceptional harmonics and their gauge overlap must still be
classified in the accompanying physical audit.

## Closed first-order operators

Choose dz=dy/rho with z(T)=0, so z ranges over (-infinity,0]. Set

    psi = rho^(3/2) h,    S=(3/2) d(log rho)/dz=(3/2)rho',
    D=d/dz-S,            D*= -d/dz-S.

The radial tensor norm is int rho^2 |h|^2 dy = int |psi|^2 dz. For each
fixed background S and its conformal derivative are bounded; their cone
limits are S->3/2 and dS/dz->0. The operators below therefore have the
ordinary half-line Sobolev domains, with no independently selectable boundary
condition at the limit-point cone.

Even sector: Dom D=H^1(-infinity,0], with free shell trace. Integration by
parts gives Dom D*=H^1 with zero shell trace. Thus D*D imposes Dpsi(0)=0,
exactly h_e'(T)=0, and DD* imposes Dirichlet on its partner chi.

Odd sector: Dom D=H^1 with zero shell trace, and Dom D*=H^1 with free trace.
Then D*D is Dirichlet, while DD* imposes

    chi'(0)+S_b chi(0)=0.

Its quadratic form has the **positive** endpoint term S_b |chi(0)|^2,
because S_b=(3/2)rho(T)a(T)>0. Reversing the orientation at this step would
invalidate the spectral conclusion.

## Partner bound and spectrum

Both tensor operators are H_+=D*D=-d_z^2+V_+, where V_+=S^2+dS/dz.
The partner is H_-=DD*=-d_z^2+V_-, with

    V_- = S^2-dS/dz
        = rho^2(3a^2/4-3a'/2)
        = 9/4 + rho^2(9p^2/16-U/8).

If U<=0 throughout the interior, V_->=9/4. The even partner's Dirichlet
form and the odd partner's positive Robin form therefore both obey

    H_- >= (9/4) I.

The polar decomposition of a densely defined closed D makes D*D and DD*
unitarily equivalent on their nonzero spectral subspaces. This compares
the full nonzero spectrum, not just normalizable eigenfunctions. It follows
that neither tensor parity has spectrum in (0,9/4) or below zero.

The even kernel is exactly psi=C rho^(3/2), the normalizable massless
graviton. Its norm is finite because T is finite and rho~y at the cone.
The odd Dirichlet condition removes this kernel. The formal partner kernel
rho^(-3/2) is not square integrable at the cone.

The smooth regular-cone expansion gives V_-=9/4+O(y^2), with
y=O(exp(z)) as z->-infinity. Its decaying bounded potential perturbation is
relatively compact with respect to the half-line Laplacian; a change of a
regular finite-endpoint condition is also a finite-rank resolvent change.
Consequently the essential spectrum begins at 9/4 in both parities. The
pointwise bound rules out subthreshold massive tensor states. It does not
assert the absence of embedded positive eigenvalues without the additional
short-range spectral argument, which is unnecessary for the lower bound.

Thus the dimensionless mass spectrum is {0} plus values >=9/4 in the even
tensor sector, and values >=9/4 in the odd tensor sector. In physical units,
all nonzero tensor masses obey m_T>=3H_b/2. This lies above the usual massive
spin-2 Higuchi value m_T^2=2H_b^2; it only rules out that mass interval for
this reduced tensor problem. A claim of a complete ghost-free theory would
require the full physical action and exceptional-mode audit.

## Exact registered-model gate

The already certified branch stays inside |eta|<=1/100. Its potential is

    U = -2 eta^6/27 - 4 eta^5/9 - eta^4/6
        + 50 eta^3/27 + 14 eta^2/9 - 2/27.

The exact coefficient bound

    U(eta) <= U(0) + sum_{j>=1} |U_j| (1/100)^j < 0

holds on this entire polynomial interval, independently of delta. The
companion `verify_tensor_partner.py` reconstructs this polynomial from W,
evaluates the upper bound rationally, verifies the tensor transformation and
partner formula, and rejects four deliberately incorrect formulas. Normal
and optimized Python outputs agree. The result therefore applies to every
0<delta<=0.002 covered by the prior static theorem, without certifying any
historical floating-point trajectory or making a uniform positive physical
mass assertion as delta->0 (H_b itself tends to zero).

The old September 27 files are preserved: their exclusion of tensor bound
states in (0,9/4) was numerical. The present analytic implication depends
on the exact branch construction, U<=0, the stated endpoint domains and the
pure-tensor junction conditions. It establishes no reheating mechanism,
cosmological observables, vector-sector theorem, nonlinear attraction or
vacuum-decay suppression.
