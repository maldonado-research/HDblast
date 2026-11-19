# The signed coupling work depends on the daughter state

On each of the same two actual first-passage backgrounds, finite-energy linear
daughter states can give either sign of the integrated coupling work. This is
a fixed-background state-dependence test. It makes no claim about a sum over
all daughter modes or about a radiation fraction.

The background seeds are +0.001 and −0.001. The brackets are the archived
first-passage times plus or minus 0.005, rather than new crossing conventions:

| Parent seed | Bracket |
|---:|---:|
| +0.001 | [0.483395411319240, 0.493395411319240] |
| −0.001 | [0.584134244832117, 0.594134244832117] |

We use the canonical Einstein-frame action and unit-lapse proper time. With
the archived potential variable u, V=3u and the homogeneous energy-density
normalization is 3h². The positive daughter kinetic weight is
Z=f(Theta_top)/f(Theta). Its constant normalization cancels from every ratio
below. Positive coupling work means energy supplied to the daughter through
the varying kinetic coupling, with work by the spatial metric accounted for
separately.

## Exact finite-dimensional optimization

For a fixed real angular harmonic with ell=2, expand the daughter in the
radial basis n=2,4,...,n_max. The original-field Hamiltonian is

$$
 E(t)=\frac12p^TM(t)^{-1}p+\frac12q^TK(t)q,
 \qquad \dot q=M^{-1}p,\quad \dot p=-Kq.
$$

M and K are positive. Differentiating only Z at fixed geometry gives

$$
 \dot E\big|_Z=\frac12q^T\dot K_Zq
             -\frac12\dot q^T\dot M_Z\dot q.
$$

Let x=(q,p)=U(t)y, and choose the initial map so that the bounce energy is
E(0)=y^Ty/2. Then the signed work over the bracket is

$$
 W_Z=\frac12 y^T A_Zy,\qquad
 A_Z=\int_{t_-}^{t_+}U^T
 \begin{pmatrix}\dot K_Z&0\\0&-M^{-1}\dot M_ZM^{-1}\end{pmatrix}
 U\,dt.
$$

Consequently the exact minimum and maximum of W_Z/E(0), at any fixed
nonzero initial energy, are the smallest and largest eigenvalues of the real
symmetric A_Z. All intermediate values are attained by real initial states.
For energy specified instead at the beginning of the bracket, replace the
identity in this Rayleigh quotient by the positive form
G_-=U(t_-)^T diag(K,M^{-1})U(t_-). The corresponding endpoints are the
generalized eigenvalues of (A_Z,G_-). The JSON supplies both conventions and
explicit bounce-energy-one states attaining the first pair of extrema.

The state is propagated continuously from the original t=0 slice. Only the
work integral starts at t_-. There is no state reset at the bracket.

| Seed | Initial and propagated band | Minimum W_Z/E(0) | Maximum W_Z/E(0) |
|---:|---:|---:|---:|
| +0.001 | n≤8 | −0.00515860921 | +0.00559455299 |
| +0.001 | n≤12 | −0.00535812906 | +0.00593608507 |
| +0.001 | n≤18 | −0.00574158143 | +0.00614097077 |
| +0.001 | n≤24 | −0.00593265030 | +0.00622520632 |
| −0.001 | n≤8 | −0.0228013528 | +0.0378496958 |
| −0.001 | n≤12 | −0.0282740139 | +0.0428531379 |
| −0.001 | n≤18 | −0.0329027033 | +0.0434722476 |
| −0.001 | n≤24 | −0.0352223159 | +0.0434783494 |

These are dimensionless ratios. For example +0.0062252 means +0.62252% of
the specified daughter initial energy, not of the parent energy. The sign
intervals are far wider than the refinement errors. A universal direction
of coupling exchange on these backgrounds is therefore excluded already
in these finite sectors. This does not determine the sign for any particular
chosen vacuum or occupied state.

## Controls and the meaning of radial convergence

Changing the full initial band changes the admissible states. Its extremal
range is not a convergence sequence for one fixed state set. We therefore
also hold the initial field and velocity support at n≤8 and propagate it in
bands n≤8,12,18,24,36. In that fixed initial subspace the optimized intervals
at propagated n≤24 are

| Seed | Minimum W_Z/E(0) | Maximum W_Z/E(0) |
|---:|---:|---:|
| +0.001 | −0.005153924892 | +0.005610030702 |
| −0.001 | −0.02329174122 | +0.03920378814 |

The comparison to n≤36 changes any of the four extrema by less than
3.52e−11. Initial support
means field and proper-time velocity support; canonical momentum contains
the inhomogeneous mass weight and need not have the same support.

The main geometry uses degree 64 and 160 projection nodes. The refined
geometry uses degree 96 and 256 nodes, tighter gravity tolerance and half
the maximum step. Daughter quadrature changes from 192 to 272 points, and
daughter ODE tolerance and maximum step are independently tightened. At
n≤18, all four extrema (bounce and bracket energy conventions) change by
less than 4.58e−14 across these controls. The initial energy-whitening error
is below 4.5e−15. A constant-Z control has exactly zero coupling work.

The independent integrated matrix ledger checks the entire state space:

$$
 A_Z+A_g=U(t_+)^TG(t_+)U(t_+)-U(t_-)^TG(t_-)U(t_-),
 \quad G=\operatorname{diag}(K,M^{-1}).
$$

The largest absolute matrix residual is below 3e−14. This checks both work
signs without selecting a favorable initial state. The derivative matrices
are computed explicitly from the spline slope, parent velocity, and metric
expansion rates; the code does not obtain coupling work by subtracting the
two endpoint energies.

Replay from the study root:

```sh
OPENBLAS_NUM_THREADS=1 python -IB checks/phase_work.py --output /tmp/PHASE_WORK.json
```

[Source](phase_work.py) · [Machine-readable values and checks](PHASE_WORK.json)

## Keeping the same maximizing initial wave

We kept the field and proper-time velocity of the maximizing state at radial cutoff 24 fixed and propagated those same initial data at cutoffs 36 and 48. We use the actual daughter energy at the beginning of the first-passage window as the denominator. At cutoff 48 the fixed wave transferred +0.00697038 of that energy on the positive-seed background and +0.04111672 on the negative-seed background. Increasing the propagation cutoff from 36 to 48 changed these ratios by less than $1.7\times10^{-11}$. Reoptimizing over the full original initial support through label 24 at propagated cutoff 48 gives maxima +0.00697041 and +0.04112516. The fixed-wave examples reach 97.65% and 97.41% of the independently evaluated continuum energy bounds. These are convergence-tested numerical examples on the prescribed geometry. [All-wavelength bounds](../bounds/NOTE.md)

## Local Ward signs and the continuum relative-energy identity

Take signature (−,+,+,+),

$$
 T^\Theta_{\mu\nu}=\partial_\mu\Theta\partial_\nu\Theta
 -g_{\mu\nu}\big[\tfrac12(\partial\Theta)^2+V\big],\qquad
 T^\chi_{\mu\nu}=Z\partial_\mu\chi\partial_\nu\chi
 -\tfrac12 g_{\mu\nu}Z(\partial\chi)^2.
$$

The field equations imply the separate local identities

$$
 \nabla_\mu T_\Theta^{\mu}{}_{\nu}
 =+\tfrac12 Z_{,\Theta}(\partial\chi)^2\partial_\nu\Theta,
 \qquad
 \nabla_\mu T_\chi^{\mu}{}_{\nu}
 =-\tfrac12 Z_{,\Theta}(\partial\chi)^2\partial_\nu\Theta.
$$

In unit lapse and zero shift let h_ij be the spatial metric, kappa_ij=hdot_ij/2,
theta=tr_h(kappa), P^ij the spatial stress, and S^i=T^{0i} the energy flux.
For the daughter

$$
 \rho_\chi=\tfrac12Z\big(\dot\chi^2+|D\chi|_h^2\big),\quad
 S_\chi^i=-Z\dot\chi D^i\chi,\quad
 Q_\chi=-\tfrac12\dot Z\big(\dot\chi^2-|D\chi|_h^2\big).
$$

The parent flux is S_Theta^i=−Thetadot D^iTheta, and Q_Theta=−Q_chi. For
each sector the energy identity is

$$
 \dot\rho+\theta\rho+D_iS^i+\kappa_{ij}P^{ij}=Q.
$$

For radial chi on ds²=−dt²+e^{2L}dr²+e^{2M}sin²r dOmega²,
S_chi^r=−Z e^{−2L}chidot chi_r. Its contribution is
sqrt(h)^{-1}partial_r[sqrt(h)S_chi^r], with the displayed positive divergence
on the left. This also fixes the outward flux sign. Regularity at both
poles makes the boundary term vanish on the complete closed slice. For
nonradial harmonics the corresponding angular divergence integrates to
zero. The complete-slice identity is therefore

$$
 \frac{dE_\chi}{dt}
 =-\int\sqrt h\,\kappa_{ij}P_\chi^{ij}\,d^3x
   +\int\sqrt h\,Q_\chi\,d^3x.
$$

The first term includes anisotropic metric work. Replacing it by 3H times
an isotropic pressure on the inhomogeneous slices would change the test.

The same identity holds for the difference of expectation values in two
states of the linear daughter on the same smooth background, using the
same stress and composite-operator renormalization. For Hadamard states
whose two-point difference is smooth, all state-independent local
subtractions cancel. In particular a coherent displacement by a smooth
classical solution gives exactly

$$
 \Delta\langle T_{\mu\nu}\rangle=T_{\mu\nu}[\chi_{cl}],\quad
 \Delta E=E[\chi_{cl}],\quad
 \Delta W_Z=W_Z[\chi_{cl}].
$$

This is a continuum statement, without a finite-mode vacuum energy. Smooth
finite initial harmonic sums remain smooth under the hyperbolic evolution;
on a compact slice their relative energy is finite. Thus the classical
state-dependence mechanism also has a finite relative-energy quantum
interpretation through coherent excitations above a common reference
state. The calculation reports converged Galerkin diagnostics for specified
initial supports, rather than a theorem that n≤24 bounds every continuum
initial state.

For a positive finite covariance the same Rayleigh interval bounds the
ratio of expectation-value quadratic forms. A difference of arbitrary
quantum states need not have positive covariance, so positivity alone does
not bound an arbitrary subtracted work/energy ratio by these extrema.
No independent claim about absolute renormalized vacuum energy follows.
