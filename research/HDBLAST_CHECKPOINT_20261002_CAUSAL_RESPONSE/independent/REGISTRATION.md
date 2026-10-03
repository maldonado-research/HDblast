# Independent forced canonical modes: prospective registration

This registration precedes every source or mode evaluation. It is a separate
implementation of the forced canonical mode equations, using no primary-route
memory routine. Execution requires the root's public freeze commit identifier.
All claims concern the scalar variance on one fixed geometry and fixed state.

## Frozen model and source

H=1, r=2, a(eta)=-1/eta, x0=2, initial eta=-6, final eta=-1.5.
The two sources are s=epsilon B(u) and s=epsilon u B(u), where epsilon=1e-4,
u=eta+4, B(u)=exp(1-1/(1-u*u)) for |u|<1 and zero otherwise.
Source IDs are `positive_B` and `signed_uB`. Observations are eta=-5.5,-4.5,-4,
-3.5,-2.5,-1.5. K is 64,128,256. The initial BD density matrix is held fixed.
The source epsilon is a linear forcing scale, not a finite-mode nonlinear test.

## Independent evolution

Write delta v=v u, with v=exp(-ik eta)/sqrt(2k). Evolve

    u'=w, w'=2ik w-s, u(eta_i)=w(eta_i)=0.

This is exactly delta v''+k^2 delta v=-s v, with delta v and delta v' initially
zero. On each time step use the exact homogeneous 2x2 propagator. Integrate
only the local step's source with an independently coded eight-point
Gauss-Legendre rule. No previously integrated global-history formula is used.
Evaluate phi(z)=expm1(2ik z)/(2ik), with its z limit if k=0; momentum nodes
are strictly positive. Preserve the complex u,w variables without imposing
the Wronskian relation after a step.

The linear Wronskian residual divided by the unperturbed W=i is
2Re(u)-Im(w)/k. Monitor its absolute maximum/epsilon over every time step
and momentum node. At observations additionally reconstruct v,delta v and
their derivatives and evaluate the physical linear Wronskian independently.

The combined momentum integrand is

    F(k)=4k Re(u)+s k^2/(k^2+M^2)^(3/2), M=a sqrt(r),
    y_K=a^2 delta Q_K/epsilon=integral_0^K F(k)dk/(8pi^2 epsilon).

Combine its terms pointwise before momentum integration. Preserve the full
mode arrays, subtraction terms and combined integrand at every observation
for inspection of high-k cancellation. The source-free pre-pulse response
must be exactly zero.

## Two settings, budget, gates

Both settings use double precision, eight time Gauss nodes, sixteen momentum
Gauss nodes per panel, and the same three K endpoints.

* Coarse: dt=1/128, momentum panel width=1/2; 576 time steps, 512 momentum
  panels, 8192 momentum nodes per source. There are 576*8=4608 distinct local
  time nodes per source and 37,748,736 local time/momentum pairs per source.
* Fine: dt=1/256, momentum panel width=1/4; 1152 time steps, 1024 momentum
  panels, 16384 momentum nodes per source. There are 1152*8=9216 distinct local
  time nodes per source and 150,994,944 local time/momentum pairs per source.

Run each source once per setting. Integrate to K=256 once and sum complete
panels up to each K endpoint. All six observations coincide with step ends.
Store only O(momentum nodes) state in memory; archived observation arrays
are written source by source. Linux ru_maxrss is interpreted in KiB. Expected peak memory below 128 MiB and wall
time below 15 minutes; stop and preserve outputs if either resource budget
is exceeded. No output-driven adjustment to dt, nodes, source or cutoffs.

For every eta,K, require |y_fine-y_coarse| <= 2e-10. Record twice that
difference as an empirical numerical-error estimate; it is not a rigorous
error bound and does not by itself exclude shared discretization error.
The independently coded finite-K primary formula provides an additional
cross-route check: require |mode_fine-primary_finiteK| <= 2e-9 plus the
primary quadrature error estimate, independently of the refinement gate.
Require both amplitude and reconstructed physical Wronskian residual/epsilon
<=1e-10. Require exact zero at eta=-5.5 for each of the two sources. These checks use
explicit exceptions and remain active under Python optimization.

Do not fold any removed-cutoff difference into the finite-K integration
estimate. The separately derived, actual source-derivative tail bound belongs
to the root's comparison with the removed-cutoff time-memory response.

## Evidence and stop condition

`MANIFEST.json` pins this registration, the implementation, and the inherited
analytic protocol. Runtime requires both `--freeze-commit` and the externally
pinned `--manifest-sha256`, and verifies the manifest identity and all file
pins before source evaluation. An edited, internally consistent manifest
therefore does not silently replace the publicly frozen one.
The execution record includes the supplied public freeze commit, runtime
versions, hashes, step/node counts, outputs, numerical differences, linear
Wronskian residuals and resource measurements. Before a run, record its active
source and setting; during it update the last step, time and Wronskian residual.
If an exception interrupts evolution, preserve the current k,u,w and linear
Wronskian residual arrays in an uncompressed failure snapshot, alongside the
exception and any already completed runs. Preserve failures. If a gate
fails, stop with a bounded diagnosis; do not retune or claim successful
calibration. No stress response, stability, heating, or particle yield is
calculated by this route.
