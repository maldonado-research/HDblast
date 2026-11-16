# Reduced nonlinear evolution, time continuity and causal support

14 November 2026 UTC. Integrated conditional theorem for the actual
normalized sixteen-field mixed system. The hypotheses retain the
authenticated smooth background, coefficient chart, selected initial
family and smooth frozen-linear construction. This is a reduced-system
result; subsidiary, harmonic and Israel-condition propagation are not
inferred from it.

## 1. Exact statement and unchanged numerical admission

Let `ell=10^-6`, `Sigma=[0,6ell] x S3`, `T=13ell/8`, with the fixed
positive product bundle and all full covariant/Frobenius norms of the
companion proofs. Use the same-state elapsed reflection
`U(s)=J_E V(-s)`, including the opposite initial elapsed velocity.
The exact reduced equation and rows are

```
U_ss-A(x,U)^i D_i U_s-B(x,U)^ij H_ij U=F(x,U,DU),
F=C(x,U)+E(x,U)DU+G(x,U)(DU,DU),
U_D=0, D_z U_N=0 at z=0;       U=0 at z=6ell.
```

There are four D and twelve N fields at the shell. The principal symbol
is scalar on all sixteen fibers; all mixed coefficients and actual
reference terms are retained. In particular `C(x,0)=0` and
`F(x,0,0)=0`. The first jet DU includes the elapsed velocity.

Take an exact Gamma-invariant `H14 x H13` Cauchy pair `(p,q)` supported
in `A0=[ell/2,5ell/2] x S3`, in the stated compact-support closure, with

```
d=||p||H14+||q||H13<=Dbar,
Dbar=10^4200 a_iter,           a_iter=10^(-10^500209).
```

For the selected canonical initial family the same restriction
`|epsilon|<=a_iter` supplies this bound, without reparameterization.
Its supported smooth approximants may have size up to2Dbar, and their
nonlinear jets are recomputed from each pair. They need not solve the
canonical constraints themselves. All convergence and reconstruction
refer to the prescribed exact initial datum.

The accepted companion arguments give a unique reduced solution with

```
partial_s^j U in C([0,T];H^(14-j)),         0<=j<=14,
X14(U;s)=sum_(j=0)^14 ||partial_s^j U(s)||H^(14-j)
                  <=C0 d exp(2Hs)<=M<R/4,
C0=10^7003210,       H=10^500212,           R=10^-12,
M=2 C0 Dbar exp(2HT)<2*10^-15.
```

Every jet is the actual derivative of the same field. All initial jets
J0 through J14 are attained in their displayed strong time topologies.
All boundary traces defined at these regularities hold; the proof does
not assign an ordinary boundary trace to the L2 fourteenth jet.

The solution is Gamma invariant, has all symmetry-allowed angular
content, and is unique in the retained lower-continuous/top-bounded
high ball. The same statements and constants hold on one copy and the
specified identical-two-copy Hilbert union with global input norms.

## 2. How the companion proofs join

The low-energy proof gives the evaluated source Lipschitz constant
`L_delta=4*10^500000`, the common kernel coefficients
`G_en=16*10^500000`, `beta=12*10^500000`, and full-slab factorial
convergence in `C H1 x C L2`. It treats the varying iterate energy
metrics by uniform norm comparison; no small-time contraction or
equation for the auxiliary seed is assumed.

The lower mixed proof first upgrades only the field/velocity pair by
ordinary spatial Rellich compactness, then uses the exact shifted PDE
recursion to get `C X13`. Row twelve gives the actual fourteenth jet
in `C H^-1`. The baseline limit proof retains every top norm by weak
continuity and lower semicontinuity, proves the nonlinear equation and
all meaningful rows, and proves approximation-independent finite-data
existence and low-energy uniqueness.

The additional top-continuity proof establishes the highest commuted
weak rows, reflects the scalar-principal system on the fourfold radial
double, and applies the closed-manifold weak linear energy lemma.
It obtains continuous energy for all canonical time/rotation words
through degree thirteen. The accepted interior normal recurrence then
recovers strong continuity of the entire X14 table. This is an added
PDE/weak-boundary argument, not an inference of strong continuity from
a weak bound alone. Its qualitative regularization constants do not
enter the numerical admission, difference or factorial-tail bounds.

For actual solutions in the common high ball, the evaluated stability
bound is

```
sup_s X1(U-V;s)<=S_low (||p-p_tilde||H1+||q-q_tilde||L2),
S_low=3 exp((91/2)*10^499994).
```

Together with the uniform high bound, this also gives qualitative
continuity into `C X13` for an initial-pair sequence converging in
H14 x H13 within this admitted class: use the same pair Rellich
upgrade and exact nonlinear row induction. In particular the original
parameter curve maps continuously into `C X13`. This is not a top-X14
solution-map continuity or a parameter differentiability theorem.

## 3. Full radial flux and its sharp retained speed ceiling

At each point set `v=U_s`, `V_i=D_i U`, and include the field U in
the first-order state `Y=(v,V_i,U)`. The positive symmetrizer is
`S=diag(1,B,1)`. In a fixed orthonormal spatial frame the normal
principal matrix has the only nonzero blocks

```
(M^z Y)_v=A^z v+B^zj V_j,
(M^z Y)_(V_i)=delta_i^z v.
```

The full symmetric normal flux is

```
Phi^z=1/2<S M^z Y,Y>
      =(A^z/2)|U_s|^2+B^zj<D_j U,U_s>.
```

Its B^zA terms are retained. Since S M^z is symmetric and S is
positive, the relevant eigenvalues are real. For a nonzero eigenvalue,
the tangential V_A components vanish and the remaining two equations
give

```
lambda^2-A^z lambda-B^zz=0.
```

The other eigenvalues are zero. Including the field's zero-order block,
the scalar-fiber characteristic polynomial is
`lambda^4(lambda^2-A^z lambda-B^zz)`. Thus off-diagonal B^zA terms
contribute to the flux but do not change these normal speeds.

The already proved margins `|A^z|<=1/19`, `B^zz<=20/19` imply

```
max|lambda|
 <=[1/19+sqrt(1/361+80/19)]/2
 =[1/19+39/19]/2=20/19=:c.
```

Consequently, with the complete energy density
`e=(|U_s|^2+B^ij<D_iU,D_jU>+|U|^2)/2`,

```
|Phi^z|<=c e,                 c=20/19<2.                       (1)
```

This is a generalized symmetric-matrix energy estimate, so it does
not discard tangential gradients or require a single angular mode.
The radial speed is in the registered z and elapsed s coordinates;
it is not a proper-time or physical-unit calibration.

## 4. Moving-domain energy proves actual finite propagation

This step applies to the exact solution whose initial pair is supported
in A0. Smooth approximants may have enlarged support, and the auxiliary
seed has enlarged heat-cutoff support. Neither of those supports is
substituted for the support of this exact initial pair.

The state and full first jet lie in the admitted patch. Put
`p0=10^7 R=10^-5`, and write
`Z^2=|U|^2+|U_s|^2+|Dsp U|^2`. The zero-state identity and the actual
coefficient bound give pointwise

```
|F(x,U,DU)|<=Cstar[|U|+(1+p0)|DU|]<=2 Cstar Z,
Cstar=10^500000.
```

The last comparison uses `1+(1+p0)^2<4`. Hence the forcing part of
energy production is at most `2Cstar Z^2<=8Cstar e`. The already
proved remaining pointwise production is at most `21Cstar e`.
Their sum is bounded by `32Cstar e`.

For the exterior domain `r_+(s)<=z<=6ell`, let
`r_+(s)=5ell/2+c s`. The original data and velocity vanish there at
s=0. The fixed outer all-Dirichlet flux is zero. Differentiating the
moving-domain energy, its moving inner face contributes

```
-Phi^z(r_+(s))-c e(r_+(s))<=0
```

by (1). Therefore its nonnegative energy E_+(s) satisfies
`E_+'<=32Cstar E_+`, with `E_+(0)=0`. Gronwall yields E_+=0.
The constructed regularity amply justifies this calculation and its
traces; alternatively a smooth radial cutoff and its limit give the
same identity. On the entire slab the moving radius remains inside
the artificial wall, since

```
r_+(T)=[5/2+(20/19)(13/8)]ell=(80/19)ell<6ell.
```

Thus the exact perturbation vanishes on the fixed outer collar
`(80ell/19,6ell] x S3` for all s, with a margin of
`6ell-80ell/19=(34/19)ell`. In particular the artificial wall is
separated from the perturbation throughout the constructed interval.

For the inner domain `0<=z<=r_-(s)`, put
`r_-(s)=ell/2-c s` while r_->0. Its fixed shell flux is zero by
the actual D/N rows, scalar principal structure and
`A^z=B^zA=0`. Its moving outer face contributes
`Phi^z-c e<=0`. The same zero-data Gronwall argument shows

```
U(s,z,omega)=0 for 0<=z<ell/2-(20/19)s,
                    0<=s<19ell/40.                           (2)
```

Continuity gives the corresponding zero shell value at the endpoint
as well. This is a lower bound on the possible arrival time, not a
claim that a nonzero signal actually arrives then.

All energy arguments include every fiber and angular mode, and add
over the retained identical copies with no new norm factor. The
zero-field background is preserved by the pointwise inverse state
chart, so these zero regions also remain unperturbed under that
registered reconstruction. No unreduced constraint-propagation
statement follows from this observation.

## 5. Exact remaining scope

This theorem proves reduced nonlinear existence, uniqueness,
individual strong top-order time continuity, evaluated low-energy
stability and actual radial support control on the unchanged admitted
interval. The extremely small parameter ceiling is a sufficient
mathematical threshold, not an observed cosmological amplitude.

It does not prove strong C X14 convergence of the original Picard or
data-approximation sequence, top-order initial-data continuity, an
X15 tame comparison, evolved parameter derivatives, harmonic/
subsidiary/Israel propagation, a proper observer tube, an observable
nonlinear remainder or empirical higher-dimensional origin. Those
are distinct later tasks. The qualitative higher-data regularity
routes must preserve the already protected low X14 ball.
