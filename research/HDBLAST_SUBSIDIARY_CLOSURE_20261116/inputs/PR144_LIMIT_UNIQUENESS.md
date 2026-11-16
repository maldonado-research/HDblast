# A nonlinear reduced-system limit and low-energy uniqueness

14 November 2026 UTC. Conditional mathematical proof. This
constructs a finite-regularity solution of the actual normalized reduced
sixteen-field mixed system from the admitted initial data. It retains
the authenticated background, coefficient-chart, smooth frozen linear
construction and initial-data premises. It does not conclude the
unreduced Einstein system, propagation of subsidiary/Israel conditions,
continuous full X14, evolved parameter differentiability or an observer
theorem. Every preceding private root and public input remains read-only.

## 1. Actual system, data, and solution class

Work on the original elapsed past slab

```
Sigma=[0,6ell] x S3,   ell=10^-6,   0<=s<=T=13ell/8,
U(s)=J_E V(-s),
E=R_v+Sym^2(R dt+T*S3)_H+R_w+(R dt+T*S3)_D.
```

The reflection acts on every tensor time slot and changes the elapsed
initial velocity sign. It represents the same original past field. Use
the fixed positive product connection, full ordered covariant derivative
slots, full sphere measure and the full Frobenius symmetric-tensor norm.
The mixed norm is the literal sum

```
X_m(U;s)=sum_(j=0)^m ||partial_s^j U(s)||H^(m-j).
```

The equation and boundary rows are

```
P_U U=F(U),
P_v=partial_s^2-A(x,v)^i D_i partial_s-B(x,v)^ij H_ij,
F(v)=C(x,v)+E^mu(x,v)D_mu v+G^mu nu(x,v)(D_mu v,D_nu v),
U_D=0, D_z U_N=0 at z=0;     U=0 at z=6ell.             (1)
```

There are sixteen independent fields, four Dirichlet and twelve Neumann
components at the shell. The Hessian is the complete product Hessian.
The normalized lower source includes the entire first-jet connection
substitution; in particular `C(x,0)=0`. Coefficients and first jets in
(1) retain the same all-state conventions as the accepted coefficient
certificate, not componentwise replacements of covariant derivatives.

Let `(p,q)` be a Gamma-invariant finite Cauchy pair in `H14 x H13`,
in the stated compact-support closure on
`A0=[ell/2,5ell/2] x S3`. Its norm is

```
d=||p||H14+||q||H13.
```

The pair may be the selected canonical initial datum after its complete
initial chart. The reduced analytic construction also applies to such
collar-supported pairs without imposing canonical constraints on their
smooth approximants. Set

```
a_iter=10^(-10^500209),    Dbar=10^4200 a_iter,
d<=Dbar,                 R=10^-12.                    (2)
```

For the selected original canonical family, (2) follows from
`|epsilon|<=a_iter` and the previously proved actual initial-pair bound.
No parameter is rescaled. The nonlinear initial recursion defines
`J_j(p,q) in H^(14-j)`, `0<=j<=14`, with `J0=p,J1=q` and

```
I14=sum_(j=0)^14 ||J_j||H^(14-j) <=10^6503000 d.         (3)
```

The guard for this recursion is `d<=10^-6503020`; the nested-power
comparison in the accepted admission proof places even `2Dbar`
strictly inside that guard. The finite jet map is continuous on the
relevant pair ball. These are the accepted temporal-seed inputs.

The conclusion below uses the class

```
U_j=partial_s^j U in C([0,T];H^(13-j))
                    intersect C_weak([0,T];H^(14-j)),  0<=j<=13,
U_14 in C([0,T];H^-1) intersect C_weak([0,T];L2),
sup_s sum_(j=0)^14 ||U_j(s)||H^(14-j) < R/4.            (4)
```

Here `H^-1=(H0^1(Sigma;E))*` with the fixed positive product norm.
It uses test sections zero on the entire spatial boundary; no boundary
condition is being imposed on U14 by this notation. All temporal jets
in (4) are distributional derivatives of the same field.

## 2. Evaluated inputs and two distinct difference identities

The accepted frozen estimates have

```
Cstar=10^500000,   C_Y=10^30,
K_u=10^500101,    K_v=10^500100,
N=10^110,        L_u=L_v=10^500208.
```

They hold for an arbitrary smooth frozen successor under the predecessor
ball `sup X14(v)<=R`; no successor-ball assumption enters their proof.
The normalized metric margins are

```
|A|<=1/19,       (931/1280)I<=B<=(20/19)I,
19/20<=Bzz<=20/19.
```

The same-point positive rho factor has been cancelled before these
ratios are estimated. For the scalar-principal auxiliary energy

```
E_v[h]=(1/2)integral (|h_s|^2+B^ij<D_i h,D_j h>+|h|^2),
Q(h)=||h||H1+||h_s||L2,
sqrt(E_v[h])<=Q(h)<=3sqrt(E_v[h]).                     (5)
```

The last inequalities follow from `Z1^2/4<=E_v<=Z1^2`,
`Z1^2=||h||H1^2+||h_s||L2^2`, and `2sqrt(2)<3`.

The companion low-energy difference proof supplies the evaluated numbers

```
L_delta=4 Cstar,    G=16 Cstar,    beta=3L_delta=12 Cstar. (6)
```

It retains the joint state/first-jet segment in F, both differentiated
slots of `G(Dv,Dv)`, and both principal coefficient differences. Its
uncommuted energy kernel is `exp(G(s-r))`.

There are two different exact uses of (6). For consecutive Picard
iterates, `W^n=U^(n+1)-U^n`, `n>=1`,

```
P_(U^n) W^n = F(U^n)-F(U^(n-1))
 +(A(U^n)-A(U^(n-1)))^i D_i U^n_s
 +(B(U^n)-B(U^(n-1)))^ij H_ij U^n.                     (7)
```

The differentiated principal factor is the preceding solution U^n.
The source is bounded in L2 by `L_delta Q(W^(n-1))`.
The seed is not a frozen solution, so (7) starts at n=1.

For two actual solutions v,w of (1), set `h=v-w`. The identity is
instead

```
P_v h = F(v)-F(w)
 +(A(v)-A(w))^i D_i w_s
 +(B(v)-B(w))^ij H_ij w.                               (8)
```

Its second-derivative factor is the second actual solution w. The
same pointwise proof gives `||right side||L2<=L_delta Q(h)`
whenever both solutions lie in the high ball R. Neither formula replaces
the indicated differentiated field by the newest unknown difference.

The admitted smooth iteration, with its actual nonlinear initial jets
and compatible auxiliary lift, has the common bound

```
X14(U^n;s)<=C0 delta exp(2Hs),
C0=10^7003210,             H=10^500212,
delta<=2Dbar,
M:=2 C0 Dbar exp(2HT) <2*10^-15<R/4.                   (9)
```

The displayed M is one common strictly interior majorant, independent
of iterate number, smoothing index and original allowed parameter.
It is not a sequence of unrelated strict estimates whose limit could
reach R/4. All union/copy norms in (3),(5),(9) are already global norms.

## 3. The limit of smooth compatible iterates

First fix a smooth supported Gamma-invariant Cauchy pair of norm
`delta<=2Dbar`. Recompute its entire nonlinear initial table through
J14, then apply the accepted heat/cutoff seed. The seed and every
successor have the same jets J0 through J14. The seed vanishes on fixed
open endpoint collars for all auxiliary times; it need not agree with
all higher nonlinear interior jets. Each subsequent smooth frozen
problem has all-order corner compatibility because its recursively
computed higher initial jets are zero on those endpoint collars.
The retained smooth frozen construction gives a smooth successor on
the whole compact coefficient interval protected by (9).

The difference proof from (7), after a single exponential weight,
gives the explicit uniform tail, for `N>=1`,

```
sup_s Q(U^m-U^N;s)
 <= (R/2) exp((G+beta)T) (beta T)^N/N!,      m>N.        (10)
```

Indeed `Q(W^0)<=R/2` by the admitted two-iterate bound, while
`Q(W^n;s)<=beta integral_0^s exp(G(s-r))Q(W^(n-1);r)dr`
for n>=1. Ordered integration gives the factorial. Summing the
tail uses `sum_(n>=N) a^n/n!<=a^N exp(a)/N!`.
For every fixed finite `beta T`, the right side tends to zero.
For example its successive ratio is at most one half once
`N+1>=2 beta T`, so this is convergence on the entire prescribed slab,
without requiring a new small-time condition. The same numerical
tail works for every admitted smooth datum and supported approximant.

Consequently the fields converge in `C H1` and their actual velocities
in `C L2`. The companion lower-mixed-limit proof then applies the
ordinary spatial interpolation separately to that Cauchy pair and its
uniform H14/H13 bounds. It gives strong `C H13 x C H12` convergence.
It does not interpolate unrelated temporal jets. The equation,
differentiated successively through temporal row eleven, gives

```
U_j^n -> U_j in C H^(13-j),             0<=j<=13.        (11)
```

In this induction the coefficient/source argument is U^n and the
successor is U^(n+1); the predecessor shift is retained. All nonlinear
products converge in the indicated spaces by the actual finite
low/high derivative allocation. This supplies the lower interface
used in the rest of this proof.

The limit solves (1). In fact `U0 in C H13`, `U1 in C H12`, and
`U2 in C H11`; its nonlinear coefficient, principal and source terms
are well-defined continuously in H11, and the uncommuted equation
passes there from the iterates. The limit is therefore an actual
reduced solution, not merely a limit of its initial jet tables.
The higher distributional derivative identities and boundary traces
are justified below. No spatial C-infinity assertion for this limit
is needed or inferred from the smoothness of the individual iterates.

## 4. The fourteenth jet in a negative space

The equation differentiated twelve times, before any spatial normal
division, gives at the iterate level

```
U_14^(n+1)= A(U^n)^i D_i U_13^(n+1)
           +B(U^n)^ij H_ij U_12^(n+1)+R_n.             (12)
```

Here R_n contains every positive coefficient-hit term and the
twelfth temporal derivative of F(U^n). In every such product each
individual actual field factor has intrinsic order at most thirteen,
and total field order is at most fourteen. A high factor is put in
L2; every other factor has order at most seven and converges in
L-infinity using the three-spatial-derivative reserve of (11).
The source's zero-state extraction and both first-jet slots are retained.
It follows, as proved in the companion negative-order handoff, that
`R_n -> R` in `C L2`.

The leading terms are not asserted to converge strongly in L2.
Instead, `U13^(n+1)->U13` in C L2 implies
`D U13^(n+1)->D U13` in C H^-1, and
`U12^(n+1)->U12` in C H1 implies
`H U12^(n+1)->H U12` in C H^-1. The composed A and B converge
in C W1-infinity by U0 convergence in C H13. Such coefficients
are continuous multipliers of `(H0^1)*`: after pairing against an
H0^1 test section, their first derivative is the only added term.
Thus (12) proves

```
U14^n -> U14 in C H^-1,      U14 in C H^-1.              (13)
```

No fifteenth time jet or order-fifteen spatial estimate is required.
Passing the identity `partial_s U13^n=U14^n` against compact time
tests in H^-1 identifies U14 with the fourteenth distributional
derivative of U. The lower derivative identities similarly follow
from (11), or by passing the fundamental theorem of calculus in the
common weaker Sobolev spaces. At s=0, (12) is precisely the nonlinear
initial row twelve. All iterates have U14^n(0)=J14, so (13) gives
U14(0)=J14 in H^-1. The L2 identification in the next section upgrades
this to equality of the L2 initial traces.

## 5. Pointwise high bounds and weak continuity at every top order

Use the following elementary weak-continuity lemma. If H is a Hilbert
space continuously and injectively embedded into a weaker Banach space
B, `f_n -> f` in C B, and `||f_n(t)||H<=M` uniformly, then
`f(t) in H`, `f_n(t) -> f(t)` weakly in H for every t, and
`f in C_weak H`. To prove it, take weakly convergent subsequences in
the bounded H ball. Their image in B has the unique already known
limit f(t), which identifies every cluster point. Apply the same
argument to arbitrary times t_k tending to t, using continuity in B.
This also proves the assertion without selecting a time-dependent
subsequence or changing the continuous representative.

Apply the lemma to every jet in (11), with
`H=H^(14-j), B=H^(13-j)`, and to (13), with
`H=L2, B=H^-1`. This proves all weak continuity assertions in (4).
The convergence is pointwise weak in the finite product of the fifteen
top Hilbert spaces. The literal sum of their norms is weakly lower
semicontinuous, so the time-dependent majorant from (9) passes at
each time, with no extra factor fifteen:

```
sum_(j=0)^14 ||U_j(s)||H^(14-j)
 <=C0 delta exp(2Hs)<=M<R/4,         every s in [0,T].    (14)
```

In particular this is stronger than a merely essential top bound, but
is weaker than continuity of the top norms. The weakly continuous
representative has U14(0)=J14 in L2 because the H^-1 trace is J14
and the embedding L2 into H^-1 is injective. For j<=13, (11) and
the same injectivity identify U_j(0)=J_j in H^(14-j). Every initial
jet through fourteen is therefore the actual trace of the same solution.

## 6. Actual-solution uniqueness and low-energy initial-data stability

Consider any two solutions v,w in the lower-continuous/high-bounded
class of (4), with the same prescribed rows, allowing the larger
common high ball R. They need not arise from the same approximation
scheme. Their difference satisfies (8). The uncommuted energy proof
is directly valid for this regularity: v,w belong to C H13,
their velocities to C H12 and their accelerations to C H11. The
coefficients have continuous first spatial/time derivatives. Every
term in the energy identity and spatial integration by parts has its
ordinary strong meaning at these orders. Thus no extension of a
smooth-only uniqueness statement by an unproved limiting argument is
being used.

At the shell, the actual state Dirichlet traces give `A^z=B^zA=0`;
the difference's orthogonal D/N traces annul the remaining flux. At
the outer wall the difference has zero Dirichlet velocity. The
principal scalar structure acts identically on all field components.
The regularized square-root energy inequality, (5),(6),(8) and
Gronwall yield

```
Q(v-w;s) <=3 exp((G+beta)s)
                 [||v(0)-w(0)||H1+||v_s(0)-w_s(0)||L2],
G+beta=28*10^500000.                                  (15)
```

Here `sqrt(E_v[h](0))<=Q(h;0)` by the same normalized upper margin,
including when the two initial pairs differ. A fully evaluated uniform
stability constant is, without expansion,

```
S_low=3 exp((91/2)*10^499994),
sup_s Q(v-w;s)<=S_low Q(v-w;0).                        (16)
```

The rational exponent follows exactly from
`28*10^500000*(13/8)*10^-6=(91/2)*10^499994`.
It is finite and independent of angular cutoff, smoothing or iteration.
Zero initial difference implies equality of the solutions on [0,T].
Equation (16) is low-energy initial-data stability in this common high
class. It is not a top-X14 Lipschitz or continuity assertion.

## 7. Finite H14/H13 data and independence of supported approximants

Choose smooth Gamma-invariant pairs `(p_k,q_k)` converging to `(p,q)`
strongly in H14 x H13, with support in the one fixed slightly enlarged
annulus

```
Aplus=[7ell/16,41ell/16] x S3
      strictly inside the seed-cutoff plateau [3ell/8,21ell/8] x S3.
```

This is ordinary compact approximation followed by the natural finite
Gamma average, an operator of norm at most one in every fixed
covariant Sobolev space. It does not assume preservation of the
nonsmooth edge of A0. For k large the pair norms d_k obey
`d_k<=2Dbar`. Discarding finitely many initial approximants has no
effect. If the datum is zero one may take all approximants zero.

Recompute each nonlinear jet table from its own pair before applying
the auxiliary lift. The accepted continuous finite jet map gives

```
J_j(p_k,q_k) -> J_j(p,q) in H^(14-j),  0<=j<=14.        (17)
```

Independent truncation of a previously fixed jet table is not used.
The locality of the nonlinear initial recursion preserves the fixed
zero boundary collars and provides all available finite corner
identities. The approximants need not solve the canonical constraints.

Sections3--5 construct from each smooth pair an actual reduced solution
U^k in (4), with the same common M. The adjective smooth here belongs
to the initial pair and its Picard iterates; an all-orders regularity
claim for the limiting U^k has not been assumed. Inequality (15)
applies to U^k,U^l by their established lower regularity. Since their
initial pairs are Cauchy even in H1 x L2, U^k is Cauchy in
`C H1 x C L2`. The same Cauchy-pair interpolation and equation-based
lower-jet induction used above, now with the exact nonlinear equation
rather than a predecessor shift, give

```
U_j^k -> U_j in C H^(13-j),             0<=j<=13.
```

The row-twelve argument gives convergence of U14^k in C H^-1.
The weak-continuity lemma and weak lower semicontinuity give the
full pointwise top bound and weak traces of Section5. Because
`d_k -> d`, the sharper actual-datum estimate is

```
X14(U;s)<=C0 d exp(2Hs) <=M<R/4.                       (18)
```

Equations (1),(17) and all initial traces pass as above. This constructs
the finite-data solution. Comparing any two supported approximation
schemes by (15), or using uniqueness after passage, shows that it is
independent of the approximants and the auxiliary seed. The large but
fixed constant S_low presents no obstruction to qualitative Cauchy
convergence. No uniform high-frequency smoothing rate is claimed.

## 8. Boundary traces, symmetry, copies and exact stopping point

All asserted rows pass in trace spaces where those traces are defined.
For j<=12, the strong lower H^(13-j) convergence passes Dirichlet
traces strongly; the highest Dirichlet jet j=13 passes weakly from
H1. Thus the shell-D and outer-all-D rows hold through j=13.
For shell Neumann rows, strong lower convergence suffices through
j=11, while j=12 passes weakly from H2. The normal derivative of
U13 in H1 and either ordinary boundary trace of U14 in L2 are not
asserted. Equivalently, all top-defined traces pass: Dirichlet requires
H^s with s>1/2, and a normal derivative requires s>3/2. Their weak
continuity follows from bounded linear trace maps applied to the weakly
continuous top jets. The original uncommuted rows are in strong,
ample trace spaces and suffice for the energy uniqueness argument.

The actual normalized equation, chart and mixed boundary domain are
equivariant under the retained natural Gamma action. A Gamma-invariant
seed gives Gamma-invariant frozen coefficients/source, and uniqueness
of each smooth frozen problem gives invariant successors. Strong lower
limits preserve this closed invariant subspace. This argument uses
the natural action on every tensor slot; no finite angular-mode
restriction is imposed on the nonlinear solution.

All arguments apply either on one copy or on the specified identical
two-copy Hilbert union, with initial norms and volume already taken on
that union. No new copy factor is inserted into (3),(9),(10),(15),(18).
Identical initial copies remain identical by the same uniqueness. No
claim for an unspecified number of copies is made.

The result is a unique reduced normalized solution in the stated high
ball, continuous at every lower mixed level X13 and weakly continuous
in every top jet. It does not by itself establish strong C X14, a
Bona--Smith comparison, an X15 tame estimate, top-norm solution-map
continuity, evolved C2 parameter dependence or the later observable
construction. It also does not prove propagation of the harmonic
subsidiary, canonical constraints or Israel conditions needed to
identify a reduced solution with the full physical system. Initial
constraints of the exact original datum remain initial hypotheses.
Finite propagation and the eventual physical outer support margin
are separate; the seed's enlarged artificial support is not used as
that conclusion. No empirical higher-dimensional origin is inferred.

## Read-only dependencies

At immutable worktree base `9b74d57b6caef415775abc069838a4bcea303831`:

- `research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/`:
  the normalized margins, complete product/rotation interfaces,
  frozen energy, interior normal recovery and numerical smooth-Picard
  admission, including the six constants used above.
- `research/HDBLAST_TEMPORAL_JETS_SEED_20261114/`:
  the nonlinear initial jets, their finite-space continuity and actual
  support, compatible auxiliary lift, Gamma projection and admissible
  supported smooth approximation.
- `research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/`
  `QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md`, Section3:
  the retained smooth frozen linear construction.
- The new same-root `LOW_ENERGY_PICARD_CONVERGENCE.md`, accepted original-proof SHA256
  `92166b209542d2ee296ca8727a8101e2ece9ae6210336f63bf35970b2dfc0b04`,
  supplies (6),(10) and its actual source/factorial proof.
- The new same-root `LOWER_MIXED_LIMIT.md`, accepted original-proof SHA256
  `bcd2f404020acfd614600b55a49c244939fe0e4390292d38fb3d64aa0b76a6cb`,
  supplies (11)--(13), including the exact shifted rows and negative-order
  handoff. Its diagonal actual-solution-sequence variant is used for
  finite-data approximation.

Authentication and finite exact replays check the declared bytes and
algebraic/numerical models; they do not replace the conditional analytic
premises just stated.
