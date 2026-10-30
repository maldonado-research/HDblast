# Exact residual-corrected natural bulk Bianchi source binding

30 October 2026. Internal finite algebra audit, assembled as an additive checkpoint. This does
not claim external specialist review. Old source files and the registered map
remain unchanged.

The checked identity is the eight-modal-row relation

```
B0(QE) + p Qphi_general e^y + 2 deltaB_q(E1) + 2 S1 d eta = 0,
B_g(T)_b = nabla^a T_ab - (1/2) nabla_b tr_g T.
```

Here QE and Qphi_general are the second parameter derivatives along the affine
path g0+epsilon q, phi0+epsilon eta. The registered natural QE is used directly;
it is not replaced by the reduced Hessian or by the principal forcing. The
first metric has seven arbitrary profiles A,B,C,D,E,L,M and eta=sY. In particular
B and C are free bulk functions. The checker permits all 70 ordinary metric
jets through total order three and all six ordinary scalar jets through order
two. No first equation, harmonic condition, shell restriction or reduced
acceleration is substituted.

The residual correction is constructed without reading the quadratic source
maps. A separate physical-frame implementation forms q, nabla q and the fully
lower connection variation

```
P_cab = (nabla_a q_bc + nabla_b q_ac - nabla_c q_ab)/2,
Ric1_ab = sum_c signature[c] (nabla_c P_cab - nabla_b P_cac),
E1_ab = Ric1_ab - (dphi0)_a (deta)_b - (deta)_a (dphi0)_b
        -(2/3)(V' eta g0_ab + V q_ab),
S1 = Box0 eta - q^{ab}(nabla^2 phi0)_ab - C1^a(dphi0)_a - V'' eta.
```

It then differentiates this E1 tensor and contracts

```
deltaB_q(E1)_b = -q^{ac} nabla_a E1_cb
                +(1/2)q^{ac} nabla_b E1_ac - C1^d E1_db.
```

The linear identity B0(E1)+p S1 e^y=0, independent explicit C1 and S1 formulas,
and E1 symmetry are checked before the quadratic identity. The correction's
full angular polynomial is reconstructed exactly in the three scalar and two
gradient modes. This is an identity of collected rational Laurent
coefficients; it does not use point samples or intervals.

There are six stored frame indices but the spacetime dimension is five:
index 0 is y, index 1 is the timelike frame direction, and indices 2--5 are
four redundant ambient tangent slots for the three-sphere. Their projector is
P_ij=delta_ij-n_i n_j. The ordinary y/t derivative chain acts on free profile
symbols. Every covariant tensor slot receives its connection term, including
slots in P. The unit-sphere polynomial quotient is n3^2=1-n0^2-n1^2-n2^2.
This constructs coherent third metric derivatives; it does not treat arbitrary
covariant third components as independent.

The background substitutions are explicit:

```
rho_y=a rho,       a_y=-rho^-2-p^2/3,       p_y=V'-4ap,
V=6rho^-2+p^2/2-6a^2,
(V')_y=V'' p,      (V'')_y=V''' p,
c_t=c theta,       theta_t=1-theta^2,       theta^2=1-c^-2.
```

The algebraic Einstein constraint is included through the displayed V
substitution. The checker independently verifies E0_yy=E0_slice=S0=0 and
D_y V=V' p. Consequently the delta-squared-B(E0) term vanishes as a
neighbourhood identity. The de Sitter quotient is stable under the stated
time derivation. Nonzero rho and c are understood in the Laurent algebra.
No potential derivative above V''' is needed for this natural identity.

The source side decodes the unchanged JSON monomials and differentiates them
by the same specified ordinary chain rule. Its exported angular tensor slots
are already (J:0,J:2,J:4,H:2,H:4,T:+,T:-); no second basis conversion or
trace reversal is applied. For lambda=l(l+2), beta=2-2lambda/3 and H=theta/rho,
the B0 rows used are

```
B_y = (K_y+N_y-3J_y)/2-L_t/rho
      +a(4K+N-3J)-3H L-lambda v/(rho c),
B_0 = L_y-(K_t+N_t+3J_t)/(2rho)
      +5a L-3H(N+J)-lambda u/(rho c),
B_a = v_y-u_t/rho+5a v-4H u
      +(-K/2+N/2-J/2+beta H_tensor)/(rho c).
```

There are three y rows, three temporal rows and two angular-gradient rows.
The two TT tensor outputs lie in the divergence/trace kernel, so success does
not independently validate those QE outputs. An independent review checks
these B0 formulas by full-frame contraction on arbitrary tensor profiles.

The homogeneous specialization

```
B0(QE) = -p Qphi_general e^y
```

requires E1=0 and its first covariant prolongation nabla E1=0, as well as
S1=0. The natural neighbourhood statement E1=0 implies the prolongation;
pointwise E1 values alone do not. Replacing Qphi_general by Qphi additionally
requires C1=0. Harmonicity alone does not remove either residual correction.
This audit proves the general residual-corrected coefficient identity and the
stated conditional implication. It does not independently establish actual
master reconstruction membership in the first Einstein/scalar equations.

`checks/run_bianchi_verification.py` pins the checker, exact Laurent arithmetic, angular
algebra and registered JSON. The checker authenticates every nonstdlib file
it executes before and after use. AST checks reject assertions, unapproved
import roots and direct dynamic-code entry points. The original geometric
producer is not imported. Authentication claims apply to these listed bytes,
not to the entire old repository or every upstream dependency.

Normal and optimized baselines run in fresh isolated stdlib-only Python
subprocesses. Their reports must be byte-identical. Thirteen actual
corruptions run in each mode, with absent output paths beforehand; each must
exit nonzero, produce a collected NONZERO witness, and leave no successful
report. The controls replace QE by QEred or forcing, halve QE, prematurely
use harmonic Qphi, omit either first-residual correction, reverse temporal
contractions, remove the half-trace factor, remove radial/FRW connection
terms, corrupt the V' chain, freeze rho radially, and reverse the scalar
source sign. Wrong radial/potential chain controls may fail an earlier linear
Bianchi precheck, which is recorded honestly in their exact witness.

Reproduce from the repository root using:

```
python -I -S -B research/HDBLAST_NATURAL_FIELD_BINDING_20261030/checks/run_bianchi_verification.py \
  --repo-root . \
  --output-dir /tmp/fresh-natural-bianchi-verification --workers 3
```

The output directory must not already exist. Every child command is
`python -I -S -B [-O] check_natural_bulk_bianchi.py --repo-root ROOT --output
OUTPUT [--mutation MUTATION]`; modes, mutations, return codes and exact failure
witnesses are recorded in VERIFICATION.json. The canonical normal/optimized
reports enumerate the eight source identities and the independent linear,
angular and background checks.

This adds an implementation binding to the existing analytic Bianchi and
subsidiary argument. It supplies no new evolution or uniqueness theorem,
nonlinear family, numerical accuracy gain or cosmological-origin result.
