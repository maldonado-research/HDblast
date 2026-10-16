# Prospective full 23-amplitude response: a finite target and its residual route

This working note is outside the accepted capsule. It specifies a nontrivial
continuous candidate construction and checks an error budget. It does **not**
provide candidate coefficients, perform their quadrature, or certify a computed
response residual. The contraction corollaries below are standard consequences
of the retained construction, not a new existence theorem or physical claim.

## Retained equation and the two-term target

Use the actual selected initial prescription and the complete 23-dimensional
quadratic-feature space, including both TT orientations. On
`D={0<=x<=3, r>=0, x+r<=8}`, the accepted equation is

    U = b + T U,
    b = U_free + K F + L G,
    T V = K(A V_x + B V_r + D0 V) + L(M P_N V|r=0).

Here `||T||_X<=q=2109/1000000` and
`||V||_X=max(sup|V|_H/3,sup|V_x|_H,sup|V_r|_H)`. H uses the
probability-normalized sphere measure, the positive metric Frobenius norm,
and the scalar square. The accepted refined bound is `||U||_X<B=1.03e33`.

The five selected Dirichlet amplitudes and 18 Robin amplitudes use
`J=P_N-P_D`, not a common even extension. Precisely,

    (K f)(x,r) = (1/2) integral_0^x integral_(r-x+s)^(r+x-s) f_ext(s,v) dv ds,
    f_ext(s,-v)=J f(s,v),
    (L g)(x,r) = -integral_0^max(x-r,0) g(s) ds,
    U_free = (U0_ext(r+x)+U0_ext(r-x))/2
           + (1/2) integral_(r-x)^(r+x) V0_ext(v) dv.

`V0=U_x(0,.)` includes the inherited mixed shift-velocity conversion and
nonzero scalar compensator. The actual first master keeps its dynamic shell
trace ODE. That already constructed solution supplies the actual F and G;
this proposal does not replace the master boundary by an ordinary Robin row.
For the second unknown the retained rows are

    P_D U|b=0,   P_N U_r|b=M P_N U|b+G,
    M=-ell rho B_y,
    G=-ell rho(-2 QK,-Qp,GN),   GN=-2 tr_gamma QK-2 QC_y.

Thus the normal harmonic row is included alongside every physical shell row.
No second tangential harmonic boundary rows are added.

A concrete exact finite target is the first correction

    N1 = b + T b.                 [two retained Neumann terms]

It retains all of `KF`, `LG`, `K(A b_x+B b_r+D0 b)`, and `L(M b|b)`.
The exact identities are

    U-N1 = T^2 U,
    N1-b-T N1 = -T^2 b.

Consequently `||U-N1||_X<q^2 B=4.58131743e27`. N1 has the exact initial
and Dirichlet rows, but its Robin defect is
`P_N(N1)_r-M P_N N1-G=-M P_N(Tb)|b`. A finite truncation therefore must
not be described as satisfying the exact Robin row. That defect is included
in the complete integral residual and tail. The fixed-point comparison space
requires C1 regularity and the inherited zero wedge `x+r<1`; it does not
require each finite iterate to solve the final Robin equation exactly.

## Why C0 source errors suffice for a finite X certificate

For continuous physical f, endpoint differentiation gives

    (Kf)_x = (1/2) integral_0^x [f_ext(s,r+x-s)+f_ext(s,r-x+s)] ds,
    (Kf)_r = (1/2) integral_0^x [f_ext(s,r+x-s)-f_ext(s,r-x+s)] ds.

Therefore `||Kf||_X<=3 sup|f|_H`. For continuous g with `g(0)=0`,
`(Lg)_x=-1_(x>r)g(x-r)`, `(Lg)_r=+1_(x>r)g(x-r)`, and
`||Lg||_X<=sup|g|_H`. These same constants apply to differences. No
derivative of actual F or G occurs in these estimates. Qualitative continuity
comes from the retained smooth construction, not from a numerical C4 bound.

The Dirichlet-reflected f may jump at `r=0`; the moving crossing has measure
zero in each line integral and still gives `Kf` in C1. More quantitatively,
for `delta=|x-x'|+|r-r'|`, `M_f=sup|f|_H`, and a uniform physical-source
spatial modulus `omega_f`, either derivative obeys

    |D Kf(x,r)-D Kf(x',r')| <= 3 omega_f(delta)+2 M_f delta.

The crossing interval has length at most `delta`; elsewhere J is fixed and
isometric. For L, continuity of its first derivatives across `x=r` requires
`g(0)=0`. Keeping the retained flat zero wedge is a sufficient stronger rule.

Given a known continuous source approximation Ftilde, full-cell interval
enclosures directly certify `sup|F-Ftilde|_H<=epsilon_F`. Hence
`||K(F-Ftilde)||_X<=3 epsilon_F`. The analogous shell cost is `epsilon_G`.
This is enough for a specified finite tolerance even when the source error
has a nonshrinking uncertainty floor. Actual-source derivative estimates
are **not essential** for this route.

There is a separate issue of arbitrarily vanishing errors: fixed C0 caps
give no effective modulus for actual F (bounded functions can oscillate
arbitrarily rapidly). A comparator only gives
`omega_F(h)<=2 epsilon_F+omega_Ftilde(h)`. Qualitative smoothness alone
does not give a certified mesh size. To push below the existing uncertainty
floor one needs tighter actual-field/comparator enclosures or an effective
continuity modulus; a quantified C4 estimate is one sufficient option, not
the only option and not a prerequisite for the finite target below.

## A continuous numerical candidate without differentiating actual F

The following is a concrete implementation specification, not implemented
coefficients in this note.

1. On a rational mesh, construct C1 piecewise-polynomial U0tilde and continuous
   piecewise-polynomial V0tilde approximating the **same** selected initial
   data. Preserve the flat corner and odd/even reflection. Initial certificates
   must bound values of U0, its ordinary r derivative, and V0 on full cells.
2. Construct continuous piecewise-polynomial Atilde, Btilde, Dtilde, Mtilde,
   Ftilde, Gtilde. Preserve all 23 amplitudes, P_N for shell values, and the
   zero wedge/corner. Shared face values must agree; unrelated selections
   from interval boxes are not continuous approximants. All errors use the
   exact angular Gram matrix and H-operator norms, not unweighted modal norms.
3. Define `btilde=Ufree_tilde+K Ftilde+L Gtilde` and
   `V=(I+Ttilde)btilde`. Integrate the piecewise polynomials exactly over
   characteristic polygons and shell intervals, splitting at every cell
   boundary and reflected branch. Products of polynomial coefficient fields
   with the C1 btilde derivatives are again continuous piecewise polynomials.
   This gives a finite C1 piecewise-polynomial candidate on D. Its coefficients
   may be rational, and the construction retains every boundary contribution.
4. Compute `Rtilde=V-btilde-Ttilde V` by the same operations. Bound its value,
   x derivative and r derivative over complete polynomial cells, for example
   using exact Bernstein bounds and the registered Gram matrix. Measure
   initial, Dirichlet and Robin defects explicitly. Only then apply the true
   residual transfer estimate below.

Exact integration avoids a quadrature modulus entirely for the approximants.
Alternatively, a validated quadrature of the explicit approximants may use
their computable polynomial derivatives/moduli. For a line integral, mesh
diameter h gives error at most `length*omega_Ftilde(h)`; split reflected
crossings first. The same endpoint formulas control K's first derivatives.
One must prove these derivatives belong to the **same continuous reconstructed
candidate**: independently tabulating an approximate value and two endpoint
derivative formulas is not a C1 certificate. Add all quadrature, geometry of
integration cells, modal conversion, rounding and reconstruction errors.

This uses no new estimate for derivatives of actual F. It does require an
actual finite interpolant and verified full-cell discrepancies; neither has
been produced here. Existing global source caps alone do not provide an
accurate initial profile or an evaluated b/Tb.

## Complete error ledger

Let the initial discrepancies in H be `epsilon_0,epsilon_r,epsilon_v`.
The d'Alembert isometry gives

    epsilon_free <= max(epsilon_0/3+epsilon_v,epsilon_r+epsilon_v),
    epsilon_b <= epsilon_free+3 epsilon_F+epsilon_G+chi_b,

where `chi_b` is any remaining certified numerical X error in btilde.
If coefficient discrepancies in H-operator norm are
`delta_A,delta_B,delta_D,delta_M`, then

    tau := 3(delta_A+delta_B+3 delta_D+delta_M) >= ||T-Ttilde||,
    ||Ttilde|| <= q+tau.

The shell trace costs `sup|V|<=3||V||_X`, which explains the last term.
Angular connection, potential, scalar coupling and all Robin blocks belong
in these operator discrepancies. With `||b-btilde||<=epsilon_b`, and a final
construction error `e_alg` relative to `(I+Ttilde)btilde`,

    ||U-V||_X <= q^2 B+(1+q)epsilon_b+tau||btilde||_X+e_alg.       (A)

This is the truncation-plus-approximation route. Independently, a certified
approximate integral residual gives

    ||R_true||_X <= ||V-btilde-Ttilde V||_X+epsilon_b+tau||V||_X,
    ||U-V||_X <= ||R_true||_X/(1-q).                              (B)

Every numerical error in the approximate residual itself belongs in its
first term. Do not add the Neumann tail from (A) again to (B). An equivalent
direct residual enclosure may bound the full true integrands
`A V_x+B V_r+D0 V+F` and `M V|b+G` against continuous comparators, with
costs three and one respectively. This keeps all sources and boundary terms.

## Nonvacuous quantitative target, still prospective

The new actual complete-cover certificate gives
`F_H<=17387513156573`. The inherited full shell join gives
`G_H<=3851400`, including GN. Thus

    3 F_H+G_H = 52162543321119.

Even continuous source approximants with only the conservative magnitude
guarantees `||Ftilde||<=F_H`, `||Gtilde||<=G_H` have
`3 epsilon_F+epsilon_G<=2(3 F_H+G_H)=104325086642238`.
This is an explicit approximation error allowance, not deletion of KF or LG.
It is far below a response budget of order `1e29`; source derivatives are
therefore not the limiting mathematical requirement for that coarse target.
This statement does not imply that a C0 global cap resolves a small TT response.

A concrete sufficient target allocation is

    epsilon_b <= 4e28,
    tau <= 4e-5,
    ||btilde||_X <= 1.04e33,
    e_alg <= 1e28.

The btilde bound follows from `||b||<= (1+q)B` and the epsilon_b target.
The total numerical part of (A) is then at most `9.168436e28<1e29`.
For example, initial/free error at most `3e28` leaves ample room for the
displayed conservative source allowance and additional b construction error.
These are targets to be attained, not already certified errors. Existing
background/initial interval widths must be audited before claiming they are
attainable without refinement.

Using the simpler total approximation target `1e29`, the full error is
`||U-V||_X<1.0458131743e29`. The sphere-mean functional has norm one and
the X value weight is three, so

    |mean_w(U)-mean_w(V)| < 3.1374395229e29.

The already certified evolved mean margins `+/-1.85e30` would then remain
separated from zero by more than `1.53625604771e30` in V on their two
rectangles. No angular point-evaluation factor belongs in this mean bound.
This supplies a nonvacuous acceptance target relative to existing features;
it is not a computed candidate, a new evolved sign proof, or a TT sign result.

## Why the coarse sign target does not resolve the forced increment

Before implementing that coarse total-response target, the scientifically
useful distinction is the forced increment. Put

    S=KF+LG,
    U_hom=(I-T)^(-1) U_free,
    W=(I-T)^(-1) S,
    U=U_hom+W.

The same full operator T, including its coupled Robin lift, acts in both
resolvents. U_hom is not the free wave. It is an auxiliary reduced-operator
comparison using the same registered second initial data and zero affine
bulk/shell sources. W has zero initial data and retains all actual affine
bulk/shell sources. Because the registered initial prescription and quadratic
subsidiary constraints themselves depend on the actual first pair, neither
summand is claimed to solve the complete physical quadratic constraint and
gauge-completion system separately.

Nevertheless the operator decomposition and its estimate are exact:

    ||W||_X <= (3 F_H+G_H)/(1-q)
             =52162543321119000000/997891 <5.228e13,
    |mean_w(W)| <=3||W||_X <1.5684e14.

The exact ratio of the mean upper bound to the existing `1.85e30` evolved
mean margin is less than `8.477e-17`. Hence those signs persist also in
U_hom. Their existing large margins do **not** resolve the forced quadratic
increment. This conclusion uses a full-source norm comparison, not a dropped
boundary term, and it does not identify an invariant backreaction observable.
It does not exclude local importance where U_hom is small, near a zero or in
a weak modal sector; an upper bound alone does not prove W is nonzero.

A `1e29` X error certificate can certify the already established total mean
signs but cannot certify resolving W at the much smaller displayed scale.
This does not prove the candidate's actual numerical error is large or that
resolving W is impossible. The more informative
next numerical target is therefore the zero-initial forced increment W,
computed separately from U_hom. Its two-term candidate `S+TS` has exact tail

    ||W-(S+TS)||_X <=q^2(3 F_H+G_H)/(1-q) <2.326e8.

For illustration, an eventual total X error below `1e10` would require a
certified approximation to `S+TS` substantially below `1e10`; the displayed
tail would then be small enough. This is an absolute target relative to an
upper bound only. It does not guarantee a nonzero or signed W: an actual
candidate feature must exceed its certified error before any such claim.
The coarse source allowance `2(3F_H+G_H)` is too wide for this smaller target.
Tighter actual-source/comparator errors are needed, with their uncertainty
floor checked. The C0-to-X method still applies and still does not logically
require derivatives of actual F.

## Work still required

The finite source-derivative objection is resolved at the level of method:
C0-to-X propagation and continuous approximants suffice. Still missing are
the coherent full-domain interpolants, actual initial/operator error bounds
at the allocated scale, characteristic cell integrations, and a certified
complete residual or (A) ledger for their concrete coefficients. Exact row
and angular normalization checks must accompany the implementation. Further
input refinement may be needed if the existing interval uncertainty floor
exceeds the target; nothing here proves that it does or does not.

No current capsule, producer, manuscript or publication artifact is modified.
The small `check_budget.py` script checks the arithmetic only and authenticates
its retained input files. It does not replay the expensive scientific inputs.
