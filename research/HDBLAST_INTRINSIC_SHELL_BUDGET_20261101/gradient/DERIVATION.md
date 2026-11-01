# Conditional intrinsic shell scalar-gradient second coefficient

1 November 2026. Portable read-only-source numerical component. The coefficient
below is fixed by the actual first field. Its computation does not select or
evaluate a nonlinear weighted completion or its scalar second monopole.

The formula below already appears as equation (3) of the 16 October source
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/shell-invariant/SHELL_INTRINSIC_SCALAR.md`.
This component certifies the frozen current enclosures and joins the actual
first field to the selected exact-family and intrinsic-clock premises; it does
not claim a new formula.

## Definition and exact second derivative

Let i_e be the physical shell embedding, h_e=i_e* g_e, and Phi_e=i_e* phi_e.
Define the Lorentzian intrinsic scalar

    Q_e = h_e^{ab} (D_a Phi_e)(D_b Phi_e).

This Q is distinct from the compact bump called Q in the selected source.
The background shell is z=0 and

    h_0 = rho_b^2 [-dt^2 + cosh(t)^2 kappa_S3], Phi_0=phi_b.

Phi_0 is constant on the entire shell. With parameter DERIVATIVES

    Phi_e=Phi_0+e eta+e^2 zeta/2+o(e^2),

we have D Phi_0=0, hence identically on the shell

    Q_0=Q_1=0,
    Q_2=2 h_0^{ab} (D_a eta)(D_b eta).

Every term involving h_1, h_2 or zeta multiplies D Phi_0. No scalar second
data or second solution enters this identity. In the common fixed-shell chart
selected in the exact-family source, eta=psi(t,0)Y. The actual scalar oracle
uses eta=v+p f. Its boundary identities impose f=0 and all its tangential time
derivatives at r=0. Thus the retained scalar `s` and ordinary-coordinate-time
derivative `st` at exact r=0 are precisely psi and psi_t. The scalar is not a
normal derivative and `st` is not scaled x differentiation: t=ell*x and the
oracle has already divided by ell.

For Y=n_0^2-n_1^2 on the unit S3, the probability-normalized angular moments are

    <Y^2>=1/6, Delta Y=-8Y, <|D Y|^2>=8<Y^2>=4/3.

For example <n_i^4>=1/8 and <n_i^2 n_j^2>=1/24 on S3 give <Y^2>=1/6;
integration by parts gives the gradient moment. Therefore

    <Q_2>(t) = [-psi_t^2+8 psi^2/cosh(t)^2]/(3 rho_b^2),
    <K_2>(t) = psi_t^2/(3 rho_b^2),
    <S_2>(t) = 8 psi^2/(3 rho_b^2 cosh(t)^2),
    <Q_2>=<S_2>-<K_2>.

K_2 and S_2 here name the positive temporal and spatial contributions to this
identity; K_2 is not an extrinsic-curvature variation. Ratios below are S_2/K_2.

The second coefficient is invariant under any smooth parameter-dependent
relabeling of the physical shell with identity background map: the pullback
correction contains only Lie derivatives of Q_0 and Q_1, which vanish.
The same observation removes first- and second-order event-clock shifts.
A smoothly varying sphere area form or normalization cannot affect the
second coefficient of the angular mean because it multiplies Q_0 or Q_1.
These statements concern perturbative relabelings around these background
spheres; no alternative background slicing is selected.

## Conditional physical identification and scope

The selected source is the closed capsule
`HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031`. Its DERIVATION.md specifies the
same Einstein-scalar model, the fixed shell chart, and derivative convention;
section 4 identifies the actual full first master/gauge pair with the exact
family first derivative using the complete harmonic initial one-jet and
mixed boundary uniqueness. Section 6 supplies C2 parameter dependence and
the physical-side shell-trace topology. Those analytic hypotheses and the
actual-solution premises remain conditional inputs. Our intervals do not
reprove the analytic theorems.

The first scalar jet oracle is the authenticated
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py`.
It encloses actual fields and derivatives with the complete inherited
differentiated C3 comparison errors, including master and gauge errors.
It is called at physical shell r=z/ell=0 EXACTLY. Neither the old bulk event
r=1/8 nor an independently chosen comparison solution is used.

Any admissible exact completion with the SAME first physical shell pullback
has this same Q_2, even when the nonlinear correction weights, exponents, or
the scalar second monopole have not been numerically selected. This proves
only the intrinsic scalar-gradient contribution to a possible shell Gauss
scalar-curvature budget. In the conventions of the selected residuals and
one-copy laws, the contracted Gauss identity would read

    R[h] = Q + 2 V(Phi) + sigma(Phi)^2/3 - sigma'(Phi)^2/4.

The remaining second derivative includes scalar second data through the last
three terms. We do not enclose it here. There is no conclusion about full
shell scalar curvature, shell expansion, radiation, force, or a nonlinear
finite-amplitude model solution. The negative angular mean is not a claim
that Q_2 is negative at every angle. No effective parameter radius is supplied.

## Separated pilot and frozen campaign

The original exploratory pilot ran three explicitly diagnostic cells and wrote
the retained `PILOT_ONLY.json`. They were the fixed target x=5/4, the preselected nearby
interval [31/25,63/50], and its central smaller interval. They all used the
actual error-enclosed oracle; no scan over other favorable target points was
performed. The broad interval already excluded zero.

After that pilot, `FREEZE.json` fixed the SAME target x=5/4 and the interval
[31/25,63/50] with 8 equal ascending cells, followed by one target call.
The pilot hash and source-pins hash are in the freeze; its own SHA-256 and the
source-pins SHA-256 are hard-coded in the final producer. There were exactly
9 actual calls. No post-freeze refinement occurred. Every complete metric
(42 entries), scalar (6 entries), and background output is retained per call.

The physical scales are ell=10^-6 and t=ell*x. The whole interval has x>1,
the initial packet's first shell-contact value, and x<3/2, the common
selected exact-family domain endpoint. Both endpoints are included.

The producer verifies both closed source capsules and every member before
any scientific import and again after all calls. The source tree is
read-only; `SOURCE_PINS.json` records the actual checkout commit and tree.
All interval decisions use exact Fraction rational endpoints with outward
rounding on denominator 2^112. Exponential and square-root enclosures are
inherited authenticated rational algorithms. Floating point appears only
in optional human-readable displays outside the final producer.

## Result

At the fixed shell target x=5/4, an outward INTEGER enclosure is

    -803004132630938124 <= <Q_2> <= -659922099878760930 < 0.

Over every point of the continuous interval x in [31/25,63/50], the hull of
the eight whole-cell enclosures gives

    -1203216564043213024 <= <Q_2> <= -408804027527669570 < 0.

The exact rational temporal, spatial, and ratio enclosures are retained in
`REPORT_NORMAL.json`. Approximate displays at the target are

    <K_2> in [6.599220998790298e17, 8.030041326311590e17],
    <S_2> in [220875.685715, 268783.029340],
    <S_2>/<K_2> in [2.75061e-13, 4.07296e-13].

The displayed S interval is rounded outward at six decimal places; K and
ratio scientific values above are illustrative, while REPORT_NORMAL.json is
the exact certificate. The continuous hull of the ratio has upper endpoint
less than 7.50780e-13. Temporal domination is very strong even after all
inherited first-field errors and whole-cell interval variation are charged.

The direct formula and the separated S_2-K_2 formula each exclude zero in
all nine calls. The independent `audit_aggregation.py` imports no scientific
or producer code: it reauthenticates the source members, verifies the fixed
coverage and full retained output inventories, and recomputes all squares,
products, divisions, ratios, and hulls with a separate exact-rational
interval implementation. `AGGREGATION_NORMAL.json` passes. That audit
does not independently reprove the inherited PDE or oracle error theorem.

The portable CLI, exact source closure, provenance normalization, normal and
optimized byte identity, strict negative controls and relocation positives are
documented in README.md. FROZEN_REFERENCE.json retains every original
scientific output unchanged. Producer and auditor use separate implementations
and authenticate their dependencies before work and after completion.
