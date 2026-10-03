# Independent action and junction review

Prepared 2026-10-02 UTC in staging outside the checkout. Scope: mathematical
review and exact algebraic controls for the next coupled-source specification.
No new physical evolution, external publication, discovery claim, or certified
bulk energy closure is asserted. The inherited checkpoint's original curved
pressure FAIL and separate K384 follow-up PASS are unchanged.

## Findings

The inherited three junction signs and their factors agree with the doubled
bulk action, two GHY boundary terms, one shell action, and outward normal
`n = exp(-B) partial_z` on the interior `z < 0`. A paired quantum source can
enter those junctions, but the present B1 explicit boundary adapter is not a
complete implementation of the new semiclassical system. Its source interface
and initial-value formulation must change or be replaced by a justified,
prospectively specified order reduction.

The local action review passes: lapse-retaining FRW variation independently
reproduces the displayed energy, pressure and scalar current. There are 21
passing exact checks and nine rejected source mutations. The calculation also
confirms two limitations that acceptance tests must address: conserved paired
double counting can evade a Ward test, and a pressure error can evade an energy
exchange test at an instant with H=0.

## Pinned inputs and action normalization

Read-only repository: `/workspace/HDblast`, HEAD
`e17a01b` (full hash recorded in `REVIEW_PROVENANCE.json`). Relevant inputs:

* Sep22 `matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md`: explicit full doubled
  action, normal, source normalization, signed flux and Weyl identities.
* Sep30 `B1_derived_chi_source/evolve_b1.py`, `shell_matter.py`, `static_w.py`,
  `b1_symbolic.py`: actual bulk model, boundary adapter, unit conversion and
  particle-source conventions.
* Oct02 `COMMON_ACTION_SMOOTH_FRW.md`: positive-reference stress/current,
  selected finite matching, and pressure derivatives.

With signature (-++++), the doubled bulk action is

    sum_a { 1/(2 kappa5^2) int_Ma sqrt(-g)[R-(grad phi)^2-2U]
            + 1/kappa5^2 int_shell sqrt(-h) K_a }
    + int_shell sqrt(-h)[-sigma(phi)/kappa5^2 + L_m].

The scalar phi is dimensionless, and physical tension is lambda=sigma/kappa5^2.
There is one shell scalar action. Field multiplicity multiplies its matched
stress/current; the mirror doubling does not supply an extra matter factor.
For the registered model, W=1-phi+phi^3/3,
U=W_phi^2/2-2W^2/3, and
sigma=2W+delta(1+c phi+d phi^2/2), with delta=.001, d=0 and
c=2/1.0357712571566784-4/3. Selecting d != 0 changes the tension model.
The actual B1 registered runs retain the distinct A1 tuned model: delta=.1,
d=-3.106933495673783. The next action must explicitly select the original
registered shell or that tuned extension; sharing a solver and a bulk potential
does not make their shell tensions identical.

Define T_munu=-2 delta_Gamma/(sqrt(-h) delta_h^munu) and
J=-delta_Gamma/(sqrt(-h) delta_phi). Variation gives

    K_munu-K h_munu = (kappa5^2/2)(-lambda h_munu+T_munu),
    -2 n.phi/kappa5^2-sigma_phi/kappa5^2-J = 0.

Writing k_s=n.A, k_0=n.B, w=n.phi, rho=T_tautau, and p=T_ii/h_ii:

    k_s = (sigma+kappa5^2 rho)/6,
    k_0 = [sigma-kappa5^2(2rho+3p)]/6,
    w   = -(sigma_phi+kappa5^2 J)/2.

All sources here are physical shell quantities. The B1 adapter instead names
its already weighted outputs `rho`, `pr`, and `J0`, which represent kappa5^2
times physical sources. Multiplying them by kappa5^2 a second time is incorrect.
Its H0-unit gas conversion is `b/r_b`, with b=kappa5^2 H0^3 and H0=1/r_b;
thus kappa5^2 rho_phys=b rho_hat/r_b, and the same conversion applies to J.
Benchmark profiles expressed in arbitrary mass units require an explicit scale
map before use in this adapter.

## Local source variation and differential order

For a common local action

    Gamma_local = int sqrt(-h)[-V(phi)+F(phi)R+alpha R^2],

with constant alpha, the covariant metric variation is

    T_local = -V h - 2[F G+(h box-D D)F]
              -2alpha[2R R_munu-(R^2/2)h+2(h box-D D)R],
    J_local = V_phi-F_phi R.

This includes local terms introduced as independent renormalized brane
couplings or moved consistently out of a matched quantum effective action.
It does not authorize adding another copy of the benchmark's subtraction
action to its already matched renormalized expectation value.

The independent check starts from `N a^3[-V+FR+alpha R^2]` with
R=6[a''/(aN^2)+a'^2/(a^2N^2)-a'N'/(aN^3)]. It varies N and a before setting
N=1. With proper-time dots, H=adot/a, R=6(Hdot+2H^2), it obtains

    rho_local = V-6F H^2-6H Fdot
                 +alpha(36Hdot^2-216H^2 Hdot-72H Hddot),
    p_local   = -V+2F(2Hdot+3H^2)+2Fddot+4H Fdot
                 +alpha(108Hdot^2+216H^2 Hdot+144H Hddot+24H'''),
    J_local   = V_phi-F_phi R.

Thus a nonconstant F couples the metric and scalar accelerations: the pressure
contains 2F_phi phiddot, and the scalar current contains
-6F_phi addot/a. This is an implicit coupled boundary system, even though
its principal derivatives are second order when alpha=0 and no additional
derivative-dependent quantum terms are retained. A constant F still contributes
an induced Einstein term, including metric acceleration in the temporal
junction. Alpha R^2 adds an a''' term to energy and an a'''' term to pressure.
At isolated H=0 the energy coefficient of the highest derivative can vanish,
so that instant does not justify a global order reduction.

The benchmark matching alpha(r)=0 specifies the extra local finite coefficient;
it does not remove derivative terms from the source prescription. A direct
check of its reference subtraction integrands gives, in conformal time with
w=sqrt(k^2+a^2r) and c=1+(k^2/3+a^2r)/w^2,

    partial p4 / partial a'''' = c[1/(32a^5w^3)+r/(64a^3w^5)],
    partial p4 / partial x''   = -c/(32a^2w^3),
    partial Q2 / partial a''   = 1/(4a^3w^3)+r/(8aw^5).

These coefficients are generically nonzero and are regular for a>0, r>0,
including physical x=0. They are coefficients of subtraction data, not a
heavy-mass approximation to the physical field. The bare modes and these
terms must be combined in the complete prescription. Dropping the high jets
independently changes it.

The existing B1 `neumann` calls `matter.totals` without H or the required metric
jets. Its state comprises bulk fields and their first coordinate-time
derivatives. It constructs normal ghost data and then computes accelerations
explicitly. A source depending on those same accelerations must instead be
solved implicitly or treated within a justified approximation. `neumann_t`
also requests rho_T, p_T and J_T. Differentiating the full R^2 pressure boundary
condition requires an a^(5) term; a solver that only needs the momentum
constraint should not gratuitously request every higher boundary derivative.
These interface observations are implementation blockers, not a proof that
the physical model has no consistent semiclassical formulation.

## Constraint and energy checks

For the quantum source Q=2S, with S=-delta_Gamma/(sqrt(-h) delta_x),

    J_Q = x_phi Q/2,
    x = m0^2 + G^2(phi_b-phi_star)^2,
    J_Q = G^2(phi_b-phi_star)Q.

The source Q denotes renormalized chi squared here; B1's `Q_hat` denotes decay
power. They have different dimensions and meanings and need distinct names in
the new source interface. The free Gaussian model has no B1 Yukawa decay or
phenomenological production ramp. Reusing `jprod` alongside exact quantum
production would add a separate force and require a separate action/ledger.

The common Ward identity is rhodot+3H(rho+p)=Jv. Substitution into the junctions
gives the exact boundary momentum relation

    M_shell/exp(2B_b) = -kappa5^2 [rhodot+3H(rho+p)-Jv]/2,

where the conformal-gauge momentum constraint is
M=-3A_Tz-3A_T(A_z-B_z)+3A_z B_T-phi_T phi_z.
This validates source compatibility, but evaluating A_Tz from the imposed
differentiated junction can conceal a boundary defect. Acceptance must also
differentiate the recorded bulk velocities independently near the shell.

The signed local flux identity is

    d(lambda+rho)/dtau + 3H(rho+p) = -2wv/kappa5^2
                                 = -2 T_bulk(n,tau).

It does not furnish a global gravitational energy or a closed four-dimensional
evolution. The inherited projected Friedmann and Weyl relations still require
the evolved bulk, including normal scalar data and its time derivative.
Computing a Weyl value by subtracting the same Friedmann source terms defines
a diagnostic; it is not independent evidence of bulk energy closure. A bulk
flux or quasi-local ledger must name its region, surfaces, time flow, boundary
terms and truncation flux, and be independently derived before certification.

The negative controls expose these cases:

* Reverse J, omit F_phi R from J, omit 2Fddot from pressure, or omit 24alpha H'''
  from pressure: the local Ward residual is generically nonzero.
* Double local metric stress with a single current: the Ward residual is Jv.
* Double both stress and current: Ward still passes. The junction coefficients
  change, so coefficient provenance must detect this model change.
* At H=0 a pressure-only error changes neither the exchange term nor its energy
  residual. It changes k_0 by -kappa5^2 delta_p/2 and requires a direct pressure
  or independent trace test.

## Prerequisites before a new coupled run

1. Pin the action, physical scale map, field multiplicity, finite couplings,
   constant positive reference r and quantum state. For every local operator,
   say whether its variation is included in the renormalized quantum source,
   in a separate brane action, or on the metric left-hand side, exactly once.
   The benchmark finite convention is an explicit possible model choice, not
   a determination of the inherited physical couplings.
2. Select full semiclassical dynamics or a controlled order reduction, including
   derivative order, branch/runaway treatment and initial data. An order
   reduction must preserve the common-action source pairing to its claimed
   truncation order; setting derivatives to zero independently does not.
3. Prepare constraint-compatible bulk plus shell data and an ultraviolet
   admissible state/history. The benchmark's exactly static past cannot be
   assumed on the actual de Sitter shell, and archive interpolation does not
   automatically supply the smooth high jets or a suitable quantum state.
4. Derive an energy ledger with independently measured bulk and boundary terms.
   Register refinement and adverse-source checks; include pressure-sensitive
   checks at turning points and static intervals. Retain zero-source and
   wrong-sign controls with the original failed model outcomes.

The action and junction specification can be completed now. A full physical
coupled evolution remains unexecuted and cannot be authorized by algebraic
checks alone. Thermalization, interaction rates and a sustained radiation era
remain further model calculations.

## Reproduction

Run `python independent_junction_checks.py` in this staging directory. It writes
`INDEPENDENT_JUNCTION_CHECKS.json`; all failures raise explicit exceptions.
The script varies the action before testing Ward, verifies the covariant R^2
components, checks junction/chain-rule identities, and rejects nine adverse
mutations. It does not run the repository's physical solvers or edit frozen
inputs. The same run under `python -O` is retained as a receipt to verify that
the checks do not depend on removable Python assertions.
