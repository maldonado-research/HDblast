# Sources and applicability

The retained mathematical sources below are linked at immutable repository
commit `c4b32c03243e2bcc04b748570e7c1bd5861703fa`. Their exact input bytes,
lengths and Git object identities are recorded in [SOURCE_PINS.json](SOURCE_PINS.json).
They are conditional research checkpoints, with no implication of external
peer review. Earlier work supplies the qualitative evolution theorem; this
checkpoint evaluates selected inputs to its application.

## Retained model and analytic sources

| Source | Role in this checkpoint |
| --- | --- |
| [Registered exact-family binding](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/REGISTERED_EXACT_FAMILY_BINDING.md), §§1–4 | Fixed background and first packet, half-action canonical normalization, selected weighted curve, affine temporal lapse and zero initial shift. Its rescaling of the generic companion's momenta is essential. |
| [Initial harmonic completion and subsidiary system](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md), §§1–5 | Exact ADM completion, relative-harmonic reduction, varying reference contraction, trace-reversed Cauchy constraints, scalar/Codazzi cancellation and conservative parameter derivative budget. Its doubled canonical normalization must be reconciled with the binding. |
| [Physical second Cauchy jet](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_WEIGHTED_REGULARITY_20261106/PHYSICAL_CAUCHY_SECOND_JET.md) | Geometric `H14/H13` input, fixed positive norm and physical density-to-velocity conversion. It is an initial derivative-at-zero estimate, not a finite-amplitude bound. |
| [Source jets](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_WEIGHTED_REGULARITY_20261106/SOURCE_JETS.md) and [finite Taylor engine](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_WEIGHTED_REGULARITY_20261106/derive_source_jets.py) | Coupled background and known first-packet derivative majorants. The present initial-gauge calculation uses those inputs without repeating the weighted inverse construction. |
| [Actual potential and background boxes](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/ACTUAL_POTENTIAL_JET_BOUNDS.json) | Exact potential coefficients and the actual wider-collar background enclosure used by the new finite derivative recurrence. |
| [Quantified boundary flow](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md), [metric/energy patch](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md) and [coordinate dictionary](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/COORDINATE_DICTIONARY.md) | Literal proper-`y` component-root flow, background subtraction, finite chart jets, physical inner patch and coercivity. The outer holomorphic flow box does not imply physical metric-patch membership. |
| [Full angular linear matrix and energy](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_ANGULAR_UNIFORM_CONTACT_20261015/ANGULAR_UNIFORM_LINEAR_CONTACT.md), §§2–5 | Four Dirichlet/twelve coupled Robin rows, full-bundle conjugation, homogeneous energy and explicit boundary-lifting cost. |
| [Finite-regularity shell evolution](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md) and [variable mixed-wave proof](https://github.com/maldonado-research/HDblast/blob/c4b32c03243e2bcc04b748570e7c1bd5861703fa/research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md) | Retained conditional `X14` evolution, compatible iteration and `C2` parameter dependence into `X10`. Their application constants and actual-family bootstrap are not supplied by the present base-energy inequality alone. |

The complete input roster also identifies the retained variable-wave audit
and full finite-angular second-order source. The reduced bulk equations
retain the spatially varying reference connection. No shell-frozen bulk
connection or scalar-only replacement of the sixteen-field matrix is used.

## Focused primary method checks

1. Jason Metcalfe and Ann Stewart, [*Almost global existence for quasilinear
   wave equations in waveguides with Neumann boundary conditions*,
   arXiv:math/0509326v1](https://arxiv.org/html/math/0509326v1),
   14 September 2005, §5 (5.1)–(5.4). The variable-coefficient estimate
   charges an `L1_t L∞_x` coefficient derivative budget and imposes a
   normal/tangential principal-tensor compatibility condition in (5.2).
   This supports checking the shell cross-trace identity at every frozen
   iterate and retaining all commutators. Its scalar waveguide result
   supplies neither this coupled gauge system nor its numerical constants.

2. H.-O. Kreiss, O. Reula, O. Sarbach and J. Winicour,
   [*The Einstein equations, boundaries and integration by parts*,
   arXiv:0707.4188v1](https://arxiv.org/html/0707.4188v1),
   27 July 2007, Theorem1 and Appendix1 (15)–(20). The explicit energies
   show lower-order stability and tangential-then-normal derivative
   recovery. Their boundary-input estimate requires a strictly positive
   Sommerfeld time coefficient. Setting it to zero does not produce an
   arbitrary `L2` boundary-input estimate for the present pure Neumann
   rows; the explicit lifting cost remains necessary.

3. Kirk D. Blazek, Christiaan Stolk and William W. Symes,
   [*A Mathematical Framework for Inverse Wave Problems in Heterogeneous
   Media*, arXiv:1211.4656v1](https://arxiv.org/html/1211.4656v1),
   20 November 2012, Theorem9 and following remark. The autonomous
   symmetric-hyperbolic coefficient-to-solution/trace map loses source
   and output regularity in its differentiability statement; the paper
   identifies those losses as sharp in its class. This supports an
   explicit parameter derivative budget. It does not obstruct the retained
   lower-space `C2` family or establish a `C3` theorem for this model.

## Closest previously retained comparisons

4. Grigorios Fournodavlos and Jacques Smulevici,
   [*The initial boundary value problem in General Relativity: the umbilic
   case*, arXiv:2104.08851v1](https://arxiv.org/html/2104.08851v1),
   18 April 2021, Theorem1 and AppendixA. The wave-gauge/umbilic strategy
   concerns smooth four-dimensional vacuum with constant umbilicity.
   Its theorem does not literally cover the present five-dimensional
   scalar shell with field-dependent tension and relative-harmonic gauge.

5. Tatsien Li, Xing Lu and Bopeng Rao,
   [*Exact boundary controllability and exact boundary synchronization for
   a coupled system of wave equations with coupled Robin boundary controls*,
   arXiv:2202.00760v1](https://arxiv.org/html/2202.00760v1),
   1 February 2022, Theorem4.1 (4.3)–(4.7),
   [DOI:10.1051/cocv/2020047](https://doi.org/10.1051/cocv/2020047).
   Their scalar exponential Robin-to-Neumann change explains a familiar
   method. The coupled weak/control theory's real symmetrizability
   assumptions do not directly cover the arbitrary nonlinear twelve-row
   boundary vector used here; the present homogeneous point-map energy
   is checked separately.

6. Zhongshan An and Michael T. Anderson,
   [*Well-posed geometric boundary data in General Relativity, II: twisted
   Dirichlet boundary data*, arXiv:2507.15567v3](https://arxiv.org/html/2507.15567v3),
   revised 30 May 2026, Theorem1.1. The vacuum theorem in arbitrary spatial
   dimension at least two uses smooth tame/Nash–Moser theory and
   derivative-free conformal-class/twisted-volume data. Those data differ
   from the present Israel/scalar Robin laws and do not supply its finite
   Sobolev constants or a proper-time observation domain.

7. Dionysios Anninos, Raúl Arias and Damián A. Galante,
   [*Umbilic Observatories*, arXiv:2610.40783v1](https://arxiv.org/html/2610.40783v1),
   29 October 2026. Its linear conformal/Dirichlet observatory analysis
   and selected nonlinear boundary solutions concern a different vacuum
   boundary problem. They do not establish the general variable-tension
   Einstein–scalar shell evolution considered here.

The focused primary passages were read on 6 November 2026; the closest
comparisons retain their previously checked exact versions. No cited
external theorem supplies the present model-specific numerical application.
The applicability differences imply neither a novelty claim nor an
ill-posedness claim. None of these references supplies the still-missing
straightened initial Sobolev budget, actual-family high-order bootstrap,
finite amplitude radius or proper-time/area observable remainder.
