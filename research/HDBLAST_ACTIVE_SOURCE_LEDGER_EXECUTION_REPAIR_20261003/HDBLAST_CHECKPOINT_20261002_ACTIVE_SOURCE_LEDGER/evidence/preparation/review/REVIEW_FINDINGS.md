# Active-source implementation review findings

These findings precede any new physical source/array/trajectory evaluation.
This document is internal review, not independent external peer review.

| Finding | Evidence and consequence | Current disposition |
| --- | --- | --- |
| Primary moment coefficient was twice the action value | The first provisional method used `E=sum(mu*k*Re(U))`, but `Rmode=2E+...`, `Pmode=2E/3+...` and invariant `2E-J`. The unrestricted direct operators require `E`, `E/3` and `E-J`. | Primary acknowledged and corrected the staged code and method before freeze. |
| Active source exposes different momentum targets | The old producer integrates pointwise rationalized finite-band baseline and full contact operators with the saved momentum rule. Analytic finite-K primitives use a continuous momentum integral. This distinction disappears in the prior source-free local contacts. | Root requires explicit signed `E_momentum`, preserved raw gaps and prospective matched-target comparisons. |
| Independent source scalar interface | `source_jet` creates `(6,)+z.shape`; scalar input makes `out[n]` a scalar and masked assignment invalid. The current planned vector source interface avoids it, but a scalar endpoint caller would fail. | Independent now normalizes to an at-least-one-dimensional array and restores the original shape. No study-source call was needed to correct this interface. |
| Generic pair-valued dense WKB route exceeds prospective resource target | Independent initially reported approximately 4223 seconds projected for the complete route. A naive full fine-grid GL24-by-GL24 kernel requires about 302 MiB, exceeding the total 256 MiB memory allowance before inputs and output. | Independent replaced repeated per-source Pair evaluation with symbolically verified moment regrouping and bounded geometry caches. The final fabricated all-case/three-control/MP80-and-100 numerical core completes at 473.506 seconds and 110420 KiB, including midpoint/end flow projections and output serialization. Its receipt clearly distinguishes executed numerical-core bytes from subsequent separately tested inexpensive bootstrap/schema guards; no actual study source/input was evaluated. |
| Arithmetic precision scope is narrower than source/phase precision | Both methods currently evaluate sources, complex phases and forced mode updates in native binary80, with MP80/100 only for weighting/contact/reduction arithmetic. | Staged methods now explicitly state this scope. Cross-route numerical comparisons use the prospectively fixed empirical `2e-7` scale; exact algebra/MP reduction/serialization checks retain `1e-12`. |
| The computed decomposition needs a reconstruction remainder | The three-term expression `D_cont=E_flow+E_momentum+E_operator` assumes exact forced propagation and exact contact integration. Fixed GL and native arithmetic instead contribute `E_reconstruction=delta(E-J)+delta(C_R^A)+M_A*J_source-Q_contact`. | Theory, primary and root agreed to retain the fourth term explicitly. Its empirical gate is `2e-7`; arithmetic closure including all four terms remains `1e-12`. |
| Initial and midpoint faults could cancel in endpoint differences | The first full primary prototype reported only the endpoint difference of stored/direct operator residuals, and omitted separate baseline and midpoint momentum gates. | Primary now reports and gates all three fixed-knot operator/contact/baseline audits, source/forcing-jet gaps, midpoint/end signed-flow projections and separate triangle envelopes. It also gates both reconstruction ingredients and compares quadrature controls at both MP contexts. |
| Direct independent CLI lacked the complete public-freeze guard | The current independent CLI initially checked the registration hash and a forty-character commit format, without authenticating the freeze receipt or member/header inventories. Resolving paths before checking symlinks also lost symlink evidence. | Independent added a stdlib-only entry wrapper, exact entry/core/helper source pins, the shared complete frozen-byte/capsule guard before numerical imports, pinned runtime, a shared complete-route timer and post-run frozen verification. The numerical core now requires wrapper authorization and validates all nineteen array shapes/dtypes/finiteness/labels before source execution. |

The current numerical methods and resource scope pass the internal design
review. The primary full fabricated route completes at 219.155 seconds and
183944 KiB; the independent completes at 473.506 seconds and 110420 KiB.
Stable-kernel normal/optimized high-phase review controls pin the current
primary and independent numerical-kernel bytes. The global schemas and exact-byte copies now match the reviewed design. An
independently authored full twelve-case fixture passes 43 normal/optimized
policy and mutation checks, and separate review policy probes pass 30 checks
in both modes. Previously accepted float/boolean metadata are now rejected;
wrapper-added provenance fields are included in the schema. The joint complete
control-universe comparison is explicit. The public remote freeze, new
scientific outcomes and standalone physical replay remain pending. No actual
physical result or inferred outcome appears here.
