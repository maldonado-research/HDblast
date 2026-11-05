# Input boundary and what each criterion establishes

The complete physical shell fields include the actual shell embedding and
their identification with a fixed reference shell. A coordinate restriction
to an unperturbed shell is insufficient. Ambient physical injectivity after
composition requires that shell identification itself to be an embedding.
Parameter derivatives hold the actual proper-clock and transported labels
fixed. A change of anchor or a parameter-dependent coordinate identification
requires its own terms and is not implicit in these notes.

| Input/result | Role | Still required |
| --- | --- | --- |
| Regular full shell fields and exact initial normal/label data | Define the local intrinsic geodesic/Jacobi ODE | Retained PDE, preparation, gauge and field regularity assumptions; evaluated common coefficient bounds |
| Fixed coefficient neighborhood and strict first-exit guard | Continue all candidate flows over the prescribed finite duration | Numerical neighborhood/state constants and a parameter interval; no injectivity follows from phase-space uniqueness alone |
| Projected parameter/Jacobi bounds and the metric gap | Quantify spacelike cuts and reference-slab retention | Actual uniform full-pullback bounds, including clock, embedding and chart factors |
| Normal proper geodesic structure, positive cuts, timelike reference time vector | Global angular graphs and strict graph-time monotonicity | The whole candidate flow and these hypotheses on it; no generic-map shortcut |
| Small ambient derivative residual | Alternative quantitative full-cylinder embedding test | A genuine uniform bound on the same product cylinder, with the sphere treated as a manifold |

The retained [nonlinear observer budget][observer] assumed a common regular,
spacelike, noncrossing tube. The present notes replace parts of that assumed
domain with explicit sufficient geometric/ODE tests. They do not establish
that the actual model passes those tests. The prior conditional uniform first
and second variation estimates still have unevaluated constants; a first
tangent or an expansion coefficient at zero does not supply their values.

The exact-source copies and hashes are recorded in [INPUT_PROVENANCE.json](INPUT_PROVENANCE.json).
Authentication is an integrity check, not a cryptographic signature, a new
proof of a retained field-equation theorem, or an evaluated nonlinear bound.
Private handoff records and editorial history are excluded from the public
source inventory.

[observer]: https://github.com/maldonado-research/HDblast/blob/d62bdda754507a29c198213946906cea6d9a9dbf/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/observer/OBSERVER_FLOW_BUDGET.md
