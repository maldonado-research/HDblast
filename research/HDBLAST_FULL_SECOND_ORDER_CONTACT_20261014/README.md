# Complete second-order coefficients through controlled shell contact

14 October 2026 · conditional Einstein–scalar mathematics · GitHub checkpoint

For the selected copy-even aligned quadrupole, the required **complete second
Taylor coefficient** can now be constructed through shell contact on a
controlled positive-cone collar. The construction includes scalar degrees
zero, two and four, the degree-two tensor sector, all physical Israel/scalar
conditions, and propagation of the gravitational constraints. It uses the
actual registered background and the inherited smooth first-order contact
solution and localized second-order initial-constraint theorem.

The [combined proof](FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md) fixes a
background-reference harmonic gauge. Four mixed metric Dirichlet conditions,
eleven physical shell conditions and one normal harmonic condition provide
the sixteen incoming conditions. The physical shell laws induce the remaining
subsidiary conditions. A finite angular characteristic construction handles
the coupled Robin matrix without assuming a symmetric boundary energy.

The scalar calculation also identifies a failed shortcut. In the fixed first
longitudinal gauge the complete degree-four angular scalar-TF residual is
\((8u^2-2v^2/3)(\operatorname{Hess}F_4)^{\rm TF}\). The unchanged first scalar
reconstruction has zero linear residual in that row and cannot accommodate
nontrivial compact displacement data. Repairing that trace relation alone
still fails a mixed equation. The full construction keeps the extra fields.

- [Independent complete scalar sources and general harmonic operator](astra-bulk/COMPLETE_L4_SCALAR_BULK.md)
- [Independent direct metric derivation and additional forced relations](sol-bulk/FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md)
- [Full physical scalar shell rows](boundary/SCALAR_SHELL_BOUNDARY_SYSTEM.md)
- [Exact harmonic boundary sources and constraint completion](boundary/GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md)
- [Independent harmonic construction and actual-note acceptance](harmonic-referee/ACTUAL_NOTE_REVIEW.md)
- [Critical Astra proof review](reviews/ASTRA_FULL_CONTACT_REVIEW.md)
- [Independent boundary review](reviews/BOUNDARY_FULL_CONTACT_REVIEW.md)
- [Primary literature, scope and attribution](literature/SCALAR_CONTACT_LITERATURE.md)

The standalone bulk/boundary notes prove their respective identities. Their
separate scope statements do not establish contact evolution on their own;
the combined proof supplies it for this selected finite-angular coefficient.
The smooth localized linear constraint correction is used here. Smoothness
is not inferred by differentiating the separate high-finite-Sobolev exact
precontact family.

A concrete first datum can use any thin annulus \(0<r_-<r_+<b\), with
\(\Delta=r_+^2-r_-^2\),
\[
 G(s)=\begin{cases}
 \exp[-\Delta^2/((s-r_-^2)(r_+^2-s))],&r_-^2<s<r_+^2,\\
 0,&\text{otherwise},
 \end{cases}
 \quad F|_{T=0}=G(R^2)R^2Y,\quad F_T|_{T=0}=\omega G(R^2)R^2Y,
\]
where \(Y=n_1^2-n_2^2\), \(\omega>0\) is fixed. Choose the annulus close
enough to the shell that the explicit contact-time inequality in the proof
holds. This is the commuting \(A=B=\operatorname{diag}(1,-1,0,0)\) family.
The proof supplies the subsequent coefficient by convergent finite-slab
characteristic iteration; no numerical spacetime evolution is claimed.

**Limits.** The coefficients satisfy the bulk and shell equations through
order \(\epsilon^2\). An exact nonlinear contact family, a uniform estimate
over all angular frequencies, convergence of an all-order formal series and
global stability remain unproved. The controlled domain is the positive-cone
collar, not an arbitrary ambient moving-ball or global spacetime slab.
The result does not establish a higher-dimensional Big Bang origin, external
novelty or APS readiness. Reviews are internal AI-assisted checks, not external
peer review. This checkpoint creates no new Zenodo version or APS submission.

## Reproduction

Use Python 3.12.14 and `sympy==1.14.0`. Authenticate the source and ZIP with
the [external frozen-source receipt](../../hdblast/continuation/2026.10.14-full-second-order-contact/FROZEN_SOURCE_RECEIPT.json)
before executing project code. The manifest pins the proof, six checkers,
saved results and four unchanged controlling source snapshots.

```sh
python -m pip install -r requirements.txt
python reproduce.py --manifest-sha256 EXTERNAL_MANIFEST_SHA256 --output-dir /tmp/hdblast-full-contact-fresh
```

The output directory must be new or empty. Each checker runs from a fresh
directory normally and under `python -O`; all scientific output must match
the saved result. The runner also requires fourteen actual corrupted-formula
executions to fail for their specified reason. It does not simulate the
spacetime or certify the analytic proof from test counts.

| Independent verifier | Positive checks per mode | Evidence scope |
|---|---:|---|
| Astra bulk | 52 | Arbitrary symbolic base jets and exact sphere decomposition |
| Sol direct metric | 69 | Exact rational local metric/Christoffel comparisons and embedded controls |
| Sol angular | 38 | Global symbolic sphere decomposition and vector selection |
| Physical boundary | 39 | Complete Taylor sources, gauge conversion and shell identities |
| Principal boundary | 31 | Frozen boundary rank and subsidiary identities |
| Harmonic referee | 53 | Boundary sources, initial gauge completion and a nonsymmetric Robin ODE |

Counts remain separate because several identities overlap. Historical
receipts retain original workspace paths; distributed bytes are pinned by
MANIFEST.json. Downloaded full-paper texts, opaque caches and private chats
are excluded. The next scientific target is a constraint-compatible boundary
estimate uniform in angular frequency, needed for an exact nonlinear contact
theorem.
