# Natural first fields, physical shell traces and bulk source compatibility

30 October 2026. Exact implementation checks against the existing analytic
Einstein–scalar construction for the retained admissible quadrupole. Original
sources, numerical candidates, error certificates and frozen manuscripts
are unchanged.

Three joins are now checked explicitly:

- The implemented longitudinal reconstruction gives the natural Einstein
  and scalar residuals as exact multiples of the master residual and its
  required radial derivative. Independent tensor calculations use different
  coordinate constructions. The cancelling 5/2 terms, variable-rho
  derivatives and nonharmonic scalar term are retained.
- The physical first scalar shell residual is an exact multiple of the
  difference between the bulk master residual and the actual dynamic trace
  residual. Metric junctions and shell-fixing gauge cancellations are
  checked, with tangential gauge values left free.
- All eight natural quadratic Bianchi source rows match independently
  constructed first-residual corrections. The registered seven first-metric
  profiles and scalar quadrupole profile have free coherent ordinary jets,
  before imposing the first equations or harmonicity. A separate full-frame
  contraction verifies the modal operator and its transverse-traceless kernel.

The sufficient analytic reconstruction already existed. This checkpoint
binds the current implementation to it; it does not claim a new PDE existence
theorem. The [combined derivation and conditional chain](theory/CONDITIONAL_CHAIN.md)
states how the joins compose.

## Reproduce

Use a complete repository checkout and Python 3.12.14. All distributed
checks use the standard library. Run from the repository root:

```sh
python -I -S -B research/HDBLAST_NATURAL_FIELD_BINDING_20261030/checks/run_checkpoint.py \
  --repo-root . --output-dir /tmp/hdblast-natural-field-fresh
```

Choose a fresh output directory. The runner authenticates this closed
capsule before and after execution, reproduces five complete canonical
scientific reports in normal and optimized Python, and checks 140 actual
corrupted subprocesses, 70 per mode. Each must fail for a mathematical
mismatch and leave no successful report. The Bianchi driver additionally
records commands, runtimes and exact nonzero witnesses in its fresh evidence.

| Canonical report | Route | Corruptions per mode |
| --- | --- | ---: |
| [NATURAL_MASTER](reference/NATURAL_MASTER.json) | Longitudinal connection, natural Einstein and direct scalar calculation; actual source joins | 22 |
| [NATURAL_INDEPENDENT](reference/NATURAL_INDEPENDENT.json) | Separate flat-slicing coordinate calculation and global quadrupole master conversion | 22 |
| [FIRST_SHELL](reference/FIRST_SHELL.json) | Physical one-copy normal/tension equations and dynamic trace identity | 12 |
| [BULK_BIANCHI](reference/BULK_BIANCHI.json) | Natural source map, independent linear residuals, coherent jets and complete angular reconstruction | 13 |
| [B0_INDEPENDENT](reference/B0_INDEPENDENT.json) | Separate connection table on arbitrary registered tensor amplitudes and first derivatives | 1 |

Counts describe the replay roster, not independent theorems. Complete report
bytes must match the references and agree between modes. Runtime evidence
with elapsed times and local command paths is generated separately.

Sparse checkouts need this capsule and the five unchanged source roots
`HDBLAST_REGULAR_CONE_EVOLUTION_20261013`,
`HDBLAST_REGULAR_SHELL_SPECTRA_20261013`,
`HDBLAST_INTERIOR_MASTER_JETS_20261015`,
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016` and
`HDBLAST_BULK_SHELL_CONSTRAINT_BINDING_20261030`, under `research/`.
The master/shell checks authenticate complete stated capsule manifests;
the Bianchi checks authenticate the listed map and helper bytes they use.
Historical numerical producers are not imported.

## Derivations and internal review

- [Natural first-master factorization](theory/NATURAL_FIRST_MASTER_BINDING.md)
- [Physical shell and actual trace equation](theory/FIRST_SHELL_MASTER_BINDING_NOTE.md)
- [Residual-corrected natural bulk Bianchi](theory/NATURAL_BULK_BIANCHI_NOTE.md)
- [Master/Bianchi composition](theory/MASTER_BIANCHI_INTEGRATION.md)
- [Independent Sol tensor derivation and producer review](review/SOL_FIRST_NATURAL_REVIEW.md)
- [Independent physical-shell review](review/INDEPENDENT_FIRST_SHELL_REVIEW.md)
- [Independent Bianchi implementation/operator review](review/BIANCHI_IMPLEMENTATION_REVIEW.md)
- [Bounded primary-literature comparison](literature/PRIMARY_METHODS_COMPARISON.md)

The independent first-master implementation uses its own arithmetic and
coordinate tensor construction. The Bianchi operator reviewer shares audited
arithmetic/sphere primitives but does not import the Bianchi producer. It
checks the operator, not every quadratic source component: the eight
covector rows cannot validate the transverse-traceless source components
annihilated by that operator. These are internal AI-assisted reviews, not
external specialist assessments.

## Assumptions and remaining work

The actual smooth background, master and gauge solutions and their local
PDE identities remain analytic inputs, together with corner regularity and
the existing subsidiary uniqueness argument. The scalar equation needs the
master residual and its derivative to vanish on a neighbourhood; an interval
containing zero or a zero value at a point is insufficient. Homogeneous
Bianchi specialization needs E1 and its first prolongation to vanish and
S1=0. Replacing the general scalar source by its harmonic formula additionally
needs C1=0. Harmonicity alone supplies neither first field equation.

No new numerical third-jet adapter, smaller response error, physical force,
exact finite-amplitude family or Big Bang origin is established. Earlier
quantitative certificates retain their assumptions and absolute/pointwise
limitations. The next implementation target is the reduced-source Bianchi
identity with full-metric-lowered QC_actual and general PC1 terms retained
before specialization. The analytic subsidiary argument already exists.

Only public research derivations, checks and references are included. No
new Zenodo version or journal submission accompanies this checkpoint.
