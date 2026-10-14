# Shell Noether identity and rotational obstruction

This capsule records an exact nonlinear shell Noether identity, an independent
second-order rotational obstruction for smooth compactly supported packets,
and the quadrupole criterion that, within the fixed radial-profile, time-scale,
and reflection-parity family, all ten Lorentz quadratic tests vanish exactly
when the two trace-free symmetric harmonic matrices commute. Its four boost
tests vanish by antipodal parity; the six rotational tests supply the
commutator condition. These are necessary
restrictions under the stated classical-family assumptions. They do not prove
nonlinear existence, sufficiency of the rotational criterion, positive physical
energy, stability, or cosmological evidence.

Read the [manuscript](manuscript/SHELL_NOETHER.pdf), its
[source](manuscript/SHELL_NOETHER.md), and [results](RESULTS.md).
The [next calculations](NEXT_CALCULATIONS.md) identify the remaining work.

The primary proofs are the [nonlinear identity](derivation/FULL_NONLINEAR_SHELL_NOETHER_IDENTITY.md),
[independent rotational argument](independent/ROTATIONAL_SECOND_ORDER_SHELL_OBSTRUCTION.md),
and [quadrupole criterion](root-audit/QUADRUPOLE_ADMISSIBILITY.md).
The [reading register](literature/PUBLIC_PRIMARY_READING_REGISTER.json) distinguishes
selected primary passages from metadata-only references.

Use Python **3.12.14**, SymPy **1.14.0**, and mpmath **1.3.0**. After verifying
the manifest and replay-source SHA256 against the separately published frozen
pins, authenticate the full capsule with the standard library:

```sh
python -I -S -B reproduce.py --verify-only --manifest-sha256 MANIFEST_SHA256
python -I -B -m pip install --only-binary=:all: -r requirements.txt
python -I -B reproduce.py --manifest-sha256 MANIFEST_SHA256 --output ../shell-noether-replay-001
```

Replace `MANIFEST_SHA256` with the externally recorded final hash. The output
directory must be fresh and outside this capsule. `--verify-only` is read-only.
The runner copies the capsule to a new location, authenticates all files and
selected prior proofs before science execution, runs all three checkers in
normal and optimized modes, compares complete result bytes and exact named
inventories, and tests failure injection and provenance corruption. Each mode
checks **81 positive algebra checks and 21 wrong-formula controls**. Optimization
mode remains explicit in the aggregate replay receipt.

The rotational public checker has only two serialization metadata fields
removed so scientific receipts are identical across modes. Its
[original](independent/original/check_rotational_constraints.py),
[exact diff](independent/SERIALIZATION_ADAPTATION.diff), and
[adaptation record](independent/SERIALIZATION_ADAPTATION.json) preserve that change.

The [prior input register](basis/SOURCE_INPUTS.json) authenticates seven selected
unchanged public proof bodies against the full prior shell-contact manifest
`fadd405ece09e04740f1ee681c6086b84918b54b5e221146c17b5229c364c52e`.
The earlier cone-evolution manifest is separately pinned to
`a873132efe11755217ddb865fd0156041b8c1ad21ae1c45310183da4fb6140ae`.
Prior scientific scripts are not replayed. Public literature notes contain no
source PDF cache; exact algebra does not mechanically certify the analytic
second-order argument or a PDE existence theorem.
