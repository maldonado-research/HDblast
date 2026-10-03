# Prospective active-source diagnostic reproduction

Retained arrays and earlier active-source failures were already published.
This is a prospective diagnostic motivated by known failure, not a blind
confirmation. At preparation, the new reference-flow/contact/direct-integral
attribution quantities are UNCOMPUTED, and its public freeze and result report
have not yet been made. The fixed CLI, numerical contracts and
output inventories are specified by `EXPERIMENT.json`, `OUTPUT_SCHEMA.json`
and the method notes. Do not run physical routes until every registered
source/input byte is publicly frozen and root has recorded independent remote
verification in `FREEZE_RECEIPT.json`.

The portable checkpoint directory is
`HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER`. Install Python 3.12.14 and
the exact package versions in `requirements-replay.txt`. Both numerical routes
use this interface:

```text
--checkpoint-root CPP --registration-sha256 HEX64 --freeze-commit FULL40 --output-dir FRESH
```

The full registration contains schema version 1 and a `files` map binding all
frozen relative paths to SHA256. It must bind `EXPERIMENT.json`, input manifest,
all four nineteen-member capsules, both route implementations, complete direct
subtraction/operator inventories, proof code, validator and reproduction code.
The full registration and later freeze receipt cannot self-reference. The
receipt binds the same externally supplied registration hash and full public
commit. SHA256 verification, duplicate JSON rejection, path guards and capsule
member/header authentication precede physical source calls and numerical array
decoding. The primary imports ordinary NumPy/mpmath libraries before byte
authentication, while keeping those physical operations guarded. The
independent entrypoint is a stdlib-only wrapper that authenticates registered
bytes before importing its numerical core.

Source/phase arithmetic is native long double. The primary generates Gaussian
nodes/weights with native binary80 Newton iteration; the independent promotes
NumPy binary64 `leggauss` constants to long double. MP80/100 contexts diagnose
final reduction arithmetic and explicitly scoped primary contacts, not
Gaussian constant accuracy or complete source/phase accuracy.

Registered reproduction entries are:

```json
{
  "routes": {
    "primary": "code/diagnostic_primary.py",
    "independent": "independent/diagnostic_independent.py"
  },
  "validator": "code/validate_active.py",
  "result_filename": "diagnostic.json",
  "validation_filename": "VALIDATION.json",
  "postfreeze_output_directories": ["outputs", "reports", "figures", "evidence"]
}
```

`EXPERIMENT.json` fixes the case/control/precision/resource contracts;
`OUTPUT_SCHEMA.json` contains the exact literal field inventories under
`source_constants`, additional derived checks under `additional_evaluator_gates`,
and numerical versus fatal failure semantics. Its literal science-field lists
alone are not the full evaluator gate inventory.
The common `--checkpoint-root` flag is mandatory for both routes. The
independent entrypoint is `independent/diagnostic_independent.py`, with its
separate numerical core and all contact/kernel modules also registered.
The validator accepts both route result paths and the outer execution receipt,
using normal and optimized Python in separate processes. Normal and optimized
substantive classifications must be identical. `LEDGER_ERROR_DEMONSTRATED`,
`NO_GATE_SCALE_ATTRIBUTION` and `CONSISTENCY_FAILURE` are distinct scientific
outcomes; all retain earlier metric FAIL. A negative scientific outcome exits
zero as a completed execution. Native long-double numerical controls above
`2e-7`, including quadrature, cross-route, operator and contact reconstruction
gaps, are scientific `CONSISTENCY_FAILURE` outcomes. The stricter `1e-12`
threshold applies only to final MP80/100 accumulation gaps, serialization and
exact definitional closure. Nonzero route exits, missing diagnostics,
resource overruns, integrity failures or schema/arithmetic failures retain all
logs/progress and prevent an accepted scientific classification.

Stored global-prefix/reset/fine-double Simpson scalars must agree across the
two readers at `1e-12`; the validator independently recomputes `DeltaR` from
stored endpoint profiles. Result JSON contains only local `F`, so those
agreements do not constitute a third replay of the complete native prefix from
`eta=-6`. Both frozen readers retain the original full-prefix long-double
arithmetic. The independent joint GL16/16-versus-GL24/24 control applies to the
full declared numerical universe at `2e-7`, alongside the two isolated controls.

The outer replay driver uses one 900-second / 262144-KiB budget per entire
route, including all twelve cases, all controls, both MP contexts, provenance
and serialization. Symbolic proofs and their receipts are checked during
preparation and bound by the full registration; the physical route clock does
not claim rerunning SymPy proofs. It covers byte authentication, decoding,
numerical evaluation and output serialization. The driver forbids
descendants/extra threads via Linux
`RLIMIT_NPROC=0`, requires a nonroot UID, monitors memory and elapsed time, and
records the actual positive integer `executor_uid` before launch alongside
signed exits and `wait4` peak RSS. It runs the other route after a
producer execution failure, skips scientific validation if either producer
failed, and rejects source mutation before, between and after commands.

The preparation resource checks recorded primary 219.155 seconds / 183944 KiB
and independent 473.506 seconds / 110420 KiB, separately budgeted. The
independent benchmark executed the unchanged numerical engines before final
entry/schema guard additions, which have separate normal/optimized guard
receipts. It did not execute every final production source byte as one whole
application. Original provenance/reading fixtures were fabricated; no physical
study source or array values were evaluated. These checks establish empirical
feasibility rather than a guaranteed production runtime. Actual complete-route
outer receipts remain authoritative.

```bash
python CPP/code/test_replay_package_guards.py
python -O CPP/code/test_replay_package_guards.py
python CPP/code/verify_public_freeze.py --checkpoint-root CPP --repository-root REPO --registration-sha256 HEX64 --freeze-commit FULL40
python CPP/code/replay_active_source.py --checkpoint-root CPP --registration-sha256 HEX64 --freeze-commit FULL40 --output-dir FRESH --plan-only
python CPP/code/replay_active_source.py --checkpoint-root CPP --registration-sha256 HEX64 --freeze-commit FULL40 --output-dir FRESH2
```

The guard tests create fabricated stdlib payloads and perform zero physical
evaluations. `--plan-only` authenticates bytes and records commands but launches
zero numerical children. The Git verifier verifies local Git object bytes,
not remote publication or prospective timing; root's remote verification is a
separate recorded step.

The round label `20261002` is retained, while the prospective public
registration date is 3 October 2026 UTC. Its freeze receipt records the actual
public publication and independent verification timestamps.

After all original outcomes are retained, deterministic packaging protects
every payload byte including the full registration, freeze receipt, inputs,
recorded outputs and reports. Only its own literal bookkeeping filenames and
Python bytecode caches are excluded. Rebuild the ZIP twice and compare bytes,
then extract and perform a fresh replay outside the checkpoint. Check every
scientific output field, retained failure and source pin; compare runtime/RSS
as observed quantities, not deterministic scientific bytes. Packaging performs
no numerical evaluation. It makes no claim of an unattended research service.
