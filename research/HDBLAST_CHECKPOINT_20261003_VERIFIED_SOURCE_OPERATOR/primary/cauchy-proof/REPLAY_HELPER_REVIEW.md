# Portable replay helper audit

Verdict: **GO for the reviewed helper's manifest/path checks, bootstrap call
ordering, two-route replay composition, and exact-comparison semantics.** No
scientific replay or production source callback was executed for this audit.

Reviewed helper:
`/workspace/hdblast-research-work/verified-integration-20261003/primary-design/replay_checkpoint.py`

Reviewed final SHA256:
`6b982c8959e046059030ac8b40ff3db4e04e53d6ed4a4c52b2d3b467ff744ad3`

This review covers the final version that rejects dangling output symlinks and
symlink output ancestors. It supersedes review of the earlier 2f8408e9... pin.
Only source files and private manufactured standard-library archives were used.

## Manifest and paths before archive code imports

The helper imports standard-library modules at top level. Its first replay
operation validates the exact externally supplied manifest SHA256. The
manifest loader rejects duplicate keys and nonfinite JSON constants. The
manifest must have exactly schema_version and files, integer schema 1, a
nonempty file dictionary, and no entry for MANIFEST.json itself. Each file
entry has exactly a nonnegative integer byte count and a lowercase 64-hex
SHA256. Required registration, execution helper, review/guard, and both
archived route output/receipt files must be present.

The checkpoint directory must be real and have no symlink ancestors. Manifest
member names are canonical relative POSIX paths: no absolute path, traversal,
dot segment, doubled separator, trailing separator, or backslash is accepted.
Every component is checked for symlinks. Members must be regular files. The
actual full file inventory must equal the manifest inventory, and every byte
count and streamed SHA256 must match. Unlisted files, member symlinks, special
files, changed members, and missing members are rejected.

Finally, the bytes of the executing helper must have the exact hash of the
archived execution/replay_checkpoint.py manifest member. Thus an externally
repinned archive with a different helper does not authorize this running
helper. All these checks occur before the first guard or review module load.

The interpreter, receipt pin, and registration pin are explicit caller inputs.
The receipt must be a regular non-symlink file matching its explicit hash.
The interpreter must be an executable file. No interpreter, credential,
custodian root, or cloud workspace path is hardcoded.

The new output directory must not exist, must be outside the checkpoint, and
must not itself be a symlink. Every output ancestor is rejected if it is a
symlink. This prevents an apparently external path from writing through an
alias into the frozen archive; the extra private control verifies that no
directory is created inside the archive in that case.

## Authentication and actual baseline before either custodian call

The ordered bootstrap is:

1. Validate the complete manifest and self-helper pin.
2. Validate explicit receipt/registration pins, interpreter, and output path.
3. Load the verified guard and authenticate the registered source/archive.
4. Load the verified completed-execution reviewer and review both archived
   actual execution receipts and route outputs.
5. Revalidate the full manifest, then create the fresh external output and
   save the archived-baseline review.
6. Revalidate the manifest before each of the two custodian invocations.

The existing completed-execution reviewer authenticates before reading the
actual receipts. Its inspected code distinguishes actual and fabricated
executions, checks route/output pins and resource conditions, checks source
chronology, and classifies the registered numerical pair. The helper relies
on that separately authenticated reviewer; it does not substitute its private
test stubs for a production review.

No subprocess/custodian call occurs on an invalid manifest, path, bootstrap
pin, or missing actual baseline receipt. Even the external output directory
is not created when actual-baseline review fails. Private import markers
confirm that malformed manifests and invalid output/checkpoint paths are
rejected before guard and reviewer imports, while the nonactual receipt case
reaches the authenticated mock review and stops before output or launch.

## Two bounded routes and preserved failure semantics

The helper launches exactly primary and independent, in that order, through
the existing execution/execute_bounded.py. It passes explicit checkpoint,
receipt, receipt hash, registration hash, interpreter, route, and fresh
external output directory. The existing custodian owns nonroot execution,
source authorization, CPU/address-space/wall limits, and authoritative route
receipts. The helper does not add a source-evaluation path.

Both fresh route results must pass the authenticated completed-execution
review. Payloads are compared by exact canonical JSON bytes; complete stable
scientific frames are also compared, using the reviewed frame definition that
omits only specified internal runtime measurements. Numerical enclosures,
models, bounds, and other stable method evidence remain in the comparison.
The full frozen manifest is checked again before the final report.

An exact replay success is distinct from scientific success. The archived and
fresh scientific status strings are both preserved and must match. A negative
or unresolved archived state can be reproduced exactly, yielding a successful
reproduction report with that same negative scientific state. It does not
upgrade the source/operator result or establish a pressure/contact certificate.
The private negative-science control verifies this distinction.

Unequal payloads or stable frames produce FAIL_EXACT_FRESH_REPLAY with exact
comparison records retained. A failed bounded launch or review retains logs,
route artifacts, commands, exit codes, and a failure report. The archived
payload files are never intentionally rewritten by the helper.

## Private fabricated evidence

The existing test_replay_fabricated.py controls pass in both normal Python and
Python -O at the final helper pin: 11 cases in each mode. They cover exact
positive and optimized replay, preservation of negative scientific status,
retained unequal values, bad manifest/receipt/registration pins, changed or
extra members, an existing output, and a missing actual archived receipt.

Eight additional private controls pass in each mode:

- Output parent symlink into the checkpoint is rejected before imports.
- Dangling output symlink is rejected before imports.
- A different archived helper, even with a repinned complete manifest, is
  rejected before imports.
- Duplicate manifest keys and extra manifest schema fields are rejected
  before imports.
- Symlink checkpoint and symlink payload members are rejected before imports.
- A nonactual baseline receipt is rejected before output creation or launch.

The checks use explicit exceptions, so their verification remains active under
-O. All archives are private temporary manufactured fixtures with standard
library guard/review/custodian stubs. Additional rejection cases launch zero
mock custodians; every control executes zero production source callbacks and
decodes zero checkpoint arrays. These controls establish helper behavior and
ordering, not the existence or validity of an actual scientific archive.

Evidence is staged in `cauchy-proof/replay-helper-check/`:
ORIGINAL_CONTROLS_NORMAL.json, ORIGINAL_CONTROLS_OPTIMIZED.json,
EXTRA_CONTROLS_NORMAL.json, EXTRA_CONTROLS_OPTIMIZED.json, and the extra private
control source check_replay_paths.py. Every receipt pins the reviewed final
helper SHA256 above.

No outstanding finding remains in this audit's stated helper scope. Root's
authenticated actual archive, source proof, runtime contract, custodian
authorization, and scientific replay remain separate prerequisites.
