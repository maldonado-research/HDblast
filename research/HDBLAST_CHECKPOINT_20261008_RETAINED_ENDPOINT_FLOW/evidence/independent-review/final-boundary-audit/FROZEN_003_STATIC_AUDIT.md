# Frozen 003 independent static boundary audit

Status: **BLOCKED — pre-authentication Python import aliases remain possible.**

Registration: `manufactured-registration-003.json`, SHA-256
`6773210daf54fd3b666e3ea15396f74799ece2f7e95bb0c71f2904784359d660`.
The registration contains 40 files. The five primary files below were read from
the implementation and verified identical to the reviewer's frozen-source-003
copy and the registration's exact size/hash pins.

This audit was static. No candidate module was imported, no original retained
array was opened, no values were decoded, and no physical source callback was
invoked. There were no implementation, repository, or remote edits.

## Blocking finding: Python source aliases execute before closure checks

`execution/bounded_launcher.py:36–41` scans for cached/native executable names,
prepends its execution directory to `sys.path`, then imports `kernel_guard`.
`execution/run_endpoint.py:7–12` similarly scans, prepends its execution directory,
then imports `registration_guard`. Complete source closure authentication occurs
later (`bounded_launcher.py:65–84`, `run_endpoint.py:25`).

The scans do not reject pure Python package aliases. An extra
`execution/kernel_guard/__init__.py` or
`execution/registration_guard/__init__.py` wins normal import resolution over
the corresponding pinned `.py` file. Its code executes before the registration
or GO is authenticated, even for the supported `python -I -B` invocation.

The pinned Python 3.12 standard-library resolver confirms this statically:
`importlib/_bootstrap_external.py:1618–1625` tests and returns a package's
`__init__` before `:1630–1640` checks a same-name module file. A package-directory
symlink is another alias form because the pre-import walks do not reject every
symlink before the first helper import.

There are also ordinary sibling-module aliases. After the launcher prepends the
execution directory, its `kernel_guard.py:8` imports `ctypes`, permitting an
unregistered `execution/ctypes.py` to run first. After the worker prepends the
execution directory, `registration_guard.py:2` imports `hashlib`, permitting an
unregistered `execution/hashlib.py` to run first; `resource` is another dependency.
Rejecting only matching helper package names does not cover this broader path.

Required correction: protect all pre-authentication import resolution, including
pure Python packages, sibling standard-library shadows and symlink aliases.
For example, load trusted bootstrap helpers by their exact file paths while
keeping candidate directories off `sys.path` until closure authentication, with
the bootstrap source trust and safe-file rules explicit. Alternatively, perform
the complete unexpected-artifact/alias check before adding any candidate import
directory. Preserve the bytecode/native protections. A new freeze is needed.

The scripts currently import standard-library modules before checking isolated
mode. Supported `-I -B` invocation protects these initial imports. An explicit
early isolated-mode check would make unsupported direct invocation fail clearly;
this is an integration hardening note, separate from the confirmed `-I` blocker.

## Other requested invariants

Subject to fixing bootstrap execution before authentication, no additional
blocker was identified in the scoped static review:

- Registration requires the mandatory semantic file set and exact complete
  noncache file closure; exact contract values are compared. Physical GO requires
  the new scope, matching registration and all file pins, explicit source/decode
  authorization, public-byte verification and independent-review success.
  Legacy scope or an omitted semantic subset cannot satisfy these checks.
- The metadata adapter validates four fixed capsule identities/order, 19 unique
  members per capsule, exactly the nine selected member names and expected
  selected types/shapes. Each codec verification authenticates an opaque whole
  snapshot and all complete member streams before headers; all four verified
  capsule objects exist before the worker constructs sources or requests its
  first selected iterator. This gives 76 complete member authentications before
  the 36 selected array iterators.
- Actual source construction calls the authorization closure before source
  algebra. Deterministic source/cell payload checks and the frozen double loop
  provide exactly 128 callbacks, followed by an explicit count guard. This path
  is separate from the manufactured polynomial provider.
- Each node's six complex saved modes are divided exactly once by the frozen
  represented epsilon. Both endpoint evaluations use normalized `u_1,w_1` at
  anchor -9/2; saved `u_2,w_2` and `u_3,w_3` are comparison values. Each endpoint
  target is evaluated directly from the same anchor, not chained from endpoint
  one. The engine performs no additional epsilon normalization.
- Final guards fix 49,152 nodes, 98,304 comparisons, 36 selected arrays, 688,152
  decoded real slots and 24 prefix cases. Fourteen signed-zero flags are
  exported per node; six observation-time slots per capsule account for the
  remaining 24 real slots. No original physical execution was attempted here.

## Exact reviewed primary file pins

| File | Bytes | SHA-256 |
|---|---:|---|
| execution/registration_guard.py | 6402 | 92b0bf0835d8caae49d89d39de3e34d8a647c4a00f224b75cf81be2ad6fd2753 |
| execution/run_endpoint.py | 2381 | fa216b8392aa405323f35d7967e0e1f0d118ab3776ed77f7d4054211a2d30c35 |
| execution/bounded_launcher.py | 25519 | 9167e2d3c75a597d8cfa041be438db1cf5463750201fd6102b5befa18bc84243 |
| execution/endpoint_worker.py | 16221 | 6572f65c45dfb8ffbc07f117af2d41c9ced4b346d845ebf5378eec5d9eb2f1c5 |
| provenance_decoder/input_spec_adapter.py | 4140 | 535983ca6a7047df997f41d807a4d0a4f95e72a4e0ffee395586e4d075b9fe61 |

## Parent's subsequent harmless dynamic controls

The parent reviewer independently reproduced three pre-authentication alias
executions under `-I -B`: the worker's `registration_guard/__init__.py`, the
launcher's `kernel_guard/__init__.py`, and the launcher's sibling `ctypes.py`.
The parent's separate worker `resource.py` control did **not** execute its marker;
help completed successfully. Accordingly, the earlier resource reference is a
static possible-dependency observation, not a confirmed bypass. The worker's
sibling `hashlib.py` possibility has not yet been dynamically tested.

The preserved parent evidence is
`../bootstrap-python-alias-frozen-003/RESULT.json`, including source pins and
logs. These results confirm the blocking finding while delimiting which aliases
were observed. This auditor did not run those controls or execute candidate code.
