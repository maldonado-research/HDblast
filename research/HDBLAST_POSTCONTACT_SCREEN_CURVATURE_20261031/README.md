# A conditional magnetic screen coefficient after shell contact

At the fixed registered event `(x,r)=(5/4,1/8)`, the selected second-parameter
magnetic screen coefficient satisfies the outward enclosure

```
h in [-0.404311163330, -0.014601815851].
```

The interval excludes zero. Its backward cone reaches the shell history
through `x-r=9/8`, after the initial packet support edge can reach the shell.
The calculation retains the original bulk source, inhomogeneous Neumann
history, complete scalar TT feedback allowance, and same-event quadratic
curvature term. It encloses the actual conditional formal TT target rather
than replacing it by a finite feedback iterate.

Read [the derivation and assumptions](DERIVATION.md),
[the exact primary report](main/POSTCONTACT_NORMAL.json), and
[the independent numerical review](numeric/PUBLIC_NUMERIC_REVIEW.md).
[The next step](NEXT_STEP.md) specifies the initial-family and observer
identifications needed before transferring this coefficient to an exact
family and a transported observation.

This is a coefficient for the specified slicing observer and screen. The
shell-source intervals include zero; retaining boundary history does not
establish boundary causation. No transported gauge-invariant observer,
detector strain, physical force, stability or Big Bang-origin conclusion
follows. Earlier conditional exact initial-data, precontact and controlled
shell-evolution theorems retain their own assumptions. Applying an exact
family theorem here requires the appropriate initial/gauge and causal-domain
identifications, including second initial jets or a sufficient TT projection
match. A shared first tangent alone is insufficient.

## Reproduce

Use a complete checkout of this repository, Python 3.12.14 and the pinned
symbolic dependencies. From the repository root:

```sh
python -m pip install -r research/HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/requirements.txt
python -I -S -B research/HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-postcontact-fresh
```

The output directory must not exist and must be outside the source checkout.
The replay authenticates the closed capsule inventory and six inherited
source capsules before scientific imports and again afterward. Selected
additional source files are separately pinned. It stages immutable copies
outside the checkout; explicit path adapters only redirect scratch wrappers.
The inherited source capsules and all scientific implementations are
preserved. No scientific network access is required after checkout and
dependency installation.

The replay must reproduce nine canonical report pairs in normal and
optimized Python, reject 68 real corrupted subprocesses with the expected
failure and no success artifact, and reproduce the primary numerical result
from a relocated code/source copy. Five distinct report-only corruptions
are also checked in both modes; they are separate from the subprocess count.
The independent arithmetic/cover code review executes 4,016 exact cases per
mode. Counts describe verification work, not numbers of new theorems.

The primary numerical campaign uses 140 inherited actual-first-jet calls.
The separate numerical implementation uses 119 calls and its own interval
arithmetic, source evaluation, cone cover, shell and screen formulas. Both
share the inherited actual-first-solution/error oracle; this is an explicit
limit on their independence. Report-only audits independently reconstruct
aggregation from authenticated interval inputs and do not reprove those
inputs. Internal AI-assisted mathematical reviews are not external
specialist review or proof-assistant verification.

The [bounded primary-literature comparison](literature/METHODS_SCOPE.md)
was written during target preparation. Its prospective wording records that
stage; the completed reports here pass the stated sign test. It distinguishes
the old exterior-only exporter from the global first-jet API used here.
Established second-order magnetic-curvature methods are not claimed as new.
The search is not exhaustive and makes no originality determination.

`MANIFEST.json` lists every capsule member except itself. Execution logs,
runtime receipts, raw third-party retrievals and private editorial history
are excluded. No new Zenodo version or journal submission accompanies this
checkpoint. Earlier capsules, manuscripts and websites are preserved.
