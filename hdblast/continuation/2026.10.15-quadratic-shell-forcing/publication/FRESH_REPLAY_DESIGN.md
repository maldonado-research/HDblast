# Fresh quadratic-forcing replay design (roster pending)

This is a preparation contract, not an approved execution roster or a claim
that any new calculation has run. Final checker paths, copied dependencies,
output schemas, golden outputs, formula controls and counts remain pending.

1. Authenticate one read of the complete source manifest against an external
   SHA256 pin before importing or executing any scientific file or its local
   verifier. Verify all regular source members, canonical paths and exact
   inventory. Authenticate the runner, verifier, requirements and roster;
   enforce the actually tested Python/dependency versions in CI.
2. Read the frozen `REPLAY_ROSTER.json`. Require exact agreement between its
   declared positive process identities, scripts, group/embedded counts and
   controls (checker, name, expected reason, mechanism) and the runner's code.
   Counting totals alone is insufficient. Refuse incomplete or unapproved
   roster declarations; no inherited positive or negative counts.
3. For each positive and corrupted process in both normal and optimized modes,
   create a distinct copied workspace outside the source. Authenticate each
   copied dependency, including any minimal retained first-variation inputs.
   Exclude golden output files from executable copies unless a declared checker
   consumes them as certificate inputs. Generated outputs start absent.
4. Invoke isolated processes with `-I -B` and `-I -B -O`. Require explicit
   checks independent of `assert`. Parse the entire expected JSON schema and
   compare all scientific fields against immutable source goldens. Compare
   complete deterministic generated JSON/CSV bytes, not selected endpoints or
   headlines. Compare normal/optimized scientific outputs; only individually
   checked execution-mode fields may be normalized.
5. Each prepared formula, input or certificate corruption runs on its own fresh
   authenticated copy. Record original and changed bytes, exact mutation recipe,
   unchanged checker hash, exact exit status and complete diagnostic. Require
   exit 1 with the specified mathematical failure reason and no stale passing
   output; unrelated exceptions or mere nonzero exits cannot count as controls.
   The receipt's `control_mechanism` equals the frozen declaration exactly.
6. Keep copied dependency/script hashes, arguments, process time, stdout/stderr,
   result bytes and all extra generated outputs in the receipt. Retain all
   result/log/CSV artifacts in public evidence; do not include disposable
   executable replay copies. Authenticate the original source again afterward.
7. The initial evidence manifest covers every retained artifact exactly once.
   Recheck the full source-plus-manifest and ZIP inventory directly from both
   Git index blobs and the final immutable commit before pushing. This catches
   ignored files that a local ZIP/source replay alone cannot detect.
8. Test the final runner on a copied prefreeze fixture, including bad external
   manifest pins, modified scientific dependencies, same-count wrong roster
   identities, stale-output defenses and each intended formula control. A final
   authenticated ZIP replay remains a separate authoritative execution.

Finite arithmetic certifies only the stated algebra/enclosure computations.
The prescribed quadratic forcing and constitutive stress source terms are not
the full `U2` response. Record unresolved metric, gauge, shell-motion and bulk
coupling terms explicitly in scientific source once the reviewers settle them.
No nonlinear closure or observational inference is supplied by a replay count.
