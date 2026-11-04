# Retained actual packet and continuous source ledger

This component binds the current canonical first packet to the retained
actual quadratic source certificate, integrates existing continuous source
boxes on the future observation cone, and records the packet's normalized
principal characteristic weights. It uses Python 3.12.14 and the standard
library only. No source/background/bump oracle or PDE is rerun.

Read [SOURCE_LEDGER.md](SOURCE_LEDGER.md) for the complete source adapter,
normalizations, continuous cell geometry and conditional energy consequence;
read [REGISTERED_PACKET_AUDIT.md](REGISTERED_PACKET_AUDIT.md) for the current
first-state binding, the principal packet diagnostic and corner constraints.

The cellwise mean bulk squared-L2 cap is about `8.52642e24`, compared with
the uniform forcing-integral ceiling about `8.50291e25`. The resulting
auxiliary zero-initial trace bound is about `2.74e13`, still looser than
the retained Duhamel bound about `5.4e12`. This is an improved source-integral
ceiling, not a numerical prediction improvement. The zero-initial reduced
piece is not alone a constraint-compatible physical second state. This
future-only certificate neither selects a past preparation nor determines
absolute proper-volume or expansion signs.

The characteristic ratio concerns precisely the normalized principal master
`w`. It is not a physical energy, Noether law or radiation flux, and does
not rule out past incoming preparation followed by potential scattering.

## One replay command and retained artifacts

Use a complete authenticated checkout containing the retained sources and a
fresh absolute output path outside both the checkout and this component:

```
python -I -S -B sources/checks.py --repo-root REPO --output OUT/SOURCES_CONTROLS.json
```

The wrapper creates `OUT/SOURCES_CONTROLS.artifacts/NORMAL.json` and
`OPTIMIZED.json`, compares both with the frozen references, and retains them
for root integration. Its one report records one canonical pair, two actual
positive subprocesses, eight actual rejected subprocesses, and zero
inprocess mutant evaluations. All children explicitly use `-I -S -B`, with
`-O` added for the optimized mode. Each produces one JSON stdout object and
empty stderr. Each rejection must return exit 2, its exact witness, and
either no output artifact or the byte-identical preexisting artifact.

The four negative cases, each in both modes, are a modified local source
pin, a modified frozen design scale, an occupied output, and a relative
output. They do not rerun inherited scientific proofs. The component rejects
symlink ancestors before resolving source/output paths. An output location
must be fresh and outside code and repository. Arithmetic and malformed
input exceptions use the same single-JSON rejection protocol.

Direct producer use is available when a separate report is needed:

```
python -I -S -B sources/source_budget.py --repo-root REPO --output OUT/NORMAL.json
python -I -S -B -O sources/source_budget.py --repo-root REPO --output OUT/OPTIMIZED.json
```

Root integration should use the wrapper once and consume its retained
canonical reports, avoiding duplicate positive producer calls.

## Exact preservation and selected dependencies

`FROZEN_DESIGN.json` records the original canonical report SHA256 and the
scientific projection digest. Compared with that original dictionary, only
the top-level `schema` and `scope` fields change and `public_provenance` is
added. Excluding exactly those three keys, every scientific value and
structure is identical. The producer verifies this projection before
writing its report. The scientific arithmetic body is unchanged.

`SOURCE_PINS.json` preserves the original 26 directly pinned inputs and its
hash. Accordingly the original report's `input_files_authenticated=26`
remains part of the unchanged scientific projection. The new
`SOURCE_CLOSURE_INPUTS.json` authenticates the union of 35 selected sources,
including packet pin declarations and the additional proof citation; the
new `public_provenance` separately records that union and the seven packet
pin comparisons. Duplicate declarations of the same first-state derivation
are explicitly retained in `PACKET_SOURCE_JOIN.json`. These selected checks
supplement the parent checkpoint's authentication of its inherited closed
capsules; they do not reprove their analytic premises.

All report provenance paths are logical roles or repository-relative names.
The old first-state/source identity remains conditional on the retained
smooth PDE and mixed uniqueness theorems. No new background state, altered
first velocity, past-time source coverage, common nonlinear parameter
window, physical energy law, matter ledger or physical scale is supplied.
