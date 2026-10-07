# Authenticate and reproduce the incoming-state checkpoint

This is an exact algebra/rational replay. It reads only the authenticated new
payload and uses no physical source callbacks, saved quantum arrays, physical
trajectory, network access, GitHub operation or Zenodo action. The first
primary receipt predating the paired-omission control is preserved under
`primary/development/initial/`; final replay uses the primary source at the
top of `primary/`.

Use **Python 3.12.14** and **SymPy 1.14.0**. The primary route requires only
the standard library; the independent route uses SymPy. In this cloud task,
the pinned environment is available at
`/workspace/hdblast-research-work/environment/bin/python`. A fresh isolated
Python environment with `python -m pip install -r requirements-replay.txt` is sufficient
for this checkpoint; application NumPy/SciPy dependencies are unnecessary.

From the repository root, use the manifest SHA256 recorded in the separately
verified review/public-freeze receipt and direct outputs to a new path outside
the checkpoint. Replace `VERIFIED_MANIFEST_SHA256` with that digest:

```sh
python -B research/HDBLAST_CHECKPOINT_20261007_INCOMING_STATE/replay.py \
  --expected-manifest-sha256 VERIFIED_MANIFEST_SHA256 \
  --output-directory /tmp/hdblast-incoming-state-fresh-replay
```

For a standalone extraction, use its `replay.py` path instead.
The output directory must not already exist. If it does, choose another fresh
path; preserve previous receipts. Outputs outside the checkout and `-B` keep
the authenticated payload unchanged.

The replay validates the exact payload inventory, every file size and SHA256,
then executes both independent routes in ordinary and optimized Python. Each
new JSON receipt must match its archived counterpart in all fields. It then
compares the ordinary/optimized scientific fields while excluding only
`python_optimization`, checks the independently implemented eight rational
coefficient families at K=64,128,256, and reauthenticates the payload.

The manifest does not hash itself. Replay requires its externally pinned hash
before reading it and checks the same hash again afterward. Obtain the digest
from authenticated bytes at the independently verified repository commit or
from its verified freeze receipt; a digest merely computed from the adjacent
downloaded manifest does not establish that provenance. The initial development receipts are
authenticated historical files; they are not used as final expected results.
The pinned inherited text documents are premises, not executable physical
inputs. Their original public paths, repository commit and identical copy
hashes are recorded in `SOURCE_PINS.json`.

Receipt success verifies the displayed conditional mathematics and stable
implementation output. It does not establish the hypotheses for the actual
incoming state or complete the unchanged twelve-case numerical certificate.
The hypothetical amplitude in the sensitivity table is not an achieved error
or a new gate. Historical metric **FAIL** and Big Bang **NOT_ESTABLISHED**
remain unchanged.
