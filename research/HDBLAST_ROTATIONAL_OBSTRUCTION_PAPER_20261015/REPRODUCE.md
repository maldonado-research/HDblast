# Authentication and exact replay

Use Python 3.12.14, SymPy 1.14.0 and mpmath 1.3.0. Install the pinned dependencies from `supplement/requirements.txt` in a suitable environment before the offline replay. Expected outputs and algebra checks are described in `supplement/README.md`.

The external manifest and archive SHA-256 values are recorded in the repository's `.github/workflows/hdblast_rotational_obstruction_paper.yml`. Obtain those pins from the reviewed repository commit. A checksum obtained only from the same untrusted download is an integrity comparison, not independent authentication.

From the capsule directory:

```sh
python -I -S -B verify_manifest.py --expected EXTERNAL_MANIFEST_SHA256
python -I -B supplement/reproduce.py --output-dir /absolute/path/to/new-replay-output
```

Replace the manifest placeholder with the workflow's full 64-character value. The output directory must not exist and must lie outside `supplement/`; the runner refuses overwrites. Its final success status is `PASS_EXACT_AGREEMENT`. Any failed or interrupted run is not a successful replay.

The manifest covers every capsule file except itself, rejecting missing/extra files, symlinks and changed bytes. The workflow first checks the manifest, verifier, runner and archive against its independent pins, checks that ZIP members match the source tree byte-for-byte, and verifies the complete manifest before installing dependencies or importing project code. It then runs the supplement from a relocated copy.

The supplement compares full captured reports, checks optimization-mode agreement with only the declared optimization metadata removed, runs deliberately injected failures, and checks that package bytes remain unchanged. These are reproducibility and finite-algebra checks. The written analytic proof, model assumptions, literature comparison and significance require separate scientific assessment.

Build the PDF in a separate copy of `manuscript/` to avoid adding LaTeX build files to the frozen capsule. Do not write replay outputs or Python bytecode into the authenticated tree.
