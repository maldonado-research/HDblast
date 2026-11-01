# Portable independent review

This component checks the explicit preparation, its five harmonic means,
the second-only conditioning gauge, and the proper scalar/volume/expansion
budget at x=5/4. All source locations are repository-relative. There is no
Git, network, installed project package, private path or source-manifest
requirement. Only the exact algebra checker requires SymPy; the numerical
reconstructions and controls use the standard library.

Each checker authenticates the declared repository inputs, canonical reports
and local essentials before third-party imports or scientific evaluation,
and again before writing. Output is stdout by default. An explicit output
must be a new file outside both the repository and this code directory.

From this directory in the assembled repository:

```sh
python algebra_checks.py --repo-root REPOSITORY --output /tmp/mean-algebra-new.json
python reconstruct_response.py --repo-root REPOSITORY --output /tmp/mean-response-new.json
python reaggregate_initial.py --repo-root REPOSITORY --output /tmp/mean-quadrature-new.json
python run_controls.py --repo-root REPOSITORY --output /tmp/mean-controls-new.json
```

Optional `--initial-report`, `--conditioned-report` and `--response-report`
paths select the same pinned bytes from another location. Their defaults are
this capsule's initial/references and response/references files. Run the
same commands with `python -O`; all scientific report bytes must agree.

References retain 22 exact algebra/Green checks; 14,910 independently derived
conditioned component inclusions; 12 exact response report matches; the
additional 994-cell I/J/primitive audit; and six exact-witness rejection
controls. The review states the shared analytic and raw interval premises.
No independent checker invokes an actual-field or background oracle.

Negative v2 and Theta2 are second changes at the same anchored proper clock.
They do not assert negative total volume, contraction, force, radiation,
a finite-amplitude parameter range, or a cosmological origin.


Revision 02 preserves the algebra, response, and quadrature reference bytes from revision 01. Control reports additionally retain the complete final witness and stdout/stderr protocol evidence; their reasons and counts are unchanged.
Its standard-library bootstrap makes local checker imports explicit under
Python isolated mode. Run the algebra checker with `python -I -B`; run the
response, quadrature, and controls with `python -I -S -B`. All support `-O`
with identical canonical scientific reports. The quadrature audit uses only
the standard library.

Failures produce exactly one JSON line on stdout with stable status, stage
and reason; stderr is empty, exit status is two, and no result artifact is
written. Controls compare the complete object and propagate isolated mode
to their child checks. The updated response reference differs only in
custody metadata, and its exact new hash is pinned in SOURCE_PINS.json.
