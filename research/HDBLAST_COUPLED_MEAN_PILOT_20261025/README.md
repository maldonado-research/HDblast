# A bounded coupled-mean response pilot

25 October 2026. Conditional mathematical and reproducibility checkpoint.

This capsule supplies an exact finite evaluator for the first reflected
five-mean feedback, together with a bounded numerical accuracy pilot.
The numerical outcome is pending the first recorded runner execution.
No accuracy pass is asserted by this draft.

The kernel reduces the half-line retarded wave integral for
f(s,r)=F(r+s)+G(r-s) to normalized one-dimensional primitives. Reflection is
applied to the complete physical source. The resulting finite piecewise
polynomial has exact rational values and first derivatives. This is a standard
integral identity and an implementation improvement; no novelty claim is made.

The native pilot keeps the previously registered rational displacement,
constructs one continuous rational P1 velocity, and forms independent
continuous P1 displacement/derivative profiles only inside the feedback.
It retains all five metric/scalar mean couplings and the shell lift.
It stops at the declared mesh limit 256 if the complete upper error bound
does not meet 1e29. A failed upper bound is not a lower bound on actual error.

Read [DERIVATION.md](DERIVATION.md) for the kernel, C1 joins, boundary conditions
and complete error identity. The exact old input-precision gate remains
frozen and is not rerun or reinterpreted.

## Replay

Use a manifest SHA256 retained outside this capsule:

```
python -I -S -B reproduce.py --manifest-sha256 EXTERNAL_PIN --output-dir /tmp/fresh-ridge-replay
```

This portable replay needs Python and Node with BigInt support. It authenticates
the complete source before and after execution, regenerates independent exact
fixtures, and checks the Python kernel in normal and optimized modes. It does
not evaluate the retained numerical response.

To run the bounded numerical pilot from a checkout containing the six complete
upstream capsules, add explicit options:

```
python -I -S -B reproduce.py --manifest-sha256 EXTERNAL_PIN \
  --output-dir /tmp/fresh-native-ridge --native --repo-root /path/to/HDblast
```

Normal and optimized numerical runs must produce identical REPORT.json,
PROFILES.json and EVALUATED_CANDIDATE.json. When frozen results are present,
every byte is also compared with those references. The native accuracy flag
is separate from software replay success. The GitHub workflow records both.

## Scientific scope

The numerical target is a continuous approximation to the coupled five-mean
coefficient response on the retained causal cone. The exact full response is
U=B+TU, with q=3409/8000000 and the inherited bound X(U)<=1.03e33.
The approximation has a complete conditional error budget for its inputs,
coefficients, finite feedback and omitted forcing/tail.

The finite candidate generally has a nonzero remaining Robin defect. It is
included in the complete error identity, not set to zero. All ordinary-time
velocities use partial_t=ell^(-1) partial_x with ell=1e-6.

This does not supply the full angular second-order spacetime, a pointwise
second-derivative Einstein residual, an intrinsic observable value/sign,
a finite-amplitude nonlinear family, or evidence for a higher-dimensional
origin of the Big Bang. Internal AI-assisted review is distinct from external
peer review and proof-assistant formalization.

The managed cloud failed to provision. The independent JavaScript controls
were executed in V8; Python and native execution use this repository's
GitHub Actions runner. Exact sources and outcomes are retained.
Earlier capsules, frozen papers and other projects are preserved.
No new Zenodo version or journal submission accompanies this checkpoint.
