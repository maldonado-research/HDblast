# Complete local evaluation of the tighter smooth candidate

This checkpoint evaluates the [tighter smooth candidate](../HDBLAST_SMOOTH_CANDIDATES_20261026/)
including its first feedback `V_s=B0_s+B1_s`, in all 23 modes at
`(x,r)=(1/1000,1499/1000)`. Earlier work had evaluated only selected initial
and free-wave components. The [local formulas](theory/LOCAL_FORMULAS.md) avoid
subtracting large global primitive intervals: integration by parts reduces the
feedback to endpoint values, one local displacement integral and exact velocity
moments. The complete rational table is in [reference/COMPLETE.json](reference/COMPLETE.json).

## Conditional scalar-mean bounds

Combining the candidate enclosure with the earlier conditional C¹ error bound
and its zero-wedge value estimate gives, in the registered normalization,

    4.17e28 < U_w:0 < 6.36e28,
    -5.288e31 < partial_r U_w:0 < -5.283e31.

The time-derivative interval contains zero. Here `w:0` is the angular mean
coefficient of the formal second scalar perturbation in the registered gauge.
The result does not fix the scalar's value at every angle, give a physical
force or gauge-invariant observable, or establish a nonlinear family. Earlier
work already obtained other scalar-mean signs; no first-discovery or originality
claim is made. The point is before shell contact and does not test the coupled
junction evolution. The broad all-mode table does not supersede earlier sharper
specialized tensor-mode certificates.

[The detailed interval calculation](theory/NUMERIC_NOTE.md) uses pinned nominal
initial inputs, a local 32-panel interval integral, proved convolution remainders
and exact smoothed P1 velocity moments. The prior error-transfer arithmetic is
recomputed from authenticated reports; its underlying analytic and actual-field
certificates remain inherited assumptions. Affine bulk/shell sources are omitted
by this candidate and are charged in the inherited error budget.

[An independent coarse method](theory/COARSE_NOTE.md) bounds feedback over the
local characteristic cone without evaluating the feedback implementation. It
uses the previously published free-wave interval and reproduces the two signs.
The detailed candidate intervals lie inside this independent coarse enclosure.
Neither calculation supplies a numerical C² error, pointwise Einstein residual,
full inhomogeneous boundary solve or observed cosmological prediction.

## Reproduction

Python 3 standard library is sufficient. From the repository root:

```sh
C=research/HDBLAST_SMOOTH_FEEDBACK_SAMPLE_20261029
python "$C/checks/check_local_identity.py"
python "$C/checks/check_local_kernel.py"
python "$C/checks/check_complete_feedback.py" --repo-root . > /tmp/feedback-complete.json
python "$C/checks/check_coarse_feedback.py" --repo-root . > /tmp/feedback-coarse.json
python "$C/checks/compare_enclosures.py" --complete /tmp/feedback-complete.json --coarse /tmp/feedback-coarse.json
```

Normal and optimized Python must agree byte for byte. The numerical programs
pin inherited files; the nominal interval producer authenticates its upstream
roots before importing scientific code and reauthenticates after calculation.
The coarse program executes a pinned earlier free-wave checker, not the new
feedback implementation. Synthetic polynomial controls independently verify
local-kernel factors and derivative signs. Internal AI-assisted reviews are not
external specialist or proof-assistant verification. [Primary literature](RELATED.md)
provides context, not a novelty determination or a theorem for this model.

Next, extend validated evaluation through shell contact and bound the full
pointwise second-derivative residual, including bulk/source, boundary and
subsidiary errors. The full second-order solution remains unfinished, and a
higher-dimensional origin of the Big Bang remains unestablished.
