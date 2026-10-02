# Original registered metric calibration: numerical procedure failure

The original `HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE` calibration has
status **FAIL**. Its primary producer stopped on a SciPy `IntegrationWarning`
that the frozen procedure treats as failure. The original source, gates,
partial output, warning traceback and exit code remain unchanged. This is a
failure to complete and validate the registered numerical calculation; it is
not evidence of physical instability or rejection of the HDBLAST hypothesis.

The public freeze is `57668b8fadd75df8738565e0bbd1eb852c1ebae8`. Its
registration SHA-256 is
`f8d6bbd17b540b46c5e4382ab21fd74e15f0b283956991a250f49994553a476d`,
and its independent manifest SHA-256 is
`4760f2a7e864dbe81ad8911797499dcd9fcc9e37a45421a607a101d35e6ce507`.
The saved replay directory is
`/workspace/hdblast-research-work/metric-original-complete-replay-001`.

## What actually ran

The replay started at 19:06:04.937813 UTC on 2 October 2026 and stopped at
19:07:01.503822 UTC. It planned 32 commands and attempted 27. The first
**26 commands passed**; all were exact algebra, synthetic guards, input
verification, derivative certification, or pure preflight. The 27th command
was the first physical primary producer and exited 1. It did not time out.
The process elapsed time was 0.565960364 seconds; its internal failure
capture recorded 0.014436731 seconds of producer work after initialization.
These are different clocks, not inconsistent runtime estimates.

The passing prerequisites include both normal and optimized primary contact
algebra, primary preflight, action/contact identities, Ward/trace algebra,
independent general WKB identities, independent contact-inventory agreement,
generic-CSE specialization agreement, actual-root prerequisite algebra, tail
algebra, synthetic validator guards, raw baseline algebra and inherited
input pins. The independent producer's preflight and the exact derivative
certificate also passed. These prerequisite passes do not make the physical
calibration pass.

The primary saved one completed source-free row:

```
source=positive_B, eta=-5.5,
continuum and finite K=64,128,256:
q=q_prime=q_second=rho=p=current=0.
```

The baseline Q0, rho0 and p0 values are retained separately. This row precedes
the pulse support `[-5,-3]` and is a trivial zero-response control. It is
not a nonzero response or independent mode validation. The active progress
marker then identified `positive_B, eta=-4.5`.

The producer failed while evaluating that row's continuum logarithmic
history integral using `scipy.integrate.quad` with `weight='alg-logb'`.
The frozen configuration was `epsabs=epsrel=1e-12`, `limit=300`, and
`warning_is_failure=true`. SciPy reported:

> The occurrence of roundoff error is detected, which prevents the requested
> tolerance from being achieved. The error may be underestimated.

The traceback reaches `continuum_q_jet -> log_history -> integrate -> quad`.
It does **not** record which of the derivative orders 1, 2 or 3 triggered
the warning. No specific derivative order or numerical error magnitude may
be inferred from this capture. No completed nonzero continuum or finite-K
row was saved.

## What did not run

There is no completed `fresh/primary/results.json`, and no independent
physical output directory. The independent forced-mode producer, normal and
optimized physical validators, summary renderer and figure renderer were
never attempted. Thus there are no original independent raw mode archives,
nonzero cross-route comparisons, numerical Ward ledgers, refinement results,
continuum-tail comparisons, or physical control results to classify.

The original manifest and replay report verify source hashes before and
after the attempt. This read-only audit additionally rechecked all **265
frozen files** and all **26 captured fresh-output hashes**, including failure,
partial results and proof captures. `ORIGINAL_FAILURE_AUDIT.json` pins the
receipts and raw log files and records these exact coverage facts.

## Scientific interpretation

The numerical implementation could not deliver a continuum quadrature that
met its declared stopping rule. Treating the warning as failure preserves
the prospective procedure. The record neither demonstrates divergence of
the renormalized response nor distinguishes a difficult integrand from
floating-point cancellation, derivative evaluation error or the quadrature
engine's roundoff estimator. It gives no evidence of physical instability,
stability, heating, a shifted-root propagator response, a completed metric
response matrix, coupled bulk/shell evolution or a cosmological discovery.

The original experiment remains failed even if a separately registered
followup later passes. Suppressing its warning, relaxing its tolerance,
silently changing arithmetic, or replacing its output would erase the
meaning of the original registration. Preserve this entire record and
publish any subsequent numerical method under a distinct checkpoint and
new prior public freeze.
