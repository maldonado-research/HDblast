# HDBLAST — checkpoint of 27 September 2026: public summary

*Ricardo Maldonado's HDBLAST research program. Prepared with AI assistance. This is independent research. It has not been peer reviewed or published in a journal.*

## The idea being tested

HDBLAST asks whether our Big Bang could have been sparked by a *higher-dimensional blast*: a gravitational event in a five-dimensional spacetime, with our universe living on a four-dimensional "shell" (a brane) inside it. The program works with one fixed, pre-registered mathematical model and tests it with Einstein's equations. It does not adjust the model to fit the answer.

## What we found this round

- **The end state is stable, but empty.** In earlier work, the shell starts in an unstable state and rolls away from it. Last week we found the static state it appears to head towards: a smoothly expanding (de Sitter) brane. This round, two independent calculations show that small disturbances of that end state die away rather than grow, in the kinds of disturbance we could test. Whether the evolution actually reaches this state, and what large disturbances do, is not yet known. *Status: numerical result, checked by an independent audit. Not a mathematical proof.*
- **No hot Big Bang yet.** Our universe had a hot, radiation-filled early phase. We computed particle production along the simulated evolution and screened six other ways to make heat. In the model as registered, none of them produces a radiation-dominated era. The main obstacle is leftover vacuum energy on the final brane: radiation would need to dominate it for about 22 e-folds, and the best available energy budget allows about 0.6. *Status: negative, under stated assumptions.*
- **Two open leads, both requiring model changes.** Adding one extra term to the shell's tension can cancel the leftover vacuum. In a pilot run (made with a larger model parameter than the registered one, to save computing time), radiation then becomes a large part of the expansion. However, "dark radiation" coming from the fifth dimension stays above what observations of the early universe allow, at least so far. A second, newly found starting configuration has not yet been followed to its end. *Status: inconclusive.*
- **Exact mathematics.** We derived the detailed structure of the stable end state exactly, to eighth order in the model's small parameter (the audit added the ninth). We also explained where its special polynomial form comes from. *Status: exact, independently re-derived.*
- **Honest corrections.** The independent audits overturned three claims made during this round, and we corrected them. They also tightened several others. All are listed in `CLAIM_REVIEW.md`.

## What this does *not* show

- It does not show that a five-dimensional event caused our Big Bang.
- It does not provide a hot Big Bang mechanism.
- It provides no observational evidence for the hypothesis.
- It claims no discovery or proof.

The broad idea of a universe born on a wall in five dimensions has substantial published prior art (for example brane-world creation, ekpyrotic collisions and "dark bubble" cosmology). What is specific to HDBLAST is its registered model and its computed results.

## Where this leaves the hypothesis

In its registered form, the five-dimensional event ends in a stable but empty expanding universe. That is self-consistent, but it is not our universe. To remain viable, the hypothesis needs an added ingredient or a model change that produces a hot phase while keeping dark radiation within observational limits. The next calculations are designed to decide this either way (`NEXT_TESTS.md`). A clear negative result would also be a useful outcome.

## Observational status

The five-dimensional model does not yet predict anything current data can test: it has no radiation era and no spectrum of primordial fluctuations. The program's separate pulsar-timing "knee" template remains untested at the level of a full likelihood analysis. Its frequency sits at 1/year, where pulsar timing is known to have a blind spot, and future analyses must account for this.

## How to check our work

Every number comes from a script with machine-readable output. `REPRODUCE.md` gives the exact commands, and `MANIFEST.sha256.json` fingerprints every file. Each part of the work has an independent audit report (`VERIFICATION.md`). Start with `00_READ_FIRST.md`.
