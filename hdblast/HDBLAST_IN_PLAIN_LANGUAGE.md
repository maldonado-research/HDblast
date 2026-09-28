# HDBLAST in plain language (September 2026)

**By Ricardo Maldonado, independent researcher.** This updates the
[July 2026 public guide](guides/HDBLAST_GENERAL_PUBLIC_GUIDE_V20.md) with the September 2026 results.

## The question

> Could a powerful event in an extra direction of space have transferred energy into our
> universe and helped create its hot beginning?

The standard Big Bang model describes, very successfully, what happened **after** the universe
was already hot and expanding: cooling, the first elements, the cosmic microwave background,
stars and galaxies. HDBLAST does not challenge any of that. It asks a speculative question about
what might have **prepared** that hot starting state.

## The picture

Imagine creatures living on a sheet of paper. They know only the directions across the sheet. If
something outside the sheet shook it, they would feel the effect without being able to see its
cause. HDBLAST pictures our universe as a similar "sheet", a **brane** or shell, inside a larger
five-dimensional space. Physicists have studied such brane-world models seriously since the
late 1990s (Randall and Sundrum, among others). Using them does not make HDBLAST true, but the
ingredients are legitimate physics.

## What was done in 2026

**1. The model was pinned down before testing it.** One precise five-dimensional model was
"registered": its equations and numbers were fixed in advance, so later results could not be
tuned to look good.

**2. Inside the model, the shell solution was shown mathematically to exist, and to be unstable.**
A computer-assisted proof showed that the registered model has a shell solution. This is a
statement about the equations, not evidence that such a shell exists in nature. A second,
computer-assisted calculation showed it sits like a ball on top of a hill: it has a way to roll
off, with the roll-off rate pinned to five digits, and numerical checks found no other such
direction. The other main kinds of wobble (gravitational-wave-like, vector, and a special class of
scalar modes) were checked exactly and have no unstable mode, though full stability has not been
proven.

**3. A simple formula predicted the five-dimensional behavior.** A closed-form
four-dimensional formula matched the full five-dimensional calculation to about one part in a
million, and it predicted where the shell would turn around **before** the big simulation
confirmed it.

**4. The "blast" was simulated, and the simplest version failed.** Full five-dimensional
computer simulations of the roll-off show **two fates**. On one side, the shell heads toward an
empty universe that keeps expanding at an accelerating rate (the simulations show it approaching
this state but cannot follow it to the end). On the other, it reverses and collapses. **Neither
produces a hot, radiation-filled universe.** This is an important negative result: the shell's
own smooth (homogeneous, classical) motion, with nothing else added, is not enough to make a hot
Big Bang.

**5. An endpoint was corrected, and matter was added.** A proposed late-time endpoint turned out
to break one of the model's boundary conditions, and a corrected solution was found. A consistent
way to add a new matter field to the shell was derived as a proposed extension. It includes an
exact energy-exchange equation, so the blast's energy could in principle be passed to particles;
no such particle production has been calculated yet.

**6. Mistakes were published.** In July 2026 an earlier claim that the calculations showed
emitted gravitational radiation was **withdrawn** after an audit found it had looked for a
wave in a direction the declared source did not actually produce. Later drafts were also corrected
where a script bug, insufficient numerical resolution or overconfident wording had produced wrong
numbers or claims. Publishing
corrections is part of doing science honestly.

## What is *not* established

- That extra dimensions exist.
- That any higher-dimensional event caused the Big Bang.
- A mechanism in this model that produces a hot, radiation-dominated early universe.
- Any observational signal. An earlier proposed pulsar-timing "knee" was screened against public
  data, but it was never derived from the five-dimensional model, and a strict registered test
  on pilot data failed.
- Peer review. The Zenodo records are public and dated, but they are not journal publications.

## What happens next

The current research round (from 27 Sept 2026) tests:

- whether the corrected late-time solution is stable;
- whether particle production during the blast can create a lasting radiation era;
- how to control the numerical errors in the long simulations;
- which **new ingredients**, still within Einstein's equations, could turn a higher-dimensional
  event into a hot Big Bang. Candidates include matter coupled to the scalar field,
  "dark radiation" from a black hole in the extra dimension, and bubble-universe constructions
  studied by other groups.

If none of these works, that result will be published too.

## How to judge this project

A fair summary: **HDBLAST is a legitimate, speculative, testable research hypothesis with an
unusually careful public record. It is not an established theory.** It would become significant
only after a working mechanism, a unique quantitative prediction, independent reproduction,
peer review, and observational confirmation. The reviews recorded so far are internal and were
run with AI assistants; experts are invited to check the work. The code and data for each step
are in [`checkpoints/`](checkpoints/).
