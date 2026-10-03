# Post hoc trajectory audit protocol

Prepared on 1 October 2026 before running the new raw-array audit. The threshold is a new diagnostic selected after examining the older summary points. It is not a preregistered replacement for September B3.

At each saved record define the linear radiation contribution Lrad = rb² sigma(phi) R/18, quadratic contribution Qrad=rb² R²/36, and the full absolute fraction

F_abs = Lrad / [Lrad + Qrad + |vac| + |Wy| + vh²/12 + |sigma_phi Y rb vh|/24 + Y² vh²/48],

where vh=v/H0, H0=1/rb and the archived transfer closure is j=Yv. All terms are normalized by H0². Compute scalar-free F0 by omitting the final three nonnegative scalar/source terms. Both require Lrad>0.

Proposed gate: F_abs>=0.9, |Wy|/Lrad<=0.03, H/H0>0, finite data and the baseline reliability policy. Report the longest contiguous sampled expanding interval as DeltaN=ln(a_end)-ln(a_start). Do not add separate islands. F0 is an upper-bound gate: failure rules out F_abs at that sample; passage is not a full result.

Reconstruct restart history by parent records strictly before the restart time followed by the continuation, with all joins disclosed. Refuse guessed files or pickle-enabled loads. Retain source hashes, record counts, sample spacing and reasons for any omitted run. Only trust intervals on each run's allowed reliable segment; unreliable cases remain labelled.

The raw report also checks radiation/kinetic/friction normalization and Friedmann closure. Any source ledger is diagnostic, not an action-derived quantum consistency certificate. This audit reads existing data and does not rerun the field equations.

An interval found between sampled records is not automatically continuum certified. Report sampled duration and resolution comparison separately; do not use undocumented interpolation to claim a precise transition.
