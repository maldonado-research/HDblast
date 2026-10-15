# Next reconstruction: a verified feasible wider state

15 October 2026 UTC. Separate feasibility and work-planning note.

The frozen waveform checkpoint uses ell=10^-10. Independent Astra and Sol
calculations now justify ell=10^-6 as a concrete next width on the **same**
retained background. This enlarges the pulse width and required collar by
10,000 in model coordinates. The new calculations are saved separately under
`planning/`; the frozen manuscript, figure, CSV, response interval and replay
remain unchanged.

The same exact shell box and outer rectangle support the strict first-exit
collar Z=8ell=8e-6. The fixed rectangle's interval argument has strict depth
threshold

\[
L_* = \frac{5669686020285051}{123672167032900272500}
\simeq4.5844478643\times10^{-5}.
\]

It certifies L<L_*. This is a sufficient threshold of that particular
rectangle and method, not the actual background's maximal existence depth.
The existing analytic branch theorem remains an input.

The independent interval evaluations additionally justify the rounded bounds
m=1/500, M=102000, B=250 and Gamma=250 on this wider collar. For the proposed
width and T=3ell, the existing characteristic/kernel estimate then gives

\[
\|H-\widehat H\|_\infty<0.000500684,\qquad
\sup_{2\le x\le3}|H(x)|<9.325\times10^{-7},
\]
\[
\|H-\widehat H\|_{L^2(0,3)}<0.000500685,
\quad\widehat H=Q'+\ell e^{-\ell x}Q.
\]

These are conditional tolerances for a proposed wider reconstruction. The
first-exit margins, coefficient boxes and majorants are exact-rational
calculations. Astra and Sol use different potential interval evaluations and
agree on the stated rounded bounds; their full intermediate boxes need not
be identical. Neither evaluation imports or executes the archived producer.

The next calculation should instantiate this width, retain a strict causal
margin for the full nonlinear correction, and export a fresh continuous
waveform enclosure with its own response integral and peak exclusions. The
narrow state's coefficient interval and peak windows cannot be transferred.
The signed physical gain still requires its declared normalization and proper
time. A controlled finite-amplitude prediction further requires an effective
nonlinear remainder and a nonzero quantified epsilon interval, including the
localized inverse, gauge and evolution estimates.

See [Sol's independent proof and exact checks](planning/wider-state-sol/WIDER_STATE_FEASIBILITY.md)
and the separately supplied Astra feasibility certificate. Their checkers use
an explicit `--source-root` pointing to the unchanged extracted
`HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015` capsule for hash-pinned
mathematical inputs; planning files are distributed outside that capsule.
The individual notes give commands and actual check counts.

This planning deduction supplies no completed wider waveform, transferred
response interval, finite-amplitude nonlinear solution, observational scale
or cosmological-origin evidence. Internal independent review and exact
arithmetic do not establish external novelty or formally verify the analytic
PDE/background hypotheses. No additional Zenodo version or APS submission is
associated with it.
