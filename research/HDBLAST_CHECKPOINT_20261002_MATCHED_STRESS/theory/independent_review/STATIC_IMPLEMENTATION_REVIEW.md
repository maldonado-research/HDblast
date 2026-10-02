# Static matched-stress implementation review

Status: no material blocker found in the five source snapshots below. This is an internal, pre-freeze static review, not an external peer review or a numerical validation. The public-freeze SHA, full registration, final validator and numerical outcome remain outside this review's signoff.

Review date: 2026-10-02 UTC. All work was read-only on the public checkout. No source/exponential samples, forced modes, response integrals or new physical numerical tests were evaluated. AST parsing and exact integer/rational polynomial algebra were allowed and used. There were 41 exact polynomial, critical-root-count, rational-bracket and contact-factor checks; reviewed implementation modules were not imported.

## Inspected snapshots

Paths are relative to `/workspace/hdblast-research-work/`.

| File | SHA-256 |
| --- | --- |
| `stress-independent/forced_stress.py` | `aae9a60e17a3586ca5a2213dce4dd8f9122de9f463df0cb82e0d1a90218d6ffc` |
| `stress-primary/code/stress_primary.py` | `8cd3270f6d38114c9ca363131aea34eb4a543652d33f404a41bbfe6f29dc3cdd` |
| `stress-primary/code/source_jet.py` | `6428f404032acf4824f93d38a4467e495a56b06434e796f5ac338f3addee0205` |
| `stress-theory/stress_tail_bounds.py` | `47c6f9d047c2fe5da9edf47ff17ab6aa3df74f30ceebaeb0723554f3128b9ad2` |
| `stress-theory/pulse-norms/pulse_derivative_budgets.py` | `fa37e7a4d0b0ee09400d17b8a9689e0bcb7c100be89bf5a55cfc755ccba120a0` |

The independent implementation's earlier inspected snapshot was `1517117f63e5095865d0852654247d5b2b72267267132f5f9d24fbf5aa3c626c`; its updated CLI/output and refinement-reporting sections were reread before recording the final hash above.

## Findings

1. **Independent minimal stress.** `forced_stress.py:140` constructs both density and pressure from the actual complex baseline mode, forced variation and their derivatives, including the opposite-sign explicit mass contacts. Pressure is not defined from density, the trace or the Ward identity. Complete variations of W2 and W4 enter both subtraction terms. The anomaly is computed from their unsimplified combination with the variance subtraction and its derivatives. The second variance derivative uses `Re(2ik*w-s)/k`, the actual forced equation; the Wronskian is monitored without projecting the evolution onto its constraint.

2. **Ward identity at finite comoving cutoff.** The independent history ledger integrates `R' = a^4 Q0,K d'/(2 epsilon) + L R - 3 L P`, where `R=a^4 delta_rho/epsilon`, `P=a^4 delta_p/epsilon`, and `L=aH`. It uses the finite-cutoff physical baseline `Q0,K`, not its continuum limit. It compares against directly evaluated stress and defines neither stress. The primary baseline and baseline derivative, density contact and pressure/anomaly contacts agree with the finite-K identities in `STATIC_STRESS_AUDIT.md`. The primary density-derivative output denotes `a^4 delta_rho'/epsilon`; it correctly subtracts `4 L R` from the derivative of the scaled density.

3. **Variance derivatives and local terms.** `stress_primary.py:79-205` includes the derivatives of the time-dependent logarithmic contact. Its finite-K primitive satisfies `A'=-L v^3` and `A''=L^2 v^3(2-3v^2)`, with `v=K/sqrt(K^2+2a^2)`. The finite baseline, its derivative, and the `q+Q0,K*f/4` current retain their required local terms. Density, pressure and anomaly normalizations agree with the direct mode route; the error estimates use absolute linear coefficients.

4. **Separate derivative and stress tails.** `stress_tail_bounds.py` carries bounds for `q`, `q'`, `q''`, the baseline, and the finite local contacts into the density, pressure and current bounds. It does not substitute a variance-only tail for a stress bound. Exact checks confirmed the factorizations of the second-derivative contact, baseline tail and J7/J9 moment tails. The pressure contact bound follows from the algebraic contact identity only; it is not the pressure producer. Interval evaluation and outward conversion are retained.

5. **Source jets and pulse derivative budgets.** The primary polynomial tables agree exactly with the independent integer recurrence for B derivatives through order five; the signed uB source follows the Leibniz rule. The derivative-budget critical polynomials have the stated counts of roots in `(0,1)`, the listed rational brackets isolate them, and the zero-root parity factors and negative-root ordering are accounted for. This supports the registered total-variation construction without evaluating a source or a norm numerically.

6. **CLI and failure evidence.** Both producer CLIs require public-freeze/provenance pins and fresh external output directories. The primary checks the full-registration hash and registered file hashes, pins runtime versions, enforces deadlines around each quadrature and records partial progress and failures. The independent route verifies its manifest and inputs, pins Python/NumPy and extended precision, monitors the registered time/memory limits, records the active run and captures available failure arrays. These are statically inspected controls; actual runtime compliance and archive integrity must be checked after execution.

## Final integration condition

`stress_tail_bounds.py:19-25` loads `pulse_derivative_budgets.py` from either the theory directory or its `pulse-norms` subdirectory. `stress_primary.py:61` names the tail driver in its minimal required manifest set but does not separately require the pulse module. The final `FULL_REGISTRATION.json` must bind the actual loaded pulse module and avoid an unregistered sibling shadowing it. Registering the complete assembled source tree, including the intended `theory/pulse-norms/pulse_derivative_budgets.py`, and confirming no sibling exists satisfies this condition. The final independent manifest must likewise reflect its updated source hash.

No claim is made about passing numerical gates, physical stability, particle yield, heating, nonlinear evolution, or a discovery beyond this fixed-geometry linear-response calibration.
