# Pure-theory preparation history

No source callback, scientific array, physical trajectory or source integral
was evaluated in any version.

- `development/first_symbolic_pass/` preserves the first successful version:
  50 checks and 12 mutations. Its constant-U state counterexample illustrated
  unrestricted amplitudes but did not preserve the linear Wronskian.
- `development/before_interface_note/` preserves the next successful version:
  55 checks and 12 mutations, normal and optimized. It replaces that
  counterexample with a Wronskian-preserving Bogoliubov state supported between
  the probes, and explicitly strengthens the integrability hypotheses.
- The derivation adds the panel-interface note and the positive
  finite-band hidden-state density formula. Its verifier is unchanged from
  the 55-check version; the initial `SYMBOLIC_RELEASE_*.json` receipts bind
  that documentation hash.
- A final review adds the fixed-plane-wave phase to the explicit normalized
  Bogoliubov coefficient: `beta=epsilon A(k)exp(-2ik a)`. Fresh
  `SYMBOLIC_RELEASE_PHASE_*.json` receipts bind the final documentation.

The former top-level `SYMBOLIC_NORMAL.json` and
`SYMBOLIC_FINAL_*.json` remain historical evidence, with their exact source and
documentation copies preserved in the above development folders. They do not
bind the final documentation bytes. Only the release-phase receipts do.
