# Retained pre-freeze test correction

The first new adapter test stopped before the full fabricated benchmark. The constant-source reference test incorrectly demanded that the L1 difference between two complex centers fit inside a proved Euclidean disk radius. For omega=1/2^39 and width1/64, the Mexp disk radius was about1.9121849353602935e−160, while the L1 center difference was about1.912184935360307e−160. Each component, with the independent reference remainder added, was correctly enclosed. The nested moment had the same norm mismatch.

The test was corrected to require containment on each real and imaginary axis, matching the exported rectangle and the historical independent reference test. The source/solver implementation and its rigorously computed radius were unchanged by that correction. This is a fabricated-test norm error, not a physical result or a relaxation of the1e−20 actual exported-L1 gate. Subsequent complete fabricated rehearsals test that actual sum of halfwidths independently. No physical source or retained array was evaluated during the failed or corrected test.

The later per-whole bookkeeping and schema additions were followed by two further complete fabricated rehearsals for the final source bytes. Earlier successful receipts are retained as development history; the final review must pin the matching final-code receipts.
