# Preserved fabricated preparation failures

No scientific or physical evaluations occurred in these attempts.

1. The initial fabricated test fixture used Arb strings of the form
   `7/11 +/- 1/1000000`; python-flint rejected that fractional-radius syntax.
   The fixture now constructs rational balls and inflates them explicitly.
2. The initial benchmark fixture passed Python Fraction directly to acb;
   python-flint rejected that type. It now constructs the exact-rational Arb
   enclosure explicitly before acb conversion.
3. Independent review found that Arb's global series cap10 silently truncated
   coefficients despite a constructor precision27. This was reproduced only
   with manufactured sources, then repaired by explicit cap27/set/restore and
   result-precision guards. Normal/optimized high-degree fabricated jet tests
   and50 independently computed exact coefficient checks pass after repair.

No old scientific source, registration, archived result or negative outcome
was altered. The new candidate remains outside the checkout and awaits a
complete prospective freeze before real source callbacks.
