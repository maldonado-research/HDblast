# Internal review record

Separate agent processes independently reviewed the mathematics, Cauchy
pullback, citations and scope. This is internal AI-assisted review, not
external human peer review or journal acceptance.

The earlier and final mathematical reviews identify their exact historical
manuscript snapshots, which are included here. The format-delta record
covers the nonconstant-background wording and equation wrapping. The final
literature review also checks the last abstract clarification: smallness
means spatial support, not perturbation amplitude. Equations and Appendix A
are unchanged by that clarification. The final build/visual receipts pin
the delivered paper.

`ORIGIN_AND_TARGETS.json` maps authoring sources to public copies. Historical
commands in adapted JSON records are evidence, not portable invocations.
The main replay covers 92 identities and 20 wrong-formula controls. The
additional appendix review has its own nine identities and four controls,
with source and saved outputs supplied here; these are not added to the
main replay count. To inspect that supplemental algebra, run
`python check_final_appendix_independently.py` from this directory in the
recorded environment. The main runner does not execute reviewer scripts.

The Cauchy review's target hash identifies the original boundary note; a
single companion-location correction in the public copy is recorded in the
origin map. No equation or proof changed in that adaptation.
