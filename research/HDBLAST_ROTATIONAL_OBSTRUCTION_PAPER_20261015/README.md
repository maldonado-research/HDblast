# Rotational linearization obstructions: paper review draft

15 October 2026. **Research draft for specialist review; not submitted.**

[Read the 14-page manuscript](manuscript/ROTATIONAL_SHELL_OBSTRUCTIONS.pdf). Editable [LaTeX](manuscript/ROTATIONAL_SHELL_OBSTRUCTIONS.tex), [primary-source comparison](literature/PRIMARY_COMPARISON.md), and [reproducible algebra supplement](supplement/README.md) are included. The companion ZIP contains exactly the same frozen files, including the manifest.

The paper calculates a necessary rotational consistency condition for two Einstein–scalar interiors joined across a timelike shell with scalar-dependent tension. Its direct constraint identity retains the shell's area and boost-angle term. At an invariant background, an arbitrary second-order bulk or shell correction cannot cancel the quadratic obstruction. A reconstructed, compactly supported linear perturbation has a strictly nonzero value on a short slab before shell contact, so it cannot be tangent to a classical parameter-C² exact family under the stated assumptions.

The draft includes a fully specified exact background with smooth interiors, fixed scalar potential and fixed shell tension. This illustrative potential and tension **differ from the registered HDBLAST model**. It demonstrates that the theorem's hypotheses are realized; it does not establish a cosmological-origin mechanism, a new exact-solution priority claim, or long-time stability.

For one fixed-profile quadrupole family, the six rotational conditions reduce to `[A,B]=0`. The paper derives the regular and singular geometry of that finite equation, including a commuting-curve velocity criterion and exact finite inverse norm. This is a specialization of standard momentum-map geometry. It supplies no inverse for the gravitational PDE and no nonlinear spacetime sufficiency theorem.

## Contribution and review status

The general quadratic-obstruction and boundary-charge mechanisms are established theory. The contribution being assessed is the explicit physical-interface calculation and smooth obstructed witness. Close primary comparisons include Harlow–Wu (2020), De Vuyst et al. (2025), Khavkine (2015), and the physical-shell Hamiltonian of Hájíček–Kijowski with its erratum. The comparison records exact inspected passages and unresolved historical-source access gaps; no exhaustive novelty determination is claimed.

The `review/` notes are independent internal AI-assisted derivations and checks of the included source, **not external human peer review**. They found no blocking defect in their recorded scope. Correct finite checks and internal review do not establish originality, significance, or acceptance. [Questions for specialist review](REVIEW_QUESTIONS.md) identify what remains.

## Reproduce and build

Follow [REPRODUCE.md](REPRODUCE.md) to authenticate the complete source inventory and replay the exact checks in ordinary and optimized Python, including deliberate failure controls. The top-level manifest inventories the runner and input-pin file as well as the scientific source and expected reports. The workflow stores fresh execution logs separately.

To rebuild the PDF, copy `manuscript/` to a writable build directory and run:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error ROTATIONAL_SHELL_OBSTRUCTIONS.tex
```

The checked build uses pdfTeX 1.40.26 / TeX Live 2025/dev (Debian), BibTeX and standard LaTeX packages. PDF rendering metadata can differ across builds; source authentication and mathematical replay are distinct from bitwise PDF rebuilding.

All previous repository checkpoints are preserved. This capsule is a focused paper draft, not a replacement of their recorded scopes or receipts. No journal submission or archival-version publication accompanies it. The manuscript describes substantive AI assistance; the human author must review its mathematics, references, interpretation and declarations before submission.
