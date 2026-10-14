# Primary sources for the short Einstein–scalar Cauchy construction

14 October 2026. Focused source audit; no timelike-boundary theorem is used.

## Hughes–Kato–Marsden: finite Sobolev Cauchy theorem

T. J. R. Hughes, T. Kato and J. E. Marsden, **“Well-posed quasi-linear second-order hyperbolic systems with applications to nonlinear elastodynamics and general relativity,”** Archive for Rational Mechanics and Analysis **63** (1977), 273–294. [DOI](https://doi.org/10.1007/BF00251584). [Primary author PDF](https://www.cds.caltech.edu/~marsden/bib/1977/02-HuKaMa1977/HuKaMa1977.pdf). Full PDF retrieved; local snapshot `cauchy-sources/HKM1977_PRIMARY.md`.

**Exact location:** §3.1, assumptions (a1)–(a3), equations (3.1)–(3.2), **Theorem III**.

The theorem treats quasilinear second-order systems in arbitrary spatial dimension `n`, with symmetric principal matrices, uniform hyperbolicity on an open coefficient domain, sufficiently smooth bounded coefficient derivatives, and an appropriate Sobolev source at zero. For data in `H^{r+1}×H^r`, it gives a common short interval for a neighborhood of initial data, uniqueness, continuity in the original data topology, and time regularity into lower Sobolev orders. It assumes `r>n/2+1`; if principal coefficients do not depend on derivatives of the unknown, `r>n/2` suffices. The theorem is reversible to a two-sided short interval.

The article's introduction explicitly warns that top-space flow maps for quasilinear hyperbolic equations need not be smooth. **It must not be cited as supplying C² parameter dependence in the unchanged top norm.** The proposed precontact result instead differentiates the equations, uses difference-quotient energy comparisons, and conservatively reserves four Sobolev orders for the actual C² flow statement. Formal first and second variations have higher energy regularity; that alone is not optimal-topology differentiability.

The primary theorem is written on `R^n`. For the present local construction, one may use its localized energy proof on a complete compact auxiliary `T^4`, or use a suitable smooth reference field and localized `R^4` extension. The compact construction avoids global nondecay assumptions on the metric and scalar potential.

## Ringström: Einstein–scalar reduction, causality and gauge propagation

Hans Ringström, **“Future stability of the Einstein-non-linear scalar field system,”** Inventiones Mathematicae **173** (2008), 123–208. [DOI](https://doi.org/10.1007/s00222-008-0117-y). [Primary author PDF](https://people.kth.se/~hansr/invgnls.pdf). Full PDF retrieved; local snapshot `cauchy-sources/RINGSTROM_PRIMARY.md`.

**Exact locations checked:**

- §1.3, equations **(12)–(13)**:
  `Ric(g)=dφ⊗dφ+[2/(n−1)]V(φ)g` and `□_gφ−V′(φ)=0`. The discussion explicitly allows `n≥3`; the five-dimensional model has `n=4`.
- §3.1 states a general smooth quasilinear wave Cauchy setup in arbitrary spatial dimension; **Theorem 5** supplies the smooth-data version.
- §§3.2–3.3 treat causality and energy uniqueness on domains of dependence; **Theorem 6** is the zero-Cauchy-data uniqueness statement on the domain of dependence. The construction uses this result for causal uniqueness, with local existence supplied separately by the Cauchy theorem.
- §4, equations **(55)–(57)**: contracted Bianchi plus the scalar equation give a homogeneous wave equation for the gauge deviation.
- §4.1 constructs harmonic-compatible initial gauge data. It expressly permits the Einstein constraints to hold only on a chosen initial subset; propagation is then restricted to its domain of dependence.
- §3.4, **Theorem 7**, supplies Cauchy-stability context for Einstein–scalar solutions, with smooth data and convergence in the indicated high norms.

The paper's global stability results impose special potential conditions `V(0)>0`, `V′(0)=0`, `V″(0)>0`. **Those assumptions and those global conclusions are not transferred to HDBLAST.** The local reduced wave system and homogeneous gauge-propagation identities require smooth coefficients, regular Lorentzian geometry and constraint-compatible initial data; the arbitrary smooth potential enters as lower-order terms. The finite-Sobolev local existence input is Theorem III of Hughes–Kato–Marsden.

## Precise hypotheses needed here

Use a fixed background-adapted generalized harmonic reference connection; the actual background must solve the same reduced system on the controlled neighborhood. Choose metric/scalar Cauchy fields in `H^s`, their velocities in `H^{s−1}`, `s≥14`, with a uniformly positive spatial metric and uniformly Lorentzian time slicing. The potential is smooth on an open neighborhood of the scalar range. The parameter curve of these data is smooth in those fixed spaces, and the local gauge construction agrees with the background on untouched collars.

Complete arbitrary reduced Cauchy extensions on a compact auxiliary `T^4` are allowed; constraints outside the physical controlled ball are not asserted. Recover the full Einstein–scalar equations only in the domain of dependence where the initial gauge vector, its normal derivative and the constraints vanish.

Choose a **uniform bound on actual perturbed causal speeds**, together with strict distance margins between the entire correction annulus, the inner protected region and the shell. A background-null-cone contact bound alone is insufficient.

No retrieved source is represented as a material thin-shell evolution theorem, or as establishing this model's evolution through contact.
