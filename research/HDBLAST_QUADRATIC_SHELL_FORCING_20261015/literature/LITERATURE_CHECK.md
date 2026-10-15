# Bounded literature check for quadratic shell forcing

Prepared 15 October 2026 UTC. The subject is the distinction between the explicit material tension contribution to a second scalar junction jet and the complete source after expressing the junction and observable in bulk coordinates. This check supplies geometric bookkeeping and source provenance. It establishes no novelty, nonlinear existence theorem, effective amplitude radius, or observational consequence.

## Result and model convention

The registered two-copy scalar law is taken in its stated outward-normal convention,
\[
 \sum_i p_i+\sigma'(\varphi_\Sigma)=0,\qquad
 p_i=(n_i\phi_i)|_\Sigma.
\]
For a smooth family, use derivative jets (this notation does not establish that such a nonlinear family exists):
\[
 \varphi_\Sigma(\varepsilon)
 =\varphi_{\Sigma,0}+\varepsilon\psi_1
   +\tfrac12\varepsilon^2\psi_2+O(\varepsilon^3),
 \qquad
 p_i(\varepsilon)=p_{i,0}+\varepsilon p_{i,1}
   +\tfrac12\varepsilon^2p_{i,2}+O(\varepsilon^3).
\]
Twice differentiating this **material law** gives the elementary chain-rule identity
\[
 \sum_i p_{i,2}+\sigma''(\varphi_{\Sigma,0})\psi_2
 =-\sigma'''(\varphi_{\Sigma,0})\psi_1^2=-4\psi_1^2.
\]
The last equality uses the registered model value \(\sigma'''=4\); the literature does not supply that model value. In a coefficient convention without the Taylor factorial, the corresponding second coefficients and factors must be converted.

This isolates the explicit nonlinear tension term. Each \(p_{i,2}\) is already the jet of a normalized-normal derivative evaluated on the moving shell. It contains more than a fixed-coordinate second scalar derivative. Converting the whole law to bulk-coordinate unknowns requires the normal, metric/lapse, shell embedding, evaluation, and material parametrization contributions in the same convention. The displayed right side therefore cannot by itself be identified with the full bulk-coordinate forcing or a measured second-order waveform. The primary sources below support this distinction; they do not calculate the missing HDBLAST terms or establish cancellation or noncancellation of the complete source.

## Direct primary evidence

**Marc Mars, _First and Second Order Perturbations of Hypersurfaces_, arXiv:gr-qc/0507005.** [Record](https://arxiv.org/abs/gr-qc/0507005), [primary text](https://arxiv.org/html/gr-qc/0507005). Canonical metadata identifies the author and first posting on 1 July 2005. The retrieved title passage identifies v1, 1 July 2005; the index's 2016 update field is not a verified new arXiv version.

The targeted read returned Section 5, Proposition 2, equations (26)–(27). With first embedding deformation \(Z_1=Q_1n+T_1\), its induced-metric variation contains the tangential Lie derivative, \(2Q_1\kappa\), and the pulled-back spacetime metric variation. Its extrinsic-curvature variation also contains normal deformation derivatives and metric/normal terms. These passages directly show why fixing the coordinate location of a shell does not permit deleting its geometric deformation contributions.

The same read returned the explanation that the second deformation vector is the covariant acceleration of the embedding curve, the end of the second-order extrinsic-curvature formula (41), and Section 7's explicit distinction between ambient-spacetime identification freedom and hypersurface identification freedom. Section 6, Proposition 3 was identified in the returned text, but the complete second-order formulas and their appendices were **not** retrieved or audited here. The paper treats non-degenerate hypersurfaces and states arbitrary dimension and signature. Its matching application is geometric; this check imports no matter-shell evolution or remainder theorem.

**Marco Bruni, Sabino Matarrese, Silvia Mollerach and Sebastiano Sonego, _Perturbations of spacetime: gauge transformations and gauge invariance at second order and beyond_, arXiv:gr-qc/9609040.** [Record](https://arxiv.org/abs/gr-qc/9609040), [primary text](https://arxiv.org/html/gr-qc/9609040). Canonical metadata identifies the authors and first posting on 13 September 1996. A specific passage-version header was not returned.

The read returned equations (1.1)–(1.3) and Section 4.2, Proposition 2, equation (4.9). In the paper's pullback direction,
\[
 \delta T^Y-\delta T^X=\mathcal L_{\xi_1}T_0,\qquad
 \delta^2T^Y-\delta^2T^X
 =2\mathcal L_{\xi_1}\delta T^X+
   (\mathcal L_{\xi_2}+\mathcal L_{\xi_1}^2)T_0.
\]
The retrieved expansion (4.7) explicitly uses the second-order factor \(1/2\). These are general tensorial transformation rules. They explain why a second-order reconstruction mixes a first-order perturbation with a first-order coordinate change. Their use here is geometric; cosmological examples in the paper provide no cosmological evidence for this model.

For example, elementary scalar composition with
\(t_\varepsilon(\tau)=\tau+\varepsilon a_1+\varepsilon^2a_2/2+\cdots\) yields
\[
 (F_\varepsilon\circ t_\varepsilon)_2
 =F_2+2a_1\dot F_1+a_1^2\ddot F_0+a_2\dot F_0.
\]
This example describes the type of clock/pullback terms to retain. It is not a computation of the actual HDBLAST clock or an assertion that any displayed background derivative is nonzero.

**Shinji Mukohyama, _Gauge-invariant gravitational perturbations of maximally symmetric spacetimes_, arXiv:hep-th/0004067.** [Record](https://arxiv.org/abs/hep-th/0004067), [primary text](https://arxiv.org/html/hep-th/0004067). Canonical metadata identifies the author and first posting on 10 April 2000. A passage-version header was not returned.

This identifier is **not** the doubly gauge-invariant junction paper suggested in the initial lead. The actual returned Section VI discussion, equations (67)–(71), nevertheless supplies a useful first-order example: changing from a Gaussian description near the brane to the generalized Regge–Wheeler gauge generally produces nonzero coordinate brane displacement. The same discussion says simplification of Israel junction conditions is a subsequent step. Section II assumes a maximally symmetric, constant-curvature bulk. We use this paper only for its demonstrated gauge/displacement caution. Its linear master equations and background assumptions are not transferred to the present matter-supported bulk.

All three entries have **canonical metadata plus targeted primary-passage** evidence, not full-paper review.

## Prior work reused and recent screen

The earlier full-contact literature note and scalar-contact note were inspected before new discovery. They already retain bulk and boundary quadratic sources, displacement/evaluation corrections, and convention conversion, using **Kudoh and Tanaka, _Second order perturbations in the radius stabilized Randall-Sundrum two branes model_, arXiv:hep-th/0112013v3**, [primary text](https://arxiv.org/html/hep-th/0112013). Their recorded locators are equations (27)–(30) and (35)–(38). That source was not reread in this check. Its static Minkowski-sliced two-brane setup and coefficient convention remain limitations. The full-contact note also mentions Mars, so the present passage check refines existing context rather than establishing priority.

The frozen inertial-waveform literature packet was also inspected. Its energy/Volterra comparisons do not furnish a second-parameter remainder for the actual Einstein–scalar shell, and they already identify intrinsic clock, embedding and physical pullback as parts of an observable estimate. Its four primary papers were not fetched again.

Two new paper-index searches returned five results each. One search had a 2020–15 October 2026 date filter. Recent candidates were screened from **index abstracts only**:

| Candidate | Retrieved scope and disposition |
|---|---|
| [_Variational formalism for generic shells in general relativity_, arXiv:2203.03049](https://arxiv.org/abs/2203.03049) | Generic-signature shell variational laws and regularization equivalence. A potentially useful variational crosscheck, with no verified second-jet or material-clock result from a body read here. |
| [_Junction Conditions for General Gravitational Theories_, arXiv:2603.04645](https://arxiv.org/abs/2603.04645) | Distributional junction conditions for curvature-invariant gravity actions. Its higher-derivative setting does not supply the present canonical Einstein–scalar second-variation source. |

Authors, exact posting dates, versions and body claims for these two recent candidates were not independently inspected. Other returns included null scalar-tensor shells, modified-gravity junctions and a cosmological effective stress-tensor study. They were not pursued. The recent screen established no directly useful replacement for the general geometric bookkeeping above. This is a limit of the returned results, not a claim that no relevant recent paper exists.

## Coverage and publication record

The Firecrawl skill and its paper-research reference were used. This check made exactly **two paper-index searches, three canonical metadata inspections, and three targeted primary-passage reads**; all eight operations returned successfully. The reads returned thirteen passages in total. There was no new web-search call, scrape, broad crawl, citation-graph expansion, or full-paper fetch. A previously recorded web-search HTTP 402 credit failure was inspected in the frozen receipt and was not retried; billing was not changed.

[LITERATURE_LEDGER.json](LITERATURE_LEDGER.json) records exact arguments, returned identifiers, metadata, evidence levels, read limits, local prior-artifact hashes and internal response hashes. It and this synthesized note are public-safe outputs. Raw responses under `_internal/` are internal reference material and must be excluded from any public packet. No raw primary text is republished by these two outputs.

The remaining model work is to compute and verify the complete quadratic coordinate source and the observable's embedding/clock/pullback conversion from the same first-order solution. The checked literature does not make those calculations optional, certify a nonlinear remainder, or establish external novelty, observational support, or a cosmological connection.

