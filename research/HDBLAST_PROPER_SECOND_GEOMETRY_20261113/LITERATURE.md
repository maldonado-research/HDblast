# Focused proper-geometry literature and method update

Scientific search through **13 November 2026 UTC**; retrieval **13 November 2026 UTC**. This focused note extends the earlier focused reviews without changing their **6 November 2026** scientific cutoff or any frozen capsule. It addresses the intrinsic Ricci/observer interface in the [proper-geometry continuation](CONTINUATION.md), rather than repeating the earlier high-order wave survey.

The usable conclusion is a direct background-connection calculation of the second intrinsic Ricci field, followed by a covariant first-observer correction and the retained Raychaudhuri transfer. The inspected recent methods support covariance and careful observer definitions. They supply neither the registered shell constants nor a uniform nonzero-amplitude remainder or physical confirmation.

## Evidence and exact coverage

| Primary source and version | Passages actually inspected | Use and limit |
| --- | --- | --- |
| Kouji Nakamura, *Second Order Gauge Invariant Perturbation Theory — Perturbative curvatures in the two-parameter case*, [arXiv:gr-qc/0410024v2](https://arxiv.org/html/gr-qc/0410024v2), **19 May 2006**. Related publication: Progress of Theoretical Physics 113 (2005), 481–511, [DOI](https://doi.org/10.1143/PTP.113.481). | §3, (38)–(44); §3.1, (45)–(53); §3.2, (81)–(83). The exact v2 HTML and relevant indexed passages were read. | Explicit connection-difference and second-curvature formulas on an arbitrary background. The primary metadata says v2 corrects typos in the published version; the mathematics here is pinned to that corrected preprint. |
| Markus B. Fröb and William C. C. Lima, *Synchronous coordinates and gauge-invariant observables in cosmological spacetimes*, [arXiv:2303.16218v2](https://arxiv.org/pdf/2303.16218v2), **4 October 2023**. Classical and Quantum Gravity 40 (2023), 215006, [DOI](https://doi.org/10.1088/1361-6382/acf98a). | Exact v2 PDF pages 1–12: §2, (1)–(23), especially (19)–(21); §3, (25)–(49). Also relevant indexed passages. The other 45 PDF pages were not inspected. | Defines proper-time free-falling observers and transported labels, and explicitly expands both coordinates and observables. Its flat-FLRW, past-infinity prescription differs from the registered spherical fixed anchor. Its perturbation coefficients use a different factorial convention. |
| Joon-Hwi Kim, *Geodesic deviation to all orders via a tangent bundle formalism*, [arXiv:2509.23600v1](https://arxiv.org/html/2509.23600v1), **28 September 2025**; primary arXiv metadata reports General Relativity and Gravitation 58, 75 (2026), [DOI](https://doi.org/10.1007/s10714-026-03580-9). | Introduction; “Geodesic Deviation in Tangent Bundle,” (1)–(8); covariant dressing, (11)–(13); Jacobi propagators, (14)–(21); all-orders GDE, (26)–(31); Summary. | Useful recent covariant transport organization. The setup explicitly assumes a **real-analytic** manifold and varies separation in one fixed metric. It does not turn a merely C2 amplitude family of varying metrics into an all-orders or quantitatively controlled amplitude expansion. Journal text and exact publication day were not separately inspected. |
| Kirill Krasnov and Adam Shaw, *First-Order Structure of the Linearised Einstein Equations: Constant-Curvature Backgrounds*, [arXiv:2610.32692v1](https://arxiv.org/html/2610.32692v1), **23 October 2026**. | §2.1 background condition (12), §3 fields and curvature (13)–(14), and Introduction/Table 1. Primary metadata and two in-body passages were read. | A recent first-order connection factorization on four-dimensional constant-curvature backgrounds with nonzero cosmological constant. The passage after (14) explicitly says the quadratic connection term is not used below. It cannot replace the second Ricci calculation. The retrieved primary record lists no journal reference. |
| Rohit Dhormare, *Scale-Invariant Non-Collapse from Raychaudhuri Comparison*, [arXiv:2610.07879v1](https://arxiv.org/html/2610.07879v1), **6 October 2026**; [DOI](https://doi.org/10.1016/j.physletb.2026.140987). | Reuses the retained §§2–7 review; fresh passages include (2)–(9), Theorem 1 (38)–(47), and the summary. Publisher PDF retrieval covered only its first **2 of 5** physical pages: normalization/shear, Raychaudhuri (3.1), volume evolution (4.1)–(4.3), and comparison footnotes. | Already present in the earlier area review. Supports the shear-contraction and pulled-back volume conventions, and distinguishes transported volumes from spatial-ball noncollapse. Primary arXiv metadata reports Physics Letters B 882, November 2026, 140987; the exact online-publication day was not established. No model constant is imported. |

For an independent sign check, Carroll’s [*Lecture Notes on General Relativity*, arXiv:gr-qc/9712019v1](https://arxiv.org/html/gr-qc/9712019v1) was read at the curvature convention (3.92)–(3.97), Ricci contraction (3.130)–(3.132), and connection/curvature variation (4.63)–(4.69). This is background geometry, not a recent development.

## Intrinsic Ricci: keep the complete quadratic connection contribution

The following formulas are an explicit derivation in the registered shell convention, supported by the standard identities above; they are not numerical estimates supplied by those papers. Work after the physical shell embedding has been included in h_e, on one fixed shell manifold. Let ∇ be the Levi-Civita connection of h0. All subscripts 1,2 below mean amplitude derivatives, with h_e=h0+e h1+(e²/2)h2+o(e²).

Define S_abd(k)=∇_a k_bd+∇_b k_ad−∇_d k_ab and the connection-difference tensor A_e=∇^(h_e)−∇. Then

```
A1^c_ab = (1/2) h0^cd S_abd(h1),
A2^c_ab = (1/2) h0^cd S_abd(h2) − h1^cd S_abd(h1),

Ric(h_e)_ab = Ric(h0)_ab
 + ∇_c A_e^c_ab − ∇_b A_e^c_ac
 + A_e^c_cd A_e^d_ab − A_e^c_bd A_e^d_ac,

Ric2_ab = ∇_c A2^c_ab − ∇_b A2^c_ac
 + 2[A1^c_cd A1^d_ab − A1^c_bd A1^d_ac].
```

Here h1^cd is raised using h0 twice. These follow directly from [Nakamura (42)–(53)](https://arxiv.org/html/gr-qc/0410024v2#S3.E42), with index order converted to the Ricci convention giving Ric(h0)=3h0/rho_b². Writing everything using the background connection avoids silently commuting derivatives and losing curvature terms. In particular A2 has an inverse-metric variation, and Ric2 has a factor-two quadratic connection term. A first-order factorization that drops that term cannot evaluate Ric2. The displayed expressions must still be bounded using the actual fixed-gauge shell traces and frame conversions.

## Why second observer variations cancel at zero

Set λ=3/rho_b², E_e=Ric(h_e)−λh_e, V_e=∂q X_e, with the **actual** proper-clock normalization h_e(V_e,V_e)=−rho_b². Since E0 vanishes as an intrinsic tensor field,

```
R_e = Ric(h_e)(V_e,V_e) = −3 + E_e(X_e)(V_e,V_e).

xi1 = ∂e X_e|0,
eta1 = ∇_(∂e) V_e|0 = ∇_(V0) xi1,

R2 = E2(V0,V0)
   + 2(∇_(xi1) E1)(V0,V0)
   + 4 E1(eta1,V0).
```

This is the background-covariant product/chain rule for the exact expression above. Terms involving xi2, eta2, or second spatial derivatives of E0 multiply E0 or its derivatives, all of which vanish. It is an algebraic simplification at e=0; E1 and the first observer corrections remain. The covariant velocity variation contains the background connection term: in a fixed chart eta1^a=∂e V_e^a|0+Γ0^a_bc xi1^b V0^c. Substituting a bare coordinate velocity derivative into the covariant formula would omit that term.

The relational-observable construction in [Fröb–Lima §2, §3](https://arxiv.org/pdf/2303.16218v2) supports keeping metric, observer prescription and transported labels together. Its coordinate-function map is inverse to the observer-position map used here, and its order-two coefficients are not second derivatives. Convert these conventions before comparing signs or factors. Its past-infinity initial conditions cannot replace the retained fixed physical anchor or the proof of vanishing initial observer variations.

A useful independent control is the toy family h_e=exp(2ce)h0 with normalized V_e=exp(−ce)V0. Its exact R_e=−3exp(−2ce), hence R2=−12c². The formula gives E2(V0,V0)=12c² and 4E1(eta1,V0)=−24c². This tests the proper-clock velocity term and factor four; it is not an admissible replacement for the registered exact family. A pure-diffeomorphism Einstein family gives a second control: its exact E_e and R2 vanish when fields and observers are transported together.

## Proper-area transfer and what the sources do not close

Using the retained intrinsic geodesic congruence, direct shear contraction and q orientation, the background is theta0=3tanh(q), Sigma0=0, R0=−3. The twice-differentiated Raychaudhuri equation is

```
theta2_q + 2tanh(q)theta2
  = −(2/3)theta1² − 2 tr(Sigma1²) − R2.
```

The integrating factor cosh(q)² and a second integration give exactly

```
L2_R(q,n) = −∫_0^q cosh(s)² R2(s,n)
                         [tanh(q)−tanh(s)] ds.
```

This confirms the sign, factor and kernel organization in the retained [quadratic proper-area note](../HDBLAST_QUADRATIC_PROPER_AREA_20261106/QUADRATIC_AREA_BOUND.md). The full-family quiet interval remains q≤ell/4. With its fixed positive Y² weight and already evaluated kernel constant, a bound Cr≥sup_q|<R2>_w| would give |<L2_R>_w|≤C_R Cr, C_R<8.614×10^−13. These numerical inputs are reused project results, not recalculated or inferred from an external comparison theorem.

The 2026 noncollapse paper assumes uniform curvature, expansion and shear bounds on a region where the congruence exists. Its own footnotes distinguish the auxiliary comparison barrier from a physical focal time. It cannot supply those missing family-wide assumptions from zero-amplitude jets. Likewise, Kim’s all-orders displacement formalism has a real-analytic setup and does not evaluate the varying-metric family’s remainder. Neither result changes the need for the actual reverse-time join, shell identification, h2 C2/h1 C3 bounds, first observer position/velocity bounds, and their numerical evaluation.

This focused check found no inspected 2025–2026 result that closes those model-specific requirements or supplies empirical evidence for the HDblast hypothesis. Related geodesic or cosmology papers are not validation of its higher-dimensional origin mechanism. The selected odd moment’s second derivative already vanishes under the retained symmetry premises; the nontrivial new target remains the full proper second scalar or weighted intrinsic-curvature/area input. Zero-amplitude cancellation cannot establish an amplitude radius, a positive denominator throughout a nonzero parameter interval, or an O(e³) remainder from only C2 dependence.

## Provenance and stopping point

Fresh searches covered second Ricci/connection variation, proper-observer covariance, recent geodesic-deviation methods and Raychaudhuri comparison. Search abstracts only guided source selection; the consequential claims above use in-body text. Primary submission histories control version dates: the paper index’s “updated” field sometimes differs by a day or by years and must not be treated as a revision date.

The arXiv HTML for Fröb–Lima returned HTTP 404; its exact v2 PDF was then read in memory. The 2026 publisher PDF was deliberately limited to two pages. Fourteen Firecrawl credits were explicitly reported across disjoint calls; fourteen research-tool calls supplied no usage receipt, so the total cost is **unknown**, not fourteen credits. No broad crawl or autonomous research job was used. Exact receipts, selected passages, date discrepancies and exclusions are in [LITERATURE_SOURCES.json](LITERATURE_SOURCES.json).

Collection stopped because the method-level questions are adequately supported. Further literature collection would not evaluate the remaining registered-field/observer constants. No exhaustive search, novelty theorem, journal acceptance for this work, or physical confirmation is asserted.
