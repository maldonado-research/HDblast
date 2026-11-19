# What a pair-production spectrum can determine

## From fermion pairs to a response function

Boyanovsky and Herring calculate the spectral density of a scalar coupled to a pair of massless Dirac fermions. Their Eqs. IV.1–IV.3 give the threshold at energy ω = q, where q is the magnitude of the comoving wave vector. Their density is (ω² − q²)/(4π) above threshold and the integration measure is dω/(2π). Combining them gives the weight used here, K_f(ω) dω = (ω² − q²) dω/(8π²), with the subscript identifying the fermion channel. This fixes the nonlocal kernel at leading order in the squared Yukawa coupling under the stated vacuum and massless assumptions. Fleming, Roura and Hu show how a set of linearly coupled oscillators realizes a causal sine kernel with a prescribed positive spectral measure. Such an oscillator representation of a two-point function supplies no higher-order fermion statistics. [Pair spectral density and normalization](https://arxiv.org/html/1904.12343#S4.E3) · [Oscillator spectral measure and dissipation, Eqs. 4–6](https://arxiv.org/html/1004.1603#S2.E4)

## Why the dynamic calibration is independent

In a subtracted dispersion relation, a real local term can change the response without changing absorption on the frequency axis. Kniehl's Eq. 8 displays the subtraction value explicitly and explains when more subtractions are needed. An experimental example comes from third-harmonic spectroscopy. Lucarini, Saarinen and Peiponen used an independently known anchor value to improve a reconstruction from a limited spectral band. In HDblast the chosen static restoring force and kinetic normalization are model inputs. Their numerical values have yet to be fixed by an independent observable or microscopic matching calculation. [Subtracted dispersion relation](https://arxiv.org/html/hep-ph/9607255#S2.E8) · [Independent optical anchor and experimental application](https://arxiv.org/html/cond-mat/0306139#S3)

For example, freeze the switching factor and strip off the coupling and external scalar normalization. At fixed q > 0 write the bath contribution at imaginary frequency as B(x) = ∫ 2ω K_f(ω)/(ω² + x) dω, where x = ν² and the integration runs from q to the cutoff Λ. Its value at zero grows quadratically with Λ, while its first derivative diverges logarithmically. Subtracting both gives the exact identity

$$
 B(x)-B(0)-xB'(0)
 = x^2\int_q^\Lambda
 \frac{2K_f(\omega)}{\omega^3(\omega^2+x)}\,d\omega .
$$

For x ≥ 0 this curvature remainder is nonnegative. The omitted continuum tail above Λ ≥ q is at most x²/(8π²Λ²). These statements refer to the bath with the sign convention shown. They concern a frozen linear response at each momentum. Deriving the corresponding stress on an evolving shell requires a time-dependent action and its variation.

## Independent checks and limits

Del Río and colleagues explicitly find separate potential, field-strength and curvature counterterms for a Yukawa-coupled scalar and Dirac field in an expanding universe. Herring and Boyanovsky, in a February 2026 revision, compute an additional energy contribution from excited fermions in a time-dependent flat-space condensate. Both identify physical terms that a static potential alone misses. A recent compact-object example shows why one must also check the large-frequency growth condition before deducing bounds on a dynamical response from positive absorption. The high-frequency assumptions and gap used there belong to that tidal-response problem. [FLRW Yukawa counterterms, Sec. V](https://arxiv.org/html/1703.00908#S5) · [Condensate energy, Eq. III.82](https://arxiv.org/html/2511.22465v2#S3.E82) · [Kehagias and Riotto, October 2026 revision, Sec. 2](https://arxiv.org/html/2607.26686v2#S2)

## Project lineage

The November 19, 2026 catalog check found 79 Zenodo record versions under the exact creator name, 69 with the author's ORCID, nine latest records under the name, and eight public GitHub repositories. The two directly relevant HDblast releases contain the spectral normalization and the cutoff-dependent gravity histories. D-Blast 3 folders 144 and 153 supply earlier scalar-pulse and shell-matter assumptions. The author's DMDE work offers a numerical check that interrogates the implemented source directly, since agreement with an integrated quantity can conceal a source error. The Unified Theory fermion-portal review gives a second example of testing an evolution with an independently implemented solver. A shell coefficient would need its own independent measurement or matching calculation. [Zenodo creator catalog](https://zenodo.org/search?q=metadata.creators.person_or_org.name%3A%22Maldonado%2C%20Ricardo%22) · [Public repositories](https://github.com/maldonado-research?tab=repositories) · [Pair kernel and energy](https://zenodo.org/records/24358059) · [Gravity histories](https://zenodo.org/records/24360096) · [DMDE source check](https://github.com/maldonado-research/dmde-research/blob/34c6347fab1018f9e4bd707c5270270bbabe5ba6/provider/docs/DMDE_v0920_SOURCE_OPERATOR_CLOSURE_GATE.md) · [Independent fermion-portal review](https://github.com/maldonado-research/Unified-Theory-of-Everything/blob/c284afaa99c9a1e1fbde92a5547892f0ff9a2ce4/checkpoint_2026_09_06/fermion_portal/independent_review.py)
