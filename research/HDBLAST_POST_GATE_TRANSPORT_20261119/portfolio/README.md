# Which archived interaction can be used?

The archived shell model gives a source for massless Dirac pairs. Its finite bath realizes the vacuum spectral density \(K_q(\omega)=(\omega^2-q^2)/(8\pi^2)\) at order \(y^2\). The interaction is proportional to a switch \(f\), and its static counterterm is proportional to \(f^2\). Both vanish when the switch closes. The retained particle occupations then redshift freely on the prescribed shell. A nonzero post-switch collision rate requires a specified additional interaction. The spectral density alone does not determine an elastic cross section. [Hamiltonian and switch](https://github.com/maldonado-research/HDblast/blob/49756c922760880a9f5da0e4da6f2d8cce36738a/research/HDBLAST_CAUSAL_YUKAWA_ENERGY_20261119/solver/METHODS.md) · [Finite-bath results and limitations](https://github.com/maldonado-research/HDblast/blob/49756c922760880a9f5da0e4da6f2d8cce36738a/research/HDBLAST_QUANTUM_MEMORY_20261119/README.md)

## An independent angular-integration benchmark

The D-Blast 3 archive contains a worked DMDE momentum-transfer cross section for nonrelativistic massive dark matter interacting through an assumed dark-photon mediator. The public February release gives a numerical benchmark. In its Born approximation, put \(u=\cos\theta\), \(v=v_{\rm rel}/c\) and \(x=m_\chi v/m_{A'}\). The two ways of computing the transport cross section are

\[
 2\pi\int_{-1}^{1}(1-u)
 \frac{\alpha_D^2 m_\chi^2\,du}
      {[m_{A'}^2+m_\chi^2v^2(1-u)/2]^2}
 =
 \frac{8\pi\alpha_D^2}{m_\chi^2v^4}
 \left[\log(1+x^2)-\frac{x^2}{1+x^2}\right].
\]

Using \(m_\chi=5.02\) GeV, \(m_{A'}=4.33\) MeV and \(\alpha_D=1.4\times10^{-4}\), direct angular quadrature and the closed form give the following values. The conversion uses \(1\,\mathrm{GeV}^{-2}=0.389379\times10^{-27}\,\mathrm{cm}^2\) and \(1\,\mathrm{GeV}=1.78266\times10^{-24}\,\mathrm g\) in natural units. The release's rounded values appear in its ZIP member ending in `docs/DMDE_Triplepoint_CosmoGate_NextLevel_RefereeCore_v0.3.md`. [Zenodo v0.4.1](https://doi.org/10.5281/zenodo.18613652)

| Relative speed (km/s) | Recalculated \(\sigma_T/m_\chi\) (cm²/g) | Published rounded value (cm²/g) |
| ---: | ---: | ---: |
| 10 | 0.76676 | 0.77 |
| 200 | 0.40613 | 0.41 |
| 1000 | 0.012589 | 0.013 |

This checks the angular measure, the transport weight \(1-u\), and unit conversion against a published example. Its masses, mediator, coupling and angular law apply to the stated dark-matter model. The current massless shell daughters have no identified DMDE charge or mediator.

The DMDE October transport audit also exercised a compiled neutrino collision routine at one 21-node test point with a stopped-muon source. It recorded an endpoint truncation that lost much of the emitted number and energy even as its grid was refined. Its relevant lesson for this source is to check daughter moments against the full kinematic interval independently of grid convergence. That neutrino fixture supplies neither a shell-particle species identification nor an elastic coefficient for it. [Transport-source audit](https://github.com/maldonado-research/dmde-research/blob/34c6347fab1018f9e4bd707c5270270bbabe5ba6/research/rounds/2026-10-01-transport-source-audit/README.md)

## Catalog coverage

On 19 November 2026 the public GitHub account listed eight repositories. HDblast main pointed at merge 49756c9 (PR 172). The HDBLAST concept record 17088132 led to version 24374531. The supplied Zenodo creator search returned nine concept series and 81 version records across four pages, including 42 in the main HDBLAST series. The ORCID query for 0009-0009-3937-6527 returned 71 versions, all among the 81. Screening covered catalog metadata, the current HDBLAST release, relevant DMDE source files, and targeted D-Blast 3 archive text on scattering and transport. No additional post-switch interaction for the emitted HDBLAST fermions was identified in those sources. [Author catalog](https://zenodo.org/search?q=metadata.creators.person_or_org.name%3A%22Maldonado%2C%20Ricardo%22) · [GitHub repositories](https://github.com/maldonado-research?tab=repositories) · [HDBLAST release](https://doi.org/10.5281/zenodo.24374531)
